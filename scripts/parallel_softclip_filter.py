#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parallel_softclip_filter.py

rDNA-specific soft-clip analysis and filtering for fungal short reads.

Optimized for:
- Fungal 45S rDNA (highly repetitive, low MAPQ expected)
- Short reads (100-300bp) after fastp preprocessing
- ITS region length variation (real biology vs artifacts)

Key Logic:
- Uses sequence composition, base quality, and read-pair consistency instead of MAPQ
- Cluster-supported clips = likely real ITS variation (KEEP)
- Low-quality isolated clips = likely sequencing artifacts (REMOVE/TRIM)

Modes:
  stats : compute detailed quality metrics and auto-threshold
  full  : remove reads with low-quality artifactual soft-clips
  trim  : trim low-quality soft-clipped ends while preserving real variation
  igv   : create visualization BAMs with quality tags

Usage examples:
  # 1) Analyze rDNA soft-clip patterns:
  parallel_softclip_filter.py -i input.bam -o outdir -m auto -t 16 --mode stats

  # 2) Remove artifactual clips:
  parallel_softclip_filter.py -i input.bam -o outdir -m auto -t 16 --mode full

  # 3) Trim artifactual clips (recommended):
  parallel_softclip_filter.py -i input.bam -o outdir -m auto -t 16 --mode trim

  # 4) IGV visualization:
  parallel_softclip_filter.py -i input.bam -o outdir -m auto -t 16 --mode igv
"""
import os
import sys
import argparse
import statistics
import itertools
from collections import Counter, defaultdict
from multiprocessing import Pool
import pysam


def parse_args():
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="rDNA-specific soft-clip filter for fungal short reads"
    )
    p.add_argument('-i','--input', required=True, help="Input BAM file")
    p.add_argument('-o','--outdir', required=True, help="Output directory")
    p.add_argument('-m','--min_softclip', default='auto',
                   help="Always uses automatic quality-based thresholds")
    p.add_argument('-t','--threads', type=int, default=1, help="Parallel threads")
    p.add_argument('--mode', choices=['stats','full','trim','igv'], default='stats',
                   help="Operation mode")
    return p.parse_args()


def assign_buckets(contigs, lengths, n_buckets):
    """Evenly assign contigs to n buckets by total reference length."""
    buckets = [[] for _ in range(n_buckets)]
    sizes = [0] * n_buckets
    for ctg, ln in sorted(zip(contigs, lengths), key=lambda x: x[1], reverse=True):
        i = min(range(n_buckets), key=lambda j: sizes[j])
        buckets[i].append(ctg)
        sizes[i] += ln
    return buckets


def calculate_sequence_complexity(seq):
    """Calculate linguistic complexity of sequence (0-1 scale)."""
    if len(seq) < 4:
        return 0.0

    # Count unique k-mers (k=3)
    kmers = set()
    for i in range(len(seq) - 2):
        kmers.add(seq[i:i+3])

    # Complexity = observed unique kmers / theoretical maximum
    max_possible = min(64, len(seq) - 2)  # 4^3 = 64 possible 3-mers
    return len(kmers) / max_possible if max_possible > 0 else 0.0


def analyze_homopolymers(seq):
    """Find longest homopolymer run in sequence."""
    if not seq:
        return 0
    max_run = max(len(list(group)) for _, group in itertools.groupby(seq))
    return max_run


def calculate_gc_content(seq):
    """Calculate GC content (0-1 scale)."""
    if not seq:
        return 0.5
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq)


def calculate_clip_quality_score(read, clip_start, clip_end, clip_seq):
    """
    Calculate quality score for a soft-clipped region (0-100 scale).

    For rDNA short reads, focuses on:
    - Sequence composition (GC, complexity, homopolymers)
    - Base quality in clipped region
    - Read pair consistency
    """
    score = 0

    # 1. Sequence composition analysis (50% weight)
    if clip_seq:
        # GC content (should be reasonable, not extreme)
        gc_content = calculate_gc_content(clip_seq)
        if 0.3 <= gc_content <= 0.7:  # Reasonable range for fungal rDNA
            score += 20
        elif 0.2 <= gc_content <= 0.8:  # Acceptable range
            score += 10

        # Sequence complexity (avoid simple repeats)
        complexity = calculate_sequence_complexity(clip_seq)
        if complexity >= 0.6:  # High complexity
            score += 20
        elif complexity >= 0.4:  # Medium complexity
            score += 10

        # Homopolymer runs (sequencing artifacts common in homopolymers)
        max_homopoly = analyze_homopolymers(clip_seq)
        if max_homopoly <= 4:  # Short homopolymers OK
            score += 10
        elif max_homopoly <= 6:  # Medium homopolymers
            score += 5

    # 2. Base quality in clipped region (30% weight)
    if read.query_qualities and clip_start < len(read.query_qualities):
        end_pos = min(clip_end, len(read.query_qualities))
        clip_qualities = read.query_qualities[clip_start:end_pos]
        if clip_qualities:
            mean_qual = statistics.mean(clip_qualities)
            if mean_qual >= 25:  # High quality
                score += 30
            elif mean_qual >= 20:  # Medium quality
                score += 20
            elif mean_qual >= 15:  # Low but acceptable
                score += 10

    # 3. Read pair information (20% weight)
    if read.is_paired:
        # Proper pair orientation
        if read.is_proper_pair:
            score += 10

        # Reasonable insert size (adjust for your library)
        if read.template_length != 0:
            insert_size = abs(read.template_length)
            if 100 <= insert_size <= 800:  # Typical short-read library range
                score += 10
            elif 50 <= insert_size <= 1200:  # Extended acceptable range
                score += 5

    return min(score, 100)  # Cap at 100


def stats_bucket(args):
    """
    Collect detailed clip statistics with quality scores.
    Returns: total_reads, clip_records, read_stats
    """
    bucket, bam_path, min_len = args
    total = 0
    clip_records = []  # (pos, length, is_left, quality_score, read_id)
    read_stats = {
        'clipped_reads': 0,
        'total_mapq': 0,
        'left_only': 0,
        'right_only': 0,
        'both_ends': 0,
        'high_quality_clips': 0,
        'low_quality_clips': 0
    }

    bam = pysam.AlignmentFile(bam_path, 'rb')

    for ctg in bucket:
        for read in bam.fetch(ctg):
            total += 1
            read_stats['total_mapq'] += read.mapping_quality

            if not read.cigartuples:
                continue

            # Analyze soft-clips
            left_clip = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
            right_clip = read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0

            if not (left_clip or right_clip):
                continue

            read_stats['clipped_reads'] += 1

            # Classify clipping pattern
            if left_clip and right_clip:
                read_stats['both_ends'] += 1
            elif left_clip:
                read_stats['left_only'] += 1
            else:
                read_stats['right_only'] += 1

            # Analyze left clip
            if left_clip >= min_len:
                clip_seq = read.query_sequence[:left_clip] if read.query_sequence else ""
                quality_score = calculate_clip_quality_score(read, 0, left_clip, clip_seq)

                clip_records.append((
                    read.reference_start, left_clip, True, quality_score, read.query_name
                ))

                if quality_score >= 60:
                    read_stats['high_quality_clips'] += 1
                else:
                    read_stats['low_quality_clips'] += 1

            # Analyze right clip
            if right_clip >= min_len and read.query_sequence:
                clip_seq = read.query_sequence[-right_clip:]
                seq_len = len(read.query_sequence)
                quality_score = calculate_clip_quality_score(
                    read, seq_len - right_clip, seq_len, clip_seq
                )

                clip_records.append((
                    read.reference_end, right_clip, False, quality_score, read.query_name
                ))

                if quality_score >= 60:
                    read_stats['high_quality_clips'] += 1
                else:
                    read_stats['low_quality_clips'] += 1

    bam.close()
    return total, clip_records, read_stats


def compute_stats(bam_path, buckets, threads, min_initial=1):
    """Aggregate statistics with rDNA-specific quality analysis."""
    tasks = [(buckets[i], bam_path, min_initial) for i in range(len(buckets))]

    with Pool(threads) as pool:
        results = pool.map(stats_bucket, tasks)

    # Aggregate results
    total_reads = sum(r[0] for r in results)
    all_clips = []
    for r in results:
        all_clips.extend(r[1])

    # Aggregate read stats
    combined_stats = defaultdict(int)
    for _, _, read_stats in results:
        for key, value in read_stats.items():
            combined_stats[key] += value

    # Analyze clustering with quality awareness
    position_clips = defaultdict(list)  # (pos, is_left) -> [quality_scores]
    quality_scores = [clip[3] for clip in all_clips]
    clip_lengths = [clip[1] for clip in all_clips]

    for pos, length, is_left, quality, read_id in all_clips:
        position_clips[(pos, is_left)].append((quality, length))

    # Define clusters: positions with ≥2 reads AND average quality ≥ 50
    clusters = {}
    high_quality_clusters = {}

    for (pos, is_left), quality_lengths in position_clips.items():
        if len(quality_lengths) >= 2:  # Cluster support
            avg_quality = statistics.mean([ql[0] for ql in quality_lengths])
            avg_length = statistics.mean([ql[1] for ql in quality_lengths])

            clusters[(pos, is_left)] = {
                'count': len(quality_lengths),
                'avg_quality': avg_quality,
                'avg_length': avg_length,
                'lengths': [ql[1] for ql in quality_lengths]
            }

            if avg_quality >= 50:  # High-quality cluster
                high_quality_clusters[(pos, is_left)] = clusters[(pos, is_left)]

    # Calculate thresholds
    if clip_lengths:
        # Use multiple approaches for threshold
        median_length = statistics.median(clip_lengths)
        pct75_length = statistics.quantiles(clip_lengths, n=4)[2] if len(clip_lengths) >= 4 else median_length

        # For rDNA, use a more conservative approach
        # Focus on quality rather than just length
        high_qual_lengths = [length for _, length, _, quality, _ in all_clips if quality >= 60]
        if high_qual_lengths:
            hq_median = statistics.median(high_qual_lengths)
            auto_thresh = max(10, min(int(hq_median * 0.8), 50))  # Conservative for short reads
        else:
            auto_thresh = max(10, int(pct75_length * 0.6))
    else:
        auto_thresh = 15  # Default for short reads

    # Count clips by category
    total_clips = len(all_clips)
    clustered_clips = sum(1 for clip in all_clips
                         if (clip[0], clip[2]) in clusters and clip[1] >= auto_thresh)
    isolated_clips = sum(1 for clip in all_clips
                        if (clip[0], clip[2]) not in clusters and clip[1] >= auto_thresh)

    return {
        'total_reads': total_reads,
        'total_clips': total_clips,
        'clipped_reads': combined_stats['clipped_reads'],
        'clusters': len(clusters),
        'high_quality_clusters': len(high_quality_clusters),
        'clustered_clips': clustered_clips,
        'isolated_clips': isolated_clips,
        'high_quality_clips': combined_stats['high_quality_clips'],
        'low_quality_clips': combined_stats['low_quality_clips'],
        'left_only': combined_stats['left_only'],
        'right_only': combined_stats['right_only'],
        'both_ends': combined_stats['both_ends'],
        'avg_mapq': combined_stats['total_mapq'] / total_reads if total_reads else 0,
        'auto_thresh': auto_thresh,
        'median_length': statistics.median(clip_lengths) if clip_lengths else 0,
        'mean_quality': statistics.mean(quality_scores) if quality_scores else 0,
        'clusters_dict': clusters,
        'high_quality_clusters_dict': high_quality_clusters
    }


def filter_bucket(args):
    """Filter/trim reads based on quality-aware clustering with detailed statistics."""
    bucket, bam_path, outdir, threshold, clusters, mode, idx = args
    out_sam = os.path.join(outdir, f'filter_{mode}_{idx}.sam')

    # Statistics counters
    stats = {
        'reads_processed': 0,
        'reads_with_clips': 0,
        'reads_removed': 0,          # For full mode
        'reads_modified': 0,         # For trim mode
        'reads_kept_unchanged': 0,
        'left_clips_analyzed': 0,
        'right_clips_analyzed': 0,
        'left_clips_removed': 0,     # Full: entire read removed, Trim: clip trimmed
        'right_clips_removed': 0,
        'left_clips_preserved': 0,   # High-quality or clustered clips kept
        'right_clips_preserved': 0,
        'total_bp_trimmed': 0,       # Only for trim mode
        'total_reads_trimmed_left': 0,
        'total_reads_trimmed_right': 0
    }

    bam = pysam.AlignmentFile(bam_path, 'rb')

    with open(out_sam, 'w') as out:
        for ctg in bucket:
            for read in bam.fetch(ctg):
                stats['reads_processed'] += 1

                if not read.cigartuples:
                    out.write(read.to_string() + "\n")
                    continue

                left_clip = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
                right_clip = read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0

                if left_clip or right_clip:
                    stats['reads_with_clips'] += 1

                # Analyze clip quality and clustering
                remove_left = False
                remove_right = False

                if left_clip >= threshold:
                    stats['left_clips_analyzed'] += 1
                    clip_seq = read.query_sequence[:left_clip] if read.query_sequence else ""
                    quality_score = calculate_clip_quality_score(read, 0, left_clip, clip_seq)

                    cluster_key = (read.reference_start, True)
                    if quality_score < 50 and cluster_key not in clusters:
                        remove_left = True
                        stats['left_clips_removed'] += 1
                    else:
                        stats['left_clips_preserved'] += 1

                if right_clip >= threshold and read.query_sequence:
                    stats['right_clips_analyzed'] += 1
                    seq_len = len(read.query_sequence)
                    clip_seq = read.query_sequence[-right_clip:]
                    quality_score = calculate_clip_quality_score(
                        read, seq_len - right_clip, seq_len, clip_seq
                    )

                    cluster_key = (read.reference_end, False)
                    if quality_score < 50 and cluster_key not in clusters:
                        remove_right = True
                        stats['right_clips_removed'] += 1
                    else:
                        stats['right_clips_preserved'] += 1

                if mode == 'full':
                    # Skip reads with low-quality clips
                    if remove_left or remove_right:
                        stats['reads_removed'] += 1
                        continue
                    else:
                        stats['reads_kept_unchanged'] += 1
                        out.write(read.to_string() + "\n")

                elif mode == 'trim':
                    # Trim only low-quality clips
                    if remove_left or remove_right:
                        stats['reads_modified'] += 1
                        seq = read.query_sequence or ''
                        qual = read.query_qualities or []
                        new_cigar = list(read.cigartuples)

                        trim_left = left_clip if remove_left else 0
                        trim_right = right_clip if remove_right else 0

                        if trim_left:
                            stats['total_reads_trimmed_left'] += 1
                            stats['total_bp_trimmed'] += trim_left
                        if trim_right:
                            stats['total_reads_trimmed_right'] += 1
                            stats['total_bp_trimmed'] += trim_right

                        if trim_left or trim_right:
                            # Update sequence and quality
                            end_pos = len(seq) - trim_right if trim_right else len(seq)
                            read.query_sequence = seq[trim_left:end_pos]
                            if qual:
                                read.query_qualities = qual[trim_left:end_pos]

                            # Update CIGAR
                            if trim_left and new_cigar[0][0] == 4:
                                new_left = left_clip - trim_left
                                if new_left > 0:
                                    new_cigar[0] = (4, new_left)
                                else:
                                    new_cigar = new_cigar[1:]

                            if trim_right and new_cigar and new_cigar[-1][0] == 4:
                                new_right = right_clip - trim_right
                                if new_right > 0:
                                    new_cigar[-1] = (4, new_right)
                                else:
                                    new_cigar = new_cigar[:-1]

                            read.cigartuples = new_cigar
                    else:
                        stats['reads_kept_unchanged'] += 1

                    out.write(read.to_string() + "\n")

    bam.close()
    return stats


def igv_bucket(args):
    """Create IGV BAM with quality and cluster tags."""
    bucket, bam_path, outdir, threshold, clusters, idx = args
    out_bam = os.path.join(outdir, f'igv_{idx}.bam')

    bam_in = pysam.AlignmentFile(bam_path, 'rb')
    bam_out = pysam.AlignmentFile(out_bam, 'wb', template=bam_in)

    for ctg in bucket:
        for read in bam_in.fetch(ctg):
            if not read.cigartuples or not any(op == 4 for op, _ in read.cigartuples):
                continue

            left_clip = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
            right_clip = read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0

            # Add quality tags
            if left_clip > 0:
                clip_seq = read.query_sequence[:left_clip] if read.query_sequence else ""
                quality_score = calculate_clip_quality_score(read, 0, left_clip, clip_seq)
                cluster_key = (read.reference_start, True)

                if quality_score >= 60:
                    read.set_tag('XL', 'high_quality')
                elif quality_score >= 30:
                    read.set_tag('XL', 'medium_quality')
                else:
                    read.set_tag('XL', 'low_quality')

                if cluster_key in clusters:
                    read.set_tag('CL', 'clustered')
                else:
                    read.set_tag('CL', 'isolated')

            if right_clip > 0:
                seq_len = len(read.query_sequence) if read.query_sequence else 0
                clip_seq = read.query_sequence[-right_clip:] if read.query_sequence else ""
                quality_score = calculate_clip_quality_score(
                    read, seq_len - right_clip, seq_len, clip_seq
                )
                cluster_key = (read.reference_end, False)

                if quality_score >= 60:
                    read.set_tag('XR', 'high_quality')
                elif quality_score >= 30:
                    read.set_tag('XR', 'medium_quality')
                else:
                    read.set_tag('XR', 'low_quality')

                if cluster_key in clusters:
                    read.set_tag('CR', 'clustered')
                else:
                    read.set_tag('CR', 'isolated')

            read.set_tag('QT', threshold)  # Quality threshold used
            bam_out.write(read)

    bam_in.close()
    bam_out.close()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Build contig buckets
    bam = pysam.AlignmentFile(args.input, 'rb')
    buckets = assign_buckets(bam.references, bam.lengths, args.threads)
    bam.close()

    # ─── STATS MODE ──────────────────────────────────────────────────────────
    if args.mode == 'stats':
        stats = compute_stats(args.input, buckets, args.threads)

        print(f"Soft-clip analysis: {stats['clipped_reads']}/{stats['total_reads']} reads ({stats['clipped_reads']/stats['total_reads']*100:.1f}%) with clips")
        print(f"Auto-threshold: {stats['auto_thresh']}bp (median: {stats['median_length']:.1f}bp)")

        sys.exit(0)

    # ─── FILTER/TRIM MODES ───────────────────────────────────────────────────
    if args.mode in ('full', 'trim'):
        stats = compute_stats(args.input, buckets, args.threads)

        threshold = stats['auto_thresh']
        print(f"Using automatic quality-based threshold: {threshold}bp")

        clusters = stats['high_quality_clusters_dict']
        print(f"Using {len(clusters)} high-quality cluster positions for filtering")

        tasks = [
            (buckets[i], args.input, args.outdir, threshold, clusters, args.mode, i)
            for i in range(len(buckets))
        ]

        with Pool(args.threads) as pool:
            bucket_stats = pool.map(filter_bucket, tasks)

        # Aggregate statistics from all buckets
        total_stats = {
            'reads_processed': sum(s['reads_processed'] for s in bucket_stats),
            'reads_with_clips': sum(s['reads_with_clips'] for s in bucket_stats),
            'reads_removed': sum(s['reads_removed'] for s in bucket_stats),
            'reads_modified': sum(s['reads_modified'] for s in bucket_stats),
            'reads_kept_unchanged': sum(s['reads_kept_unchanged'] for s in bucket_stats),
            'left_clips_analyzed': sum(s['left_clips_analyzed'] for s in bucket_stats),
            'right_clips_analyzed': sum(s['right_clips_analyzed'] for s in bucket_stats),
            'left_clips_removed': sum(s['left_clips_removed'] for s in bucket_stats),
            'right_clips_removed': sum(s['right_clips_removed'] for s in bucket_stats),
            'left_clips_preserved': sum(s['left_clips_preserved'] for s in bucket_stats),
            'right_clips_preserved': sum(s['right_clips_preserved'] for s in bucket_stats),
            'total_bp_trimmed': sum(s['total_bp_trimmed'] for s in bucket_stats),
            'total_reads_trimmed_left': sum(s['total_reads_trimmed_left'] for s in bucket_stats),
            'total_reads_trimmed_right': sum(s['total_reads_trimmed_right'] for s in bucket_stats)
        }

        # Report detailed statistics
        if args.mode == 'full':
            print(f"Filtered: {total_stats['reads_removed']:,}/{total_stats['reads_processed']:,} reads ({total_stats['reads_removed']/total_stats['reads_processed']*100:.1f}%) removed")
        elif args.mode == 'trim':
            print(f"Trimmed: {total_stats['reads_modified']:,}/{total_stats['reads_processed']:,} reads ({total_stats['reads_modified']/total_stats['reads_processed']*100:.1f}%) modified")

        print(f"Preservation rate: {(total_stats['left_clips_preserved'] + total_stats['right_clips_preserved'])/(total_stats['left_clips_analyzed'] + total_stats['right_clips_analyzed'])*100:.1f}%")
        total_preserved = total_stats['left_clips_preserved'] + total_stats['right_clips_preserved']
        total_analyzed = total_stats['left_clips_analyzed'] + total_stats['right_clips_analyzed']
        if total_analyzed > 0:
            print(f"  Overall preservation rate: {total_preserved}/{total_analyzed} ({total_preserved/total_analyzed*100:.1f}%)")

        print(f"\nFiltered {args.mode} SAM files created in {args.outdir}")
        sys.exit(0)

    # ─── IGV MODE ─────────────────────────────────────────────────────────────
    if args.mode == 'igv':
        stats = compute_stats(args.input, buckets, args.threads)

        threshold = stats['auto_thresh']
        clusters = stats['high_quality_clusters_dict']

        tasks = [
            (buckets[i], args.input, args.outdir, threshold, clusters, i)
            for i in range(len(buckets))
        ]

        with Pool(args.threads) as pool:
            pool.map(igv_bucket, tasks)

        print(f"IGV BAM files created with quality tags:")
        print(f"  XL/XR: left/right clip quality (high/medium/low)")
        print(f"  CL/CR: left/right cluster status (clustered/isolated)")
        print(f"  QT: quality threshold used ({threshold})")
        sys.exit(0)


if __name__ == '__main__':
    main()
