#!/usr/bin/env python3
"""
parallel_softclip_filter.py

Parallel soft-clip analysis and filtering with cluster support + percentile cutoff.

Modes:
  stats : compute and report detailed stats (total segments, cluster count,
          distribution stats, count≥threshold, percent≥threshold, auto-threshold)
  igv   : extract reads with ≥threshold soft-clips into per-bucket BAMs for IGV
  full  : remove reads whose soft-clips meet the chosen threshold at supported positions
  trim  : trim soft-clipped ends meeting the chosen threshold at supported positions


Usage examples:
  # 1) Stats with auto threshold:
  parallel_softclip_filter.py -i input.bam -o outdir -m auto -t 16 --mode stats

  # 2) Full removal at numeric threshold:
  parallel_softclip_filter.py -i input.bam -o outdir -m 108 -t 16 --mode full

  # 3) Trim with auto threshold:
  parallel_softclip_filter.py -i input.bam -o outdir -m auto -t 16 --mode trim

  # 4) Build IGV BAM chunks:
  parallel_softclip_filter.py -i input.bam -o outdir -m 30 -t 16 --mode igv
"""
import os
import sys
import argparse
import statistics
from collections import Counter
from multiprocessing import Pool
import pysam


def parse_args():
    """
    Parse command-line arguments:
      -i/--input       : input BAM file path
      -o/--outdir      : output directory for SAM/BAM chunks or stats
      -m/--min_softclip: 'auto' or integer soft-clip length cutoff
      -t/--threads     : number of parallel threads
      --mode           : one of 'stats','full','trim','igv'
    """
    p = argparse.ArgumentParser(description="Soft-clip filter/trimmer/stats/IGV export")
    p.add_argument('-i','--input', required=True, help="Input BAM file")
    p.add_argument('-o','--outdir', required=True, help="Output directory for chunks or stats")
    p.add_argument('-m','--min_softclip', required=True,
                   help="'auto' or integer soft-clip length cutoff")
    p.add_argument('-t','--threads', type=int, default=1, help="Parallel threads")
    p.add_argument('--mode', choices=['stats','full','trim','igv'], default='stats',
                   help="Operation mode: stats|full|trim|igv")
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


def stats_bucket(args):
    """
    Collect clip positions and lengths for stats. Returns:
      total reads,
      list of (pos,length) for all clips ≥ min_len,
      left-only count, right-only count, both-ends count,
      sum MAPQ all reads, sum MAPQ clipped reads
    """
    bucket, bam_path, min_len = args
    total = 0
    left_only = right_only = both = 0
    mapq_all = mapq_clip = 0
    records = []
    bam = pysam.AlignmentFile(bam_path, 'rb')
    for ctg in bucket:
        for r in bam.fetch(ctg):
            total += 1
            mapq_all += r.mapping_quality
            if not r.cigartuples:
                continue
            left = r.cigartuples[0][1] if r.cigartuples[0][0] == 4 else 0
            right = r.cigartuples[-1][1] if r.cigartuples[-1][0] == 4 else 0
            if left or right:
                if left >= min_len:
                    records.append((r.reference_start, left))
                    mapq_clip += r.mapping_quality
                if right >= min_len:
                    records.append((r.reference_end, right))
                    mapq_clip += r.mapping_quality
                if left and right:
                    both += 1
                elif left:
                    left_only += 1
                else:
                    right_only += 1
    bam.close()
    return total, records, left_only, right_only, both, mapq_all, mapq_clip


def compute_stats(bam_path, buckets, threads, min_initial=1):
    """
    Aggregate stats: cluster support + percentile, plus detailed counts.
    Returns a dict with all metrics including 'clusters_set'.
    """
    tasks = [(buckets[i], bam_path, min_initial) for i in range(len(buckets))]
    with Pool(threads) as p:
        results = p.map(stats_bucket, tasks)
    total = sum(r[0] for r in results)
    all_records = [rec for r in results for rec in r[1]]
    left_only = sum(r[2] for r in results)
    right_only = sum(r[3] for r in results)
    both = sum(r[4] for r in results)
    mapq_all = sum(r[5] for r in results)
    mapq_clip = sum(r[6] for r in results)

    # cluster support: positions with ≥2 reads
    pos_counts = Counter(pos for pos, _ in all_records)
    clusters = {pos for pos, c in pos_counts.items() if c >= 2}
    cluster_lengths = [ln for pos, ln in all_records if pos in clusters]

    if cluster_lengths:
        mean_clip = statistics.mean(cluster_lengths)
        median_clip = statistics.median(cluster_lengths)
        pct90 = statistics.quantiles(cluster_lengths, n=100)[89]
    else:
        mean_clip = median_clip = pct90 = 0

    auto_thresh = int(pct90)
    soft_count = sum(1 for ln in cluster_lengths if ln >= auto_thresh)
    soft_pct = soft_count / total * 100 if total else 0
    avg_mapq_all = mapq_all / total if total else 0
    avg_mapq_clip = mapq_clip / soft_count if soft_count else 0

    return {
        'total': total,
        'clusters': len(clusters),
        'mean_clip': mean_clip,
        'median_clip': median_clip,
        'pct90': pct90,
        'left_only': left_only,
        'right_only': right_only,
        'both': both,
        'auto_thresh': auto_thresh,
        'soft_count': soft_count,
        'soft_pct': soft_pct,
        'avg_mapq_all': avg_mapq_all,
        'avg_mapq_clip': avg_mapq_clip,
        'clusters_set': clusters,
        'cluster_lengths': cluster_lengths
    }


def filter_bucket(args):
    """Filter or trim reads based on threshold + support clusters."""
    bucket, bam_path, outdir, threshold, clusters, mode, idx = args
    out_sam = os.path.join(outdir, f'filter_{mode}_{idx}.sam')
    bam = pysam.AlignmentFile(bam_path, 'rb')
    with open(out_sam, 'w') as out:
        for ctg in bucket:
            for r in bam.fetch(ctg):
                keep = True
                if r.cigartuples:
                    left = r.cigartuples[0][1] if r.cigartuples[0][0] == 4 else 0
                    right = r.cigartuples[-1][1] if r.cigartuples[-1][0] == 4 else 0
                    cond = ((left >= threshold and r.reference_start in clusters) or
                            (right >= threshold and r.reference_end in clusters))
                    if cond:
                        if mode == 'full':
                            keep = False
                        else:  # trim
                            seq = r.query_sequence or ''
                            qual = r.query_qualities or []
                            r.query_sequence = seq[left:len(seq)-right]
                            r.query_qualities = qual[left:len(qual)-right]
                            r.cigartuples = [(op, ln) for op, ln in r.cigartuples if op != 4]
                if keep:
                    out.write(r.to_string() + "\n")
    bam.close()


def igv_bucket(args):
    """
    Write a BAM of all soft-clipped reads for IGV per bucket.
    Outputs one BAM: igv_{idx}.bam
    """
    bucket, bam_path, outdir, min_clip, idx = args
    out_bam = os.path.join(outdir, f'igv_{idx}.bam')
    bam_in = pysam.AlignmentFile(bam_path, 'rb')
    bam_out = pysam.AlignmentFile(out_bam,   'wb', template=bam_in)
    for ctg in bucket:
        for r in bam_in.fetch(ctg):
            # include any read with at least one soft-clip, regardless of length
            if r.cigartuples and any(op == 4 for op, ln in r.cigartuples):
                bam_out.write(r)
    bam_in.close()
    bam_out.close()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # build contig buckets once
    bam = pysam.AlignmentFile(args.input, 'rb')
    buckets = assign_buckets(bam.references, bam.lengths, args.threads)
    bam.close()

    # ─── stats mode ──────────────────────────────────────────────────────────
    if args.mode == 'stats':
        stats = compute_stats(args.input, buckets, args.threads)
        print(f"Total mapped segments: {stats['total']}")
        print(f"Soft-clipped segments ≥{stats['auto_thresh']}: {stats['soft_count']} ({stats['soft_pct']:.2f}%)")
        print(f"Mean clip length: {stats['mean_clip']:.1f}")
        print(f"Median clip length: {stats['median_clip']}")
        print(f"90th pct clip length: {stats['pct90']}")
        print(f"Clipped ends (left/right/both): {stats['left_only']}/{stats['right_only']}/{stats['both']}")
        print(f"Avg MAPQ (all/soft): {stats['avg_mapq_all']:.1f}/{stats['avg_mapq_clip']:.1f}")
        sys.exit(0)

    # ─── full/trim modes ─────────────────────────────────────────────────────
    if args.mode in ('full', 'trim'):
        stats = compute_stats(args.input, buckets, args.threads)
        # choose threshold: auto or user-provided
        if str(args.min_softclip).lower() == 'auto':
            threshold = stats['auto_thresh']
            print(f"Auto threshold used: {threshold}")
        else:
            threshold = int(args.min_softclip)
            count_user = sum(1 for ln in stats['cluster_lengths'] if ln >= threshold)
            pct_user = count_user / stats['total'] * 100 if stats['total'] else 0
            print(f"Soft-clipped segments ≥{threshold}: {count_user} ({pct_user:.2f}%)")
            print(f"Auto MIN_SOFTCLIP: {stats['auto_thresh']}")
        clusters = stats['clusters_set']

        # dispatch parallel filtering/trimming
        tasks = [
            (buckets[i], args.input, args.outdir, threshold, clusters, args.mode, i)
            for i in range(len(buckets))
        ]
        with Pool(args.threads) as pool:
            pool.map(filter_bucket, tasks)
        sys.exit(0)

    # ─── igv mode ────────────────────────────────────────────────────────────
    if args.mode == 'igv':
        tasks = [
            (buckets[i], args.input, args.outdir, args.min_softclip, i)
            for i in range(len(buckets))
        ]
        with Pool(args.threads) as pool:
            pool.map(igv_bucket, tasks)
        sys.exit(0)

    # unknown mode
    raise ValueError(f"Unknown mode: {args.mode}")

if __name__ == '__main__':
    main()
