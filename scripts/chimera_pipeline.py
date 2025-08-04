#!/usr/bin/env python3
def main():
    import os
    import re
    import math
    import numpy as np
    from Bio import AlignIO, Phylo
    import matplotlib.pyplot as plt

    raw_aln_path    = os.environ["PHYLO_RAW"]
    head_trim       = int(os.environ["TRIM_HEAD"])
    aln_path        = os.environ["PHYLO_ALN"]
    tree_file       = os.environ["TREEFILE"]
    contigs_fasta   = os.environ["CONTIGS_FASTA"]
    genus_ref_fasta = os.environ["GENUS_REF_FASTA"]
    out_tsv         = os.environ["CHIMERA_REPORT"]
    plot_path       = os.environ["CHIMERA_PLOT"]
    sample_id       = os.environ["sample"]
    itsx_dir        = os.environ["ITSX_DIR"]
    USE_ALIGNMENT = os.environ.get("CHIMERA_USE_ALIGNMENT", "false").lower() == "true"


    # ─── Helper: normalize any ID to alphanumeric+underscores ──────────────x
    def clean_id(seq_id):
        return re.sub(r"[^0-9a-zA-Z]+", "_", seq_id)

    # ———————————————————————————————————————————
    # Load ITSx–hit contig IDs so we can override filters
    itsx_ids = set()
    itsx_dir = os.environ.get("ITSX_DIR", "")
    sample  = os.environ.get("sample", "")
    for suffix in ("large_nhmmer.positions.txt", "small.positions.txt"):
        pf = os.path.join(itsx_dir, f"{sample}.{suffix}")
        if os.path.exists(pf):
            with open(pf) as fh:
                for line in fh:
                    contig = line.split()[0]    # ID is the first token
                    itsx_ids.add(contig)
    # ———————————————————————————————————————————
    # ─── Load the raw MSA and build contig→raw‐MSA maps ────────────────────
    raw = AlignIO.read(raw_aln_path, "fasta")
    contig2raw = {}
    for rec in raw:
        sid     = clean_id(rec.id)
        mapping = {}
        cpos    = 0
        for idx, base in enumerate(str(rec.seq)):
            if base != "-":
                cpos += 1
                mapping[cpos] = idx + 1       # raw‐MSA col (1‐based)
        contig2raw[sid] = mapping

    # ———————————————————————————————————————————
    # Load ITSx–hit contig IDs so we can override filters
    itsx_ids = set()
    itsx_dir = os.environ.get("ITSX_DIR", "")
    sample  = os.environ.get("sample", "")
    for suffix in ("large_nhmmer.positions.txt", "small.positions.txt"):
        pf = os.path.join(itsx_dir, f"{sample}.{suffix}")
        if os.path.exists(pf):
            with open(pf) as fh:
                for line in fh:
                    contig = line.split()[0]    # ID is the first token
                    itsx_ids.add(contig)
    # ———————————————————————————————————————————


    # ─── Load & clean tree ───────────────────────────────────────
    if not USE_ALIGNMENT:
        tree = Phylo.read(tree_file, "newick")
        for leaf in tree.get_terminals():
            leaf.name = re.sub(r'^_?R_', '', clean_id(leaf.name))
        tree_leaves = {leaf.name for leaf in tree.get_terminals()}
    else:
        tree = None
        # in alignment mode, we won’t use tree_leaves — all filtering will be via alignment_dict
        tree_leaves = set()


    # ─── Load & clean alignment ──────────────────────────────────
    alignment = AlignIO.read(aln_path, "fasta")
    for rec in alignment:
        rec.id   = rec.name = re.sub(r'^_?R_', '', clean_id(rec.id))
    alignment_dict = {rec.id: rec for rec in alignment}

    # ─── Build contig→alignment maps ───────────────────────────────────────
    # For each contig ID we’ll record, for each non-gap index in the contig
    # (1-based), what column (1-based) in the alignment that corresponds to.
    contig2aln = {}
    for sid, rec in alignment_dict.items():
        seq = str(rec.seq)
        mapping = {}
        contig_pos = 0
        # aln_idx is 0-based; we store col = aln_idx+1
        for aln_idx, base in enumerate(seq):
            if base != '-':
                contig_pos += 1
                mapping[contig_pos] = aln_idx + 1
        contig2aln[sid] = mapping


    # ─── Determine which ITSx positions file to use ────────────────────────
    small_file = os.path.join(itsx_dir, f"{sample_id}.small.positions.txt")
    large_file = os.path.join(itsx_dir, f"{sample_id}.large_nhmmer.positions.txt")
    if os.path.exists(large_file):
        pos_file = large_file
    elif os.path.exists(small_file):
        pos_file = small_file
    else:
        pos_file = None

    # ─── Parse region positions from ITSx ──────────────────────────────────
    regions = {}
    if pos_file:
        print(f"[DBG] Using ITSx positions file: {pos_file}")
        with open(pos_file) as pf:
            for line in pf:
                parts = line.strip().split('\t')
                raw_id = parts[0]
                node = clean_id(parts[0])
                regmap = {}
                for fld in parts[2:]:
                    if ':' not in fld:
                        continue
                    name, coords = fld.split(':',1)
                    coords = coords.strip()
                    m = re.match(r'^(\d+)\s*-\s*(\d+)$', coords)
                    if not m:
                        # “No end”, “Not found”, etc. → skipped
                        continue
                    start, end = map(int, m.groups())
                    regmap[name] = (start, end)
                if regmap:
                  regions[node] = regmap



    # ─── Colors set on “Okabe–Ito” color‐blind‐safe set ───────────────────────────────
    region_colors = {
        'SSU':  '#E69F00',  # orange
        'LSU':  '#D55E00',  # vermillion
        'ITS1': '#56B4E9',  # sky-blue
        'ITS2': '#009E73',  # bluish-green
        '5.8S': '#F0E442',  # yellow
    }


    # ─── Read & filter contigs by depth ≥60 (but always include ITSx hits) ────
    # load ITSx-hit contig IDs
    itsx_ids = set()
    itsx_dir = os.environ["ITSX_DIR"]
    sample   = os.environ["sample"]
    for suffix in ("large_nhmmer.positions.txt", "small.positions.txt"):
        pf = os.path.join(itsx_dir, f"{sample}.{suffix}")
        if os.path.exists(pf):
            with open(pf) as fh:
                for line in fh:
                    itsx_ids.add(line.split()[0])

    with open(contigs_fasta) as f:
        headers = [line[1:].strip() for line in f if line.startswith(">")]

    sample_contig_ids = []
    for hdr in headers:
        cid = hdr.split()[0]

        # 1) always include ITSx-hit contigs
        if cid in itsx_ids:
            sample_contig_ids.append(cid)
            continue

        # 2) otherwise apply your existing depth/multi filter
        if cid.startswith("NODE_"):
            m = re.search(r"cov_([0-9]+(?:\.[0-9]+)?)$", cid)
            if m and float(m.group(1)) >= 60:
                sample_contig_ids.append(cid)
        else:
            info = dict(p.split("=",1) for p in hdr.split()[1:] if "=" in p)
            if float(info.get("multi",0.0)) >= 60:
                sample_contig_ids.append(cid)


    # DIAGNOSTIC PRINTS
    print("=== Diagnostic: Contig name correspondence ===")
    print("Raw sample_contig_ids:")
    for cid in sample_contig_ids:
        print("  ", cid)
    print("Cleaned sample_contig_ids:")
    for cid in sample_contig_ids:
        print("  ", clean_id(cid))
    print("Alignment IDs:")
    for aid in alignment_dict.keys():
        print("  ", aid)
    print("Tree leaf IDs:")
    for leaf in tree_leaves:
        print("  ", leaf)
    print("=== End diagnostic ===\n")

    sample_seqs = [
    clean_id(cid) for cid in sample_contig_ids
    if clean_id(cid) in alignment_dict
       and (USE_ALIGNMENT or clean_id(cid) in tree_leaves)
    ]

    # ─── Build genus‐level reference list ──────────────────────────────────
    with open(genus_ref_fasta) as f:
        raw_refs = [line[1:].strip().split()[0] for line in f if line.startswith(">")]

    if USE_ALIGNMENT:
        # in alignment mode, allow any ref in the trimmed alignment
        ref_seq_ids = [clean_id(rid) for rid in raw_refs if clean_id(rid) in alignment_dict]
    else:
        ref_seq_ids = [
            clean_id(rid) for rid in raw_refs
            if clean_id(rid) in alignment_dict and clean_id(rid) in tree_leaves
        ]

    # ─── Helpers ───────────────────────────────────────────────────────────
    def short_id(seq_id):
        parts = [p for p in seq_id.split("_") if p]
        return "_".join(parts[:2]) if len(parts)>=2 else (parts[0] if parts else seq_id)

    # ─── Nearest-reference function ─────────────────────────────────────
    if USE_ALIGNMENT:
        # alignment-based p-distance fallback
        def get_closest_ref(sid):
            rec_q = alignment_dict.get(sid)
            if not rec_q:
                return None
            def p_distance(a, b):
                s1, s2 = str(a.seq), str(b.seq)
                diffs = sum(1 for x, y in zip(s1, s2)
                            if x != '-' and y != '-' and x != y)
                comps = sum(1 for x, y in zip(s1, s2)
                            if x != '-' and y != '-')
                return diffs / comps if comps else 0.0
            # compute all distances and pick minimal
            dlist = [(p_distance(rec_q, alignment_dict[rid]), rid)
                     for rid in ref_seq_ids]
            return min(dlist, key=lambda t: t[0])[1] if dlist else None
    else:
        # original tree-based logic
        def get_closest_ref(sid):
            dists = [(rid, tree.distance(sid, rid)) for rid in ref_seq_ids]
            return min(dists, key=lambda x: x[1])[0] if dists else None


    def kimura2p(s1, s2):
        pairs = [(a.upper(),b.upper()) for a,b in zip(s1,s2) if a!='-' and b!='-']
        L = len(pairs)
        if L==0: return float("nan")
        ti = sum(1 for a,b in pairs if (a,b) in (('A','G'),('G','A'),('C','T'),('T','C')))
        tv = sum(1 for a,b in pairs if a!=b and (a,b) not in (('A','G'),('G','A'),('C','T'),('T','C')))
        P,Q = ti/L, tv/L
        try:
            return -0.5*math.log(1-2*P-Q) - 0.25*math.log(1-2*Q)
        except ValueError:
            return float("nan")

    def sliding_k2p(seq1_str, seq2_str, win=50, step=10):
        """
        Compute a sliding-window Kimura-2P profile between sequences sid and rid
        in the trimmed alignment. Only compute K2P if at least half the window
        consists of ungapped site-pairs.
        """
        L = len(seq1_str)

        min_sites = win // 2    # require at least half the window in non-gaps

        prof = []
        for i in range(0, L - win + 1, step):
            w1 = seq1_str[i : i+win]
            w2 = seq2_str[i : i+win]

            # drop columns where both have gaps
            filtered      = [(a,b) for a,b in zip(w1, w2) if not (a=='-' and b=='-')]

            # if filtered is empty, that window was 100% gaps in both => skip it entirely
            if not filtered:
                continue

            # keep only fully non-gap pairs
            pairs_nominal = [(a,b) for a,b in filtered if a!='-' and b!='-']

            # only compute K2P if we have enough sites
            if len(pairs_nominal) < min_sites:
                d = float('nan')
            else:
                s1_clean = ''.join(a for a,b in pairs_nominal)
                s2_clean = ''.join(b for a,b in pairs_nominal)
                d = kimura2p(s1_clean, s2_clean)

            gap_frac = 1 - len(pairs_nominal)/win
            center   = i + win//2
            prof.append((center, d, gap_frac))

        return prof

    def detect_chimera_with_stats(profile, alpha=0.01):
        """
        Detect chimeras using statistical validation and optimal thresholds.
        Returns (is_chimera, breakpoint_start, breakpoint_end, p_value, adaptive_threshold)
        """
        import numpy as np
        from scipy import stats

        centers, distances, gaps = zip(*profile)
        distances = np.array(distances)

        # Remove NaN values for statistical analysis
        valid_idx = ~np.isnan(distances)
        if np.sum(valid_idx) < 10:  # Need at least 10 valid points
            return False, None, None, 1.0, None

        valid_distances = distances[valid_idx]
        valid_centers = np.array(centers)[valid_idx]

        # Calculate baseline statistics for adaptive threshold
        baseline_mean = np.median(valid_distances)
        baseline_std = np.std(valid_distances)

        # Adaptive threshold based on data distribution
        adaptive_thresh = baseline_mean + 2.5 * baseline_std
        min_thresh = 0.05  # Minimum biological threshold
        adaptive_thresh = max(adaptive_thresh, min_thresh)

        # Find potential breakpoints using change point detection
        breakpoints = []
        min_stable_span = 2  # Minimum 2 windows for stability

        for i in range(1, len(valid_distances) - min_stable_span):
            current_dist = valid_distances[i]
            prev_dist = valid_distances[i-1]

            # Check for significant jump
            jump_magnitude = abs(current_dist - prev_dist)

            if jump_magnitude > adaptive_thresh:
                # Verify stability after jump
                stable = True
                post_jump_mean = np.mean(valid_distances[i:i+min_stable_span])

                for k in range(min_stable_span):
                    if i+k < len(valid_distances):
                        if abs(valid_distances[i+k] - post_jump_mean) > baseline_std:
                            stable = False
                            break

                if stable:
                    # Statistical test for significance
                    before_segment = valid_distances[max(0, i-3):i]
                    after_segment = valid_distances[i:min(len(valid_distances), i+3)]

                    if len(before_segment) >= 2 and len(after_segment) >= 2:
                        # Mann-Whitney U test for different distributions
                        try:
                            _, p_value = stats.mannwhitneyu(before_segment, after_segment,
                                                          alternative='two-sided')

                            if p_value < alpha:
                                # Find breakpoint extent
                                breakpoint_start = valid_centers[i-1]
                                breakpoint_end = valid_centers[min(i+min_stable_span-1, len(valid_centers)-1)]

                                breakpoints.append({
                                    'start': breakpoint_start,
                                    'end': breakpoint_end,
                                    'p_value': p_value,
                                    'jump_magnitude': jump_magnitude,
                                    'threshold': adaptive_thresh
                                })
                        except ValueError:
                            continue

        if breakpoints:
            # Return the most significant breakpoint
            best_breakpoint = min(breakpoints, key=lambda x: x['p_value'])
            return (True,
                   best_breakpoint['start'],
                   best_breakpoint['end'],
                   best_breakpoint['p_value'],
                   best_breakpoint['threshold'])
        else:
            return False, None, None, 1.0, adaptive_thresh

    # ─── Run chimera detection, plot, and write report ─────────────────────
    n = len(sample_seqs)
    if n==0:
        print("[WARN] No contigs passed filters – exiting.")
        exit(0)

    ncols = 2
    nrows = (n+ncols-1)//ncols
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=(12,4*nrows), squeeze=False)

    with open(out_tsv, "w") as out:
        out.write("contig_id\tclosest_ref\tglobal_k2p\tis_chimera\tbreakpoint_start\tbreakpoint_end\tp_value\tadaptive_threshold\tavg_gap_frac\n")
        for idx, sid in enumerate(sample_seqs):
            ref = get_closest_ref(sid)
            if not ref:
                continue

            # ─── Build the two‐sequence, gap-pruned alignment ───────────────────────────
            from Bio.Align import MultipleSeqAlignment
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord

            # extract only our node and its closest reference
            pair_aln = MultipleSeqAlignment([
                alignment_dict[sid],
                alignment_dict[ref]
            ])
            # drop any column both sequences gap out
            cols_to_keep = [
                i for i in range(pair_aln.get_alignment_length())
                if any(rec.seq[i] != '-' for rec in pair_aln)
            ]
            # rebuild the pairwise MSA with only those columns
            filtered = []
            for rec in pair_aln:
                new_seq = "".join(rec.seq[i] for i in cols_to_keep)
                filtered.append(SeqRecord(Seq(new_seq), id=rec.id, description=""))
            pair_aln = MultipleSeqAlignment(filtered)

            # ─── DEBUG #1: show pair_aln length ───────────────────────────────
            #print(f"[PLOT DEBUG] {sid}: pair_aln length = {pair_aln.get_alignment_length()}")


            seq1_pw = str(pair_aln[0].seq)
            seq2_pw = str(pair_aln[1].seq)

            # 1) global K2P distance on the gap-pruned pairwise alignment
            global_k2p = kimura2p(seq1_pw, seq2_pw)


            # 2) sliding-window K2P *and* gap‐fraction on the pairwise alignment
            profile = sliding_k2p(seq1_pw, seq2_pw, win=50, step=10)

            # Enhanced chimera detection with statistics
            chim, bp_start, bp_end, p_val, adapt_thresh = detect_chimera_with_stats(profile)
            avg_gap = sum(g for _, _, g in profile) / len(profile) if profile else 0

            centers_pw, dists_pw, gaps_pw = zip(*profile)

            # Format breakpoint positions
            bp_start_str = f"{bp_start:.0f}" if bp_start is not None else "NA"
            bp_end_str = f"{bp_end:.0f}" if bp_end is not None else "NA"
            p_val_str = f"{p_val:.2e}" if p_val < 1.0 else "NA"
            thresh_str = f"{adapt_thresh:.4f}" if adapt_thresh is not None else "NA"

            # write report line
            out.write(f"{sid}\t{ref}\t{global_k2p:.6f}\t{chim}\t{bp_start_str}\t{bp_end_str}\t{p_val_str}\t{thresh_str}\t{avg_gap:.3f}\n")

            # prepare plot panel
            row, col = divmod(idx, ncols)
            ax = axes[row][col]
            if not profile:
                ax.set_visible(False)
                continue

            # ─── DEBUG #2: confirm axes selection ───────────────────────────────
            #print(f"[PLOT DEBUG] plotting on axes[{row}][{col}] for {sid}")


            # 3) shade ITSx regions on the *pairwise* alignment (using regions_pair)
            # ─── 1) Build a map from trimmed‐MSA columns → pairwise‐MSA columns ──────────────────
            # `cols_to_keep` is a list of 0-based indices in the trimmed alignment that survive when you remove columns where both sequences are gaps.
            # We need to know, for any given trimmed‐alignment index, what its new index will be in the gap-pruned two-sequence alignment (pair_aln).
            trimmed_to_pair = {
                trimmed_idx: pair_idx
                for pair_idx, trimmed_idx in enumerate(cols_to_keep)
            }


            # ─── 2) Build contig2pair: map from raw contig coordinate → column in pair_aln ─────────
            pair_mapping = {}
            for cpos, aln_col1 in contig2aln[sid].items():
                # convert trimmed-aln 1-based → 0-based
                t_idx = aln_col1 - 1
                if t_idx in trimmed_to_pair:
                    pair_mapping[cpos] = trimmed_to_pair[t_idx]


            # ─── 3) Recalculate each ITSx region in the pairwise alignment’s coordinates ─────────
            regions_pair = {}
            for name, (cstart, cend) in regions.get(sid, {}).items():
                # Lookup the new start/end columns for this region
                a = pair_mapping.get(cstart)   # 0-based start in pair_aln
                b = pair_mapping.get(cend)     # 0-based end   in pair_aln

                # Only include if both ends still map into the pairwise alignment
                if a is not None and b is not None:
                    regions_pair[name] = (a, b)


            # ─── 4) Shade each ITSx region as a semi-transparent span on the axes ────────────────
            seen = set()
            for name, (a, b) in regions_pair.items():
                # Only label the first occurrence of each region in the legend
                lbl = name if name not in seen else None

                # Draw the colored rectangle from column a → b
                ax.axvspan(
                    a,      # left edge (0-based)
                    b,      # right edge
                    color=region_colors.get(name, 'gray'),
                    alpha=0.3,
                    lw=0,
                    label=lbl
                )
                seen.add(name)


            # 5) gap-fraction shading (behind everything)
            ax.fill_between(
                centers_pw,
                0,
                gaps_pw,
                step='mid',
                alpha=0.2,
                color='gray',
                label='pairwise window gap fraction',
                zorder=1
            )


            # 6) adaptive threshold line
            if adapt_thresh is not None:
                baseline_median = np.median([d for d in dists_pw if not np.isnan(d)])
                ax.axhline(
                    baseline_median + adapt_thresh,
                    color='tab:red', linestyle=':',
                    label=f'adaptive threshold ({adapt_thresh:.3f})',
                    zorder=2
                )

            # 7) global K2P line
            ax.axhline(
                global_k2p,
                color='tab:orange', linestyle='--',
                label='global K2P',
                zorder=2
            )

            # 8) sliding K2P curve
            ax.plot(
                centers_pw,
                dists_pw,
                color='tab:blue',
                label='sliding K2P (non-gap)',
                zorder=3
            )

            # 9) highlight breakpoint if detected
            if bp_start is not None and bp_end is not None:
                ax.axvspan(bp_start, bp_end, color='red', alpha=0.3,
                          label=f'chimeric region (p={p_val:.2e})', zorder=4)

            # ───────── crop the x‐axis to the contig’s coverage ─────────────
            # This uses contig2aln, which maps each contig‐position → 1‐based trimmed column.
            cols = list(contig2aln[sid].values())

            # crop to the span of our two-sequence alignment
            ax.set_xlim(0, pair_aln.get_alignment_length() - 1)

            # 10) finalize panel
            title_text = f"{short_id(sid)} vs {short_id(ref)}\n"
            if chim:
                title_text += f"CHIMERA (p={p_val:.2e}, pos {bp_start_str}-{bp_end_str})"
            else:
                title_text += "clean"

            ax.set_title(title_text, color=('red' if chim else 'darkgreen'))
            ax.set_xlabel("Position (bp)")
            ax.set_ylabel("K2P distance")
            ax.legend()

    # drop any empty subplots
    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        fig.delaxes(axes[r][c])

    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()


    # ─── Build HTML report ───────────────────────────────────────────────────
    sample_name = os.path.basename(out_tsv).replace("_chimera_report.tsv","")
    html_path   = os.path.join(os.path.dirname(out_tsv),
                               f"{sample_name}_chimera_report.html")
    with open(out_tsv) as tsv, open(html_path,"w") as html:
        html.write(f"<h2>Chimera report for {sample_name}</h2>\n")
        html.write(f'<img src="{os.path.basename(plot_path)}" '
                   'style="max-width:100%;"><br>\n')
        html.write("<table border=1 cellpadding=4>"
                   "<tr><th>contig_id</th><th>closest_ref</th>"
                   "<th>global_k2p</th><th>is_chimera</th>"
                   "<th>breakpoint_start</th><th>breakpoint_end</th>"
                   "<th>p_value</th><th>adaptive_threshold</th>"
                   "<th>avg_gap_frac</th></tr>\n")
        next(tsv)
        for line in tsv:
            parts = line.strip().split("\t")
            contig, ref, dist, is_chim, bp_start, bp_end, p_val, thresh, avg_gap = parts
            color = "red" if is_chim=="True" else "darkgreen"
            html.write(
              f"<tr><td>{contig}</td><td>{ref}</td>"
              f"<td>{dist}</td>"
              f"<td style='color:{color};font-weight:bold'>{is_chim}</td>"
              f"<td>{bp_start}</td><td>{bp_end}</td>"
              f"<td>{p_val}</td><td>{thresh}</td>"
              f"<td>{avg_gap}</td></tr>\n"
            )
        html.write("</table>\n")

if __name__=="__main__":
    main()
