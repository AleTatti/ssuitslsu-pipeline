#!/usr/bin/env python3
def main():
    import os
    import re
    import math
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

    # ─── Helper: normalize any ID to alphanumeric+underscores ──────────────x
    def clean_id(seq_id):
        return re.sub(r"[^0-9a-zA-Z]+", "_", seq_id)


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

    # ─── Convert raw‐MSA cols → trimmed‐MSA cols by subtracting head_trim ───
    contig2trim = {}
    for sid, mapping in contig2raw.items():
        mapping_trim = {}
        for cpos, rawcol in mapping.items():
            if rawcol > head_trim:
                mapping_trim[cpos] = rawcol - head_trim
        contig2trim[sid] = mapping_trim

    # ─── Load & clean tree ───────────────────────────────────────
    tree = Phylo.read(tree_file, "newick")
    for leaf in tree.get_terminals():
        leaf.name = re.sub(r'^_?R_', '', clean_id(leaf.name))
    tree_leaves = {leaf.name for leaf in tree.get_terminals()}

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


    # colour map for shading
    region_colors = {
        'SSU':   'lightblue',
        'ITS1':  'lightgreen',
        '5.8S':  'khaki',
        'ITS2':  'lightgreen',
        'LSU':   'lightblue',
    }

    # ─── Read & filter contigs by depth ≥100 ───────────────────────────────
    with open(contigs_fasta) as f:
        headers = [line[1:].strip() for line in f if line.startswith(">")]

    sample_contig_ids = []
    for hdr in headers:
        parts = hdr.split()
        cid   = parts[0]
        if cid.startswith("NODE_"):
            m = re.search(r"cov_([0-9]+(?:\.[0-9]+)?)$", cid)
            if m and float(m.group(1)) >= 100:
                sample_contig_ids.append(cid)
        else:
            info = dict(p.split("=",1) for p in parts[1:] if "=" in p)
            if float(info.get("multi",0.0)) >= 100:
                sample_contig_ids.append(cid)

    sample_seqs = [
        clean_id(cid) for cid in sample_contig_ids
        if clean_id(cid) in alignment_dict and clean_id(cid) in tree_leaves
    ]

    # ─── Build genus‐level reference list ──────────────────────────────────
    with open(genus_ref_fasta) as f:
        raw_refs = [line[1:].strip().split()[0] for line in f if line.startswith(">")]
    ref_seq_ids = [
        clean_id(rid) for rid in raw_refs
        if clean_id(rid) in alignment_dict and clean_id(rid) in tree_leaves
    ]

    # ─── Helpers ───────────────────────────────────────────────────────────
    def short_id(seq_id):
        parts = [p for p in seq_id.split("_") if p]
        return "_".join(parts[:2]) if len(parts)>=2 else (parts[0] if parts else seq_id)

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

    def sliding_k2p(sid, rid, win=100, step=20):
        """
        Compute a sliding-window Kimura-2P profile between sequences sid and rid
        in the trimmed alignment. Only compute K2P if at least half the window
        consists of ungapped site-pairs.
        """
        seq1 = alignment_dict[sid].seq
        seq2 = alignment_dict[rid].seq
        L = len(seq1)

        min_sites = win // 2    # require at least half the window in non-gaps

        prof = []
        for i in range(0, L - win + 1, step):
            w1 = seq1[i : i+win]
            w2 = seq2[i : i+win]

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

    def detect_chimera(profile, jump_thresh=0.15, stable_span=3):
        pd = [p for _,p,_ in profile]
        for j in range(1, len(pd)-stable_span):
            if abs(pd[j]-pd[j-1])>jump_thresh and all(
               abs(pd[j+k]-pd[j])<0.05 for k in range(1,stable_span+1)):
                return True
        return False

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
        out.write("contig_id\tclosest_ref\tglobal_k2p\tis_chimera\tavg_gap_frac\n")
        for idx, sid in enumerate(sample_seqs):
            ref = get_closest_ref(sid)
            if not ref:
                continue

            # 1) global K2P distance
            seq1 = str(alignment_dict[sid].seq)
            seq2 = str(alignment_dict[ref].seq)
            global_k2p = kimura2p(seq1, seq2)

            # 2) sliding windows (only skips double-gappy windows, keeps deletions)
            profile = sliding_k2p(sid, ref)
            chim    = detect_chimera(profile)
            avg_gap = sum(g for _, _, g in profile) / len(profile) if profile else 0

            ## DBG
            centers = [c for c,_,_ in profile]
            print(f"[DEBUG] {sid} sliding‐window centers (len={len(centers)}): {centers[:10]} …")

            # write report line
            out.write(f"{sid}\t{ref}\t{global_k2p:.6f}\t{chim}\t{avg_gap:.3f}\n")

            # prepare plot panel
            row, col = divmod(idx, ncols)
            ax = axes[row][col]
            if not profile:
                ax.set_visible(False)
                continue

            # 3) shade ITSx regions in the *trimmed* alignment directly, clamped to bounds
            seen_regions = set()
            L_trim = alignment.get_alignment_length()
            # pre-compute min/max mapped column in case we need to clamp
            mapping_vals = contig2aln.get(sid, {}).values()
            min_col = min(mapping_vals) if mapping_vals else 1
            max_col = max(mapping_vals) if mapping_vals else L_trim

            for name, (cstart, cend) in regions.get(sid, {}).items():
                mapping = contig2aln.get(sid, {})

                # get start; if missing (region begins before trimmed‐alignment), clamp to first mapped col
                start_col = mapping.get(cstart, min_col)
                # get end; if missing (region extends beyond trimmed‐alignment), clamp to last mapped col
                end_col   = mapping.get(cend,   max_col)

                # convert to 0-based alignment coords and clamp within [0, L_trim-1]
                aln_start = max(0, start_col - 1)
                aln_end   = min(L_trim - 1, end_col - 1)
                if aln_end < aln_start:
                    # nothing to draw
                    continue
                print(f"[DEBUG] {sid} {name}: contig {cstart}–{cend} → align {aln_start}–{aln_end}")
                lbl = name if name not in seen_regions else None
                ax.axvspan(
                    aln_start, aln_end,
                    color=region_colors.get(name, 'gray'),
                    alpha=0.3, lw=0,
                    label=lbl
                )
                seen_regions.add(name)


            # 4) unpack sliding window data
            x, y, gaps = zip(*profile)

            # 5) gap-fraction shading (behind everything)
            ax.fill_between(
                x, 0, gaps,
                step='mid',
                alpha=0.2,
                color='gray',
                label='gap fraction',
                zorder=1
            )

            # 6) global K2P line
            ax.axhline(
                global_k2p,
                color='tab:orange', linestyle='--',
                label='global K2P',
                zorder=2
            )

            # 7) sliding K2P curve
            ax.plot(
                x, y,
                color='tab:blue',
                label='sliding K2P',
                zorder=3
            )

            # ───────── crop the x‐axis to the contig’s coverage ─────────────
            # This uses contig2aln, which maps each contig‐position → 1‐based trimmed column.
            cols = list(contig2aln[sid].values())
            ax.set_xlim(min(cols), max(cols))

            # 8) finalize panel
            ax.set_title(
                f"{short_id(sid)} vs {short_id(ref)}\n"
                f"{'chimera' if chim else 'clean'}",
                color=('red' if chim else 'darkgreen')
            )
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
                   "<th>avg_gap_frac</th></tr>\n")
        next(tsv)
        for line in tsv:
            contig, ref, dist, is_chim, avg_gap = line.strip().split("\t")
            color = "red" if is_chim=="True" else "darkgreen"
            html.write(
              f"<tr><td>{contig}</td><td>{ref}</td>"
              f"<td>{dist}</td>"
              f"<td style='color:{color};font-weight:bold'>{is_chim}</td>"
              f"<td>{avg_gap}</td></tr>\n"
            )
        html.write("</table>\n")

if __name__=="__main__":
    main()
