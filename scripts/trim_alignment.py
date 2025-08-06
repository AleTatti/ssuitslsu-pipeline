#!/usr/bin/env python3
import os, re, csv
from Bio import AlignIO
import sys
from Bio.Align import MultipleSeqAlignment

# Helper to normalize any ID to alphanumeric+underscores
def clean_id(seq_id):
    return re.sub(r"[^0-9a-zA-Z]+", "_", seq_id)

# Environment variables expected:
#  ITSX_DIR, sample, PHYLO_RAW, trimmed_aln, PHYLO_DIR

# 0) parse sample-level ITSx position files into regions dict
regions = {}
small_file = os.path.join(os.environ["ITSX_DIR"], f"{os.environ['sample']}.small.positions.txt")
large_file = os.path.join(os.environ["ITSX_DIR"], f"{os.environ['sample']}.large_nhmmer.positions.txt")
pos_file = large_file if os.path.exists(large_file) else (small_file if os.path.exists(small_file) else None)
if pos_file:
    with open(pos_file) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            sid_raw = parts[0]
            sid = clean_id(sid_raw)
            regmap = {}
            for fld in parts[2:]:
                if ':' not in fld:
                    continue
                name, coords = fld.split(':', 1)
                m = re.match(r"^(\d+)\s*-\s*(\d+)$", coords.strip())
                if not m:
                    continue
                start, end = map(int, m.groups())
                regmap[name] = (start, end)
            if regmap:
                regions[sid] = regmap

# 1) load raw MSA and compute occupancy
raw = AlignIO.read(os.environ["PHYLO_RAW"], "fasta")
n = len(raw)
L = raw.get_alignment_length()
occ = [sum(rec.seq[i] != "-" for rec in raw) for i in range(L)]
half = n / 2.0
start_occ = next(i for i, o in enumerate(occ) if o >= half)
end_occ = L - next(i for i, o in enumerate(reversed(occ)) if o >= half)

# 2) build contig→raw map for every sequence (1-based columns)
contig2raw = {}
for rec in raw:
    rid = re.sub(r'^_?R_', '', rec.id)
    sid = clean_id(rid)
    mapping = {}
    ungapped = 0
    for col, nt in enumerate(str(rec.seq), start=1):
        if nt != "-":
            ungapped += 1
            mapping[ungapped] = col
    contig2raw[sid] = mapping

# 3) extract ITSx start/end columns from contig→raw mapping
ssu_cols = []
lsu_cols = []
for sid, mapping in contig2raw.items():
    if sid not in regions:
        continue
    for region, (cs, ce) in regions[sid].items():
        if region.upper().startswith('SSU') and cs in mapping:
            ssu_cols.append(mapping[cs])
        if region.upper().startswith('LSU') and ce in mapping:
            lsu_cols.append(mapping[ce])

# 4) Force the trim window to include every SSU start and LSU end
ssu_zeros = [col - 1 for col in ssu_cols]
lsu_excls = [col for col in lsu_cols]
trim_start = min([start_occ] + ssu_zeros) if (ssu_zeros or lsu_excls) else start_occ
trim_end = max([end_occ] + lsu_excls) if (ssu_zeros or lsu_excls) else end_occ
trim_start = max(0, trim_start)
trim_end = min(L, trim_end)



# 5) slice & drop gap-only sequences
trimmed = raw[:, trim_start:trim_end]
removed = [rec.id for rec in trimmed if all(nt == '-' for nt in rec.seq)]
filtered = MultipleSeqAlignment(rec for rec in trimmed if any(nt != '-' for nt in rec.seq))
if removed:
    sys.stderr.write(
        "Sequences removed because they were only gaps in the trimmed window: "
        + ", ".join(removed)
        + "\n"
    )
if len(filtered) == 0:
    sys.stderr.write("[ERROR] All sequences are gap-only after trimming → exiting\n")
    sys.exit(1)
AlignIO.write(filtered, os.environ["trimmed_aln"], "fasta")

# 6) rebuild contig→trim mapping in 0-based trimmed coordinates
contig2trim = {}
for sid, mapping in contig2raw.items():
    newmap = {}
    for cpos, rawcol in mapping.items():
        idx_trim = rawcol - 1 - trim_start
        if 0 <= idx_trim < filtered.get_alignment_length():
            newmap[cpos] = idx_trim
    contig2trim[sid] = newmap
    if sid in regions:
        for region, (cs, ce) in regions[sid].items():
            _ = newmap.get(cs)
            _ = newmap.get(ce)

# 7) Write raw vs trimmed ITSx coords to TSV
out_tsv = os.path.join(
    os.environ["PHYLO_DIR"],
    f"{os.environ['sample']}_ITSx_coords_aln_trimmed.tsv"
)
with open(out_tsv, "w") as out:
    out.write("contig_id\tregion\traw_start\traw_end\ttrimmed_start\ttrimmed_end\n")
    for sid, regmap in regions.items():
        raw_map = contig2raw.get(sid, {})
        trim_map = contig2trim.get(sid, {})
        for region, (cs, ce) in regmap.items():
            out.write(f"{sid}\t{region}\t{raw_map.get(cs, '')}\t{raw_map.get(ce, '')}\t{trim_map.get(cs, '')}\t{trim_map.get(ce, '')}\n")

# 8) report how many columns trimmed off each end
print(trim_start, L - trim_end)
