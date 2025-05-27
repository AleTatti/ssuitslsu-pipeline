#!/usr/bin/env bash

set -euo pipefail

# ──────────────────────────────────────────────────────────────────────────────
# Version checks
echo "➤ Bash version: $BASH_VERSION"
conda --version >/dev/null 2>&1 && echo "➤ Conda version: $(conda --version | cut -d' ' -f2)" \
                             || echo "⚠️  conda not found in PATH"
# ──────────────────────────────────────────────────────────────────────────────


export PYTHONUNBUFFERED=1

# build your log‐filename
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG="ssuitslsu_${TIMESTAMP}.log"
# redirect all stdout+stderr into tee, which writes to both console and log
exec > >(tee -a "$LOG") 2>&1

t_align_start=0
t_align_end=0
t_phylo_start=0
t_phylo_end=0

# ──────────────────────────────────────────────────────────────────────────────
# Usage/help message
usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  -n, --no-trim           Skip fastp trimming
  -m, --max-memory MB     Override memory (GB) for SPAdes
  -t, --threads N         Override number of threads for all steps
  -c, --config FILE       Pipeline config YAML (default: config/pipeline.yaml)
      --assembler [spades|megahit]
      --mapq N            Override mapping quality filter (samtools -q)
      --outdir DIR        Override output directory
  -h, --help              Show this help message and exit
EOF
}
# ──────────────────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────────────────
# Initialize CLI-override vars with defaults
SKIP_TRIM_CLI=false
MEM_CLI=""
THREADS_CLI=""
CONFIG_CLI=""
ASSEMBLER_CLI=""
MAPQ_CLI=""
OUTDIR_CLI=""

# ──────────────────────────────────────────────────────────────────────────────

# Parse CLI args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)         usage; exit 0 ;;
    -c|--config)       CONFIG_CLI="$2"; shift 2 ;;
    -n|--no-trim)      SKIP_TRIM_CLI=true; shift ;;
    -m|--max-memory)   MEM_CLI="$2"; shift 2 ;;
    -t|--threads)      THREADS_CLI="$2"; shift 2 ;;
    --mapq)            MAPQ_CLI="$2"; shift 2 ;;
    --assembler)       ASSEMBLER_CLI="$2"; shift 2 ;;
    --outdir)          OUTDIR_CLI="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

# ──────────────────────────────────────────────────────────────────────────────
# Initialize Conda in non-interactive shells
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
# ──────────────────────────────────────────────────────────────────────────────

# 1) Activate utils environment for shyaml
conda activate ssuitslsu-utils

# 2) Load config
CONFIG="$(dirname "$0")/../config/pipeline.yaml"
READS_DIR=$(shyaml get-value reads_dir           < "$CONFIG")
TAXO_XLSX=$(shyaml get-value taxonomy_file       < "$CONFIG")
TAXO_SHEET=$(shyaml get-value taxonomy_sheet     < "$CONFIG")
REF_FASTA=$(shyaml get-value ref_fasta           < "$CONFIG")
THREADS=$(shyaml get-value threads               < "$CONFIG")
MEM_GB=$(shyaml get-value mem_gb                 < "$CONFIG")
OUTDIR=$(shyaml get-value outdir                 < "$CONFIG")
auto_subsample=$(shyaml get-value auto_subsample < "$CONFIG")
MAPQ=$(shyaml get-value mapq < "$CONFIG")
max_coverage=$(shyaml get-value max_coverage < "$CONFIG")
target_coverage=$(shyaml get-value target_coverage < "$CONFIG")

# ──────────────────────────────────────────────────────────────────────────────
# Apply CLI overrides (if provided)
[[ -n "$CONFIG_CLI"     ]] && CONFIG="$CONFIG_CLI"
[[ "$SKIP_TRIM_CLI" == true ]] && SKIP_TRIM=true
[[ -n "$MEM_CLI"        ]] && MEM="$MEM_CLI"
[[ -n "$THREADS_CLI"    ]] && THREADS="$THREADS_CLI"
[[ -n "$ASSEMBLER_CLI"  ]] && ASSEMBLER="$ASSEMBLER_CLI"
[[ -n "$MAPQ_CLI"       ]] && MAPQ="$MAPQ_CLI"
[[ -n "$OUTDIR_CLI"     ]] && OUTDIR="$OUTDIR_CLI"
# ──────────────────────────────────────────────────────────────────────────────

# fall back to defaults if not set in pipeline.yaml
auto_subsample=${auto_subsample:-true}
max_coverage=${max_coverage:-100}
target_coverage=${target_coverage:-90}
MAPQ=${MAPQ:-0}

# Strip stray quotes (so paths with spaces work)
for var in READS_DIR TAXO_XLSX REF_FASTA; do
  val="${!var}"
  val="${val#\"}"; val="${val%\"}"
  val="${val#\'}"; val="${val%\'}"
  printf -v "$var" '%s' "$val"
done

echo "Starting ssuitslsu pipeline at $(date)"

# ─── Download & build filtered reference if missing ─────────────────────────
DATA_DIR="$(dirname "$0")/../data"
mkdir -p "$DATA_DIR"

DB_NAME="General_EUK_longread_v1.9.4"
ZIP_URL="https://sisu.ut.ee/wp-content/uploads/sites/643/${DB_NAME}.zip"
ZIP_FILE="$DATA_DIR/${DB_NAME}.zip"
RAW_FASTA="$DATA_DIR/${DB_NAME}.fasta"
FILTERED_FASTA="$DATA_DIR/${DB_NAME}_nomito_Fungi.fasta"

if [[ ! -f "$FILTERED_FASTA" ]]; then
  echo "[$(date)] Downloading Eukaryome database…"
  curl -L "$ZIP_URL" -o "$ZIP_FILE"

  echo "[$(date)] Extracting FASTA from ZIP…"
  unzip -p "$ZIP_FILE" "${DB_NAME}.fasta" > "$RAW_FASTA"

  echo "[$(date)] Filtering to k__Fungi, excluding mitochondrial sequences…"
  awk '
    /^>/ { keep = ($0 ~ /k__Fungi/ && tolower($0) !~ /mitochondri/) }
    keep
  ' "$RAW_FASTA" > "$FILTERED_FASTA"

  echo "[$(date)] Cleaning up…"
  rm -f "$ZIP_FILE" "$RAW_FASTA"
fi

# choose user-provided reference if valid, otherwise use our filtered download
if [[ -n "$REF_FASTA" && -f "$REF_FASTA" ]]; then
  echo "[$(date)] Using user-provided reference: $REF_FASTA"
else
  REF_FASTA="$FILTERED_FASTA"
  echo "[$(date)] Using auto-downloaded reference: $REF_FASTA"
fi


# Apply CLI overrides
if [[ "$SKIP_TRIM_CLI" == true ]]; then
  SKIP_TRIM=true
else
  SKIP_TRIM=$(shyaml get-value skip_trimming < "$CONFIG")
fi
[[ -n "$MEM_CLI"     ]] && MEM_GB="$MEM_CLI"
[[ -n "$THREADS_CLI" ]] && THREADS="$THREADS_CLI"

# ─── Assembler default & validation ────────────────────────────────────────
ASSEMBLER=${ASSEMBLER:-$(shyaml get-value assembler < "$CONFIG")}

if [[ "$ASSEMBLER" != "spades" && "$ASSEMBLER" != "megahit" ]]; then
  echo "ERROR: --assembler must be 'spades' or 'megahit' (you passed '$ASSEMBLER')" >&2
  exit 1
fi

# ─── 3) Validate config ─────────────────────────────────────────────
die(){ echo "ERROR: $*" >&2; exit 1; }

# ensure reads / taxonomy / reference exist
[[ -d "$READS_DIR" ]]  || die "reads_dir not found: $READS_DIR"
[[ -f "$TAXO_XLSX" ]]  || die "taxonomy_file not found: $TAXO_XLSX"
[[ -f "$REF_FASTA" ]]  || die "ref_fasta not found: $REF_FASTA"


echo "Running with:"
echo "  skip_trimming    = $SKIP_TRIM"
echo "  threads          = $THREADS"
echo "  mem_gb           = $MEM_GB"
echo "  auto_subsample   = $auto_subsample"
echo "  max_coverage     = ${max_coverage}×"
echo "  target_coverage  = ${target_coverage}×"
echo "  mapping quality treshold = ${MAPQ}"
echo

mkdir -p "$OUTDIR"
TAX_CSV="$OUTDIR/$(basename "${TAXO_XLSX%.*}").csv"

# 4) Taxonomy conversion (auto-detect sheet & header) #only once
if [[ ! -f "$TAX_CSV" || "$TAXO_XLSX" -nt "$TAX_CSV" ]]; then
  conda activate ssuitslsu-taxo
  echo "[`date`] Converting taxonomy → CSV"
  PYTHON_EXE="$CONDA_PREFIX/bin/python"
  [[ ! -x "$PYTHON_EXE" ]] && PYTHON_EXE=python

  "$PYTHON_EXE" <<PYCODE
import sys, os, pandas as pd

tax_file   = "${TAXO_XLSX}"
out_csv    = "${TAX_CSV}"
requested  = "${TAXO_SHEET}"
wanted     = [
  'Sample ID','Phylum','Class','Order',
  'Family','Subfamily','Tribe',
  'Genus','Species','Subspecies'
]
ext = os.path.splitext(tax_file)[1].lower()

if ext in ('.xlsx','.xls'):
    xlsx   = pd.ExcelFile(tax_file, engine='openpyxl')
    sheets = xlsx.sheet_names
    # exact match, else substring match, else first sheet
    sheet = next((s for s in sheets if s.lower()==requested.lower()), None)
    if sheet is None:
        sheet = next((s for s in sheets if requested.lower() in s.lower()), sheets[0])
        sys.stderr.write("Using sheet '%s'\\n" % sheet)
    # read without header to find the header row where Sample ID appears in column A
    raw = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    hdr = next((i for i,v in raw.iloc[:,0].astype(str).str.strip().items() if v=='Sample ID'), None)
    if hdr is None:
        raise ValueError("Cannot find header row with 'Sample ID'")
    df = pd.read_excel(xlsx, sheet_name=sheet, header=hdr)
elif ext == '.csv':
    df = pd.read_csv(tax_file)
elif ext in ('.tsv','.txt'):
    df = pd.read_csv(tax_file, sep='\\t')
else:
    raise ValueError("Unsupported taxonomy extension: %s" % ext)

missing = [c for c in wanted if c not in df.columns]
if missing:
    raise ValueError("Missing columns: %s" % missing)

df[wanted].to_csv(out_csv, index=False)
PYCODE

else
  echo "[`date`] Skipping taxonomy conversion (up to date)"
fi

# 5) Define read‐suffix arrays
R1_SUFFIXES=(
  _1.fastq.gz _1.fastq _1.fq.gz _1.fq
  _R1.fastq.gz _R1.fastq _R1.fq.gz _R1.fq
  _1.fastq.bz2 _1.fq.bz2 _R1.fastq.bz2 _R1.fq.bz2
)
R2_SUFFIXES=(
  _2.fastq.gz _2.fastq _2.fq.gz _2.fq
  _R2.fastq.gz _R2.fastq _R2.fq.gz _R2.fq
  _2.fastq.bz2 _2.fq.bz2 _R2.fastq.bz2 _R2.fq.bz2
)

# 6) Main loop
for fp in "$READS_DIR"/*; do
  for suf in "${R1_SUFFIXES[@]}"; do
    if [[ "$fp" == *"$suf" ]]; then
      sample="${fp##*/}"; sample="${sample%$suf}"; R1="$fp"

      # find R2
      R2=""
      for sfx in "${R2_SUFFIXES[@]}"; do
        [[ -e "$READS_DIR/${sample}${sfx}" ]] && { R2="$READS_DIR/${sample}${sfx}"; break; }
      done
      [[ -z "$R2" ]] && echo "⚠️ No R2 for $sample; skipping" >&2 && continue 2

      SAMPLE_DIR="$OUTDIR/$sample"; mkdir -p "$SAMPLE_DIR"
      start_ts=$(date +"%Y-%m-%d %H:%M:%S")
      echo "[`date`] ===> Processing $sample"; t0=$(date +%s)

      # ─── Reference Selection ────────────────────────────────────────────────
      t_ref_start=$(date +%s)
      taxo_line=$(grep -m1 "^${sample}," "$TAX_CSV" || true)
      [[ -z "$taxo_line" ]] && echo "!! No taxonomy for $sample; skipping" >&2 && continue 2
      IFS=, read -r _ phylum cls ord fam subfam tribe gen spp subspp <<<"$taxo_line"
      for v in phylum cls ord fam subfam tribe gen spp subspp; do
        eval "$v=\$(echo \"\${$v}\" | sed -e 's/^[^_]*__//' -e 's/^ *//;s/ *\$//')"
      done
      [[ "$spp"    == *" "* ]] && spp="${spp##* }"
      [[ "$subspp" == *" "* ]] && subspp="${subspp##* }"
      tags=()
      [[ -n "$spp"    ]] && tags+=( "g__${gen};s__${spp}" )
      [[ -n "$subspp" ]] && tags+=( "s__${subspp}" )
      [[ -n "$gen"    ]] && tags+=( "g__${gen}" )
      [[ -n "$tribe"  ]] && tags+=( "t__${tribe}" )
      [[ -n "$subfam" ]] && tags+=( "sf__${subfam}" )
      [[ -n "$fam"    ]] && tags+=( "f__${fam}" )
      [[ -n "$ord"    ]] && tags+=( "o__${ord}" )
      [[ -n "$cls"    ]] && tags+=( "c__${cls}" )
      [[ -n "$phylum" ]] && tags+=( "p__${phylum}" )

      SAMPLE_REF=""
      for tag in "${tags[@]}"; do
        safe="${tag//;/__}"
        out="$SAMPLE_DIR/${sample}_ref_${safe}.fasta"
        if [[ "$tag" == *";"* ]]; then
          gtag="${tag%%;*}"; stag="${tag##*;}"
          awk -v g="$gtag" -v s="$stag" 'BEGIN{RS=">";ORS=""} $0~g&&$0~s{print ">"$0}' "$REF_FASTA" > "$out"
        else
          awk -v t="$tag" 'BEGIN{RS=">";ORS=""} $0~("(^|;)"t"(;|$)"){print ">"$0}' "$REF_FASTA" > "$out"
        fi
        if [[ -s "$out" ]]; then
          echo "    → using reference for level '$tag'"
          SAMPLE_REF="$out"
          break
        else
          rm -f "$out"
        fi
      done
      [[ -z "$SAMPLE_REF" ]] && echo "!! No reference found; skipping" >&2 && continue 2
      t_ref_end=$(date +%s)

      # ─── Trimming ────────────────────────────────────────────────
      t_trim_start=$(date +%s)

      # remember raw inputs
      RAW_R1="$R1"
      RAW_R2="$R2"

      out_fw="$SAMPLE_DIR/${sample}_trimmed_1P.fastq.gz"
      out_rev="$SAMPLE_DIR/${sample}_trimmed_2P.fastq.gz"
      un_fw="$SAMPLE_DIR/${sample}_trimmed_1U.fastq.gz"
      un_rev="$SAMPLE_DIR/${sample}_trimmed_2U.fastq.gz"

      if [[ "$SKIP_TRIM" == "true" ]]; then
        echo "[`date`] skip_trimming; continuing with raw FASTQs"
        # leave R1/R2 pointing at original RAW_R1/RAW_R2
        R1="$RAW_R1"
        R2="$RAW_R2"

      else
        # only run fastp if trimmed outputs do *not* exist
        if [[ ! -f "$out_fw" ]] || [[ ! -f "$out_rev" ]]; then
          conda activate ssuitslsu-fastp
          echo "[`date`] Trimming reads for $sample with fastp"
          fastp \
            -i "$R1" -I "$R2" \
            -o "$out_fw" -O "$out_rev" \
            --unpaired1 "$un_fw" --unpaired2 "$un_rev" \
            --detect_adapter_for_pe \
            --length_required 70 \
            --thread "$THREADS" \
            --html "$SAMPLE_DIR/${sample}_fastp.html" \
            --json "$SAMPLE_DIR/${sample}_fastp.json"

          # after trimming, point R1/R2 at the new files
          R1="$out_fw"
          R2="$out_rev"

        else
          echo "[`date`] Skipping fastp (outputs exist)"
          # point R1/R2 at existing trimmed files
          R1="$out_fw"
          R2="$out_rev"
        fi
      fi

      t_trim_end=$(date +%s)


      # ─── Mapping + auto‐subsampling ──────────────────────────────────────────────
      t_map_start=$(date +%s)
      BAM="$SAMPLE_DIR/${sample}.sorted.bam"
      sub1="$SAMPLE_DIR/${sample}.sub_1P.fastq.gz"
      sub2="$SAMPLE_DIR/${sample}.sub_2P.fastq.gz"

      conda activate ssuitslsu-mapping
      # ensure we have a python binary to run our inline script
      PYTHON_EXE="${CONDA_PREFIX}/bin/python"
      [[ ! -x "$PYTHON_EXE" ]] && PYTHON_EXE=python

      # 1) initial mapping if needed
      if [[ ! -f "$BAM" ]]; then
        # build the index if needed
        [[ ! -f "${SAMPLE_REF}.bwt.2bit.64" ]] && bwa-mem2 index "$SAMPLE_REF"

        echo "[$(date)] Mapping reads (MAPQ ≥ $MAPQ)"
        bwa-mem2 mem -v0 -t "$THREADS" "$SAMPLE_REF" "$R1" "$R2" \
          2> "$SAMPLE_DIR/mapping.log" \
        | samtools view -b -@ "$THREADS" -F 4 -q "$MAPQ" \
        | samtools sort -@ "$THREADS" -o "$BAM"

        samtools index "$BAM"
      else
        echo "[$(date)] Skipping mapping (BAM exists)"
        [[ ! -f "${BAM}.bai" ]] && samtools index "$BAM"
      fi

      # fail early if no alignments
      [[ ! -s "$BAM" ]] && echo "!! Empty BAM; skipping sample" && continue 2

      # 2) coverage check + optional capping
      subbam="${SAMPLE_DIR}/${sample}.depthcapped.bam"
      if [[ ! -f "${subbam}" ]]; then
        # Always compute and display mean coverage for reporting
        mean_cov=$(samtools depth -a "${BAM}" \
                   | awk '{sum+=$3; cnt++} END {print (cnt? sum/cnt : 0)}')
        printf "[%s] Mean coverage: %.1f×\n" "$(date)" "${mean_cov}"

        # Check if auto_subsample is enabled
        if [[ "$auto_subsample" == "true" ]]; then
          if (( $(echo "${mean_cov} > ${max_coverage}" | bc -l) )); then
            echo "[$(date)] auto_subsample=true: Coverage ≥ ${max_coverage}×, downsampling to ${target_coverage}…"

"${PYTHON_EXE}" <<PYCODE
import pysam, random

IN_BAM  = "${BAM}"
OUT_BAM = "${subbam}"
TARGET  = ${target_coverage}

# 1) Build a contig‐aware depth map
depth = {}
with pysam.AlignmentFile(IN_BAM, "rb") as bam:
    for read in bam.fetch():
        rname = read.reference_name
        for pos in range(read.reference_start, read.reference_end):
            depth[(rname, pos)] = depth.get((rname, pos), 0) + 1

# 2) Decide keep‐probability per read
keep_names = set()
with pysam.AlignmentFile(IN_BAM, "rb") as bam:
    for read in bam.fetch():
        rname = read.reference_name
        # for each base, compute local “keep” ratio = min(1, TARGET/depth)
        ratios = [
            min(1.0, TARGET / depth[(rname, p)])
            for p in range(read.reference_start, read.reference_end)
        ]
        # average ratio across the read = probability to KEEP it
        p_keep = sum(ratios) / len(ratios)
        if random.random() < p_keep:
            keep_names.add(read.query_name)

# 3) Write out exactly those reads
with pysam.AlignmentFile(IN_BAM, "rb") as bam, \
     pysam.AlignmentFile(OUT_BAM, "wb", template=bam) as out:
    for read in bam.fetch():
        if read.query_name in keep_names:
            out.write(read)
PYCODE

        else
          echo "[$(date)] auto_subsample=true: Coverage ≤ ${max_coverage}×; no downsampling needed (copying original)."
          cp "${BAM}" "${subbam}"
        fi
      else
        echo "[$(date)] auto_subsample=false: Skipping coverage-based downsampling (copying original)."
        cp "${BAM}" "${subbam}"
      fi

      samtools index "${subbam}"
      BAM="${subbam}"
    else
      echo "[$(date)] Depth-capped BAM exists; skipping coverage analysis."
      BAM="${subbam}"
    fi


      # ─── Extract mapped reads ─────────────────────────────────────────────────
      P1="$SAMPLE_DIR/${sample}.mapped_1.fastq.gz"
      P2="$SAMPLE_DIR/${sample}.mapped_2.fastq.gz"
      U0="$SAMPLE_DIR/${sample}.mapped_0.fastq.gz"

      # re-extract if missing, empty, or older than the BAM
      if [[ ! -s "$P1" || ! -s "$P2" || "$BAM" -nt "$P1" || "$BAM" -nt "$P2" ]]; then
        echo "[$(date)] Extracting properly-paired primary alignments"
        rm -f "$P1" "$P2" "$U0"
        
        # Extract only properly paired reads to avoid mate-pair mismatches
        samtools view -h -@ "$THREADS" \
          -b \
          -f 1 -f 2 \
          -F 4 -F 8 -F 256 -F 2048 \
          -q "$MAPQ" \
          "$BAM" |
        samtools collate -uO -@ "$THREADS" - |
        samtools fastq -@ "$THREADS" \
          -1 "$P1" \
          -2 "$P2" \
          -0 /dev/null \
          -s /dev/null

        # Create empty singletons file for compatibility
        touch "$U0"

        # sanity check: ensure mate-pairs match
        p1=$(zgrep -c '^@' "$P1")
        p2=$(zgrep -c '^@' "$P2")
        if [[ "$p1" -ne "$p2" ]]; then
          echo "[ERROR] mate-pair counts mismatch: R1=$p1 vs R2=$p2" >&2
          exit 1
        fi
        echo "[$(date)] Extracted $p1 paired reads each in P1/P2"

      else
        echo "[$(date)] Mapped FASTQs up-to-date; skipping extraction"
      fi

      t_map_end=$(date +%s)


      # ─── Assembly ────────────────────────────────────────────────────────
      t_asm_start=$(date +%s)
      ASM_DIR="$SAMPLE_DIR/assembly"
      mkdir -p "$ASM_DIR"

      if [[ "$ASSEMBLER" == "spades" ]]; then
        SP_DIR="$ASM_DIR/spades"
        CONTIG="$SP_DIR/contigs.fasta"

        if [[ ! -f "$CONTIG" ]]; then
          conda activate ssuitslsu-spades
          echo "[`date`] Assembling with SPAdes"

          # 1) clear any old run
          rm -rf "$SP_DIR"
          mkdir -p "$SP_DIR"

          # 2) count true singletons and decide on the -s flag
          if [[ -f "$U0" ]]; then
            n0=$(zgrep -c '^@' "$U0" || echo 0)
          else
            n0=0
          fi

          if (( n0 > 0 )); then
            echo "→ Including $n0 singletons"
            SING_OPTS=(-s "$U0")
          else
            echo "→ No singletons; not adding -s"
            SING_OPTS=()
          fi

          # 3) build and run the SPAdes command
          cmd=(
            spades.py --careful
              -1 "$P1" -2 "$P2"
              "${SING_OPTS[@]}"
              -t "$THREADS" -m "$MEM_GB"
              -o "$SP_DIR"
          )

          "${cmd[@]}" 2>&1 | tee "$SP_DIR/spades.log"
        else
          echo "[`date`] Skipping SPAdes assembly (found $CONTIG)"
        fi

      elif [[ "$ASSEMBLER" == "megahit" ]]; then
        MH_DIR="$ASM_DIR/megahit"
        CONTIG="$MH_DIR/final.contigs.fa"

        if [[ ! -f "$CONTIG" ]]; then
          conda activate ssuitslsu-megahit
          echo "[`date`] Assembling with MEGAHIT"

          # — delete U0 if it contains no FASTQ records
          if [[ -f "$U0" ]]; then
              if ! zgrep -q '^@' "$U0"; then
                  echo "→ Removing empty singleton"
                  rm -f "$U0"
              else
                  echo "→ Including $(zgrep -c '^@' "$U0") singletons"
              fi
          fi

          # ─── Clear any old output
          rm -rf "$MH_DIR"

          # build and run Megahit 
          cmd=(megahit \
            -1 "$P1" \
            -2 "$P2" \
            --out-dir "$MH_DIR" \
            --num-cpu-threads "$THREADS" \
            --min-contig-len 400 \
            --memory "$MEM_GB" \
            --verbose \
            --keep-tmp-files
          )
          # include singletons if present
          [[ -f "$U0" ]] && cmd+=( -r "$U0" )

          # run and capture the log
          "${cmd[@]}" 2>&1 | tee "$MH_DIR/megahit.log"
        else
          echo "[`date`] Skipping MEGAHIT assembly (exists)"
        fi

      else
        die "Unknown assembler: $ASSEMBLER (must be spades or megahit)"
      fi

      t_asm_end=$(date +%s)


      # ─── Assembly stats (contigs ≥1 kb) + 45S coverage ────────────────────────────────────
      CONTIG_DIR=$(dirname "$CONTIG")
      STATS_MINLEN=1000
      FILTERED_CONTIGS="${CONTIG_DIR}/contigs.${STATS_MINLEN}bp.fasta"
      ASSEMBLY_STATS="${CONTIG_DIR}/assembly_stats.${STATS_MINLEN}bp.txt"

      if [[ -s "$ASSEMBLY_STATS" ]]; then
        echo "[$(date)] Skipping assembly stats (found): $ASSEMBLY_STATS"
      else
        echo "[$(date)] Computing assembly stats (contigs ≥ ${STATS_MINLEN} bp)"
        conda activate ssuitslsu-itsx

        # 1) extract only contigs ≥1 kb
        seqkit seq -m${STATS_MINLEN} "$CONTIG" -o "$FILTERED_CONTIGS"

        # 2) if none survive, skip
        if [[ ! -s "$FILTERED_CONTIGS" ]]; then
          echo "[$(date)] No contigs ≥${STATS_MINLEN} bp; skipping stats."
        else
          # 3) grab the one data line from seqkit stats (tabular, all stats)
          stats_data=$(seqkit stats -a -T "$FILTERED_CONTIGS" | tail -n1)

          NUM_SEQS=$(echo "$stats_data" | cut -f4)
          SUM_LEN=$( echo "$stats_data" | cut -f5)
          MIN_LEN=$( echo "$stats_data" | cut -f6)
          AVG_LEN=$( echo "$stats_data" | cut -f7)
          MAX_LEN=$( echo "$stats_data" | cut -f8)
          N50=$(     echo "$stats_data" | cut -f13)
          L50=$(     echo "$stats_data" | cut -f14)
          GC_PCT=$(  echo "$stats_data" | cut -f18 | tr -d '%')

          # 4) compute mean coverage across the SSU+ITS+LSU reference(s)
          mean_45S_cov=$(samtools depth -a "$BAM" \
                           | awk '{sum+=$3; cnt++} END{ if(cnt) printf("%.1f", sum/cnt); else print "0.0" }')

          {
            echo "CONTIGS-${STATS_MINLEN}BP:    $NUM_SEQS"
            echo "ASSEMBLY_LEN-${STATS_MINLEN}BP:    $SUM_LEN"
            echo "LARGEST_CONTIG:   $MAX_LEN"
            echo "N50-${STATS_MINLEN}BP:           $N50"
            echo "L50-${STATS_MINLEN}BP:           $L50"
            echo "GC-${STATS_MINLEN}BP:            ${GC_PCT}%"
            echo "MEAN_COV-45S:     ${mean_45S_cov}×"
          } | tee "$ASSEMBLY_STATS"
        fi
      fi

      echo


      # ─── ITS extraction ─────────────────────────────────────────────────────────
      t_its_start=$(date +%s)
      ITSX_DIR="${SAMPLE_DIR}/itsx"
      mkdir -p "$ITSX_DIR"

      # define “done” markers so we don’t rerun completed branches
      SMALL_DONE="$ITSX_DIR/${sample}.small.done"
      LARGE_DONE="$ITSX_DIR/${sample}.large_nhmmer.done"

      # if both are done, skip entire ITSx step
      if [[ -f "$SMALL_DONE" && -f "$LARGE_DONE" ]]; then
        echo "[`date`] Skipping ITSx on $sample (already completed both small & large runs)"
      else
        echo "[`date`] Preparing input for ITSx on $sample"
        conda activate ssuitslsu-itsx

        # split contigs by length
        SMALL_FASTA="$ITSX_DIR/${sample}.small.fasta"
        LARGE_FASTA="$ITSX_DIR/${sample}.large.fasta"
        : > "$SMALL_FASTA"
        : > "$LARGE_FASTA"

        awk -v L=100000 -v out_small="$SMALL_FASTA" -v out_large="$LARGE_FASTA" '
          BEGIN { RS=">"; ORS="" }
          NR>1 {
            header = $1
            seq = ""
            for (i=2; i<=NF; i++) seq = seq $i
            if (length(seq) > L) {
              printf(">%s\n%s\n", header, seq) > out_large
            } else {
              printf(">%s\n%s\n", header, seq) > out_small
            }
          }
        ' "$CONTIG"

        # 1) run ITSx on the small contigs (default hmmsearch)
        if [[ -s "$SMALL_FASTA" && ! -f "$SMALL_DONE" ]]; then
          echo "[`date`] Running ITSx (hmmsearch) on small contigs for $sample"
          ITSx -i "$SMALL_FASTA" \
               -o "${ITSX_DIR}/${sample}.small" \
               --preserve T \
               --save_regions all \
               --anchor HMM \
               --region all \
               --cpu "$THREADS" \
            && touch "$SMALL_DONE"
        else
          if [[ -f "$SMALL_DONE" ]]; then
            echo "[`date`] Skipping small-contig ITSx (already done)"
          else
            echo "[`date`] No small contigs to process"
          fi
        fi

        # 2) run ITSx on the large contigs (nhmmer)
        if [[ -s "$LARGE_FASTA" && ! -f "$LARGE_DONE" ]]; then
          echo "[`date`] Running ITSx (nhmmer) on large contigs for $sample"
          ITSx -i "$LARGE_FASTA" \
               -o "${ITSX_DIR}/${sample}.large_nhmmer" \
               --preserve T \
               --save_regions all \
               --anchor HMM \
               --region all \
               --nhmmer T \
               --cpu "$THREADS" \
            && touch "$LARGE_DONE"
        else
          if [[ -f "$LARGE_DONE" ]]; then
            echo "[`date`] Skipping large-contig ITSx (already done)"
          else
            echo "[`date`] No large contigs (>100 kb) to process"
          fi
        fi
      fi

      t_its_end=$(date +%s)

      # ─── Phylogenetic analysis ────────────────────────────────────────────────
      PHYLO_DIR="$SAMPLE_DIR/phylogeny"
      mkdir -p "$PHYLO_DIR"

      # 1) Gather sample sequences: all ITSx regions + contigs ≥400 bp
      SAMPLE_SEQ_FASTA="$PHYLO_DIR/${sample}_sequences.fasta"
      : > "$SAMPLE_SEQ_FASTA"  # truncate/create

      # a) ITSx outputs
      for size in small large_nhmmer; do
        for region in full ITS1 ITS2 SSU LSU 5_8S; do
          src="$ITSX_DIR/${sample}.${size}.${region}.fasta"
          if [[ -s "$src" ]]; then
            # log what we’re adding
            printf "[%s] Adding ITSx %-4s (%s) → %s\n" "$(date)" "$region" "$size" "$src"
            # prefix headers with sample, size, and region
            awk -v smp="$sample" -v sz="$size" -v rgn="$region" '
              /^>/ { sub(/^>/, ">" smp "_" sz "_" rgn "_"); print; next }
              { print }
            ' "$src" >> "$SAMPLE_SEQ_FASTA"
          fi
        done
      done

      # b) assembled contigs ≥400 bp
      printf "[%s] Adding contigs ≥400 bp from %s\n" "$(date)" "$CONTIG"
      conda activate ssuitslsu-itsx
      seqkit seq -m400 "$CONTIG" >> "$SAMPLE_SEQ_FASTA"

      # 2) Build genus‐level reference
      REF_GEN_FASTA="$PHYLO_DIR/${gen}_ref.fasta"
      printf "[%s] Building genus reference for g__%s\n" "$(date)" "$gen"
      awk -v tag="g__${gen}" '
        BEGIN { RS=">"; ORS="" }
        $0 ~ ("(^|;)" tag "(;|$)") { print ">" $0 }
      ' "$REF_FASTA" > "$REF_GEN_FASTA"

      # 3) Combine sample + reference into one FASTA
      PHYLO_FASTA="$PHYLO_DIR/${gen}_phylo.fasta"
      cat "$REF_GEN_FASTA" "$SAMPLE_SEQ_FASTA" > "$PHYLO_FASTA"
      printf "[%s] Combined phylo FASTA: %s (%d refs + %d sample seqs)\n\n" \
        "$(date)" "$PHYLO_FASTA" \
        "$(grep -c '^>' "$REF_GEN_FASTA")" \
        "$(grep -c '^>' "$SAMPLE_SEQ_FASTA")"

      # 4) MAFFT alignment
      PHYLO_ALN="$PHYLO_DIR/${gen}.aln"
      if [[ -s "$PHYLO_ALN" ]]; then
        printf "[%s] Skipping MAFFT (alignment exists): %s\n\n" "$(date)" "$PHYLO_ALN"
      else
        t_align_start=$(date +%s)
        conda activate ssuitslsu-mafft
        printf "[%s] Running MAFFT → %s\n" "$(date)" "$PHYLO_ALN"
        mafft --thread "$THREADS" --auto --adjustdirection "$PHYLO_FASTA" > "$PHYLO_ALN"
        t_align_end=$(date +%s)
        printf "\n"
      fi

      # 5) IQ-TREE ML inference
      PHYLO_PREFIX="$PHYLO_DIR/${gen}"
      TREEFILE="${PHYLO_PREFIX}.treefile"
      if [[ -f "${PHYLO_PREFIX}.ckp.gz" ]]; then
        printf "[%s] Skipping IQ-TREE (tree exists): %s\n\n" "$(date)" "$TREEFILE"
      else
        t_phylo_start=$(date +%s)
        conda activate ssuitslsu-iqtree
        printf "[%s] Running IQ-TREE → prefix %s\n" "$(date)" "$PHYLO_PREFIX"
        iqtree \
          -s "$PHYLO_ALN" \
          -m MFP \
          -nt AUTO \
          -ntmax "${THREADS}" \
          -mem "${MEM_GB}G" \
          -nt AUTO \
          -bb 1000 \
          -pre "$PHYLO_PREFIX"
        t_phylo_end=$(date +%s)
        printf "\n"
      fi

      # ─── Timing & report ───────────────────────────────────────────────
      t1=$(date +%s)
      end_ts=$(date +"%Y-%m-%d %H:%M:%S")
      elapsed=$((t1 - t0))

      ref_dur=$((t_ref_end - t_ref_start))
      trim_dur=$((t_trim_end - t_trim_start))
      map_dur=$((t_map_end  - t_map_start))
      asm_dur=$((t_asm_end  - t_asm_start))
      its_dur=$((t_its_end  - t_its_start))
      align_dur=$((t_align_end  - t_align_start))
      phylo_dur=$((t_phylo_end  - t_phylo_start))

      printf "→ Sample: %s\n  Started: %s\n  Finished: %s\n  Elapsed: %02dh:%02dm:%02ds\n\n" \
        "$sample" "$start_ts" "$end_ts" \
        $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

      printf "    Taxonomical Rank Selection:   %02dm%02ds\n" $((ref_dur/60)) $((ref_dur%60))
      printf "    Trimming:   %02dm%02ds\n" $((trim_dur/60)) $((trim_dur%60))
      printf "    Mapping:    %02dm%02ds\n" $((map_dur/60)) $((map_dur%60))
      printf "    Assembly:   %02dm%02ds\n" $((asm_dur/60)) $((asm_dur%60))
      printf "    ITSx:       %02dm%02ds\n" $((its_dur/60)) $((its_dur%60))
      printf "    Alignment:  %02dm%02ds\n" $((align_dur/60)) $((align_dur%60))
      printf "    Phylogeny:  %02dm%02ds\n\n" $((phylo_dur/60)) $((phylo_dur%60))


      break
    fi
  done
done

echo "[`date`] ALL SAMPLES COMPLETE."
