#!/usr/bin/env bash

set -euo pipefail

# ──────────────────────────────────────────────────────────────────────────────
# Version checks
echo "➤ Bash version: $BASH_VERSION"
conda --version >/dev/null 2>&1 && echo "➤ Conda version: $(conda --version | cut -d' ' -f2)" \
                             || echo "⚠️  conda not found in PATH"
# ──────────────────────────────────────────────────────────────────────────────

export TMPDIR="/tmp"
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
  -n, --no-trim             Skip fastp trimming
  -m, --max-memory MB       Override memory (GB) for SPAdes
  -t, --threads N           Override number of threads for all steps
      --assembler [spades|megahit]
                            Choose assembler (spades or megahit)
      --mapq N              Override mapping quality filter (samtools -q)
      --outdir DIR          Override output directory
      --reads-dir DIR       Specify directory containing raw read files
      --taxonomy-file FILE  Path to taxonomy file
      --taxonomy-sheet SHEET
                            Sheet name within the taxonomy file
      --ref-fasta FILE      Reference FASTA for mapping/indexing
      --filter-softclip     Enable filtering of soft‐clipped reads (uses automatic quality-based thresholds)
      --softclip-mode MODE    Override soft-clip filter mode (full|trim)
      --auto-subsample [true|false] Override auto_subsample behavior
      --max-cov N             Override max_coverage threshold
      --target-cov N          Override target_coverage threshold
      --skip-phylogeny        Skip IQ-TREE phylogenetic inference

  -h, --help                Show this help message and exit
EOF
}

# ──────────────────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────────────────
# Initialize CLI-override vars with defaults
SKIP_TRIM_CLI=false
SKIP_TRIM=false
MEM=""
THREADS_CLI=""
ASSEMBLER_CLI=""
MAPQ_CLI=""
FILTER_SOFTCLIP_CLI=""
OUTDIR_CLI=""
READS_DIR_CLI=""
TAXONOMY_FILE_CLI=""
TAXONOMY_SHEET_CLI=""
REF_FASTA_CLI=""
SKIP_PHYLO_CLI=false

# ──────────────────────────────────────────────────────────────────────────────

# Parse CLI args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)         usage; exit 0 ;;
    --reads-dir)       READS_DIR_CLI="$2"; shift 2 ;;
    --taxonomy-file)   TAXONOMY_FILE_CLI="$2"; shift 2 ;;
    --taxonomy-sheet)  TAXONOMY_SHEET_CLI="$2"; shift 2 ;;
    --ref-fasta)       REF_FASTA_CLI="$2"; shift 2 ;;
    -n|--no-trim)      SKIP_TRIM_CLI=true; shift ;;
    -m|--max-memory)   MEM_CLI="$2"; shift 2 ;;
    -t|--threads)      THREADS_CLI="$2"; shift 2 ;;
    --mapq)            MAPQ_CLI="$2"; shift 2 ;;
    --filter-softclip) FILTER_SOFTCLIP_CLI=true; shift ;;
    --assembler)       ASSEMBLER_CLI="$2"; shift 2 ;;
    --outdir)          OUTDIR_CLI="$2"; shift 2 ;;
    --softclip-mode)   SOFTCLIP_MODE_CLI="$2";     shift 2 ;;
    --auto-subsample)  AUTO_SUBSAMPLE_CLI="$2";    shift 2 ;;
    --max-cov)         MAX_COV_CLI="$2";          shift 2 ;;
    --target-cov)      TARGET_COV_CLI="$2";       shift 2 ;;
    --skip-phylogeny)  SKIP_PHYLO_CLI=true;        shift ;;

    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

# ──────────────────────────────────────────────────────────────────────────────
# Initialize Conda in non-interactive shells
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
# ──────────────────────────────────────────────────────────────────────────────

# 1) Activate utils environment for shyaml
set +u
conda activate ssuitslsu-utils
set -u

# 2) Load config
CONFIG="$(dirname "$0")/../config/pipeline.yaml"

READS_DIR=$(shyaml get-value reads_dir              < "$CONFIG")
TAXONOMY_FILE=$(shyaml get-value taxonomy_file      < "$CONFIG")
TAXONOMY_SHEET=$(shyaml get-value taxonomy_sheet    < "$CONFIG")
REF_FASTA=$(shyaml get-value ref_fasta              < "$CONFIG")

# trimming
SKIP_TRIM_RAW=$(shyaml get-value skip_trimming < "$CONFIG")

# mapping/filtering
MAPQ=$(shyaml get-value mapq                        < "$CONFIG")
FILTER_SOFTCLIP_RAW=$(shyaml get-value filter_softclipped_reads < "$CONFIG")
SOFTCLIP_MODE=$(shyaml get-value softclip_filter_mode < "$CONFIG")

# coverage
AUTO_SUBSAMPLE_RAW=$(shyaml get-value auto_subsample    < "$CONFIG")

# Optional coverage values - only load if they exist in config
MAX_COV=$(shyaml get-value max_coverage 2>/dev/null < "$CONFIG" || echo "")
TARGET_COV=$(shyaml get-value target_coverage 2>/dev/null < "$CONFIG" || echo "")

# Track if user actually specified these values (for mode selection)
USER_SPECIFIED_COVERAGE=false
if [[ -n "$MAX_COV" && -n "$TARGET_COV" ]]; then
    USER_SPECIFIED_COVERAGE=true
fi


# assembly
ASSEMBLER=$(shyaml get-value assembler              < "$CONFIG")

#phylogeny
SKIP_PHYLO_RAW=$(shyaml get-value skip_phylogeny < "$CONFIG")

# resources
THREADS=$(shyaml get-value threads                  < "$CONFIG")
MEM=$(shyaml get-value mem_gb                       < "$CONFIG")

# output
OUTDIR=$(shyaml get-value outdir                    < "$CONFIG")


SKIP_TRIM="${SKIP_TRIM_RAW,,}"
ASSEMBLER="${ASSEMBLER,,}"
FILTER_SOFTCLIP="${FILTER_SOFTCLIP_RAW,,}"
SOFTCLIP_MODE="${SOFTCLIP_MODE,,}"
AUTO_SUBSAMPLE="${AUTO_SUBSAMPLE_RAW,,}"
SKIP_PHYLO="${SKIP_PHYLO_RAW,,}"

# ─── Apply CLI overrides (if provided) ────────────────────────────────────────
[[ -n "${READS_DIR_CLI:-}"      ]] && READS_DIR="$READS_DIR_CLI"
[[ -n "${TAXONOMY_FILE_CLI:-}"  ]] && TAXONOMY_FILE="$TAXONOMY_FILE_CLI"
[[ -n "${TAXONOMY_SHEET_CLI:-}" ]] && TAXONOMY_SHEET="$TAXONOMY_SHEET_CLI"
[[ -n "${REF_FASTA_CLI:-}"      ]] && REF_FASTA="$REF_FASTA_CLI"

[[ "${SKIP_TRIM_CLI:-}" == "true" ]] && SKIP_TRIM=true

[[ -n "${MAPQ_CLI:-}"           ]] && MAPQ="$MAPQ_CLI"
[[ -n "${FILTER_SOFTCLIP_CLI:-}" ]] && FILTER_SOFTCLIP="$FILTER_SOFTCLIP_CLI"
[[ -n "${SOFTCLIP_MODE_CLI:-}"  ]] && SOFTCLIP_MODE="$SOFTCLIP_MODE_CLI"

[[ -n "${AUTO_SUBSAMPLE_CLI:-}" ]] && AUTO_SUBSAMPLE="$AUTO_SUBSAMPLE_CLI"
# Only override if CLI values are provided
[[ -n "${MAX_COV_CLI:-}" ]] && MAX_COV="$MAX_COV_CLI"
[[ -n "${TARGET_COV_CLI:-}" ]] && TARGET_COV="$TARGET_COV_CLI"

# Update user coverage detection after CLI overrides
USER_SPECIFIED_COVERAGE=false
if [[ -n "$MAX_COV" && -n "$TARGET_COV" ]]; then
    USER_SPECIFIED_COVERAGE=true
fi

[[ -n "${ASSEMBLER_CLI:-}"      ]] && ASSEMBLER="$ASSEMBLER_CLI"

[[ -n "${THREADS_CLI:-}"        ]] && THREADS="$THREADS_CLI"
[[ -n "${MEM_CLI:-}"            ]] && MEM="$MEM_CLI"

[[ -n "${OUTDIR_CLI:-}"         ]] && OUTDIR="$OUTDIR_CLI"

[[ "${SKIP_PHYLO_CLI:-}" == "true" ]] && SKIP_PHYLO=true


# ──────────────────────────────────────────────────────────────────────────────

# fall back to defaults if not set in pipeline.yaml
MAPQ=${MAPQ:-0}
AUTO_SUBSAMPLE=${AUTO_SUBSAMPLE:-true}

# Strip stray quotes (so paths with spaces work)
for var in READS_DIR TAXONOMY_FILE REF_FASTA; do
  val="${!var:-}"
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

REFERENCE_REMOVE_FILE="$DATA_DIR/reference_to_be_removed.txt"

# Check if we need to rebuild the filtered reference:
# - If filtered FASTA doesn't exist, OR
# - If reference_to_be_removed.txt is newer than the filtered FASTA
if [[ ! -f "$FILTERED_FASTA" ]]; then
  echo "[$(date)] Filtered reference database not found, building it…"
elif [[ -f "$REFERENCE_REMOVE_FILE" && "$REFERENCE_REMOVE_FILE" -nt "$FILTERED_FASTA" ]]; then
  echo "[$(date)] Reference removal file has been updated, rebuilding filtered database…"
fi

# Check if user provided a custom reference database
if [[ -n "$REF_FASTA" && -f "$REF_FASTA" ]]; then
  echo "[$(date)] Using user-provided reference: $REF_FASTA"
else
  # Only download and process the eukaryome database if no custom reference is provided
  REF_FASTA="$FILTERED_FASTA"
  echo "[$(date)] Using auto-downloaded reference: $REF_FASTA"

  if [[ ! -f "$FILTERED_FASTA" ]] || [[ -f "$REFERENCE_REMOVE_FILE" && "$REFERENCE_REMOVE_FILE" -nt "$FILTERED_FASTA" ]]; then
    # Download and extract the database
    echo "[$(date)] Downloading Eukaryome database…"
    curl -L "$ZIP_URL" -o "$ZIP_FILE"

    echo "[$(date)] Extracting FASTA from ZIP…"
    unzip -p "$ZIP_FILE" "${DB_NAME}.fasta" > "$RAW_FASTA"

    echo "[$(date)] Filtering to k__Fungi, excluding mitochondrial sequences…"
    awk '
      /^>/ { keep = ($0 ~ /k__Fungi/ && tolower($0) !~ /mitochondri/) }
      keep
    ' "$RAW_FASTA" > "$FILTERED_FASTA"

    # ─── Remove all "(Fungi)" from headers ────────────────────────────────────
    echo "[$(date)] Stripping \"(Fungi)\" suffix from genus names…"
    sed -i'' -e 's/(Fungi)//g' "$FILTERED_FASTA"
    # ───────────────────────────────────────────────────────────────────────────

    echo "[$(date)] Pruning out unwanted sequences listed in /data/references_to_be_removed.txt…"
    awk '
      BEGIN {
        # read IDs (one per line, without the ">") into remove[]
        while ((getline id < "'"$DATA_DIR"'/reference_to_be_removed.txt") > 0) {
          # Extract just the ID part (before the first semicolon) from remove list
          semicolon_pos = index(id, ";")
          if (semicolon_pos > 0)
            id = substr(id, 1, semicolon_pos - 1)
          remove[id] = 1
        }
      }
      /^>/ {
        seq_id = substr($0, 2)
        # Remove underscore prefix if present
        if (substr(seq_id, 1, 1) == "_")
          seq_id = substr(seq_id, 2)
        # Extract just the ID part (before the first semicolon)
        semicolon_pos = index(seq_id, ";")
        if (semicolon_pos > 0)
          seq_id = substr(seq_id, 1, semicolon_pos - 1)
        skip = (seq_id in remove)
        if (!skip) print
        next
      }
      { if (!skip) print }
    ' "$FILTERED_FASTA" > "${FILTERED_FASTA%.fasta}_pruned.fasta"

    mv "${FILTERED_FASTA%.fasta}_pruned.fasta" "$FILTERED_FASTA"

    echo "[$(date)] Cleaning up temporary files…"
    rm -f "$ZIP_FILE" "$RAW_FASTA"
  fi
fi

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
[[ -f "$TAXONOMY_FILE" ]] || die "taxonomy file not found: $TAXONOMY_FILE"
[[ -f "$REF_FASTA" ]]  || die "ref_fasta not found: $REF_FASTA"


echo "Running with:"
echo "  skip_trimming    = $SKIP_TRIM"
echo "  threads          = $THREADS"
echo "  mem_gb           = $MEM"
echo "  max_coverage     = ${MAX_COV}×"
echo "  target_coverage  = ${TARGET_COV}×"
echo "  mapping quality treshold = ${MAPQ}"
echo "  filter_softclipped_reads = $FILTER_SOFTCLIP"
if [[ $FILTER_SOFTCLIP == true ]]; then
  echo "  softclip_filter_mode     = $SOFTCLIP_MODE (automatic quality-based thresholds)"
fi

echo

mkdir -p "$OUTDIR"
TAX_CSV="$OUTDIR/$(basename "${TAXONOMY_FILE%.*}").csv"

# 4) Taxonomy conversion (auto-detect sheet & header) #only once
if [[ ! -f "$TAX_CSV" || "$TAXONOMY_FILE" -nt "$TAX_CSV" ]]; then
  set +u
  conda activate ssuitslsu-taxo
  set -u
  echo "[`date`] Converting taxonomy → CSV"
 PYTHON_EXE="$CONDA_PREFIX/bin/python"
  [[ ! -x "$PYTHON_EXE" ]] && PYTHON_EXE=python

  "$PYTHON_EXE" <<PYCODE
import sys, os, pandas as pd

tax_file   = "${TAXONOMY_FILE}"
out_csv    = "${TAX_CSV}"
requested  = "${TAXONOMY_SHEET}"
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
  _1.fastq.bz2 _1.fq.bz2 _R1.fastq.bz2 _R1.fq.bz2 _R1_001.fastq.gz
)
R2_SUFFIXES=(
  _2.fastq.gz _2.fastq _2.fq.gz _2.fq
  _R2.fastq.gz _R2.fastq _R2.fq.gz _R2.fq
  _2.fastq.bz2 _2.fq.bz2 _R2.fastq.bz2 _R2.fq.bz2 _R2_001.fastq.gz
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

      # Save original stdout/stderr
      exec 3>&1 4>&2

      # Create individual log file for this sample
      SAMPLE_LOG="$SAMPLE_DIR/${sample}.log"
      exec 1> >(tee -a "$SAMPLE_LOG")
      exec 2> >(tee -a "$SAMPLE_LOG" >&2)

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

        # Create reference file
        if [[ "$tag" == *";"* ]]; then
          gtag="${tag%%;*}"; stag="${tag##*;}"
          awk -v g="$gtag" -v s="$stag" 'BEGIN{RS=">";ORS=""} $0~g&&$0~s{print ">"$0}' "$REF_FASTA" > "$out"
        else
          awk -v t="$tag" 'BEGIN{RS=">";ORS=""} $0~("(^|;)"t"(;|$)"){print ">"$0}' "$REF_FASTA" > "$out"
        fi

        # Check if file has content (not just empty or whitespace)
        if [[ -s "$out" ]] && [[ $(grep -c '^>' "$out" 2>/dev/null || echo 0) -gt 0 ]]; then
          echo "    → using reference for level '$tag'"
          SAMPLE_REF="$out"
          PHYLO_REF_TAG="$tag"
          break
        else
          # Remove empty or invalid file
          rm -f "$out"
        fi
      done
      [[ -z "$SAMPLE_REF" ]] && echo "!! No reference found; skipping" >&2 && continue 2
      t_ref_end=$(date +%s)

      # ─── Trimming ────────────────────────────────────────────────
      t_trim_start=$(date +%s)

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
          set +u
          conda activate ssuitslsu-fastp
          set -u
          echo "[`date`] Trimming reads for $sample with fastp"
          fastp \
            -i "$R1" -I "$R2" \
            -o "$out_fw" -O "$out_rev" \
            --unpaired1 "$un_fw" --unpaired2 "$un_rev" \
            --detect_adapter_for_pe \
            --length_required 70 \
            --thread "$THREADS" \
            --html /dev/null \
            --json /dev/null

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

      set +u
      CONDA_BASE=$(conda info --base)
      source "$CONDA_BASE/etc/profile.d/conda.sh"
      conda activate ssuitslsu-minimap
      set -u

      PYTHON_EXE="${CONDA_PREFIX}/bin/python"
      [[ ! -x "$PYTHON_EXE" ]] && PYTHON_EXE=python

      # 1) initial mapping if needed
      if [[ ! -f "$BAM" ]]; then
        # build the minimap2 index if needed
        [[ ! -f "${SAMPLE_REF}.mmi" ]] && minimap2 -d "${SAMPLE_REF}.mmi" "$SAMPLE_REF"

        LOG="$SAMPLE_DIR/mapping.log"

        # remove any stale BAM/index
        rm -f "$BAM" "${BAM}.bai" "$LOG"

        echo "[$(date)] Mapping paired and unpaired reads (MAPQ ≥ $MAPQ)"
        (
          # redirect FD-4 to the log
          exec 4> "$LOG"

          samtools merge -u -@ "$THREADS" - \
            <( minimap2 -t "$THREADS" -ax sr "${SAMPLE_REF}.mmi" "$R1" "$R2" 2>&4 \
                | samtools view -hb -@ "$THREADS" -q "$MAPQ" - ) \
            <( minimap2 -t "$THREADS" -ax sr "${SAMPLE_REF}.mmi" "$un_fw"    2>&4 \
                | samtools view -hb -@ "$THREADS" -q "$MAPQ" - ) \
            <( minimap2 -t "$THREADS" -ax sr "${SAMPLE_REF}.mmi" "$un_rev"   2>&4 \
                | samtools view -hb -@ "$THREADS" -q "$MAPQ" - ) \
          | samtools sort -@ "$THREADS" -o "$BAM"

          # close FD-4
          exec 4>&-
        )

        echo "[$(date)] Indexing final BAM"
        samtools index "$BAM"

      else
        echo "[$(date)] Skipping mapping (BAM exists)"
        [[ ! -f "${BAM}.bai" ]] && samtools index "$BAM"
      fi

      # ─── Ensure FASTA index for IGV ──────────────────────────────────────────────
      # before any mapping steps, index the reference if needed
      if [[ ! -f "${SAMPLE_REF}.fai" ]]; then
        echo "[$(date)] Generating FASTA index for ${SAMPLE_REF}"
        samtools faidx "${SAMPLE_REF}"
      fi
      # ────────────────────────────────────────────────────────────────────────────


      # fail early if no alignments
      [[ ! -s "$BAM" ]] && echo "!! Empty BAM; skipping sample" && continue 2


       # ─── Soft-clip statistics (auto or user) ────────────────────────────────

      set +u
      CONDA_BASE=$(conda info --base)
      source "$CONDA_BASE/etc/profile.d/conda.sh"
      conda activate ssuitslsu-mapping
      set -u
      PYTHON_EXE="${CONDA_PREFIX}/bin/python"

      SC_DIR="$SAMPLE_DIR/chimera/softclipped_reads"
      mkdir -p "$SC_DIR"   # ensure output directory exists
      STAT_FILE="$SC_DIR/${sample}_softclip_stats.txt"

      # Generate stats file for record-keeping (if not already present)
      if [[ ! -f "$STAT_FILE" ]]; then
        echo "[$(date)] Generating soft-clip statistics for record-keeping…"
        echo "[$(date)] Computing soft-clip stats for $sample" > "$STAT_FILE"

        # run the Python stats helper and append its output
        "$PYTHON_EXE" scripts/parallel_softclip_filter.py \
          -i "$BAM" \
          -o "$SC_DIR" \
          -m auto \
          -t "$THREADS" \
          --mode stats \
        | tee -a "$STAT_FILE"
      else
        echo "[$(date)] Soft-clip stats already present; skipping stats generation"
      fi

      # ─── Build IGV BAM of soft-clipped reads ─────────────────
      SC_DIR="$SAMPLE_DIR/chimera/softclipped_reads"
      mkdir -p "$SC_DIR"
      SC_OUT="$SC_DIR/${sample}.softclipped.sorted.bam"

      if [[ ! -f "$SC_OUT" ]]; then
        echo "[$(date)] Extracting soft-clipped reads into per-bucket BAMs…"
        "$PYTHON_EXE" scripts/parallel_softclip_filter.py \
          -i "$BAM" \
          -o "$SC_DIR" \
          -m auto \
          -t "$THREADS" \
          --mode igv

        echo "[$(date)] Merging and sorting IGV BAM…"
        # merge all the bucket BAMs, then one sort+index
        samtools merge -@ "$THREADS" -u -f - "$SC_DIR"/igv_*.bam \
          | samtools sort -@ "$THREADS" -o "$SC_OUT" -
        samtools index "$SC_OUT"

        echo "[$(date)] Cleaning up temporary IGV BAM chunks…"
        rm "$SC_DIR"/igv_*.bam
      fi


      # ─── Only perform soft-clip filtering/trimming if the flag is enabled ────
      if [[ $FILTER_SOFTCLIP == true ]]; then
        echo "[$(date)] filter_softclipped_reads = $FILTER_SOFTCLIP"
        echo "[$(date)] softclip_filter_mode     = $SOFTCLIP_MODE (automatic quality-based thresholds)"

        # check if we already produced a fixmate’d BAM
        CS_BAM="$SAMPLE_DIR/${sample}.fixmate.bam"
        if [[ -f "$CS_BAM" ]]; then
          echo "[$(date)] $CS_BAM already exists; skipping soft-clip and mate-fixing."
          BAM="$CS_BAM"
        else
          # ─── Perform the chosen soft-clip operation ────────────────────────
          if [[ "$SOFTCLIP_MODE" == "full" || "$SOFTCLIP_MODE" == "trim" ]]; then
            if [[ "$SOFTCLIP_MODE" == "full" ]]; then
              echo "[$(date)] Removing reads with low-quality artifactual soft-clips…"
            else
              echo "[$(date)] Trimming low-quality artifactual soft-clipped ends…"
            fi

            PROCESSED_BAM="$SAMPLE_DIR/${sample}.filtered.softclip.bam"

            # 1) Extract and save the SAM header
            samtools view -@ "$THREADS" -H "$BAM" > "$SC_DIR/header.sam"

            # 2) Launch Python filter/trimmer in parallel
            "$PYTHON_EXE" scripts/parallel_softclip_filter.py \
              -i "$BAM" \
              -o "$SC_DIR" \
              -m auto \
              -t "$THREADS" \
              --mode "$SOFTCLIP_MODE"

            # 3) Merge header + filtered SAM fragments, convert to BAM
            cat "$SC_DIR/header.sam" "$SC_DIR"/filter_"$SOFTCLIP_MODE"_*.sam \
              | samtools view -@ "$THREADS" -bS -o "$PROCESSED_BAM" -

            rm "$SC_DIR/header.sam" "$SC_DIR"/filter_"$SOFTCLIP_MODE"_*.sam

            BAM="$PROCESSED_BAM"
            echo "[$(date)] Soft-clip step done; downstream BAM is $BAM"
          fi

          # ─── Fix pairing/info flags and sort/index ───────────────────────────
          echo "[$(date)] Fixing mate flags on ${BAM} …"
          NS_BAM="$SAMPLE_DIR/${sample}.fix.ns.bam"
          FM_BAM="$SAMPLE_DIR/${sample}.fixmate.ns.bam"

          # name-sort for fixmate
          samtools sort -n -@ "$THREADS" -o "$NS_BAM" "$BAM"
          # fix mate-pair tags
          samtools fixmate -m -@ "$THREADS" "$NS_BAM" "$FM_BAM"
          # drop unpaired, coordinate-sort, and index final
          samtools view -@ "$THREADS" -b -F 4 "$FM_BAM" \
            | samtools sort -@ "$THREADS" -o "$CS_BAM"
          samtools index "$CS_BAM"

          rm -f "$NS_BAM" "$FM_BAM"

          BAM="$CS_BAM"
          echo "[$(date)] Finished fixmate; downstream BAM is $BAM"
        fi  # end existing-CS_BAM check

        echo "[$(date)] Soft-clip analysis ends"
      fi  # end FILTER_SOFTCLIP guard


      # ─── Coverage check + intelligent subsampling ────────────────────────────────
      # Two-tier strategy:
      # - Conservative by default: gentle reduction, assembly-safe, fast
      # - Aggressive on demand: per-base precision, slower but preserves low-coverage regions
      subbam="${SAMPLE_DIR}/${sample}.depthcapped.bam"
      SUBSAMPLING_OCCURRED=false

      if [[ ! -f "${subbam}" ]]; then
        # Always compute and display mean coverage for reporting
        mean_cov=$(samtools depth -a "${BAM}" \
                   | awk '{sum+=$3; cnt++} END {print (cnt? sum/cnt : 0)}')
        printf "[%s] Mean coverage: %.1f×\n" "$(date)" "${mean_cov}"

        # Define conservative defaults optimized for assembly
        CONSERVATIVE_MAX=1000     # Only subsample if coverage is extreme
        CONSERVATIVE_TARGET=500   # Gentle reduction, still assembly-friendly

        # Use user values if provided, otherwise use conservative defaults
        EFFECTIVE_MAX_COV=${MAX_COV:-$CONSERVATIVE_MAX}
        EFFECTIVE_TARGET_COV=${TARGET_COV:-$CONSERVATIVE_TARGET}

        # Determine if subsampling is needed
        if [[ "$AUTO_SUBSAMPLE" == "true" && $(echo "$mean_cov > $EFFECTIVE_MAX_COV" | bc -l) -eq 1 ]]; then

          # Check if user specified custom values (triggers per-base mode)
          if [[ "$USER_SPECIFIED_COVERAGE" == "true" ]]; then
            echo "[$(date)] 🎯 User-specified thresholds detected: MAX=${EFFECTIVE_MAX_COV}×, TARGET=${EFFECTIVE_TARGET_COV}×"
            echo "[$(date)] 🔬 Using assembly-optimized window-based subsampling (preserves assembly continuity)"

            # ═══ AGGRESSIVE MODE: Per-base subsampling ═══
            set +u
            conda activate ssuitslsu-mapping
            set -u
            PYTHON_EXE="$CONDA_PREFIX/bin/python"
            export BAM
            export subbam
            export EFFECTIVE_TARGET_COV

            "${PYTHON_EXE}" << 'PYCODE'
import pysam, numpy as np, random
import sys
import os
from collections import defaultdict

print("🔬 Assembly-optimized window-based subsampling activated", file=sys.stderr)
print("⚠️  This ensures continuous coverage for optimal assembly", file=sys.stderr)

# Get parameters from environment
IN_BAM = os.environ['BAM']
OUT_BAM = os.environ['subbam']
TARGET = float(os.environ['EFFECTIVE_TARGET_COV'])

# Assembly-optimized parameters for short reads
WINDOW_SIZE = 25          # 25bp windows for fine-grained control with short reads
MIN_COVERAGE = max(5, TARGET * 0.2)   # Minimum 5× or 20% of target
OVERLAP_BONUS = 2.0       # Higher bonus for reads spanning multiple windows

# Set random seed for reproducibility
random.seed(42)

print(f"🎯 Target coverage: {TARGET}×, Min coverage: {MIN_COVERAGE}×, Window size: {WINDOW_SIZE}bp", file=sys.stderr)

reads_to_keep = set()
total_reads = 0

with pysam.AlignmentFile(IN_BAM, "rb") as bam:
    total_contigs = len(bam.references)
    print(f"📊 Processing {total_contigs} reference sequences...", file=sys.stderr)

    for contig_idx, (rname, length) in enumerate(zip(bam.references, bam.lengths)):
        print(f"   🧬 Analyzing contig {contig_idx+1}/{total_contigs}: {rname} ({length:,}bp)", file=sys.stderr)

        # Step 1: Get all reads for this contig
        contig_reads = []
        for read in bam.fetch(rname):
            if not read.is_unmapped and read.reference_start < read.reference_end:
                contig_reads.append(read)

        total_reads += len(contig_reads)
        print(f"      📖 Found {len(contig_reads):,} mapped reads", file=sys.stderr)

        if not contig_reads:
            continue

        # Step 2: Calculate current coverage per window
        num_windows = (length + WINDOW_SIZE - 1) // WINDOW_SIZE
        window_coverage = np.zeros(num_windows)
        window_reads = defaultdict(list)

        for read in contig_reads:
            start_win = read.reference_start // WINDOW_SIZE
            end_win = min((read.reference_end - 1) // WINDOW_SIZE, num_windows - 1)

            # Count read contribution to each window it spans
            for win_idx in range(start_win, end_win + 1):
                window_coverage[win_idx] += 1
                window_reads[win_idx].append(read)

        initial_mean_cov = window_coverage.mean()
        print(f"      📈 Initial mean coverage: {initial_mean_cov:.1f}×", file=sys.stderr)

        # Step 3: Window-based intelligent subsampling
        selected_reads_this_contig = set()

        # Phase 1: Ensure minimum coverage in all windows
        for win_idx in range(num_windows):
            win_reads = window_reads[win_idx]
            current_cov = window_coverage[win_idx]

            if current_cov <= MIN_COVERAGE:
                # Keep all reads in low-coverage windows
                for read in win_reads:
                    selected_reads_this_contig.add(read.query_name)
            else:
                # Subsample high-coverage windows but maintain minimum
                target_reads = max(int(MIN_COVERAGE), int(current_cov * TARGET / initial_mean_cov))
                target_reads = min(target_reads, len(win_reads))

                # Prioritize reads that span multiple windows (better for assembly)
                read_scores = []
                for read in win_reads:
                    span_windows = ((read.reference_end - 1) // WINDOW_SIZE) - (read.reference_start // WINDOW_SIZE) + 1
                    score = span_windows * OVERLAP_BONUS + random.random()  # Add randomness for tie-breaking
                    read_scores.append((score, read))

                # Select top-scoring reads
                read_scores.sort(reverse=True)
                for i in range(target_reads):
                    selected_reads_this_contig.add(read_scores[i][1].query_name)

        # Phase 2: Fill remaining quota with best spanning reads
        current_selection_size = len(selected_reads_this_contig)
        total_target = int(len(contig_reads) * TARGET / initial_mean_cov)

        if current_selection_size < total_target:
            remaining_quota = total_target - current_selection_size

            # Get unselected reads, prioritize long spans
            unselected_reads = [r for r in contig_reads if r.query_name not in selected_reads_this_contig]
            span_scores = []

            for read in unselected_reads:
                span_length = read.reference_end - read.reference_start
                span_windows = ((read.reference_end - 1) // WINDOW_SIZE) - (read.reference_start // WINDOW_SIZE) + 1
                score = span_length * span_windows + random.random()
                span_scores.append((score, read))

            span_scores.sort(reverse=True)
            for i in range(min(remaining_quota, len(span_scores))):
                selected_reads_this_contig.add(span_scores[i][1].query_name)

        reads_to_keep.update(selected_reads_this_contig)

        # Calculate final coverage for this contig
        final_coverage = np.zeros(num_windows)
        for read in contig_reads:
            if read.query_name in selected_reads_this_contig:
                start_win = read.reference_start // WINDOW_SIZE
                end_win = min((read.reference_end - 1) // WINDOW_SIZE, num_windows - 1)
                for win_idx in range(start_win, end_win + 1):
                    final_coverage[win_idx] += 1

        final_mean = final_coverage.mean()
        final_min = final_coverage.min()
        final_max = final_coverage.max()
        kept_pct = (len(selected_reads_this_contig) / len(contig_reads)) * 100

        print(f"      ✅ Selected {len(selected_reads_this_contig):,}/{len(contig_reads):,} reads ({kept_pct:.1f}%)", file=sys.stderr)
        print(f"      📊 Final coverage: mean={final_mean:.1f}×, min={final_min:.0f}×, max={final_max:.0f}×", file=sys.stderr)

# Step 4: Write selected reads to output
print(f"🔄 Writing selected reads to output...", file=sys.stderr)
reads_processed = 0
reads_kept = 0

with pysam.AlignmentFile(IN_BAM, "rb") as bam, \
     pysam.AlignmentFile(OUT_BAM, "wb", template=bam) as out:

    for read in bam.fetch():
        reads_processed += 1

        if reads_processed % 100000 == 0:
            kept_pct = (reads_kept/reads_processed)*100 if reads_processed > 0 else 0
            print(f"   📝 Written {reads_processed:,} reads, kept {kept_pct:.1f}%", file=sys.stderr)

        if read.is_unmapped or read.query_name in reads_to_keep:
            out.write(read)
            reads_kept += 1

reduction_pct = (1.0 - reads_kept/reads_processed) * 100 if reads_processed > 0 else 0

print(f"✅ Assembly-optimized subsampling complete:", file=sys.stderr)
print(f"   📊 Processed: {reads_processed:,} reads", file=sys.stderr)
print(f"   ✂️  Kept: {reads_kept:,} reads ({100-reduction_pct:.1f}%)", file=sys.stderr)
print(f"   🗑️  Removed: {reads_processed-reads_kept:,} reads ({reduction_pct:.1f}%)", file=sys.stderr)
print(f"   🎯 Assembly continuity preserved with minimum {MIN_COVERAGE}× coverage", file=sys.stderr)
PYCODE

            # Index the subsampled BAM file
            echo "[$(date)] 🔧 Indexing subsampled BAM file..."
            samtools index "${subbam}"

            # Verify final coverage achieved
            echo "[$(date)] 🔍 Verifying final coverage..."
            final_cov=$(samtools depth -a "${subbam}" \
                       | awk '{sum+=$3; cnt++} END {print (cnt? sum/cnt : 0)}')
            printf "[%s] ✅ Final coverage: %.1f× (target: %s×)\n" "$(date)" "${final_cov}" "${EFFECTIVE_TARGET_COV}"

            SUBSAMPLING_OCCURRED=true
            echo "[$(date)] ✅ Aggressive per-base subsampling completed: ${mean_cov}× → ${EFFECTIVE_TARGET_COV}×"

          else
            echo "[$(date)] 🛡️  Conservative mode: gentle subsampling for assembly safety"
            echo "[$(date)] 📉 Reducing ${mean_cov}× → ${EFFECTIVE_TARGET_COV}× (preserving assembly quality)"

            # ═══ CONSERVATIVE MODE: Simple but safe ═══
            FRAC=$(awk -v t="$EFFECTIVE_TARGET_COV" -v c="$mean_cov" \
              'BEGIN { frac = t/c; print (frac > 1.0 ? 1.0 : frac) }')

            echo "[$(date)] 📊 Gentle subsampling fraction: $FRAC"

            # Fast uniform subsampling with proper random seed
            samtools view -h -@ "$THREADS" -s "1234${FRAC}" "$BAM" \
              | samtools sort -@ "$THREADS" -o "$subbam" -T "$SAMPLE_DIR/subsample_tmp"

            SUBSAMPLING_OCCURRED=true
            echo "[$(date)] ✅ Conservative subsampling completed"
          fi

        else
          echo "[$(date)] 🎯 Coverage acceptable (${mean_cov}× ≤ ${EFFECTIVE_MAX_COV}×), no subsampling needed"
          cp "${BAM}" "${subbam}"
          SUBSAMPLING_OCCURRED=false
        fi

      else
        echo "[$(date)] Depth-capped BAM exists; skipping coverage analysis."
        # Check if this was a subsampled file (rough heuristic)
        original_size=$(stat -f%z "$BAM" 2>/dev/null || stat -c%s "$BAM" 2>/dev/null || echo "0")
        subsampled_size=$(stat -f%z "$subbam" 2>/dev/null || stat -c%s "$subbam" 2>/dev/null || echo "0")

        if [[ "$subsampled_size" -lt "$original_size" ]]; then
          SUBSAMPLING_OCCURRED=true
          echo "[$(date)] Existing file appears to be subsampled (smaller than original)"
        else
          SUBSAMPLING_OCCURRED=false
        fi
      fi

      # ─── Conditional fixmate based on what operations occurred ────────────────────
      NEEDS_FIXMATE=false

      # Check if subsampling occurred
      if [[ "$SUBSAMPLING_OCCURRED" == "true" ]]; then
        echo "[$(date)] 🔧 Fixmate needed: subsampling occurred"
        NEEDS_FIXMATE=true
      fi

      # Apply fixmate only if needed
      if [[ "$NEEDS_FIXMATE" == "true" ]]; then
        echo "[$(date)] 🔧 Running fixmate to correct mate-pair flags after subsampling..."

        FIXMATE_NS="${SAMPLE_DIR}/${sample}.depthcapped.namesorted.bam"
        FIXMATE_NS_OUT="${SAMPLE_DIR}/${sample}.depthcapped.fixmate.namesorted.bam"
        FIXMATE_CS="${SAMPLE_DIR}/${sample}.depthcapped.fixmate.coordsorted.bam"

        if [[ ! -f "$FIXMATE_CS" ]]; then
          # 1) Name-sort
          echo "[$(date)] Name-sorting for fixmate..."
          samtools sort -n -@ "$THREADS" -o "$FIXMATE_NS" "$subbam"

          # 2) fixmate
          echo "[$(date)] Applying fixmate corrections..."
          samtools fixmate -m -@ "$THREADS" "$FIXMATE_NS" "$FIXMATE_NS_OUT"

          # 3) coord-sort + filter unmapped + index
          echo "[$(date)] Coordinate-sorting and indexing..."
          samtools view -bh -@ "$THREADS" -F 4 "$FIXMATE_NS_OUT" \
            | samtools sort -@ "$THREADS" -o "$FIXMATE_CS"
          samtools index "$FIXMATE_CS"

          # Clean up intermediate files
          rm -f "$FIXMATE_NS" "$FIXMATE_NS_OUT"

          BAM="$FIXMATE_CS"
          echo "[$(date)] ✅ Fixmate completed, using corrected BAM"
        else
          echo "[$(date)] Fixmate files already exist, using existing corrected BAM"
          BAM="$FIXMATE_CS"
        fi
      else
        echo "[$(date)] 🎯 No fixmate needed, using BAM as-is"
        BAM="$subbam"
      fi



      # ─── Extract mapped reads with conditional filtering ────────────────────────
      P1="$SAMPLE_DIR/${sample}.mapped_1.fastq.gz"
      P2="$SAMPLE_DIR/${sample}.mapped_2.fastq.gz"
      U0="$SAMPLE_DIR/${sample}.mapped_0.fastq.gz"

      # re-extract if missing, empty, or older than the BAM
      if [[ ! -s "$P1" || ! -s "$P2" || "$BAM" -nt "$P1" || "$BAM" -nt "$P2" ]]; then
        echo "[$(date)] Checking read loss from strict filtering"

        # Total mapped reads (removing the unmapped)
        total=$(samtools view -@ "$THREADS" -c -F 4 "$BAM")

        # Reads that would pass the strict filter
        filtered=$(samtools view -@ "$THREADS" -c -F 2316 -f 2 -q "$MAPQ" "$BAM")  # 2316 = 4(unmapped) + 8(unpaired) + 256(secondary aln) + 2048(supplementary aln)



        if (( total > 0 )); then
          pct_kept=$(awk -v t="$total" -v f="$filtered" 'BEGIN { printf("%.2f", (f/t)*100) }')
          echo "[$(date)] Strict filtering keeps $pct_kept% of mapped reads"
        else
          echo "[$(date)] No mapped reads (total=0); skipping strict-filter percentage calculation"
          pct_kept=0
        fi

        # if we lose >20% of the reads we go back to relaxed mapping
        FILTER=""
        ORPHAN_OPT=""
        if (( $(echo "$pct_kept > 80.0" | bc -l) )); then
          echo "[$(date)] Applying strict filtering (properly paired, primary, high MAPQ)"
          FILTER="-F 2316 -f 2 -q $MAPQ"
          ORPHAN_OPT="-s $U0"
        else
          echo "[$(date)] Too much data loss, relaxing filtering to include unpaired and non-primary reads"
          FILTER="-F 4 -q $MAPQ"
          ORPHAN_OPT="-0 $U0"
        fi

        rm -f "$P1" "$P2" "$U0"

        # 1) write out a temporary BAM with the chosen filter
        samtools view -bh -@ "$THREADS" -b $FILTER "$BAM" -o "${SAMPLE_DIR}/${sample}.filtered.tmp.bam"

        # 2) name-sort that BAM (so fixmate can see both mates together)
        samtools sort -n -@ "$THREADS" \
          -o "${SAMPLE_DIR}/${sample}.filtered.namesorted.bam" \
          "${SAMPLE_DIR}/${sample}.filtered.tmp.bam"

        # 3) run fixmate to correct/mask the pairing flags on any orphans
        samtools fixmate -m -@ "$THREADS" \
          "${SAMPLE_DIR}/${sample}.filtered.namesorted.bam" \
          "${SAMPLE_DIR}/${sample}filtered.fixmate.ns.bam"
        samtools sort -@ "$THREADS" \
          -o "${SAMPLE_DIR}/${sample}.filtered.fixmate.bam" \
          "${SAMPLE_DIR}/${sample}filtered.fixmate.ns.bam"


        FINAL_BAM="${SAMPLE_DIR}/${sample}.filtered.fixmate.bam"

        if [[ "$FILTER" == *"-f 2"* ]]; then
          # ─── Strict mode: paired+proper only, orphans as singleton ───────────────
          samtools collate -uO -@ "$THREADS" "$FINAL_BAM" |
          samtools fastq -@ "$THREADS" \
            -1 "$P1" \
            -2 "$P2" \
            $ORPHAN_OPT \
            -0 /dev/null  # swallow any stray reads

        else
          # ─── Relaxed mode: split primary pairs vs everything else ───────────────

          # A) Primary, properly-paired, mapped → P1/P2
          samtools view -bh -@ "$THREADS" \
            -f 1 -f 2 \
            -F 4 \
            -F 256 -F 2048 \
            "$FINAL_BAM" \
          | samtools collate -uO -@ "$THREADS" /dev/stdin \
          | samtools fastq -@ "$THREADS" \
              -1 "$P1" \
              -2 "$P2" \
              -0 /dev/null

          # B) All other mapped reads → singleton file U0
          samtools view -h -@ "$THREADS" -F 4 "$FINAL_BAM" \
          | samtools collate -uO -@ "$THREADS" /dev/stdin \
          | samtools fastq -@ "$THREADS" \
               $ORPHAN_OPT
        fi

        # sanity check: ensure mate-pair counts match
        p1=$(zgrep -c '^@' "$P1" 2>/dev/null || echo 0)
        p2=$(zgrep -c '^@' "$P2" 2>/dev/null || echo 0)
        if [[ "$p1" -ne "$p2" ]]; then
          echo "[ERROR] mate-pair counts mismatch: R1=$p1 vs R2=$p2" >&2
          exit 1
        fi
        echo "[$(date)] Extracted $p1 paired reads each in P1/P2"

      else
        echo "[$(date)] Mapped FASTQs up-to-date; skipping extraction"
      fi

      # ─── Bail out if filtering produced no reads ────────────────────────────────
      p1=$(zgrep -c '^@' "$P1" 2>/dev/null || echo 0)
      p1=$(echo "$p1" | tr -d '\n\r\t ' | grep -o '^[0-9]*$' || echo 0)
      p2=$(zgrep -c '^@' "$P2" 2>/dev/null || echo 0)
      p2=$(echo "$p2" | tr -d '\n\r\t ' | grep -o '^[0-9]*$' || echo 0)

      if [[ $p1 -eq 0 ]] && [[ $p2 -eq 0 ]]; then
        echo "[`date`] No reads after filtering for $sample (R1=$p1, R2=$p2); skipping this sample"
        continue
      fi

      t_map_end=$(date +%s)



      # ─── Assembly ────────────────────────────────────────────────────────
      t_asm_start=$(date +%s)
      ASM_DIR="$SAMPLE_DIR/assembly"
      mkdir -p "$ASM_DIR"

      if [[ "$ASSEMBLER" == "spades" ]]; then
        SP_DIR="$ASM_DIR/spades"
        CONTIG="$SP_DIR/contigs.fasta"

        if [[ ! -s "$CONTIG" ]]; then
          set +u
          conda activate ssuitslsu-spades
          set -u
          echo "[`date`] Assembling with SPAdes"

          # 1) clear any old run
          rm -rf "$SP_DIR"
          mkdir -p "$SP_DIR"

          # 2) count true singletons and decide on the -s flag

          if [[ -f "$U0" ]]; then
            # capture only zgrep’s single “0” (or a positive count)
            n0=$(zgrep -c '^@' "$U0" 2>/dev/null || true)
            # if zgrep failed entirely, or printed nothing, default to 0
            n0=${n0:-0}
          else
            n0=0
          fi

          if [ "$n0" -gt 0 ]; then
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
              -t "$THREADS" -m "$MEM"
              -o "$SP_DIR"
          )

          set +e
          "${cmd[@]}" 2>&1 | tee "$SP_DIR/spades.log"
          set -e

          # Check for contigs; if none, retry with relaxed parameters
          CONTIG="$SP_DIR/contigs.fasta"
          if [[ ! -s "$CONTIG" ]]; then
            printf "[%s] No contigs for %s from SPAdes → retrying with relaxed parameters\n" \
              "$(date)" "$sample"

            # Clean out the old SPAdes directory
            rm -rf "$SP_DIR"
            mkdir -p "$SP_DIR"

            # Retry SPAdes with lower coverage cutoff and smaller min-contig
            printf "[DEBUG] retry command: spades.py -t %s -m %s --cov-cutoff auto --min-contig 100 -k 21,33,55,77 -o %s -1 %s -2 %s\n" \
            "$THREADS" "$MEM" "$SP_DIR" "$P1" "$P2" >&2
            rm -rf "$SP_DIR"
            mkdir -p "$SP_DIR"
            set +e
            spades.py \
              -t "$THREADS" \
              -m "$MEM" \
              --cov-cutoff auto \
              -k 21,33,55,77 \
              -o "$SP_DIR" \
              -1 "$P1" \
              -2 "$P2" \
            2>&1 | tee "$SP_DIR/spades.retry.log"
            set -e

          fi

          # If still no contigs, skip to the next FASTQ pair
          if [[ ! -s "$CONTIG" ]]; then
            printf "[%s] Still no contigs for %s after retry; skipping this sample\n" \
              "$(date)" "$sample"
            continue
          fi

        else
          echo "[`date`] Skipping SPAdes assembly (found $CONTIG)"
        fi

      elif [[ "$ASSEMBLER" == "megahit" ]]; then
        MH_DIR="$ASM_DIR/megahit"
        CONTIG="$MH_DIR/final.contigs.fa"

        if [[ ! -s "$CONTIG" ]]; then
          set +u
          conda activate ssuitslsu-megahit
          set -u
          echo "[`date`] Assembling with MEGAHIT"

          # — delete U0 if it contains no FASTQ records
          if [[ -f "$U0" ]]; then
              if ! zgrep -q '^@' "$U0"; then
                  echo "→ Removing empty singleton"
                  rm -f "$U0"
              else
                  count0=$(zgrep -c '^@' "$U0" 2>/dev/null || echo 0)
                  echo "→ Including $count0 singletons"
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
            --memory "$MEM" \
            --verbose \
            --keep-tmp-files
          )
          # include singletons if present
          [[ -f "$U0" ]] && cmd+=( -r "$U0" )

          # run and capture the log
          set +e
          "${cmd[@]}" 2>&1 | tee "$MH_DIR/megahit.log"
          set -e

          # Check for contigs; if none, skip to the next FASTQ pair
          CONTIGS="$MH_DIR/final.contigs.fa"
          if [[ ! -s "$CONTIGS" ]]; then
            printf "[%s] No MEGAHIT contigs for %s; skipping this sample\n" \
              "$(date)" "$sample"
            continue
          fi

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
        set +u
        conda activate ssuitslsu-itsx
        set -u

        # 1) extract contigs ≥1 kb
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
          set +u
          conda activate ssuitslsu-mapping
          set -u
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
        set +u
        conda activate ssuitslsu-itsx
        set -u

        # split contigs by length
        SMALL_FASTA="$ITSX_DIR/${sample}.small.fasta"
        LARGE_FASTA="$ITSX_DIR/${sample}.large.fasta"
        : > "$SMALL_FASTA"
        : > "$LARGE_FASTA"

        awk -v L=100000 \
          -v out_small="$SMALL_FASTA" \
          -v out_large="$LARGE_FASTA" \
          'BEGIN { RS=">"; ORS=""; FS="\n" }
          NR>1 {
            header = $1
            seq = ""
            for (i=2; i<=NF; i++) seq = seq $i
            out = (length(seq)>L ? out_large : out_small)
            printf(">%s\n%s\n", header, seq) > out
          }
        ' "$CONTIG"

        # 1) run ITSx on the small contigs (default hmmsearch)
        SMALL_PREFIX="${ITSX_DIR}/${sample}.small"
        SMALL_DONE="${SMALL_PREFIX}.done"
        RECHECK_DONE="${SMALL_PREFIX}.recheck.done"

        # (ri)run principale solo se non fatto
        if [[ -s "$SMALL_FASTA" && ! -f "$SMALL_DONE" ]]; then
          echo "[`date`] Running ITSx (hmmsearch) on small contigs for $sample"
          ITSx -i "$SMALL_FASTA" \
               -o "$SMALL_PREFIX" \
               --preserve T \
               --save_regions all \
               --anchor HMM \
               --region all \
               --cpu "$THREADS" \
            && touch "$SMALL_DONE"
        elif [[ -f "$SMALL_DONE" ]]; then
          echo "[`date`] Skipping small-contig ITSx (already done)"
        else
          echo "[`date`] No small contigs to process"
        fi

        # fallback: re-scan just no_detection
        NODET="${SMALL_PREFIX}_no_detections.fasta"
        if [[ -s "$NODET" && ! -f "$RECHECK_DONE" ]]; then
          echo "[`date`] Re-running no-detections with nhmmer for $sample"
          ITSx -i "$NODET" \
               -o "${ITSX_DIR}/${sample}.small.recheck" \
               --preserve T \
               --save_regions all \
               --anchor HMM \
               --region all \
               --nhmmer T \
               --cpu "$THREADS" \
            && touch "$RECHECK_DONE"

          # append eventuali hit nuovi nei file small.*.fasta
          for region in full ITS1 ITS2 SSU LSU 5_8S; do
            recheck_fa="${ITSX_DIR}/${sample}.small.recheck.${region}.fasta"
            main_fa="${SMALL_PREFIX}.${region}.fasta"
            if [[ -s "$recheck_fa" ]]; then
              echo "[`date`] Adding recheck ${region} → ${main_fa}"
              cat "$recheck_fa" >> "$main_fa"
            fi
          done
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
      CONTIGS_FASTA="$PHYLO_DIR/${sample}_contigs.fasta"
      ITSX_FASTA="$PHYLO_DIR/${sample}_ITSx.fasta"

      if [[ -s "$SAMPLE_SEQ_FASTA" ]]; then
        echo "[$(date)] Sample FASTA already exists → skipping ITSx+contigs merge: $SAMPLE_SEQ_FASTA"
      else

        # 2) Truncate/create each file
        : > "$ITSX_FASTA"
        : > "$CONTIGS_FASTA"
        : > "$SAMPLE_SEQ_FASTA"

        # 3a) Collect ITSx outputs and count by region type
        ITS_COUNT=0
        SSU_COUNT=0
        LSU_COUNT=0
        OTHER_COUNT=0

        for size in small large_nhmmer; do
          for region in full ITS1 ITS2 SSU LSU 5_8S; do
            src="$ITSX_DIR/${sample}.${size}.${region}.fasta"
            if [[ -s "$src" ]]; then
              seq_count=$(grep -c "^>" "$src" 2>/dev/null || echo "0")
              printf "[%s] Adding ITSx %-4s (%s): %d sequences → %s\n" "$(date)" "$region" "$size" "$seq_count" "$src"

              # Count by region type
              if [[ "$region" =~ ^(full|ITS1|ITS2|5_8S)$ ]]; then
                ITS_COUNT=$((ITS_COUNT + seq_count))
              elif [[ "$region" == "SSU" ]]; then
                SSU_COUNT=$((SSU_COUNT + seq_count))
              elif [[ "$region" == "LSU" ]]; then
                LSU_COUNT=$((LSU_COUNT + seq_count))
              else
                OTHER_COUNT=$((OTHER_COUNT + seq_count))
              fi

              awk -v smp="$sample" -v sz="$size" -v rgn="$region" '
                /^>/ { sub(/^>/, ">" smp "_" sz "_" rgn "_"); print; next }
                { print }
              ' "$src" >> "$ITSX_FASTA"
            fi
          done
        done

        # Report what we found
        echo "[$(date)] ITSx region summary: $ITS_COUNT ITS, $SSU_COUNT SSU, $LSU_COUNT LSU, $OTHER_COUNT other"
        echo "[$(date)] DEBUG: Starting ITS concatenation section" >&2

        # ─── 3a.5) ITS concatenation fallback (only ITS1+5.8S+ITS2) ────────────────
        FULL_ITS_COUNT=0
        SMALL_REGIONS_FOUND=0

        # Count full ITS sequences
        for full_file in "$ITSX_DIR"/${sample}.*.full.fasta; do
          if [[ -s "$full_file" ]]; then
            file_count=$(grep -c "^>" "$full_file" 2>/dev/null || echo "0")
            FULL_ITS_COUNT=$((FULL_ITS_COUNT + file_count))
          fi
        done
        echo "[$(date)] DEBUG: Found $FULL_ITS_COUNT full ITS sequences" >&2

        CONCATENATED_ITS_FASTA="$PHYLO_DIR/${sample}_concatenated_ITS.fasta"

        # Count small ITS regions for potential concatenation
        for region in ITS1 5_8S ITS2; do
          if [[ -s "$ITSX_DIR/${sample}.small.${region}.fasta" ]]; then
            SMALL_REGIONS_FOUND=$((SMALL_REGIONS_FOUND + 1))
          fi
        done
        echo "[$(date)] DEBUG: Found $SMALL_REGIONS_FOUND small regions" >&2

        # Apply concatenation logic
        if [[ "$FULL_ITS_COUNT" -eq 0 && "$SMALL_REGIONS_FOUND" -ge 2 ]]; then
          echo "[$(date)] No full ITS sequences; attempting concatenation fallback ($SMALL_REGIONS_FOUND small regions found)"

          set +u
          conda activate ssuitslsu-utils
          set -u

          set +e
          python3 "$(dirname "$0")/its_concatenation_fallback.py" \
            "$ITSX_DIR" "$sample" "$CONCATENATED_ITS_FASTA"
          PY_EXIT=$?
          set -e

          if [[ $PY_EXIT -eq 0 && -s "$CONCATENATED_ITS_FASTA" ]]; then
            echo "[$(date)] Successfully concatenated ITS sequences, adding to collection"
            cat "$CONCATENATED_ITS_FASTA" >> "$ITSX_FASTA"
          else
            echo "[$(date)] ITS concatenation failed or produced no sequences"
          fi
        elif [[ "$FULL_ITS_COUNT" -eq 0 && "$SMALL_REGIONS_FOUND" -lt 2 ]]; then
          if [[ "$SSU_COUNT" -gt 0 || "$LSU_COUNT" -gt 0 ]]; then
            echo "[$(date)] No ITS sequences, but found $SSU_COUNT SSU and $LSU_COUNT LSU sequences - proceeding with alignment"
          else
            echo "[$(date)] No ITS, SSU, or LSU sequences found - alignment will use contigs only"
          fi
        else
          echo "[$(date)] Found $FULL_ITS_COUNT full ITS sequences - skipping concatenation"
        fi
        echo "[$(date)] DEBUG: Finished ITS concatenation logic" >&2

        echo "[$(date)] Proceeding with sequence selection and alignment"
        echo "[$(date)] DEBUG: Starting contig selection" >&2


        # 3b) Extract contigs ≥300 bp & cov>40 (SPAdes cov_ or MEGAHIT multi), plus always include ITSx-hit contigs

        # ITSx guard: only do ITSx-based selection if ITSx found hits
        LARGE_POS="$ITSX_DIR/${sample}.large_nhmmer.positions.txt"
        SMALL_POS="$ITSX_DIR/${sample}.small.positions.txt"
        if [[ -s "$LARGE_POS" ]]; then
          POS_FILE="$LARGE_POS"
          echo "[$(date)] DEBUG: Using large positions file" >&2
        elif [[ -s "$SMALL_POS" ]]; then
          POS_FILE="$SMALL_POS"
          echo "[$(date)] DEBUG: Using small positions file" >&2
        else
          POS_FILE=""
          echo "[$(date)] DEBUG: No position files found" >&2
        fi
        combined_pos=$(mktemp)
        if [[ -z "$POS_FILE" ]]; then
          printf "[%s] No ITSx hits for %s; including *all* contigs\n" "$(date)" "$sample"
          # fallback: just copy every contig through
          cp "$CONTIG" "$CONTIGS_FASTA"
          echo "[$(date)] DEBUG: Copied all contigs" >&2
        else
          echo "[$(date)] DEBUG: Processing ITSx hits and filtering contigs" >&2
          printf "[%s] Selecting contigs for alignment (ITSx hits + length/coverage filter) from %s → %s\n" \
            "$(date)" "$CONTIG" "$CONTIGS_FASTA"

          for pf in \
            "$ITSX_DIR/${sample}.large_nhmmer.positions.txt" \
            "$ITSX_DIR/${sample}.small.positions.txt"; do
            if [[ -f "$pf" ]]; then
              echo "[DEBUG]  Adding IDs from $pf" >&2
              cat "$pf" >> "$combined_pos"
            else
              echo "[DEBUG]  Missing positions file: $pf" >&2
            fi
          done
          head -n5 "$combined_pos" >&2

          awk -v posfile="$combined_pos" '
            BEGIN {
              # 1) read ITSx IDs line-by-line
              RS = "\n"; OFS = ""
              if (posfile != "") {
                while ((getline fline < posfile) > 0) {
                  # trim after first space or tab
                  sub(/[ \t].*$/, "", fline)
                  itsx[fline] = 1
                }
                close(posfile)
              }
              # debug how many unique IDs we loaded
              count = 0
              for (k in itsx) count++
              print "[DEBUG] AWK: loaded " count " ITSx IDs" > "/dev/stderr"
              # now switch to FASTA record mode
              RS = ">"; ORS = ""
            }
            NR > 1 {
              # split header + sequence
              split($0, P, "\n")
              hdr = P[1]
              seq = ""
              for (i = 2; i <= length(P); i++) seq = seq P[i]

              # extract contig ID (first token)
              split(hdr, H, /[ \t]+/)
              id = H[1]

              # keep if ITSx hit
              keep = (id in itsx)

              # otherwise apply length>=300 & (cov>40 or multi>40)
              if (!keep && length(seq) >= 300) {
                cov = -1; multi = -1
                if (hdr ~ /cov_[0-9]/) {
                  split(hdr, A, "cov_"); split(A[2], B, /[^0-9.]/)
                  cov = B[1]
                }
                if (hdr ~ /multi=[0-9]/) {
                  split(hdr, C, "multi="); split(C[2], D, /[^0-9.]/)
                  multi = D[1]
                }
                if (cov > 40 || multi > 40) keep = 1
              }

              if (keep) print ">" $0
            }
          ' "$CONTIG" > "$CONTIGS_FASTA"
        fi

        grep '^>' "$CONTIGS_FASTA" | head -n5 >&2

        rm -f "$combined_pos"

        # 4) Merge ITSx + contigs into the final sample FASTA
        cat "$CONTIGS_FASTA" "$ITSX_FASTA" > "$SAMPLE_SEQ_FASTA"

        # Report sequence counts for debugging
        ITSX_COUNT=$(grep -c "^>" "$ITSX_FASTA" 2>/dev/null || echo "0")
        CONTIGS_COUNT=$(grep -c "^>" "$CONTIGS_FASTA" 2>/dev/null || echo "0")
        TOTAL_COUNT=$(grep -c "^>" "$SAMPLE_SEQ_FASTA" 2>/dev/null || echo "0")

        echo "[$(date)] Sequence summary for $sample:"
        echo "[$(date)]   - ITSx sequences: $ITSX_COUNT ($ITS_COUNT ITS, $SSU_COUNT SSU, $LSU_COUNT LSU)"
        echo "[$(date)]   - Contigs: $CONTIGS_COUNT"
        echo "[$(date)]   - Total sequences: $TOTAL_COUNT"

        # Check if we have any sequences at all
        if [[ "$TOTAL_COUNT" -eq 0 ]]; then
          echo "[$(date)] ⚠️  No sequences found for $sample - creating minimal placeholder sequence"
          # Create a minimal sequence so alignment can still proceed
          echo ">${sample}_no_sequences_found" > "$SAMPLE_SEQ_FASTA"
          echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> "$SAMPLE_SEQ_FASTA"
        else
          echo "[$(date)] ✅ Proceeding with alignment using available sequences"
        fi
      fi

      # 5a) Build reference at the same taxonomic level used for mapping
      safe_phylo_tag="${PHYLO_REF_TAG//;/__}"
      REF_PHYLO_FASTA="$PHYLO_DIR/${safe_phylo_tag}_ref.fasta"

      if [[ ! -s "$REF_PHYLO_FASTA" ]]; then
        printf "[%s] Building phylogeny reference for %s → %s\n" "$(date)" "${PHYLO_REF_TAG}" "$REF_PHYLO_FASTA"
        # extract all entries matching this taxonomic level
        if [[ "$PHYLO_REF_TAG" == *";"* ]]; then
          # Handle combined tags like "g__Genus;s__species"
          gtag="${PHYLO_REF_TAG%%;*}"; stag="${PHYLO_REF_TAG##*;}"
          awk -v g="$gtag" -v s="$stag" 'BEGIN{RS=">";ORS=""} $0~g&&$0~s{print ">"$0}' "$REF_FASTA" > "$REF_PHYLO_FASTA"
        else
          # Handle single tags like "g__Genus", "f__Family", etc.
          awk -v tag="$PHYLO_REF_TAG" '
            BEGIN { RS=">"; ORS="" }
            $0 ~ ("(^|;)" tag "(;|$)") { print ">" $0 }
          ' "$REF_FASTA" > "$REF_PHYLO_FASTA"
        fi

        # count and prune if too many
        num=$(grep -c '^>' "$REF_PHYLO_FASTA")
        if (( num > 100 )); then
          echo "[$(date)] $num sequences in phylogeny reference; pruning to max 100…"

          # a) drop all "s__unclassified" entries
          awk '
            BEGIN { RS=">"; ORS="" }
            NR>1 && $0 !~ /s__unclassified/ { print ">" $0 }
          ' "$REF_PHYLO_FASTA" > "${REF_PHYLO_FASTA%.fasta}_pruned.fasta"
          mv "${REF_PHYLO_FASTA%.fasta}_pruned.fasta" "$REF_PHYLO_FASTA"

          # 5b) if still >100, keep only one per species (based on the "s__" field)
          num=$(grep -c '^>' "$REF_PHYLO_FASTA")
          if (( num > 80 )); then
            awk '
              BEGIN {
                RS=">"; ORS=""
                kept=0
              }
              NR>1 && kept<80 {
                # split header by ";" to find the species tag
                n = split($0, fields, ";")
                species = ""
                for (i=1; i<=n; i++) {
                  if (fields[i] ~ /^s__/) {
                    species = substr(fields[i], 4)
                    break
                  }
                }
                if (species == "") {
                  species = "__no_species__"
                }
                if (!(species in seen)) {
                  seen[species] = 1
                  print ">" $0
                  kept++
                }
              }
            ' "$REF_PHYLO_FASTA" > "${REF_PHYLO_FASTA%.fasta}_dedup.fasta"
            mv "${REF_PHYLO_FASTA%.fasta}_dedup.fasta" "$REF_PHYLO_FASTA"
          fi

          echo "[$(date)] Reference count after pruning: $(grep -c '^>' "$REF_PHYLO_FASTA")"
        fi

      else
        printf "[%s] Skipping phylogeny reference for %s (already exists): %s\n" \
          "$(date)" "$PHYLO_REF_TAG" "$REF_PHYLO_FASTA"
      fi

      # 5c) Combine sample + reference into one FASTA
      # Only rebuild the .fasta if it doesn't exist or either input is newer
      PHYLO_FASTA="$PHYLO_DIR/${safe_phylo_tag}_phylo.fasta"
      if [[ ! -s "$PHYLO_FASTA" ]] \
         || [[ "$REF_PHYLO_FASTA"   -nt "$PHYLO_FASTA" ]] \
         || [[ "$SAMPLE_SEQ_FASTA" -nt "$PHYLO_FASTA" ]]; then

        echo "[$(date)] (Re)building combined FASTA → $PHYLO_FASTA"
        {
          cat "$REF_PHYLO_FASTA"
          printf "\n"
          cat "$SAMPLE_SEQ_FASTA"
        } > "$PHYLO_FASTA"
      fi


      REF_COUNT=$(grep -c '^>' "$REF_PHYLO_FASTA" 2>/dev/null || echo "0")
      SAMPLE_COUNT=$(grep -c '^>' "$SAMPLE_SEQ_FASTA" 2>/dev/null || echo "0")
      TOTAL_PHYLO_COUNT=$(grep -c '^>' "$PHYLO_FASTA" 2>/dev/null || echo "0")

      printf "[%s] Combined phylo FASTA: %s (%d refs + %d sample seqs = %d total)\n\n" \
        "$(date)" "$PHYLO_FASTA" "$REF_COUNT" "$SAMPLE_COUNT" "$TOTAL_PHYLO_COUNT"

      # Ensure we have enough sequences for meaningful alignment
      if [[ "$TOTAL_PHYLO_COUNT" -lt 2 ]]; then
        echo "[$(date)] ⚠️  Only $TOTAL_PHYLO_COUNT sequences total - alignment may not be meaningful"
        echo "[$(date)] Proceeding anyway to maintain pipeline consistency"
      fi

      # 6a) MAFFT alignment: run if missing OR FASTA is newer than the .aln
      PHYLO_ALN="$PHYLO_DIR/${safe_phylo_tag}.aln"
      if [[ ! -s "$PHYLO_ALN" ]] || [[ "$PHYLO_FASTA" -nt "$PHYLO_ALN" ]]; then
        # either no alignment yet, or the FASTA changed—rebuild it
        t_align_start=$(date +%s)
        set +u; conda activate ssuitslsu-mafft; set -u
        printf "[%s] Running MAFFT → %s\n" "$(date)" "$PHYLO_ALN"
        mafft --thread "$THREADS" --auto --adjustdirection "$PHYLO_FASTA" > "$PHYLO_ALN"
        t_align_end=$(date +%s)
        printf "\n"
      else
        printf "[%s] Skipping MAFFT; alignment is up to date: %s\n\n" \
          "$(date)" "$PHYLO_ALN"
        t_align_end=$t_align_start
      fi

      trimmed_aln="${PHYLO_ALN%.aln}.trimmed.aln"
      if [[ ! -s "$trimmed_aln" || "$PHYLO_ALN" -nt "$trimmed_aln" ]]; then
        echo "[$(date)] Trimming alignment → $trimmed_aln"
        # 6b) Trimming the alignment
        echo "Trimming the alignment..."
        source "$(conda info --base)/etc/profile.d/conda.sh"
        set +u
        conda activate ssuitslsu-chimera
        set -u
        PYTHON_EXE="$CONDA_PREFIX/bin/python"
        trimmed_aln="${PHYLO_ALN%.aln}.trimmed.aln"
        export trimmed_aln="${PHYLO_ALN%.aln}.trimmed.aln"
        export PHYLO_RAW="$PHYLO_ALN"
        export PHYLO_ALN="$trimmed_aln"
        export PHYLO_DIR="$PHYLO_DIR"
        export ITSX_DIR="$ITSX_DIR"
        export sample="$sample"

        read START TAIL < <( "$PYTHON_EXE" scripts/trim_alignment.py )

        export TRIM_HEAD=$START
        export TRIM_TAIL=$TAIL

        # sanity check & override
        if [[ ! -s "$trimmed_aln" ]]; then
          echo "[ERROR] trimmed file not created: $trimmed_aln"
          exit 1
        fi

        TRIMMED_ALN="$trimmed_aln"
      else
        echo "[$(date)] Skipping trim; up to date: $trimmed_aln"
      fi
      export PHYLO_ALN="$trimmed_aln"
      export PHYLO_RAW="$PHYLO_ALN"
      export TRIM_HEAD=0
      export TRIM_TAIL=0

      # 7) IQ-TREE ML inference (conditionally skipped)
      if [[ "$SKIP_PHYLO" == "true" ]]; then
        echo "[$(date)] skip_phylogeny=true: skipping IQ-TREE phylogenetic inference"
        TREEFILE=""
        export CHIMERA_USE_ALIGNMENT=true
      else
        PHYLO_PREFIX="$PHYLO_DIR/${safe_phylo_tag}"
        TREEFILE="${PHYLO_PREFIX}.treefile"
        export CHIMERA_USE_ALIGNMENT=false

        if [[ -f "${PHYLO_PREFIX}.ckp.gz" ]]; then
          printf "[%s] Skipping IQ-TREE (tree exists): %s\n\n" "$(date)" "$TREEFILE"
        else
          t_phylo_start=$(date +%s)
          set +u; conda activate ssuitslsu-iqtree; set -u
          printf "[%s] Running IQ-TREE → prefix %s\n" "$(date)" "$PHYLO_PREFIX"
          iqtree \
            -s "$TRIMMED_ALN" \
            -m MFP \
            -nt AUTO \
            -ntmax "${THREADS}" \
            -mem "${MEM}G" \
            -bb 1000 \
            -pre "$PHYLO_PREFIX"
          t_phylo_end=$(date +%s)
          printf "\n"
        fi
      fi


      # ─── Chimera Detection via Sliding p-distance ────────────────────────────────
      echo "[`date`] Starting chimera detection for $sample..."
      t_chimera_start=$(date +%s)
      set +u
      conda activate ssuitslsu-chimera
      set -u
      PYTHON_EXE="$CONDA_PREFIX/bin/python"

      export CHIMERA_DIR="${SAMPLE_DIR}/chimera"
      export CHIMERA_USE_ALIGNMENT="$CHIMERA_USE_ALIGNMENT"
      export PHYLO_RAW="$PHYLO_RAW"
      export TRIM_HEAD="$TRIM_HEAD"
      export PHYLO_ALN="$PHYLO_ALN"
      export TREEFILE="$TREEFILE"
      export CONTIGS_FASTA="$CONTIGS_FASTA"
      export GENUS_REF_FASTA="$REF_PHYLO_FASTA"
      export CHIMERA_REPORT="${CHIMERA_DIR}/${sample}_chimera_report.tsv"
      export CHIMERA_PLOT="${CHIMERA_DIR}/${sample}_chimera_plot.png"
      export sample="$sample"
      export ITSX_DIR="$ITSX_DIR"

      mkdir -p "$CHIMERA_DIR"

      "${PYTHON_EXE}" scripts/chimera_pipeline.py


      t_chimera_end=$(date +%s)
      echo "[`date`] Chimera detection completed for $sample in $((t_chimera_end - t_chimera_start)) seconds."


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
      chimera_dur=$((t_chimera_end  - t_chimera_start))

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
      printf "    Chimera detection:  %02dm%02ds\n\n" $((chimera_dur/60)) $((chimera_dur%60))

      # Restore original stdout/stderr
      exec 1>&3 2>&4
      exec 3>&- 4>&-

      break
    fi
  done
done

echo "[$(date)] ALL SAMPLES COMPLETE."
