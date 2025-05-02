#!/usr/bin/env bash
set -euo pipefail

# ──────────────────────────────────────────────────────────────────────────────
# Usage/help message
usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  -n, --no-trim           Skip Trimmomatic trimming (use existing *_trimmed FASTQs)
  -m, --max-memory MB     Override memory (GB) for SPAdes
  -t, --threads N         Override number of threads for all steps
  -h, --help              Show this help message and exit
EOF
}
# ──────────────────────────────────────────────────────────────────────────────

# Default CLI overrides
SKIP_TRIM_CLI=false
MEM_CLI=""
THREADS_CLI=""

# Parse CLI args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--no-trim)      SKIP_TRIM_CLI=true; shift ;;
    -m|--max-memory)   MEM_CLI="$2"; shift 2 ;;
    -t|--threads)      THREADS_CLI="$2"; shift 2 ;;
    -h|--help)         usage; exit 0 ;;
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
TRIM_ADAPTERS=$(shyaml get-value trimmomatic_adapters < "$CONFIG")

# Strip stray quotes (so paths with spaces work)
for var in READS_DIR TAXO_XLSX REF_FASTA TRIM_ADAPTERS; do
  val="${!var}"
  val="${val#\"}"; val="${val%\"}"
  val="${val#\'}"; val="${val%\'}"
  printf -v "$var" '%s' "$val"
done

 ─── Download & build filtered reference if missing ─────────────────────────
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

# 3) Validate config
die(){ echo "ERROR: $*" >&2; exit 1; }
[[ -d "$READS_DIR"     ]] || die "reads_dir not found: $READS_DIR"
[[ -f "$TAXO_XLSX"     ]] || die "taxonomy_file not found: $TAXO_XLSX"
[[ -f "$REF_FASTA"     ]] || die "ref_fasta not found: $REF_FASTA"
[[ -f "$TRIM_ADAPTERS" ]] || die "trimmomatic_adapters not found: $TRIM_ADAPTERS"

echo "Running with:"
echo "  skip_trimming = $SKIP_TRIM"
echo "  threads       = $THREADS"
echo "  mem_gb        = $MEM_GB"
echo "  adapters      = $TRIM_ADAPTERS"
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

      # Reference selection
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

      # Trimming
      t_trim_start=$(date +%s)

      out_fw="$SAMPLE_DIR/${sample}_trimmed_1P.fastq.gz"
      out_rev="$SAMPLE_DIR/${sample}_trimmed_2P.fastq.gz"
      un_fw="$SAMPLE_DIR/${sample}_trimmed_1U.fastq.gz"
      un_rev="$SAMPLE_DIR/${sample}_trimmed_2U.fastq.gz"

      if [[ "$SKIP_TRIM" == "true" ]]; then
        echo "[`date`] skip_trimming; checking for existing trimmed files"
        if [[ -f "$out_fw" ]] && [[ -f "$out_rev" ]]; then
          echo "[`date`] Found trimmed files; using them"
        else
          die "skip_trimming is true but missing: $out_fw or $out_rev"
        fi

      else
        # decide if we need to run Trimmomatic
        need_trim=false
        if [[ ! -f "$out_fw" ]]; then
          need_trim=true
        elif [[ ! -f "$out_rev" ]]; then
          need_trim=true
        fi

        if [[ "$need_trim" == "true" ]]; then
          conda activate ssuitslsu-trimmomatic
          echo "[`date`] Trimming reads for $sample"
          trimmomatic PE \
            -threads "$THREADS" -phred33 \
            "$R1" "$R2" \
            "$out_fw" "$un_fw" \
            "$out_rev" "$un_rev" \
            ILLUMINACLIP:"${TRIM_ADAPTERS}:2:30:10" \
            SLIDINGWINDOW:4:15 MINLEN:70 HEADCROP:10
        else
          echo "[`date`] Skipping trimming (outputs exist)"
        fi
      fi
      t_trim_end=$(date +%s)
      R1="$out_fw"
      R2="$out_rev"


      # Mapping
      t_map_start=$(date +%s)
      BAM="$SAMPLE_DIR/${sample}.sorted.bam"
      if [[ ! -f "$BAM" ]]; then
        conda activate ssuitslsu-mapping
        [[ ! -f "$SAMPLE_REF.bwt.2bit.64" ]] && bwa-mem2 index "$SAMPLE_REF"
        echo "[`date`] Mapping reads"
        bwa-mem2 mem -v0 -t "$THREADS" "$SAMPLE_REF" "$R1" "$R2" \
          2> "$SAMPLE_DIR/mapping.log" \
          | samtools view -bS - \
          | samtools sort -@ "$THREADS" -o "$BAM"
      else
        echo "[`date`] Skipping mapping (exists)"
      fi
      [[ ! -s "$BAM" ]] && echo "!! Empty BAM; skipping" && continue 2

      # Extract mapped
      P1="$SAMPLE_DIR/${sample}.mapped_1.fastq.gz"
      P2="$SAMPLE_DIR/${sample}.mapped_2.fastq.gz"
      U0="$SAMPLE_DIR/${sample}.mapped_0.fastq.gz"
      if [[ ! -f "$P1" ]]; then
        echo "[`date`] Extracting mapped reads"
        samtools fastq -@ "$THREADS" -1 "$P1" -2 "$P2" -0 "$U0" "$BAM"
      else
        echo "[`date`] Skipping extraction (exists)"
      fi
      t_map_end=$(date +%s)

      # Assembly
      t_asm_start=$(date +%s)
      SP_DIR="$SAMPLE_DIR/spades"
      CONTIG="$SP_DIR/contigs.fasta"
      if [[ ! -f "$CONTIG" ]]; then
        conda activate ssuitslsu-spades
        echo "[`date`] Assembling"
        mkdir -p "$SP_DIR"
        if [[ -f "$U0" ]]; then
          if ! zgrep -q '^@' "$U0"; then
            echo "→ Removing empty singleton"
            rm -f "$U0"
          fi
        fi
        cmd=(spades.py --only-assembler -1 "$P1" -2 "$P2" -t "$THREADS" -m "$MEM_GB" -o "$SP_DIR")
        [[ -f "$U0" ]] && cmd+=( -s "$U0" )
        "${cmd[@]}" 2>&1 | tee "$SP_DIR/spades.log"
      else
        echo "[`date`] Skipping assembly (exists)"
      fi
      t_asm_end=$(date +%s)

      # ITS extraction
      t_its_start=$(date +%s)
      ITSX_DIR="$SAMPLE_DIR/itsx"; mkdir -p "$ITSX_DIR"
      FULL_OUT="$ITSX_DIR/${sample}.full.fasta"

      if [[ -s "$FULL_OUT" ]]; then
        echo "[`date`] Skipping ITSx (already have full region: $FULL_OUT)"
      else
        echo "[`date`] Running ITSx for $sample"
        conda activate ssuitslsu-itsx
        seqkit seq -g -m200 -M100000 "$CONTIG" -o "$ITSX_DIR/filtered.fasta"
        ITSx -i "$ITSX_DIR/filtered.fasta" \
             -o "$ITSX_DIR/${sample}" \
             --preserve T \
             --only_full T \
             --save_regions all \
             -t F \
             --region ITS \
             --cpu "$THREADS"
      fi
      t_its_end=$(date +%s)

      # Timing & report
      t1=$(date +%s)
      end_ts=$(date +"%Y-%m-%d %H:%M:%S")
      elapsed=$((t1 - t0))
      trim_dur=$((t_trim_end - t_trim_start))
      map_dur=$((t_map_end - t_map_start))
      asm_dur=$((t_asm_end - t_asm_start))
      its_dur=$((t_its_end - t_its_start))

      printf "→ Sample: %s\n  Started: %s\n  Finished: %s\n  Elapsed: %02dh:%02dm:%02ds\n\n" \
             "$sample" "$start_ts" "$end_ts" \
             $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
      printf "    Trimming: %02dm%02ds\n" \
         $((trim_dur/60)) $((trim_dur%60))
      printf "    Mapping:  %02dm%02ds\n" \
             $((map_dur/60)) $((map_dur%60))
      printf "    Assembly: %02dm%02ds\n" \
             $((asm_dur/60)) $((asm_dur%60))
      printf "    ITSx:     %02dm%02ds\n\n" \
             $((its_dur/60)) $((its_dur%60))

      break
    fi
  done
done

echo "[`date`] ALL SAMPLES COMPLETE."
