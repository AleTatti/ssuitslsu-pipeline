#!/usr/bin/env bash
set -euo pipefail

# ──────────────────────────────────────────────────────────────────────────────
# Usage/help message
usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  -n, --no-trim           Skip Trimmomatic trimming (use existing *_trimmed fastqs)
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
    -n|--no-trim)
      SKIP_TRIM_CLI=true
      shift
      ;;
    -m|--max-memory)
      if [[ -n "${2-}" && ! "$2" =~ ^- ]]; then
        MEM_CLI="$2"
        shift 2
      else
        echo "ERROR: --max-memory requires an argument" >&2
        exit 1
      fi
      ;;
    -t|--threads)
      if [[ -n "${2-}" && ! "$2" =~ ^- ]]; then
        THREADS_CLI="$2"
        shift 2
      else
        echo "ERROR: --threads requires an argument" >&2
        exit 1
      fi
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# ──────────────────────────────────────────────────────────────────────────────
# Initialize Conda in non-interactive shells
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
# ──────────────────────────────────────────────────────────────────────────────

# 1) Activate utils for shyaml
conda activate ssuitslsu-utils

# 2) Load config
CONFIG="$(dirname "$0")/../config/pipeline.yaml"
READS_DIR=$(shyaml get-value reads_dir           < "$CONFIG")
TAXO_XLSX=$(shyaml get-value taxonomy_file       < "$CONFIG")
TAXO_SHEET=$(shyaml get-value taxonomy_sheet      < "$CONFIG")
REF_FASTA=$(shyaml get-value ref_fasta           < "$CONFIG")
THREADS=$(shyaml get-value threads               < "$CONFIG")
MEM_GB=$(shyaml get-value mem_gb                 < "$CONFIG")
OUTDIR=$(shyaml get-value outdir                 < "$CONFIG")
TRIM_ADAPTERS=$(shyaml get-value trimmomatic_adapters < "$CONFIG")
TAX_CSV="${OUTDIR}/$(basename "${TAXO_XLSX%.*}").csv"

# Apply CLI overrides
if [[ "$SKIP_TRIM_CLI" == "true" ]]; then
  SKIP_TRIM=true
else
  SKIP_TRIM=$(shyaml get-value skip_trimming < "$CONFIG")
fi
[[ -n "$MEM_CLI" ]]    && MEM_GB="$MEM_CLI"
[[ -n "$THREADS_CLI" ]] && THREADS="$THREADS_CLI"

# 3) Validate config
die(){ echo "ERROR: $*" >&2; exit 1; }
[[ -d "$READS_DIR" ]]  || die "reads_dir not found: $READS_DIR"
[[ -f "$TAXO_XLSX" ]]  || die "taxonomy_file not found: $TAXO_XLSX"
[[ -f "$REF_FASTA" ]]  || die "ref_fasta not found: $REF_FASTA"
[[ -f "$TRIM_ADAPTERS" ]] || die "trimmomatic_adapters not found: $TRIM_ADAPTERS"

echo "Running with:"
echo "  skip_trimming = $SKIP_TRIM"
echo "  threads       = $THREADS"
echo "  mem_gb        = $MEM_GB"
echo "  adapters      = $TRIM_ADAPTERS"
echo

mkdir -p "$OUTDIR"

# 4) Taxonomy conversion (Excel/CSV/TSV → unified CSV)
if [[ ! -f "$TAX_CSV" || "$TAXO_XLSX" -nt "$TAX_CSV" ]]; then
  conda activate ssuitslsu-taxo
  echo "[`date`] Converting taxonomy file '$TAXO_XLSX' → CSV"
  python <<PYCODE
import sys, os
import pandas as pd

tax_file = "${TAXO_XLSX}"
out_csv  = "${TAX_CSV}"
requested = "${TAXO_SHEET}"
ext = os.path.splitext(tax_file)[1].lower()

if ext in [".xlsx", ".xls"]:
    xlsx = pd.ExcelFile(tax_file, engine="openpyxl")
    sheet = requested if requested in xlsx.sheet_names else xlsx.sheet_names[0]
    if sheet != requested:
        print(f"Warning: sheet '{requested}' not found; using '{sheet}'", file=sys.stderr)
    df = pd.read_excel(xlsx, sheet_name=sheet)
elif ext == ".csv":
    df = pd.read_csv(tax_file)
elif ext in [".tsv", ".txt"]:
    df = pd.read_csv(tax_file, sep="\t")
else:
    raise ValueError(f"Unsupported taxonomy file extension: {ext}")

wanted = [
    'Sample ID','Phylum','Class','Order',
    'Family','Subfamily','Tribe',
    'Genus','Species','Subspecies'
]
missing = [c for c in wanted if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

df = df[wanted]
df.to_csv(out_csv, index=False)
PYCODE
else
  echo "[`date`] Skipping taxonomy conversion (up to date)"
fi

# 5) Define read‐suffix patterns
R1_SUFFIXES=( _1.fastq.gz _1.fastq _1.fq.gz _1.fq _R1.fastq.gz _R1.fastq _R1.fq.gz _R1.fq _1.fastq.bz2 _1.fq.bz2 _R1.fastq.bz2 _R1.fq.bz2 )
R2_SUFFIXES=( _2.fastq.gz _2.fastq _2.fq.gz _2.fq _R2.fastq.gz _R2.fastq _R2.fq.gz _R2.fq _2.fastq.bz2 _2.fq.bz2 _R2.fastq.bz2 _R2.fq.bz2 )

# 6) Main loop
for filepath in "${READS_DIR}"/*; do
  for suf in "${R1_SUFFIXES[@]}"; do
    if [[ "$filepath" == *"$suf" ]]; then
      sample=$(basename "$filepath" "$suf")
      R1="$filepath"

      # find matching R2
      R2=""
      for suf2 in "${R2_SUFFIXES[@]}"; do
        cand="${READS_DIR}/${sample}${suf2}"
        [[ -e "$cand" ]] && { R2="$cand"; break; }
      done
      if [[ -z "$R2" ]]; then
        echo "⚠️ No R2 for $sample; skipping"
        continue 2
      fi

      SAMPLE_DIR="$OUTDIR/$sample"
      mkdir -p "$SAMPLE_DIR"

      # timing
      start_ts=$(date +"%Y-%m-%d %H:%M:%S")
      echo "[`date`] ===> Processing $sample"
      t0=$(date +%s)

      # --- Reference selection & hierarchical extraction ---
      taxo_line=$(grep -m1 "^${sample}," "$TAX_CSV") || taxo_line=""
      if [[ -z "$taxo_line" ]]; then
        echo "!! No taxonomy for $sample; skipping"
        continue 2
      fi
      IFS=, read -r _ phylum cls ord fam subfam tribe gen spp subspp <<<"$taxo_line"
      for var in phylum cls ord fam subfam tribe gen spp subspp; do
        eval "$var=\$(echo \"\${$var}\" | sed -e 's/^[^_]*__//' -e 's/^ *//;s/ *\$//')"
      done
      [[ "$spp" == *" "* ]] && spp="${spp##* }"
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
        safe_tag="${tag//;/__}"
        out="$SAMPLE_DIR/${sample}_ref_${safe_tag}.fasta"
        if [[ "$tag" == *";"* ]]; then
          gtag="${tag%%;*}"
          stag="${tag##*;}"
          awk -v g="$gtag" -v s="$stag" 'BEGIN{RS=">";ORS=""} $0~g && $0~s{print ">"$0}' "$REF_FASTA" > "$out"
          level="$gtag;$stag"
        else
          awk -v tag="$tag" 'BEGIN{RS=">";ORS=""} $0~("(^|;)"tag"(;|$)"){print ">"$0}' "$REF_FASTA" > "$out"
          level="$tag"
        fi
        if [[ -s "$out" ]]; then
          echo "    → using reference for level '$level'"
          SAMPLE_REF="$out"
          break
        else
          echo "    → no sequences for level '$level'; trying next"
        fi
      done
      if [[ -z "$SAMPLE_REF" ]]; then
        echo "!! ERROR: no reference extracted; skipping $sample"
        continue 2
      fi

      # --- Trimming ---
      out_fw="$SAMPLE_DIR/${sample}_trimmed_1P.fastq.gz"
      out_rev="$SAMPLE_DIR/${sample}_trimmed_2P.fastq.gz"
      un_fw="$SAMPLE_DIR/${sample}_trimmed_1U.fastq.gz"
      un_rev="$SAMPLE_DIR/${sample}_trimmed_2U.fastq.gz"
      if [[ "$SKIP_TRIM" == "true" ]]; then
        echo "[`date`] skip_trimming=true; using existing trimmed files"
        out_fw="${READS_DIR}/${sample}_trimmed_1P.fastq.gz"
        out_rev="${READS_DIR}/${sample}_trimmed_2P.fastq.gz"
        [[ -f "$out_fw" && -f "$out_rev" ]] || die "skip_trimming but missing trimmed files"
      else
        if [[ ! -f "$out_fw" || ! -f "$out_rev" ]]; then
          conda activate ssuitslsu-trimmomatic
          echo "[`date`] Trimming reads"
          trimmomatic PE -threads "$THREADS" -phred33 \
            "$R1" "$R2" \
            "$out_fw" "$un_fw" \
            "$out_rev" "$un_rev" \
            ILLUMINACLIP:"${TRIM_ADAPTERS}:2:30:10" \
            SLIDINGWINDOW:4:15 MINLEN:70 HEADCROP:10
        else
          echo "[`date`] Skipping trimming (exists)"
        fi
      fi
      R1="$out_fw"; R2="$out_rev"

      # --- Mapping & BAM creation ---
      BAM="$SAMPLE_DIR/${sample}.sorted.bam"
      if [[ ! -f "$BAM" ]]; then
        conda activate ssuitslsu-mapping
        if [[ ! -f "${SAMPLE_REF}.bwt.2bit.64" ]]; then
          echo "[`date`] Indexing reference"
          bwa-mem2 index "$SAMPLE_REF"
        fi
        echo "[`date`] Mapping reads"
        set +o pipefail
        bwa-mem2 mem -v 0 -t "$THREADS" "$SAMPLE_REF" "$R1" "$R2" \
          2> "$SAMPLE_DIR/mapping.log" \
          | samtools view -bS - \
          | samtools sort -@ "$THREADS" -o "$BAM"
        status=$?
        set -o pipefail
        (( status != 0 )) && echo "!! Warning: mapping exited with $status"
      else
        echo "[`date`] Skipping mapping (exists)"
      fi
      if [[ ! -s "$BAM" ]]; then
        echo "!! No alignments for $sample; skipping extraction & assembly"
        continue 2
      fi

      # --- Extract mapped reads ---
      P1="$SAMPLE_DIR/${sample}.mapped_1.fastq.gz"
      P2="$SAMPLE_DIR/${sample}.mapped_2.fastq.gz"
      U0="$SAMPLE_DIR/${sample}.mapped_0.fastq.gz"
      if [[ ! -f "$P1" || ! -f "$P2" || ! -f "$U0" ]]; then
        echo "[`date`] Extracting mapped reads"
        samtools fastq -@ "$THREADS" -1 "$P1" -2 "$P2" -0 "$U0" "$BAM"
      else
        echo "[`date`] Skipping extraction (exists)"
      fi

      # --- Assembly with SPAdes ---
      SP_DIR="$SAMPLE_DIR/spades"
      CONTIGS="$SP_DIR/contigs.fasta"
      if [[ ! -f "$CONTIGS" ]]; then
        conda activate ssuitslsu-spades
        echo "[`date`] Assembling"
        mkdir -p "$SP_DIR"
        [[ -f "$U0" && ! zgrep -q '^@' "$U0" ]] && { echo "→ Removing empty singleton"; rm -f "$U0"; }
        cmd=( spades.py --only-assembler --isolate -1 "$P1" -2 "$P2" -t "$THREADS" -m "$MEM_GB" -o "$SP_DIR" )
        [[ -f "$U0" ]] && cmd+=( -s "$U0" )
        "${cmd[@]}" 2>&1 | tee "$SP_DIR/spades.log"
      else
        echo "[`date`] Skipping assembly (exists)"
      fi

      # --- ITSx extraction ---
      ITSX_DIR="$SAMPLE_DIR/itsx"
      mkdir -p "$ITSX_DIR"
      MARKER="$ITSX_DIR/${sample}.ITSx.fungi_ITS.fasta"
      if [[ -s "$MARKER" ]]; then
        echo "[`date`] Skipping ITSx (already done)"
      else
        # drop >100kb contigs
        FILTERED="$SP_DIR/contigs.ITSx.fasta"
        seqkit seq -m 1 -M 100000 "$SP_DIR/contigs.fasta" -o "$FILTERED"
        conda activate ssuitslsu-itsx
        echo "[`date`] Extracting ITSx for $sample"
        ITSx \
          -i "$FILTERED" \
          -o "$ITSX_DIR/${sample}" \
          --preserve T \
          --only_full T \
          --save_regions all \
          --taxon Fungi \
          --region all \
          --cpu "$THREADS"
        echo "[`date`] ITSx outputs in $ITSX_DIR"
      fi

      # --- Timing & report ---
      t1=$(date +%s)
      end_ts=$(date +"%Y-%m-%d %H:%M:%S")
      elapsed=$((t1 - t0))
      printf "→ Sample: %s\n  Started: %s\n  Finished: %s\n  Elapsed: %02dh:%02dm:%02ds\n\n" \
        "$sample" "$start_ts" "$end_ts" \
        $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

      break
    fi
  done
done

echo "[`date`] ALL SAMPLES COMPLETE."
