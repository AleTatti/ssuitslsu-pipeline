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
  --assembler [spades|megahit]
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
    --assembler) ASSEMBLER="$2"; shift 2 ;;
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
TRIMMOMATIC_ADAPTERS=$(shyaml get-value trimmomatic_adapters < "$CONFIG")
# ensure it's never unset or empty
TRIMMOMATIC_ADAPTERS=${TRIMMOMATIC_ADAPTERS:-auto}

# new coverage thresholds (×)
max_coverage=$(shyaml get-value max_coverage   < "$CONFIG")
target_coverage=$(shyaml get-value target_coverage < "$CONFIG")
# fall back to defaults if not set in pipeline.yaml
max_coverage=${max_coverage:-100}
target_coverage=${target_coverage:-90}


# Strip stray quotes (so paths with spaces work)
for var in READS_DIR TAXO_XLSX REF_FASTA TRIMMOMATIC_ADAPTERS; do
  val="${!var}"
  val="${val#\"}"; val="${val%\"}"
  val="${val#\'}"; val="${val%\'}"
  printf -v "$var" '%s' "$val"
done

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

# trimmomatic adapters if not "auto"
if [[ -n "$TRIMMOMATIC_ADAPTERS" && "$TRIMMOMATIC_ADAPTERS" != "auto" ]]; then
  [[ -f "$TRIMMOMATIC_ADAPTERS" ]] \
    || die "trimmomatic_adapters not found: $TRIMMOMATIC_ADAPTERS"
fi

# new: coverage thresholds (for auto‐subsampling)
max_coverage=$(shyaml get-value max_coverage   < "$CONFIG")
target_coverage=$(shyaml get-value target_coverage< "$CONFIG")
# must be numeric and >0
[[ "$max_coverage" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "max_coverage must be numeric: $max_coverage"
[[ "$target_coverage" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "target_coverage must be numeric: $target_coverage"

echo "Running with:"
echo "  skip_trimming    = $SKIP_TRIM"
echo "  threads          = $THREADS"
echo "  mem_gb           = $MEM_GB"
echo "  adapters         = $TRIMMOMATIC_ADAPTERS"
echo "  max_coverage     = ${max_coverage}×"
echo "  target_coverage  = ${target_coverage}×"
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

      # ─── Trimming ────────────────────────────────────────────────
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
          # resolve adapter path
          if [[ "$TRIMMOMATIC_ADAPTERS" == "auto" || -z "$TRIMMOMATIC_ADAPTERS" ]]; then
            # this is where the adapters live in the activated env
            ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq2-PE.fa"
          else
            ADAPTERS="$TRIMMOMATIC_ADAPTERS"
          fi
          echo "[`date`] Trimming reads for $sample"
          trimmomatic PE \
            -threads "$THREADS" -phred33 \
            "$R1" "$R2" \
            "$out_fw" "$un_fw" \
            "$out_rev" "$un_rev" \
            ILLUMINACLIP:"$ADAPTERS:2:30:10" \
            SLIDINGWINDOW:4:15 MINLEN:70 HEADCROP:10
        else
          echo "[`date`] Skipping trimming (outputs exist)"
        fi
      fi
      t_trim_end=$(date +%s)
      R1="$out_fw"
      R2="$out_rev"

      init_p1=$(zgrep -c '^@' "$R1")
      init_p2=$(zgrep -c '^@' "$R2")
      echo "Initial trimmed: $init_p1 reads in R1, $init_p2 reads in R2"


      # ─── Mapping + auto‐subsampling ──────────────────────────────────────────────
      t_map_start=$(date +%s)
      BAM="$SAMPLE_DIR/${sample}.sorted.bam"
      sub1="$SAMPLE_DIR/${sample}.sub_1P.fastq.gz"
      sub2="$SAMPLE_DIR/${sample}.sub_2P.fastq.gz"

      conda activate ssuitslsu-mapping

      # 1) initial mapping if needed
      if [[ ! -f "$BAM" ]]; then
        # build the index if needed
        [[ ! -f "${SAMPLE_REF}.bwt.2bit.64" ]] && bwa-mem2 index "$SAMPLE_REF"

        echo "[$(date)] Mapping reads (MAPQ ≥ 30, properly paired only)"
        bwa-mem2 mem -v0 -t "$THREADS" "$SAMPLE_REF" "$R1" "$R2" \
          2> "$SAMPLE_DIR/mapping.log" \
        | samtools view -b -q 30 -f 2 -F 4 -@ "$THREADS" \
        | samtools sort -@ "$THREADS" -o "$BAM"

        samtools index "$BAM"
      else
        echo "[$(date)] Skipping mapping (BAM exists)"
        [[ ! -f "${BAM}.bai" ]] && samtools index "$BAM"
      fi

      # fail early if no alignments
      [[ ! -s "$BAM" ]] && echo "!! Empty BAM; skipping sample" && continue 2

      # 2) coverage check + optional subsampling
      if [[ ! -f "$sub1" && ! -f "$sub2" ]]; then
        mean_cov=$(samtools depth -a "$BAM" \
                   | awk '{sum+=$3; cnt++} END {print (cnt? sum/cnt : 0)}')
        printf "[%s] Mean coverage: %.1f×\n" "$(date)" "$mean_cov"

        if (( $(echo "$mean_cov > $max_coverage" | bc -l) )); then
          printf "[%s] Coverage > %d×, subsampling to ~%d×\n" \
                 "$(date)" "$max_coverage" "$target_coverage"

          frac=$(awk -v t=$target_coverage -v m=$mean_cov 'BEGIN{print t/m}')

          seqtk sample -s100 "$R1" "$frac" | gzip > "$sub1"
          seqtk sample -s100 "$R2" "$frac" | gzip > "$sub2"

          # swap in subsampled reads
          R1="$sub1"
          R2="$sub2"

          echo "[$(date)] Remapping subsampled reads"
          rm -f "$BAM" "${BAM}.bai"
          bwa-mem2 mem -v0 -t "$THREADS" "$SAMPLE_REF" "$R1" "$R2" \
            2> "$SAMPLE_DIR/mapping.log" \
          | samtools view -b -q 30 -f 2 -F 4 -@ "$THREADS" \
          | samtools sort -@ "$THREADS" -o "$BAM"
          samtools index "$BAM"

          printf "[%s] Subsampled mapping complete\n" "$(date)"
        fi
      else
        echo "[$(date)] Subsample already done; skipping coverage check"
      fi

      # ─── Extract mapped reads ─────────────────────────────────────────────────
      P1="$SAMPLE_DIR/${sample}.mapped_1.fastq.gz"
      P2="$SAMPLE_DIR/${sample}.mapped_2.fastq.gz"
      U0="$SAMPLE_DIR/${sample}.mapped_0.fastq.gz"

      # re-extract if missing, empty, or older than the BAM
      if [[ ! -s "$P1" || "$BAM" -nt "$P1" ]]; then
        echo "[$(date)] Extracting mapped reads"
        rm -f "$P1" "$P2" "$U0"
        samtools fastq -@ "$THREADS" \
        -1 "$P1" -2 "$P2" -0 "$U0" "$BAM"
      else
        echo "[$(date)] Mapped FASTQs up-to-date; skipping extraction"
      fi

      # Count *mapped* reads
      map_p1=$(zgrep -c '^@' "$P1")
      map_p2=$(zgrep -c '^@' "$P2")
      echo "Mapped reads:    $map_p1 reads in R1, $map_p2 reads in R2"

      #mapping rates
      rate1=$(awk "BEGIN { printf \"%.2f\", ($map_p1/$init_p1)*100 }")
      rate2=$(awk "BEGIN { printf \"%.2f\", ($map_p2/$init_p2)*100 }")
      echo "Mapping rate:    R1: $rate1%   R2: $rate2%"

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
          mkdir -p "$SP_DIR"

        # only pass -s if U0 actually has reads
        if [[ -s "$U0" ]]; then
          echo "→ Including $(zgrep -c '^@' "$U0") singletons"
        else
          echo "→ No singletons; not adding -s flag"
          rm -f "$U0"
        fi
        cmd=(spades.py --careful \
          -1 "$P1" -2 "$P2" \
          -t "$THREADS" -m "$MEM_GB" \
          -o "$SP_DIR")
          [[ -s "$U0" ]] && cmd+=( -s "$U0" )
          "${cmd[@]}" 2>&1 | tee "$SP_DIR/spades.log"
        else
          echo "[`date`] Skipping SPAdes assembly (exists)"
        fi

      elif [[ "$ASSEMBLER" == "megahit" ]]; then
        MH_DIR="$ASM_DIR/megahit"
        CONTIG="$MH_DIR/final.contigs.fa"

        if [[ ! -f "$CONTIG" ]]; then
          conda activate ssuitslsu-megahit
          echo "[`date`] Assembling with MEGAHIT"

          # — remove empty singleton just like before
          if [[ -f "$U0" ]] && ! zgrep -q '^@' "$U0"; then
            echo "→ Removing empty singleton"
            rm -f "$U0"
          fi

          # ─── Clear any old output
          rm -rf "$MH_DIR"

          # build and run Megahit 
          cmd=(megahit \
            -1 "$P1" \
            -2 "$P2" \
            --out-dir "$MH_DIR" \
            --num-cpu-threads "$THREADS" \
            --min-contig-len 1000 \
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
      ITSX_DIR="$SAMPLE_DIR/itsx"; mkdir -p "$ITSX_DIR"
      FULL_OUT="$ITSX_DIR/${sample}.full.fasta"

      if [[ -s "$FULL_OUT" ]]; then
        echo "[`date`] Skipping ITSx (found full region: $FULL_OUT)"
      else
        echo "[`date`] Preparing input for ITSx on $sample"
        conda activate ssuitslsu-itsx

        # 1) chunk any contig >100 kb into ≤100 kb pieces
        CHUNKS="$ITSX_DIR/contigs.chunks.fasta"
        awk -v L=100000 '
          BEGIN { RS=">"; ORS="" }
          NR>1 {
            # separate header from sequence
            header = $1
            seq = ""; 
            # join all the rest of the lines into one long seq
            for(i=2; i<=NF; i++) seq = seq $i
            len = length(seq)

            if (len > L) {
              # how many pieces?
              k = int((len - 1)/L) + 1
              for (i = 0; i < k; i++) {
                start = i*L + 1
                printf(">%s_%03d\n%s\n", header, i+1, substr(seq, start, L))
              }
            } else {
              # just print the contig unchanged
              printf(">%s\n%s\n", header, seq)
            }
          }
        ' "$CONTIG" > "$CHUNKS"

        # 2) filter out anything <200 bp
        FILTERED="$ITSX_DIR/filtered.fasta"
        seqkit seq -m200 "$CHUNKS" -o "$FILTERED"

        # 3) if no fragments ≥200 bp, skip ITSx; otherwise run it
        if [[ ! -s "$FILTERED" ]]; then
          echo "[`date`] No fragments ≥200 bp after chunking; skipping ITSx for $sample"
          touch "$ITSX_DIR/${sample}_no_detections.txt"
        else
          echo "[`date`] Running ITSx for $sample"
          ITSx -i "$FILTERED" \
               -o "$sample" \
               --preserve T \
               --only_full T \
               --save_regions all \
               -t F \
               --region all \
               --anchor HMM \
               --cpu "$THREADS"
        fi
      fi
      #               
      t_its_end=$(date +%s)


      # ─── Phylogenetic analysis ────────────────────────────────────────────────
      PHYLO_DIR="$SAMPLE_DIR/phylogeny"
      mkdir -p "$PHYLO_DIR"

      # 0) pick the first ITSx output that exists
      ITS_SRC=""
      for region in full ITS1 ITS2; do
        candidate="$ITSX_DIR/${sample}.${region}.fasta"
        if [[ -s "$candidate" ]]; then
          printf "[%s] Using ITSx %s for phylogeny\n" "$(date)" "$region"
          ITS_SRC="$candidate"
          break
        fi
      done

      if [[ -z "$ITS_SRC" ]]; then
        printf "[%s] No ITSx output (full, ITS1 or ITS2); skipping phylogeny for %s\n\n" \
               "$(date)" "$sample"
      else
        # extract the original contig name
        orig_node=$(grep -m1 '^>' "$ITS_SRC" | sed 's/^>//;s/ .*//')
        orig_contig_fasta="$PHYLO_DIR/${sample}_${orig_node}.fasta"
        awk -v node="$orig_node" '
          BEGIN { RS=">"; ORS="" }
          $1 == node { print ">" $0 }
        ' "$CONTIG" > "$orig_contig_fasta"

        # 1) Build genus‐level FASTA
        GEN_FASTA="$PHYLO_DIR/${gen}.fasta"
        printf "[%s] Building genus FASTA: %s\n" "$(date)" "$GEN_FASTA"
        awk -v tag="g__${gen}" '
          BEGIN { RS=">"; ORS="" }
          $0 ~ ("(^|;)" tag "(;|$)") { print ">" $0 }
        ' "$REF_FASTA" > "$GEN_FASTA"
        # append our contig with sample tag
        seqkit replace \
          -p '^>(.+)' -r ">${sample}_\1" \
          "$orig_contig_fasta" >> "$GEN_FASTA"

        # 2) MAFFT alignment
        PHYLO_ALN="$PHYLO_DIR/${gen}.aln"
        if [[ -s "$PHYLO_ALN" ]]; then
          printf "[%s] Skipping MAFFT (alignment exists): %s\n\n" "$(date)" "$PHYLO_ALN"
        else
          t_align_start=$(date +%s)
          conda activate ssuitslsu-mafft
          printf "[%s] Running MAFFT → %s\n" "$(date)" "$PHYLO_ALN"
          mafft --thread "$THREADS" --auto "$GEN_FASTA" > "$PHYLO_ALN"
          t_align_end=$(date +%s)
          printf "\n"
        fi

        # 3) IQ-TREE ML inference
        PHYLO_PREFIX="$PHYLO_DIR/${gen}"
        TREEFILE="${PHYLO_PREFIX}.treefile"
        if [[ -s "$TREEFILE" ]]; then
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
            -mem "${MEM_GB}" \
            -nt AUTO \
            -bb 1000 \
            -pre "$PHYLO_PREFIX"
          t_phylo_end=$(date +%s)
          printf "\n"
        fi
      fi


      # ─── Timing & report ───────────────────────────────────────────────
      t1=$(date +%s)
      end_ts=$(date +"%Y-%m-%d %H:%M:%S")
      elapsed=$((t1 - t0))

      trim_dur=$((t_trim_end - t_trim_start))
      map_dur=$((t_map_end  - t_map_start))
      asm_dur=$((t_asm_end  - t_asm_start))
      its_dur=$((t_its_end  - t_its_start))
      align_dur=$((t_align_end  - t_align_start))
      phylo_dur=$((t_phylo_end  - t_phylo_start))

      printf "→ Sample: %s\n  Started: %s\n  Finished: %s\n  Elapsed: %02dh:%02dm:%02ds\n\n" \
        "$sample" "$start_ts" "$end_ts" \
        $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

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
