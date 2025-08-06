# ssuitslsu-pipeline

> ⚠️ **Warning:** This pipeline is under active development.  
> If you've cloned it previously, please run `git pull origin main` (or re‐clone) to get the latest changes before use.

**SSU + ITS + LSU extraction**
Extracts ribosomal small subunit (SSU), internal transcribed spacer (ITS), and large subunit (LSU) regions from fungal genome‐skimming data, assembles contigs, annotates ITS regions, builds per‐genus alignments and phylogenetic trees, and includes sliding‐window chimera detection.

---

## 🔧 Prerequisites

* **Git** (to clone this repo)
* **Conda** (or [Mamba](https://mamba.readthedocs.io/) for faster installs)
* Unix‐like environment (macOS, Linux, or WSL on Windows)

---

## ⚙️ Installation

1. **Clone the repository**

   ```bash
   git clone https://github.com/AleTatti/ssuitslsu-pipeline.git
   cd ssuitslsu-pipeline
   ```

2. **Set up Conda environments**
   
   **Option A: Automated setup (Recommended)**
   ```bash
   ./scripts/setup_environments.sh
   ```
   
   **Option B: Manual setup**
   ```bash
   conda env create -f envs/utils.yaml
   conda env create -f envs/fastp.yaml
   conda env create -f envs/mapping.yaml
   conda env create -f envs/itsx.yaml
   conda env create -f envs/chimera.yaml
   conda env create -f envs/mafft.yaml
   conda env create -f envs/iqtree.yaml
   conda env create -f envs/spades.yaml
   conda env create -f envs/megahit.yaml
   conda env create -f envs/taxo.yaml
   ```

3. **Keep environments updated**
   
   In addition to updating the repository, make sure to regularly update the Conda environments to keep dependencies current. You can do this by re-running:
   
   ```bash
   ./scripts/setup_environments.sh
   ```
   
   Or update individual environments:
   ```bash
   conda env update -f envs/utils.yaml --prune
   conda env update -f envs/mapping.yaml --prune
   # ... etc for other environments
   ```

---

## 📝 Configuration

Edit `config/pipeline.yaml` to point at your data and tweak parameters:

```yaml
reads_dir: "/path/to/raw_reads"
# Taxonomy file (Excel or CSV in BOLD format)
taxonomy_file: "/path/to/your_taxonomy.xlsx"
taxonomy_sheet: "Taxonomy"
# Reference FASTA (optional; auto-downloaded if empty)
ref_fasta: "/path/to/custom_reference.fasta"
# Trimming
skip_trimming: false  # Set to true to skip fastp and use existing *_trimmed FASTQs
# Coverage-based downsampling
auto_subsample: true         # Enable automatic coverage check and downsampling
max_coverage: 100            # Threshold (×) above which to downsample
target_coverage: 90          # Target (×) coverage after subsampling
# Mapping quality filtering
mapq: 0                      # Minimum MAPQ to retain alignments (samtools -q)
# Soft-clipped read filtering
filter_softclipped_reads: false    # Enable filtering of soft-clipped reads
min_softclip_bases: 10             # 'auto' for automatic or integer minimum number of soft-clipped bases to trigger filtering
softclip_filter_mode: trim         # Options: "full" (remove entire read) or "trim" (trim softclipped ends only)
# Assembly
assembler: spades           # Options: "spades" or "megahit"
# Computational resources that controls CPU and RAM for mapping, assembly, alignment, and tree building.
threads: 8                  # Number of threads to use
mem_gb: 16                  # Total RAM in GB
# Output
outdir: "results"
# Phylogenetic analysis
skip_phylogeny: true  # if true, skip IQ-TREE phylogenetic inference for faster run
```

---

## 🚀 Running the pipeline

1. **Prepare your configuration** (edit `config/pipeline.yaml`)

2. **Run the pipeline:**
   ```bash
   ./scripts/run_pipeline.sh
   ```

### Pipeline Steps

The pipeline will automatically:

1. **Environment Setup**: Activate required conda environments
2. **Reference Download**: Download & filter the fungal SSU+ITS+LSU reference database
3. **Taxonomy Conversion**: Convert your taxonomy sheet (BOLD format) into CSV
4. **Quality Control**: Trim reads with fastp (or skip if already trimmed)
5. **Read Mapping**: Map reads to best-matching reference with coverage-based downsampling
6. **Assembly**: Assemble contigs using SPAdes or MEGAHIT
7. **ITS Extraction**: Extract SSU, ITS, and LSU regions with ITSx
8. **Phylogenetic Analysis**: Build per-genus alignments and ML trees (MAFFT + IQ-TREE)
9. **Chimera Detection**: Detect chimeras using sliding p-distance analysis
10. **Reporting**: Generate comprehensive reports and timings

All output is logged to `ssuitslsu_YYYYMMDD_HHMMSS.log`

---


## 🔧 CLI / Usage Options

Usage: `./scripts/run_pipeline.sh [OPTIONS]`

| Option | Description |
|--------|-------------|
| `-n`, `--no-trim` | Skip fastp trimming (use existing FASTQs) |
| `-m`, `--max-memory MB` | Override memory (GB) for SPAdes/MEGAHIT |
| `-t`, `--threads N` | Override number of threads for all steps |
| `--assembler [spades\|megahit]` | Choose assembler ("spades" or "megahit") |
| `--mapq N` | Override mapping-quality filter (samtools -q) |
| `--filter-softclip` | Enable soft-clip filtering of BAM |
| `--min-softclip N` | Override minimum soft-clip bases |
| `--softclip-mode MODE` | Override soft-clip filter mode (`full` \| `trim`) |
| `--auto-subsample [true\|false]` | Override auto_subsample behavior |
| `--max-cov N` | Override `max_coverage` threshold |
| `--target-cov N` | Override `target_coverage` threshold |
| `--reads-dir DIR` | Specify directory containing raw read files |
| `--taxonomy-file FILE` | Path to taxonomy file (Excel .xlsx or .csv) |
| `--taxonomy-sheet SHEET` | Sheet name within taxonomy file |
| `--skip_phylogeny [true\|false]` | Override skip_phylogeny behavior |
| `--ref-fasta FILE` | Reference FASTA for mapping/indexing |
| `--outdir DIR` | Override output directory |
| `-h`, `--help` | Show this help message and exit |

---

## 📂 Output structure

```
results/
└── <sample>/
    ├── <sample>_trimmed_1P.fastq.gz   # after fastp (if run)
    ├── <sample>_trimmed_2P.fastq.gz
    ├── <sample>.sorted.bam             # mapped & subsampled (if applicable)
    ├── <sample>.mapped_1.fastq.gz      # re-extracted paired mapped reads
    ├── <sample>.mapped_2.fastq.gz
    ├── chimera/
    │   ├── softclipped_reads/
    │   │   ├── <sample>_softclip_stats.txt
    │   │   ├── <sample>.softclipped.sam
    │   │   ├── <sample>.softclipped.sorted.bam
    │   ├── <sample>_chimera_report.tsv       # chimera detection report
    │   ├── <sample>_chimera_plot.png         # pandas/Matplotlib figure
    │   └── <sample>_chimera_report.html      # HTML summary
    ├── assembly/
    │   ├── spades/        # contigs.fasta, scaffolds.fasta, spades.log  (if spades)
    │   ├── megahit/       # final.contigs.fa, megahit.log             (if megahit)
    │   └── assembly_stats.1000bp.txt   # stats for contigs ≥1 kb
    ├── itsx/
    │   ├── <sample>.small.full.fasta         # from ITSx (small contigs)
    │   ├── <sample>.small.SSU.fasta          # …
    │   ├── <sample>.small.ITS1.fasta
    │   ├── <sample>.small.5_8S.fasta
    │   ├── <sample>.small.ITS2.fasta
    │   ├── <sample>.small.LSU.fasta
    │   ├── <sample>.large_nhmmer.full.fasta  # from ITSx (large contigs)
    │   ├── <sample>.large_nhmmer.SSU.fasta
    │   ├── <sample>.large_nhmmer.ITS1.fasta
    │   ├── <sample>.large_nhmmer.5_8S.fasta
    │   ├── <sample>.large_nhmmer.ITS2.fasta
    │   └── <sample>.large_nhmmer.LSU.fasta
    └── phylogeny/
        ├── <genus>_contigs.fasta   # merged contigs ≥400 bp & ITSx sequences
        ├── <genus>_ITSx.fasta
        ├── <genus>_sequences.fasta  # sample + reference
        ├── <genus>_ref.fasta        # per-genus reference pulled from 45S DB
        ├── <genus>_phylo.fasta      # concatenated FASTA (ref+sample)
        ├── <genus>.aln              # MAFFT alignment (raw)
        ├── <genus>.trimmed.aln      # trimmed alignment (head/tail removed)
        ├── <genus>.treefile         # IQ-TREE ML tree
        ├── <genus>.log              # IQ-TREE run log
        └── <sample>_ITSx_coords_aln_trimmed.tsv  # mapping of raw→trimmed columns
```

---

## 🛠️ Troubleshooting

### Common Issues

**"Could not find conda environment" error:**
```bash
# Run the setup script to create all environments
./scripts/setup_environments.sh
```
**Environment activation issues:**
The pipeline includes robust environment checking that handles various conda configurations. If you still encounter issues, try recreating the problematic environment:
```bash
conda env remove -n ssuitslsu-[environment-name]
conda env create -f envs/[environment-name].yaml
```

---

## 📄 License

This pipeline is released under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## ❓ Questions or Issues

Please open an issue on GitHub if you run into any problems or have feature requests:
[https://github.com/AleTatti/ssuitslsu-pipeline/issues](https://github.com/AleTatti/ssuitslsu-pipeline/issues)


