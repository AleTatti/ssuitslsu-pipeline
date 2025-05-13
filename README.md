# ssuitslsu-pipeline

**SSU + ITS + LSU extraction**
Extracts ribosomal small subunit (SSU), internal transcribed spacer (ITS), and large subunit (LSU) regions from fungal genome‐skimming data, assembles contigs, annotates ITS regions, and builds per‐genus alignments and phylogenetic trees.

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

2. **Create the Conda environments**
   This pipeline uses separate environments to isolate dependencies:

   ```bash
   conda env create -f envs/utils.yaml
   conda env create -f envs/taxo.yaml
   conda env create -f envs/trimmomatic.yaml
   conda env create -f envs/mapping.yaml
   conda env create -f envs/spades.yaml
   conda env create -f envs/megahit.yaml
   conda env create -f envs/itsx.yaml
   conda env create -f envs/mafft.yaml
   conda env create -f envs/iqtree.yaml
   ```

---

## 📝 Configuration

Edit `config/pipeline.yaml` to point at your data and tweak parameters:

```yaml
reads_dir: "/path/to/raw_reads"
taxonomy_file: "/path/to/your_taxonomy.xlsx"  # BOLD format
taxonomy_sheet: "Taxonomy"
# ref_fasta is auto-downloaded if unset
# ref_fasta: "/path/to/custom_reference.fasta"
trimmomatic_adapters: auto  # auto‐located in the trimmomatic env
tskip_trimming: false # true to use existing *_trimmed FASTQs
auto_subsample: true       # always check coverage and subsample if needed
max_coverage: 100          # threshold (×) above which to downsample
target_coverage: 90        # aim (×) when subsampling
assembler: spades           # or "megahit"
threads: 8
mem_gb: 16
outdir: "results"
```

* **skip\_trimming**: set to `true` to use existing `*_trimmed` FASTQs.
* **threads / mem\_gb**: controls CPU and RAM for mapping, assembly, alignment, and tree building.

---

## 🚀 Running the pipeline

Once your configuration is set, run:

```bash
bash scripts/run_pipeline.sh
```

This will automatically:

1. Download & filter the Eukaryome SSU+ITS+LSU reference (fungi-only).
2. Convert your taxonomy spreadsheet (BOLD format) into a CSV.
3. Trim reads with Trimmomatic (or skip if already trimmed).
4. Map to the best‐matching reference and extract mapped reads.
5. Assemble contigs with SPAdes or Megahit.
6. Computes stats for contigs ≥ 1 kb
7. Extract ITS regions with ITSx.
8. Build per‐genus 45S FASTA sets, align with MAFFT, and infer ML trees with IQ-TREE.
9. Provide per‐sample and per‐step timing breakdown.

---

## 📂 Output structure

```
results/
└── <sample>/
    ├── <sample>_trimmed_1P.fastq.gz
    ├── <sample>_trimmed_2P.fastq.gz
    ├── <sample>.sorted.bam
    ├── <sample>.mapped_1.fastq.gz
    ├── <sample>.mapped_2.fastq.gz
    ├── assembly/
    │   ├── spades/      # contigs.fasta, scaffolds.fasta, spades.log
    │   ├── megahit/     # final.contigs.fa, megahit.log
    │   └── assembly_stats.1000bp.txt
    ├── itsx/
    │   ├── <sample>.ITS1.fasta
    │   ├── <sample>.5_8S.fasta
    │   ├── <sample>.ITS2.fasta
    │   └── <sample>.full.fasta
    └── phylogeny/
        ├── <genus>.fasta       # FASTA of reference + sample contig
        ├── <genus>.aln         # MAFFT alignment
        ├── <genus>.treefile    # IQ-TREE tree
        └── <genus>.log         # IQ-TREE run log
```

---

## 📄 License

This pipeline is released under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## ❓ Questions or Issues

Please open an issue on GitHub if you run into any problems or have feature requests:
[https://github.com/AleTatti/ssuitslsu-pipeline/issues](https://github.com/AleTatti/ssuitslsu-pipeline/issues)
