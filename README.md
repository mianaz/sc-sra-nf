# sc-sra-nf

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A525.04-blue.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Automated Nextflow pipeline for processing single-cell RNA-seq data from NCBI SRA. Downloads data, detects chemistry, and quantifies gene expression automatically.

## Features

- Automatic format detection (FASTQ, BAM, or SRA)
- Chemistry identification (10x v2/v3, Seqwell S3, GeXscope)
- Flexible quantification (CellRanger or STARsolo)
- 10x barcode whitelists included
- Optional RNA velocity analysis

## Quick Start

```bash
# Clone repository
git clone https://github.com/mianaz/sc-sra-nf.git
cd sc-sra-nf

# Quick test with included test dataset (1k PBMCs, ~15-45 min)
cd test
./run_test.sh /path/to/star_index

# Or download your own dataset from GEO
# Visit: https://www.ncbi.nlm.nih.gov/geo/
# Search for your study → "SRA Run Selector" → Download "Metadata"

# Run pipeline on your data
nextflow run main.nf \
  --sra_table SraRunTable.csv \
  --star_index /path/to/star_index \
  --outdir results
```

## Installation

### Requirements

**Required:**
- [Nextflow](https://www.nextflow.io/) ≥ 25.04
- Python ≥ 3.8
- `aria2c` or `wget`

**Optional (depending on data type):**
- [STAR](https://github.com/alexdobin/STAR) ≥ 2.7.10 (for STARsolo)
- [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) ≥ 8.0 (for 10x data)
- [SRA Toolkit](https://github.com/ncbi/sra-tools) (if direct FASTQ/BAM unavailable)
- Docker (optional, for containerized CellRanger)

### Setup

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Install aria2 (Ubuntu/Debian)
sudo apt-get install aria2

# Install aria2 (macOS)
brew install aria2
```

### Reference Data

Barcode whitelists are already included in `resources/whitelists/`.

For quantification, download a reference genome:

**STAR index (recommended):**
```bash
# Download genome and annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
gunzip *.gz

# Build index (requires 32GB RAM)
STAR --runMode genomeGenerate \
  --genomeDir ./star_index \
  --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile gencode.v44.primary_assembly.annotation.gtf \
  --runThreadN 8
```

**CellRanger reference:**
```bash
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
```

## Usage

### Basic Usage

**With STARsolo:**
```bash
nextflow run main.nf \
  --sra_table SraRunTable.csv \
  --star_index ./star_index \
  --outdir results
```

**With CellRanger:**
```bash
nextflow run main.nf \
  --sra_table SraRunTable.csv \
  --cellranger_ref ./refdata-gex-GRCh38-2024-A \
  --outdir results
```

**With RNA velocity:**
```bash
nextflow run main.nf \
  --sra_table SraRunTable.csv \
  --cellranger_ref ./refdata-gex-GRCh38-2024-A \
  --run_velocyto true \
  --velocyto_gtf genes.gtf \
  --velocyto_rmsk repeatmasker.gtf \
  --cellranger_create_bam true \
  --outdir results
```

### Parameters

```nextflow
--sra_table              SraRunTable.csv from GEO (required)
--star_index             Path to STAR index (optional)
--cellranger_ref         Path to CellRanger reference (optional)
--outdir                 Output directory (default: results)
--keep_fastq             Keep FASTQ files (default: false)
--run_velocyto           Run RNA velocity analysis (default: false)
--velocyto_gtf           GTF file for Velocyto
--velocyto_rmsk          RepeatMasker GTF for Velocyto
```

## Output

```
results/
├── GSM1234567/
│   ├── quantification/
│   │   ├── outs/
│   │   │   ├── filtered_feature_bc_matrix/
│   │   │   ├── metrics_summary.csv
│   │   │   └── web_summary.html
│   │   └── velocyto/           (if enabled)
│   ├── fastq/                   (if --keep_fastq)
│   └── sample_summary.txt
└── summary/
    └── all_samples_summary.txt
```

## Troubleshooting

**Memory errors:**
Increase memory in `nextflow.config`:
```nextflow
params.cellranger_localmem_gb = 128
```

**SRA download fails (file too large):**
Edit `modules/download.nf`, add `--max-size 100G` to prefetch command.

**Nextflow version errors:**
Ensure you have Nextflow ≥ 25.04 installed.

## Citation

If you use this pipeline, please cite:

```bibtex
@software{sc_sra_nf,
  title = {sc-sra-nf: Automated single-cell RNA-seq processing from SRA},
  year = {2025},
  url = {https://github.com/mianaz/sc-sra-nf},
  version = {1.0}
}
```

## License

MIT License - see [LICENSE](LICENSE) file.

## Support

Issues and questions: [GitHub Issues](https://github.com/mianaz/sc-sra-nf/issues)
