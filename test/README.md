# Test Dataset

This directory contains a minimal test case for the sc-sra-nf pipeline.

## Test Dataset: 10x PBMC 1k cells

We use a small publicly available 10x Chromium PBMC dataset for testing.

**Dataset**: 1k PBMCs from a healthy donor
**GEO Accession**: GSM3490286
**SRA Run**: SRR8206317
**BioProject**: PRJNA515520
**Technology**: 10x Chromium Single Cell 3'
**Sequencer**: Illumina NovaSeq 6000
**Reads**: ~9M spots (544 MB)
**Expected cells**: ~1,000
**Estimated runtime**: 15-45 minutes (depending on system and download speed)

## Quick Test

```bash
# Navigate to test directory
cd test

# Run automated test script (recommended)
./run_test.sh /path/to/star_index
# OR
./run_test.sh /path/to/cellranger_ref

# Manual run with STARsolo (requires STAR index)
nextflow run ../main.nf \
  --sra_table SraRunTable.csv \
  --star_index /path/to/star_index \
  --outdir test_results

# Manual run with CellRanger (requires CellRanger reference)
nextflow run ../main.nf \
  --sra_table SraRunTable.csv \
  --cellranger_ref /path/to/cellranger_ref \
  --outdir test_results
```

## Expected Results

- **Format Detection**: Pipeline will detect available format (FASTQ, BAM, or SRA)
- **Chemistry**: 10x Chromium (auto-detected from read structure)
- **Quantification**: ~800-1,200 cells detected
- **Runtime**: ~15-45 minutes total (including download)
- **Output**: Expression matrix in `test_results/GSM3490286/quantification/outs/`

## Files Included

- `SraRunTable.csv`: Metadata table for SRR8206317 (1k PBMC sample)
- `run_test.sh`: Automated test script
- `README.md`: This file

## Troubleshooting

**Download fails:**
- Check internet connection
- Verify SRA accession is accessible: `prefetch --max-size 1G SRR8206317`

**Out of memory:**
- Reduce memory requirements in test run:
  ```bash
  nextflow run ../main.nf --sra_table SraRunTable.csv \
    --star_index /path/to/star_index \
    --outdir test_results \
    --star_threads 4  # Reduce threads if needed
  ```

**No quantifier available:**
- Ensure you have either STAR index or CellRanger reference installed
- See main README.md for reference setup instructions
