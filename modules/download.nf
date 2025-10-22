// Download processes for FASTQ and SRA data

process DOWNLOAD_FASTQ_DIRECT {
    tag { sample }
    cpus 4
    time '12h'
    input:
    tuple val(sample), path(runs_tsv), path(format_json)
    path bin_dir
    output:
    tuple val(sample), path('downloads'), emit: out
    script:
    """
    set -euo pipefail
    mkdir -p downloads

    # Parse JSON to get FASTQ URLs
    python3 -c "
import json
import sys
with open('${format_json}', 'r') as f:
    data = json.load(f)
    for run_id, info in data['runs'].items():
        if info['urls']:
            for url in info['urls']:
                print(f'{run_id}\\t{url}')
    " > fastq_urls.tsv

    # Download FASTQs using aria2c for speed
    while IFS=\$'\\t' read -r run_id url; do
        echo "Downloading FASTQ for \$run_id from \$url" >&2

        filename=\$(basename "\$url")
        output_file="downloads/\${run_id}_\${filename}"

        if command -v aria2c >/dev/null 2>&1; then
            aria2c -x 16 -s 16 -j 1 --max-tries=5 --retry-wait=3 \\
                   --allow-overwrite=true --auto-file-renaming=false \\
                   -d downloads -o "\${run_id}_\${filename}" "\$url"
        elif command -v wget >/dev/null 2>&1; then
            wget -O "\$output_file" "\$url"
        else
            curl -L -o "\$output_file" "\$url"
        fi
    done < fastq_urls.tsv

    echo "FASTQ download completed for ${sample}" >&2
    """
}

process DOWNLOAD_RUNS {
    tag { "${sample}" }
    input:
    tuple val(sample), path(runs_tsv), path(format_json)
    path bin_dir
    output:
    tuple val(sample), path('downloads'), emit: out
    script:
    """
    set -euo pipefail
    mkdir -p downloads

    # Get SRA file URLs from SDL API
    ${bin_dir}/get_sra_urls.py ${runs_tsv} sra_urls.tsv

    # Download and convert SRA files
    tail -n +2 sra_urls.tsv | tr -d '\\r' | while IFS=\$'\\t' read -r SRR url; do
      echo "Processing \$SRR" >&2
      [[ -z "\$SRR" ]] && continue

      # Download SRA file - use aria2c if URL available, otherwise prefetch
      if [[ -n "\$url" && "\$url" != "N/A" ]]; then
        echo "  Downloading \$SRR via aria2c from \$url" >&2
        if command -v aria2c >/dev/null 2>&1; then
          aria2c -x 16 -s 16 -j 1 --max-tries=5 --retry-wait=3 \\
                 --allow-overwrite=true --auto-file-renaming=false \\
                 -d downloads -o "\${SRR}.sra" "\$url"
        else
          echo "  aria2c not available, falling back to wget" >&2
          wget -O "downloads/\${SRR}.sra" "\$url"
        fi
      else
        echo "  Downloading \$SRR via prefetch" >&2
        prefetch --max-size 200G --output-directory downloads "\$SRR"
        if [[ -f "downloads/\$SRR/\$SRR.sra" ]]; then
          mv "downloads/\$SRR/\$SRR.sra" "downloads/\$SRR.sra"
          rmdir "downloads/\$SRR"
        fi
      fi

      # Convert SRA to FASTQ
      echo "  Converting \$SRR to FASTQ" >&2
      fasterq-dump --threads ${params.download_threads} --split-files --include-technical -O downloads "downloads/\$SRR.sra"

      # Compress FASTQ files
      for fq in downloads/\${SRR}_*.fastq; do
        [[ -f "\$fq" ]] || continue
        gzip -f "\$fq"
      done

      # Clean up SRA file
      rm -f "downloads/\$SRR.sra"
    done
    """
}
