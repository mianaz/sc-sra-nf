// BAM to FASTQ conversion processes

process EXTRACT_BAM_URLS {
    tag { "${sample}" }
    cpus 1
    input:
    tuple val(sample), path(runs_tsv), path(format_json)
    path bin_dir
    output:
    tuple val(sample), path('bam_urls.tsv'), emit: out
    script:
    """
    set -euo pipefail

    # First check if URLs are already in the format JSON
    URLS_IN_JSON=\$(python3 -c "
import json
with open('${format_json}', 'r') as f:
    data = json.load(f)
    has_urls = any(info.get('urls', []) for info in data.get('runs', {}).values())
    print('yes' if has_urls else 'no')
")

    if [[ "\$URLS_IN_JSON" == "yes" ]]; then
        # Extract URLs from JSON
        python3 -c "
import json
with open('${format_json}', 'r') as f:
    data = json.load(f)
    with open('bam_urls.tsv', 'w') as out:
        out.write('RunAccession\\tBamUrl\\n')
        for run_id, info in data['runs'].items():
            if info.get('urls'):
                for url in info['urls']:
                    if url.endswith('.bam') or url.endswith('.bam.1'):
                        out.write(f'{run_id}\\t{url}\\n')
                        break  # One URL per run
"
    else
        # URLs not in JSON - use extract_bam_urls.py to find them
        echo "Extracting BAM URLs using extract_bam_urls.py..." >&2
        ${bin_dir}/extract_bam_urls.py ${runs_tsv} bam_urls.tsv ${sample}
    fi
    """
}

process DOWNLOAD_BAMS {
    tag { "${sample}" }
    cpus 4
    time '24h'
    input:
    tuple val(sample), path(bam_urls_tsv)
    output:
    tuple val(sample), path('bams'), emit: out
    script:
    """
    set -euo pipefail
    mkdir -p bams

    # Strip carriage returns and process TSV
    tail -n +2 '${bam_urls_tsv}' | tr -d '\\r' | while IFS=\$'\\t' read -r run_id bam_url; do
      if [[ -n "\$bam_url" && "\$bam_url" != "N/A" ]]; then
        echo "Downloading BAM for \$run_id from \$bam_url" >&2
        if command -v aria2c >/dev/null 2>&1; then
          aria2c -x 16 -s 16 -j 1 -d bams -o "\${run_id}.bam" "\$bam_url"
        else
          wget -O "bams/\${run_id}.bam" "\$bam_url"
        fi
      fi
    done

    # Check if any BAM files were actually downloaded
    bam_count=\$(find bams -name "*.bam" -type f 2>/dev/null | wc -l)
    if [[ \$bam_count -eq 0 ]]; then
      echo "ERROR: No BAM files were downloaded. All URLs were N/A or invalid." >&2
      echo "This sample should use SRA Tools instead of direct BAM download." >&2
      exit 1
    fi

    echo "Successfully downloaded \$bam_count BAM file(s)" >&2
    """
}

process BAMTOFASTQ {
    tag { sample }
    container 'cumulusprod/cellranger:9.0.1'
    cpus 8
    memory '16 GB'
    time '12h'
    input:
    tuple val(sample), path(bam_dir)
    output:
    tuple val(sample), path('fastq'), emit: out
    script:
    """
    set -euo pipefail
    mkdir -p fastq
    for bam_file in ${bam_dir}/*.bam; do
      [[ -f "\$bam_file" ]] || continue
      run_id=\$(basename "\$bam_file" .bam)
      echo "Converting \$bam_file to FASTQ" >&2
      cellranger bamtofastq --nthreads=${task.cpus} "\$bam_file" "fastq/\${run_id}"
    done
    """
}
