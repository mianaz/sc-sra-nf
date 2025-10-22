// Utility processes for pipeline setup and validation

process PRECHECK {
    tag "precheck"
    executor 'local'
    input:
    val debug
    path bin_dir
    output:
    path 'precheck.ok'
    script:
    """
    set -euo pipefail
    echo "Running comprehensive precheck..." >&2

    # Set environment variables for precheck script
    export DEBUG=${debug}
    export USE_CONTAINER="${params.use_container}"
    export CONTAINER_ENGINE="${params.container_engine}"
    export QUANTIFIER="${params.quantifier}"

    # Run enhanced precheck script
    ${bin_dir}/precheck.sh \${DEBUG}

    touch precheck.ok
    """
}

process CLEAN_TABLE {
    tag "clean_table"
    input:
    path sra_table
    path bin_dir
    output:
    path 'sample_run_map.tsv', emit: map
    path 'run_metadata.tsv',  emit: meta
    script:
    """
    set -euo pipefail
    ${bin_dir}/parse_sra_table.py ${sra_table}
    """
}

process VALIDATE_SAMPLE {
    tag { sample }
    input:
    tuple val(sample), path(runs_tsv)
    output:
    tuple val(sample), path("${sample}.runs.tsv"), emit: out
    script:
    """
    set -euo pipefail
    echo "Validating sample ${sample}" >&2

    # Only copy if source and destination are different
    if [[ "${runs_tsv}" != "${sample}.runs.tsv" ]]; then
        cp "${runs_tsv}" "${sample}.runs.tsv"
    else
        echo "File already has correct name, skipping copy" >&2
    fi
    """
}

process RENAME_VALIDATE {
    tag { sample }
    cpus 1
    memory '4 GB'
    time '2h'
    input:
    tuple val(sample), path(run_dir)
    path bin_dir
    output:
    tuple val(sample), path("fastq/${sample}"), emit: out
    script:
    """
    set -euo pipefail
    mkdir -p fastq/${sample}

    # Parse runs from directory structure
    # Handle both flat structure (SRA) and nested structure (bamtofastq)
    runs_list=""
    for run_path in ${run_dir}/*/; do
        if [ -d "\$run_path" ]; then
            run_id=\$(basename "\$run_path")
            # Skip the output directory we just created
            if [ "\$run_id" == "${sample}" ]; then
                continue
            fi
            # Check if this directory contains FASTQ files or subdirectories with FASTQ files
            fastq_count=\$(find "\$run_path" -maxdepth 2 -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
            if [ "\$fastq_count" -gt 0 ]; then
                if [ -n "\$runs_list" ]; then
                    runs_list="\${runs_list},\${run_id}"
                else
                    runs_list="\$run_id"
                fi
            fi
        fi
    done

    # Process each run
    for run_path in ${run_dir}/*/; do
        if [ -d "\$run_path" ]; then
            run_id=\$(basename "\$run_path")
            # Skip the output directory
            if [ "\$run_id" == "${sample}" ]; then
                continue
            fi
            # Check if this directory contains FASTQ files
            fastq_count=\$(find "\$run_path" -maxdepth 2 -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
            if [ "\$fastq_count" -gt 0 ]; then
                echo "Processing run \$run_id for sample ${sample}" >&2

                # Detect if this is bamtofastq output (nested structure with many files)
                # bamtofastq creates subdirectories like: SRR_ID/flowcell_name/*.fastq.gz
                flat_fastq_count=\$(find "\$run_path" -maxdepth 1 -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
                nested_fastq_count=\$(find "\$run_path" -mindepth 2 -maxdepth 2 -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)

                if [ "\$flat_fastq_count" -gt 0 ]; then
                    # Flat structure - use standard script
                    echo "Detected flat FASTQ structure (SRA-style)" >&2
                    ${bin_dir}/rename_validate_fastq.py \
                        "\$run_id" \
                        "${sample}" \
                        "\$run_path" \
                        "fastq/${sample}" \
                        "\$runs_list"
                elif [ "\$nested_fastq_count" -gt 0 ]; then
                    # Nested structure - use BAM script
                    echo "Detected nested FASTQ structure (bamtofastq-style) with \$fastq_count files" >&2
                    ${bin_dir}/rename_validate_bam_fastq.py \
                        "${sample}" \
                        "\$run_path" \
                        "fastq/${sample}"
                else
                    echo "ERROR: FASTQs found but cannot determine structure" >&2
                    exit 1
                fi
            fi
        fi
    done

    echo "Renamed and validated FASTQs for ${sample}" >&2
    """
}

process PUBLISH_RESULTS {
    tag { sample }
    publishDir "${params.outdir}", mode: 'copy'
    input:
    tuple val(sample), path(sample_dir)
    output:
    path "${sample}", emit: out
    script:
    """
    cp -r ${sample_dir} ${sample}
    """
}

process AGGREGATE_SUMMARY {
    tag "pipeline_summary"
    publishDir params.outdir, mode: 'copy'
    input:
    path summaries
    output:
    path 'pipeline_summary.txt'
    script:
    """
    echo "Pipeline Summary - \$(date)" > pipeline_summary.txt
    echo "=================================" >> pipeline_summary.txt
    echo "" >> pipeline_summary.txt
    for summary in ${summaries}; do
        if [ -f "\$summary" ]; then
            cat "\$summary" >> pipeline_summary.txt
            echo "" >> pipeline_summary.txt
            echo "---" >> pipeline_summary.txt
            echo "" >> pipeline_summary.txt
        fi
    done
    """
}
