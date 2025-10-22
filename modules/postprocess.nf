// Post-processing: Velocyto and sample finalization

process VELOCYTO {
    tag { sample }
    cpus 12
    memory '32 GB'
    time '24h'
    conda params.velocyto_conda ?: 'velocyto'
    input:
    tuple val(sample), path(cellranger_dir)
    output:
    tuple val(sample), path("${sample}_velocyto"), emit: out
    when:
    params.run_velocyto && params.velocyto_gtf && params.velocyto_rmsk
    script:
    """
    set -euo pipefail
    mkdir -p ${sample}_velocyto

    echo "Running velocyto for ${sample}" >&2

    velocyto run10x -v \\
        -m ${params.velocyto_rmsk} \\
        ${cellranger_dir} \\
        ${params.velocyto_gtf} \\
        -@ ${task.cpus} \\
        --samtools-memory 10000 \\
        -t uint32

    # Move output loom file
    if [[ -f "${cellranger_dir}/velocyto/${sample}.loom" ]]; then
        mv "${cellranger_dir}/velocyto/${sample}.loom" "${sample}_velocyto/"
    fi
    """
}

process FINALIZE_SAMPLE {
    tag { sample }
    input:
    tuple val(sample), path(fq_dir, stageAs: 'fastqs/*'), path(count_dir, stageAs: 'quantification/*')
    output:
    tuple val(sample), path("output/${sample}"), emit: out
    path "output/${sample}_summary.txt", emit: summary
    script:
    """
    set -euo pipefail

    # Determine if this is CellRanger or STARsolo output
    if [[ -d "${count_dir}/outs" ]]; then
        # CellRanger or STARsolo with compatible output
        OUTPUT_TYPE="standard"
        outs_dir="${count_dir}/outs"
    elif [[ -d "${count_dir}/Gene" ]]; then
        # STARsolo native output
        OUTPUT_TYPE="starsolo"
        outs_dir="${count_dir}/Gene"
    else
        OUTPUT_TYPE="unknown"
        outs_dir=""
    fi

    mkdir -p output/${sample}
    cp -r ${count_dir}/* output/${sample}/

    # Keep or remove FASTQs
    if [[ "${params.keep_fastq}" == "true" ]]; then
      cp -r "${fq_dir}" "output/${sample}/fastqs"
      find "output/${sample}/fastqs" -type f -name "*.fastq" ! -name "*.gz" -exec gzip {} \\;
    fi

    # Clean up junk files
    for junk in SC_RNA_COUNTER_CS _cmdline _filelist _finalstate _invocation _jobmode _log _mrosource _sitecheck _tmpdir _vdrkill; do
      rm -rf "output/${sample}/\$junk" 2>/dev/null || true
    done

    # Extract metrics
    cells="N/A"
    reads="N/A"
    tool="\$OUTPUT_TYPE"

    if [[ "\$OUTPUT_TYPE" == "standard" ]] && [[ -f "\$outs_dir/metrics_summary.csv" ]]; then
      cells=\$(grep "Estimated Number of Cells" "\$outs_dir/metrics_summary.csv" | cut -d',' -f2 | tr -d '"' 2>/dev/null || echo "N/A")
      reads=\$(grep "Number of Reads" "\$outs_dir/metrics_summary.csv" | cut -d',' -f2 | tr -d '"' 2>/dev/null || echo "N/A")
      tool="CellRanger"
    elif [[ "\$OUTPUT_TYPE" == "starsolo" ]] && [[ -f "\$outs_dir/Summary.csv" ]]; then
      cells=\$(grep "Estimated Number of Cells" "\$outs_dir/Summary.csv" | cut -d',' -f2 2>/dev/null || echo "N/A")
      reads=\$(grep "Number of Reads" "\$outs_dir/Summary.csv" | cut -d',' -f2 2>/dev/null || echo "N/A")
      tool="STARsolo"
    fi

    cat <<EOF > "output/${sample}_summary.txt"
Sample: ${sample}
Status: SUCCESS
Tool: \$tool
Cells: \$cells
Reads: \$reads
Completed: \$(date -u +"%Y-%m-%dT%H:%M:%SZ")
EOF

    # Also copy to sample directory for convenience
    cp "output/${sample}_summary.txt" "output/${sample}/sample_summary.txt"
    """
}

process FINALIZE_SAMPLE_WITH_VELOCYTO {
    tag { sample }
    input:
    tuple val(sample), path(fq_dir, stageAs: 'fastqs/*'), path(count_dir, stageAs: 'quantification/*'), path(velocyto_dir, stageAs: 'velocyto_input/*')
    output:
    tuple val(sample), path("output/${sample}"), emit: out
    path "output/${sample}_summary.txt", emit: summary
    script:
    """
    set -euo pipefail

    mkdir -p output/${sample}
    cp -r ${count_dir}/* output/${sample}/

    # Add velocyto results
    if [[ -d "${velocyto_dir}" ]]; then
      cp -r "${velocyto_dir}" "output/${sample}/velocyto"
    fi

    # Keep or remove FASTQs
    if [[ "${params.keep_fastq}" == "true" ]]; then
      cp -r "${fq_dir}" "output/${sample}/fastqs"
      find "output/${sample}/fastqs" -type f -name "*.fastq" ! -name "*.gz" -exec gzip {} \\;
    fi

    # Clean up junk
    for junk in SC_RNA_COUNTER_CS _cmdline _filelist _finalstate _invocation _jobmode _log _mrosource _sitecheck _tmpdir _vdrkill; do
      rm -rf "output/${sample}/\$junk" 2>/dev/null || true
    done

    # Extract metrics
    cells=\$(grep "Estimated Number of Cells" "output/${sample}/outs/metrics_summary.csv" | cut -d',' -f2 | tr -d '"' 2>/dev/null || echo "N/A")
    reads=\$(grep "Number of Reads" "output/${sample}/outs/metrics_summary.csv" | cut -d',' -f2 | tr -d '"' 2>/dev/null || echo "N/A")

    cat <<EOF > "output/${sample}_summary.txt"
Sample: ${sample}
Status: SUCCESS
Tool: CellRanger
Cells: \$cells
Reads: \$reads
Velocyto: Completed
Completed: \$(date -u +"%Y-%m-%dT%H:%M:%SZ")
EOF

    # Also copy to sample directory for convenience
    cp "output/${sample}_summary.txt" "output/${sample}/sample_summary.txt"
    """
}
