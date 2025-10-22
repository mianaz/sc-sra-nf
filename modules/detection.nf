// Detection processes for data format and chemistry identification

process DETECT_DATA_FORMAT {
    tag { sample }
    cpus 1
    input:
    tuple val(sample), path(runs_tsv)
    path bin_dir
    output:
    tuple val(sample), path("${sample}.runs.tsv"), path("${sample}_format.json"), emit: out
    script:
    """
    set -euo pipefail
    ${bin_dir}/detect_data_format.py ${runs_tsv} ${sample}_format.json ${sample}

    # Only copy if source and destination are different
    if [[ "${runs_tsv}" != "${sample}.runs.tsv" ]]; then
        cp ${runs_tsv} ${sample}.runs.tsv
    else
        echo "File already has correct name, skipping copy" >&2
    fi
    """
}

process DETECT_CHEMISTRY {
    tag { sample }
    cpus 1
    input:
    tuple val(sample), path(fq_dir)
    path bin_dir
    output:
    tuple val(sample), path(fq_dir), path("${sample}_chemistry.json"), emit: out
    script:
    """
    set -euo pipefail
    ${bin_dir}/detect_chemistry.py ${fq_dir} ${sample}_chemistry.json ${sample}

    # Log chemistry detection result
    CHEMISTRY=\$(python3 -c "
import json
with open('${sample}_chemistry.json', 'r') as f:
    data = json.load(f)
    print(data['chemistry'] + ' (' + data['platform'] + ') -> ' + data['quantifier'])
    ")

    echo "Sample ${sample}: \$CHEMISTRY" >&2
    """
}
