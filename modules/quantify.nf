// Quantification processes - CellRanger, STARsolo, Alevin, Kallisto

process CELLRANGER_COUNT {
    tag { sample }
    cpus params.cellranger_localcores
    memory "${params.cellranger_localmem_gb} GB"
    time '48h'
    input:
    tuple val(sample), path(fq_dir)
    output:
    tuple val(sample), path("cellranger/${sample}"), emit: out
    when:
    params.cellranger_ref != null
    script:
    """
    set -euo pipefail
    mkdir -p cellranger

    # Check if cellranger is available (only if not using container)
    if [[ "${params.use_container}" != "true" ]]; then
        if ! command -v cellranger &> /dev/null; then
            echo "ERROR: cellranger not found in PATH" >&2
            echo "Please ensure cellranger is installed or loaded via: module load cellranger" >&2
            echo "Or set use_container=true to use Docker/Singularity" >&2
            exit 1
        fi
        echo "Using cellranger: \$(which cellranger)" >&2
        cellranger --version >&2
    fi

    echo "Running CellRanger for ${sample} (10x chemistry detected)" >&2

    cellranger count \\
      --id=${sample}_cr \\
      --transcriptome=${params.cellranger_ref} \\
      --fastqs=${fq_dir} \\
      --sample=${sample} \\
      --localcores=${params.cellranger_localcores} \\
      --localmem=${params.cellranger_localmem_gb} \\
      --create-bam=${params.cellranger_create_bam ? 'true' : 'false'} \\
      ${params.cellranger_nosecondary ? '--nosecondary' : ''} \\
      ${params.cellranger_disable_ui ? '--disable-ui' : ''} \\
      ${params.cellranger_expect_cells ? '--expect-cells=' + params.cellranger_expect_cells : ''}

    mv '${sample}_cr' cellranger/${sample}
    """
}

process STARSOLO_COUNT_AUTO {
    tag { sample }
    cpus { params.star_threads }
    memory params.memory_per_sample
    time '48h'
    container params.use_container ? 'quay.io/biocontainers/star:2.7.11a--h0033a41_0' : null

    // Use compiled STAR from ~/bin (2.7.11b_alpha) instead of Homebrew
    beforeScript 'export PATH="$HOME/bin:$PATH" && echo "Using STAR: \$(which STAR)" && STAR --version'

    input:
    tuple val(sample), path(fq_dir), path(chemistry_json)
    output:
    tuple val(sample), path("starsolo/${sample}"), emit: out
    when:
    params.star_index != null
    script:
    """
    set -euo pipefail
    mkdir -p starsolo

    # Parse chemistry detection results
    CHEMISTRY=\$(python3 -c "
import json
with open('${chemistry_json}', 'r') as f:
    data = json.load(f)
    print(data.get('chemistry', 'unknown'))
    ")

    PLATFORM=\$(python3 -c "
import json
with open('${chemistry_json}', 'r') as f:
    data = json.load(f)
    print(data.get('platform', 'unknown'))
    ")

    echo "Detected chemistry: \$CHEMISTRY (platform: \$PLATFORM)" >&2

    # Find FASTQ files (with -L to follow symlinks)
    R1_FILES=\$(find -L ${fq_dir} -name "*_R1_*.fastq.gz" -o -name "*_1.fastq.gz" | sort | paste -s -d ',' -)
    R2_FILES=\$(find -L ${fq_dir} -name "*_R2_*.fastq.gz" -o -name "*_2.fastq.gz" | sort | paste -s -d ',' -)

    if [[ -z "\$R1_FILES" ]] || [[ -z "\$R2_FILES" ]]; then
      echo "ERROR: Missing R1 or R2 FASTQ files for ${sample}" >&2
      exit 1
    fi

    # Set STARsolo parameters based on chemistry
    case "\$CHEMISTRY" in
      "10x_v2"|"10x_v3"|"10x"*)
        echo "Using 10x Genomics parameters for \$CHEMISTRY" >&2
        CB_START=1
        CB_LEN=16
        UMI_START=17
        UMI_LEN=12
        BARCODE_MATE=0
        SOLO_TYPE="CB_UMI_Simple"
        WHITELIST="None"
        ;;

      "seqwell_s3")
        echo "Using Seqwell S^3 parameters" >&2
        CB_START=1
        CB_LEN=12
        UMI_START=13
        UMI_LEN=8
        BARCODE_MATE=1
        SOLO_TYPE="CB_UMI_Simple"
        WHITELIST="${params.seqwell_whitelist}"
        if [[ ! -f "\$WHITELIST" ]]; then
          echo "WARNING: Seqwell whitelist not found, using None" >&2
          WHITELIST="None"
        fi
        ;;

      "gexscope")
        echo "Using GeXscope parameters" >&2
        CB_START=1
        CB_LEN=16
        UMI_START=17
        UMI_LEN=12
        BARCODE_MATE=1
        SOLO_TYPE="CB_UMI_Simple"
        WHITELIST="None"
        ;;

      "custom"|"unknown")
        echo "Using generic single-cell parameters" >&2
        CB_START=\$(python3 -c "
import json
with open('${chemistry_json}', 'r') as f:
    data = json.load(f)
    params = data.get('starsolo_params', {})
    print(params.get('soloCBstart', 1))
        ")
        CB_LEN=\$(python3 -c "
import json
with open('${chemistry_json}', 'r') as f:
    data = json.load(f)
    params = data.get('starsolo_params', {})
    print(params.get('soloCBlen', 16))
        ")
        UMI_START=\$(python3 -c "
import json
with open('${chemistry_json}', 'r') as f:
    data = json.load(f)
    params = data.get('starsolo_params', {})
    print(params.get('soloUMIstart', 17))
        ")
        UMI_LEN=\$(python3 -c "
import json
with open('${chemistry_json}', 'r') as f:
    data = json.load(f)
    params = data.get('starsolo_params', {})
    print(params.get('soloUMIlen', 10))
        ")
        # For 10x platform with custom parameters, barcode is in separate R1 read (mate 0)
        # For non-10x, barcode may be embedded in cDNA read (mate 1) - check platform
        if [[ "\$PLATFORM" == "10x" ]]; then
          BARCODE_MATE=0
          # Select whitelist based on UMI length: v2 (UMI=10bp) vs v3 (UMI=12bp)
          if [[ "\$UMI_LEN" == "10" ]]; then
            WHITELIST="${params.whitelist_10x_v2}"
            echo "Using 10x v2 whitelist (737K)" >&2
          elif [[ "\$UMI_LEN" == "12" ]]; then
            WHITELIST="${params.whitelist_10x_v3}"
            echo "Using 10x v3 whitelist (3M)" >&2
          else
            WHITELIST="None"
            echo "WARNING: Unknown 10x UMI length (\$UMI_LEN), no whitelist used" >&2
          fi
        else
          BARCODE_MATE=1
          WHITELIST="None"
        fi
        SOLO_TYPE="CB_UMI_Simple"
        ;;

      *)
        echo "ERROR: Unsupported chemistry \$CHEMISTRY" >&2
        exit 1
        ;;
    esac

    echo "STARsolo parameters:" >&2
    echo "  CB: start=\$CB_START, len=\$CB_LEN" >&2
    echo "  UMI: start=\$UMI_START, len=\$UMI_LEN" >&2
    echo "  Type: \$SOLO_TYPE" >&2
    echo "  Whitelist: \$WHITELIST" >&2

    # Build STAR command
    STAR_CMD="STAR --runThreadN ${params.star_threads} \\
      --genomeDir ${params.star_index} \\
      --readFilesIn \$R2_FILES \$R1_FILES \\
      --readFilesCommand gunzip -c \\
      --outFileNamePrefix starsolo/${sample}_ \\
      --outSAMtype BAM SortedByCoordinate \\
      --outSAMattributes CB UB GX GN \\
      --soloType \$SOLO_TYPE \\
      --soloCBstart \$CB_START \\
      --soloCBlen \$CB_LEN \\
      --soloUMIstart \$UMI_START \\
      --soloUMIlen \$UMI_LEN \\
      --soloStrand Forward \\
      --soloFeatures Gene GeneFull \\
      --soloUMIdedup 1MM_All \\
      --soloCellReadStats Standard"

    # Add GTF file if provided (required for STARsolo transcriptome features)
    if [[ -n "${params.star_gtf ?: ''}" ]]; then
      STAR_CMD="\$STAR_CMD --sjdbGTFfile ${params.star_gtf}"
    fi

    # Add whitelist if available
    if [[ "\$WHITELIST" != "None" ]]; then
      STAR_CMD="\$STAR_CMD --soloCBwhitelist \$WHITELIST"
    else
      STAR_CMD="\$STAR_CMD --soloCBwhitelist None"
    fi

    # Add barcode read if specified
    if [[ "\$BARCODE_MATE" == "1" ]]; then
      STAR_CMD="\$STAR_CMD --soloBarcodeMate 1"
    fi

    # Add additional parameters for non-10x platforms (Seqwell, GeXscope, etc.)
    # Check platform, not chemistry, to handle custom 10x configurations correctly
    if [[ "\$PLATFORM" != "10x" ]]; then
      STAR_CMD="\$STAR_CMD \\
        --soloUMIfiltering MultiGeneUMI_CR \\
        --soloMultiMappers EM \\
        --limitOutSJcollapsed 10000000 \\
        --limitIObufferSize 300000000"
    fi

    echo "Running STAR command:" >&2
    echo "\$STAR_CMD" >&2

    # Execute STAR
    eval \$STAR_CMD

    # Move and organize output
    mv starsolo/${sample}_Solo.out starsolo/${sample}

    # Create summary compatible with CellRanger format
    mkdir -p starsolo/${sample}/outs
    if [[ -d "starsolo/${sample}/Gene/filtered" ]]; then
      ln -s ../Gene/filtered starsolo/${sample}/outs/filtered_feature_bc_matrix
    fi
    if [[ -d "starsolo/${sample}/Gene/raw" ]]; then
      ln -s ../Gene/raw starsolo/${sample}/outs/raw_feature_bc_matrix
    fi

    echo "STARsolo completed for ${sample} (chemistry: \$CHEMISTRY)" >&2
    """
}

process STARSOLO_COUNT {
    tag { sample }
    cpus { params.star_threads }
    memory params.memory_per_sample
    time '48h'
    container params.use_container ? 'quay.io/biocontainers/star:2.7.11a--h0033a41_0' : null
    input:
    tuple val(sample), path(fq_dir)
    output:
    tuple val(sample), path("starsolo/${sample}"), emit: out
    when:
    params.star_index != null
    shell:
    '''
    set -euo pipefail
    mkdir -p starsolo

    # Find FASTQ files
    R1_FILES=$(find !{fq_dir} -name "*_R1_*.fastq.gz" | sort | paste -sd,)
    R2_FILES=$(find !{fq_dir} -name "*_R2_*.fastq.gz" | sort | paste -sd,)

    if [[ -z "$R1_FILES" ]] || [[ -z "$R2_FILES" ]]; then
      echo "ERROR: Missing R1 or R2 FASTQ files for !{sample}" >&2
      exit 1
    fi

    # Run STARsolo with standard 10x v3 parameters
    STAR \
      --runThreadN !{params.star_threads} \
      --genomeDir !{params.star_index} \
      --readFilesIn $R2_FILES $R1_FILES \
      --readFilesCommand zcat \
      --outFileNamePrefix starsolo/!{sample}_ \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes CB UB GX GN \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist !{params.star_soloCBwhitelist} \
      --soloCBstart 1 \
      --soloCBlen 16 \
      --soloUMIstart 17 \
      --soloUMIlen 12 \
      --soloStrand Forward \
      --soloFeatures Gene GeneFull \
      --soloUMIdedup 1MM_All

    mv starsolo/!{sample}_Solo.out starsolo/!{sample}
    '''
}

process ALEVIN_COUNT {
    tag { sample }
    cpus { params.alevin_threads }
    memory params.memory_per_sample
    time '48h'
    container params.use_container ? 'quay.io/biocontainers/alevin-fry:0.9.0--h9f5acd7_0' : null
    input:
    tuple val(sample), path(fq_dir)
    output:
    tuple val(sample), path("alevin/${sample}"), emit: out
    when:
    params.salmon_index != null
    shell:
    '''
    set -euo pipefail
    mkdir -p alevin

    # Find FASTQ files
    R1_FILES=$(find !{fq_dir} -name "*_R1_*.fastq.gz" | sort | paste -sd" ")
    R2_FILES=$(find !{fq_dir} -name "*_R2_*.fastq.gz" | sort | paste -sd" ")

    if [[ -z "$R1_FILES" ]] || [[ -z "$R2_FILES" ]]; then
      echo "ERROR: Missing R1 or R2 FASTQ files for !{sample}" >&2
      exit 1
    fi

    # Run salmon alevin
    salmon alevin \
      -l ISR \
      -i !{params.salmon_index} \
      -1 $R1_FILES \
      -2 $R2_FILES \
      -o alevin/!{sample}_alevin \
      -p !{params.alevin_threads} \
      --chromiumV3 \
      --dumpFeatures

    # Convert to alevin-fry format
    alevin-fry generate-permit-list \
      -i alevin/!{sample}_alevin \
      -d fw \
      -o alevin/!{sample}_quant \
      -k

    alevin-fry collate \
      -i alevin/!{sample}_quant \
      -r alevin/!{sample}_alevin \
      -t !{params.alevin_threads}

    alevin-fry quant \
      -i alevin/!{sample}_quant \
      -o alevin/!{sample} \
      -t !{params.alevin_threads} \
      -m !{params.salmon_index}/t2g.tsv \
      -r cr-like \
      --use-mtx

    rm -rf alevin/!{sample}_alevin alevin/!{sample}_quant
    '''
}

process KALLISTO_COUNT {
    tag { sample }
    cpus { params.kallisto_threads }
    memory params.memory_per_sample
    time '48h'
    container params.use_container ? 'quay.io/biocontainers/kb-python:0.28.2--pyhdfd78af_0' : null
    input:
    tuple val(sample), path(fq_dir)
    output:
    tuple val(sample), path("kallisto/${sample}"), emit: out
    when:
    params.kallisto_index != null && params.kallisto_t2g != null
    shell:
    '''
    set -euo pipefail
    mkdir -p kallisto

    # Find FASTQ files
    R1_FILES=$(find !{fq_dir} -name "*_R1_*.fastq.gz" | sort | paste -sd" ")
    R2_FILES=$(find !{fq_dir} -name "*_R2_*.fastq.gz" | sort | paste -sd" ")

    if [[ -z "$R1_FILES" ]] || [[ -z "$R2_FILES" ]]; then
      echo "ERROR: Missing R1 or R2 FASTQ files for !{sample}" >&2
      exit 1
    fi

    # Run kallisto bus
    kallisto bus \
      -i !{params.kallisto_index} \
      -o kallisto/!{sample}_bus \
      -x !{params.kallisto_technology} \
      -t !{params.kallisto_threads} \
      $R1_FILES $R2_FILES

    # Sort and correct barcodes
    cd kallisto/!{sample}_bus
    mkdir -p sort_tmp
    bustools sort \
      -t !{params.kallisto_threads} \
      -T sort_tmp \
      -o output.sorted.bus \
      output.bus

    bustools correct \
      -w !{params.kallisto_t2g}/10xv3_whitelist.txt \
      -o output.corrected.bus \
      output.sorted.bus

    # Count
    bustools count \
      -o ../!{sample}/counts \
      -g !{params.kallisto_t2g} \
      -e matrix.ec \
      -t transcripts.txt \
      --genecounts \
      output.corrected.bus

    cd ../..
    rm -rf kallisto/!{sample}_bus
    '''
}
