#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.help = false

// Import modules
include { PRECHECK; CLEAN_TABLE; VALIDATE_SAMPLE; RENAME_VALIDATE; PUBLISH_RESULTS; AGGREGATE_SUMMARY } from './modules/utils'
include { DETECT_DATA_FORMAT; DETECT_CHEMISTRY } from './modules/detection'
include { DOWNLOAD_FASTQ_DIRECT } from './modules/download'
include { CELLRANGER_COUNT; STARSOLO_COUNT_AUTO } from './modules/quantify'
include { VELOCYTO; FINALIZE_SAMPLE; FINALIZE_SAMPLE_WITH_VELOCYTO } from './modules/postprocess'

// Import workflows
include { MAKE_SAMPLE_GROUPS } from './workflows/helpers'
include { PROCESS_BAM_ROUTE; PROCESS_SRA_ROUTE; PROCESS_MIXED_ROUTE } from './workflows/routes'

// Setup input channels
Channel
    .fromPath(params.sra_table, checkIfExists: true)
    .ifEmpty { error "Missing --sra_table= path to SraRunTable.csv" }
    .set { sra_table_ch }

workflow {
    bin_ch = Channel.of(file("bin"))

    PRECHECK(params.debug ? 1 : 0, bin_ch)

    CLEAN_TABLE(sra_table_ch, bin_ch)

    MAKE_SAMPLE_GROUPS(CLEAN_TABLE.out.map)

    def validated = VALIDATE_SAMPLE(MAKE_SAMPLE_GROUPS.out)

    // Detect available data formats for each sample
    def format_info = DETECT_DATA_FORMAT(validated.out, bin_ch)

    // Branch based on detected format
    def branched = format_info.out.branch { sample, runs_tsv, format_json ->
        def json = new groovy.json.JsonSlurper().parse(format_json)
        def strategy = json.strategy

        fastq_direct: strategy == 'FASTQ_DIRECT'
            return tuple(sample, runs_tsv, format_json)
        bam_direct: strategy == 'BAM_DIRECT'
            return tuple(sample, runs_tsv, format_json)
        sra_tools: strategy == 'SRA_TOOLS'
            return tuple(sample, runs_tsv, format_json)
        mixed: strategy == 'MIXED'
            return tuple(sample, runs_tsv, format_json)
    }

    // Process each strategy
    def fastq_from_direct = DOWNLOAD_FASTQ_DIRECT(branched.fastq_direct, bin_ch)
    def fastq_from_bam = PROCESS_BAM_ROUTE(branched.bam_direct, bin_ch)
    def fastq_from_sra = PROCESS_SRA_ROUTE(branched.sra_tools, bin_ch)
    def fastq_from_mixed = PROCESS_MIXED_ROUTE(branched.mixed, bin_ch)

    // Merge all FASTQ outputs
    def all_fastqs = fastq_from_direct.mix(
        fastq_from_bam,
        fastq_from_sra,
        fastq_from_mixed
    )

    // Rename and validate all FASTQs
    def renamed = RENAME_VALIDATE(all_fastqs, bin_ch)

    // Detect chemistry for each sample
    def chemistry_detected = DETECT_CHEMISTRY(renamed.out, bin_ch)

    // Branch based on detected chemistry
    def chemistry_branched = chemistry_detected.out.branch { sample, fq_dir, chemistry_json ->
        def json = new groovy.json.JsonSlurper().parse(chemistry_json)
        def quantifier = json.quantifier
        def chemistry = json.chemistry

        cellranger: quantifier == 'cellranger' && json.confidence == 'high'
            return tuple(sample, fq_dir)
        starsolo: quantifier == 'starsolo' || json.confidence != 'high'
            return tuple(sample, fq_dir, chemistry_json)
    }

    // Process with appropriate quantifier
    def cellranger_results = Channel.empty()
    def starsolo_results = Channel.empty()

    // Run CellRanger for 10x samples
    if (params.cellranger_ref != null) {
        cellranger_results = CELLRANGER_COUNT(chemistry_branched.cellranger)
    }

    // Run STARsolo for non-10x samples
    if (params.star_index != null) {
        starsolo_results = STARSOLO_COUNT_AUTO(chemistry_branched.starsolo)
    }

    // Merge quantification results
    def counted = cellranger_results.mix(starsolo_results)

    // Optional Velocyto processing (only for samples with BAM files)
    def with_velocyto = Channel.empty()
    def without_velocyto = Channel.empty()

    if (params.run_velocyto && params.cellranger_create_bam) {
        // Branch based on whether sample used CellRanger (has BAMs)
        def for_velocyto = counted.branch { sample, count_dir ->
            def has_bam = file("${count_dir}/outs/possorted_genome_bam.bam").exists()
            with_bam: has_bam
                return tuple(sample, count_dir)
            without_bam: !has_bam
                return tuple(sample, count_dir)
        }

        // Run Velocyto on samples with BAMs
        def velocyto_results = VELOCYTO(for_velocyto.with_bam)

        // Combine results
        with_velocyto = renamed.out
            .join(for_velocyto.with_bam, by: 0)
            .join(velocyto_results.out, by: 0)
            .map { sample, fq_dir, count_dir, velocyto_dir ->
                tuple(sample, fq_dir, count_dir, velocyto_dir)
            }

        without_velocyto = renamed.out
            .join(for_velocyto.without_bam, by: 0)
            .map { sample, fq_dir, count_dir ->
                tuple(sample, fq_dir, count_dir)
            }
    } else {
        without_velocyto = renamed.out.join(counted, by: 0)
    }

    // Finalize samples
    def finalized_with_out = Channel.empty()
    def finalized_with_summary = Channel.empty()
    def finalized_without_out = Channel.empty()
    def finalized_without_summary = Channel.empty()

    if (params.run_velocyto && params.cellranger_create_bam) {
        def finalized_with = FINALIZE_SAMPLE_WITH_VELOCYTO(with_velocyto)
        finalized_with_out = finalized_with.out
        finalized_with_summary = finalized_with.summary
    }

    def finalized_without = FINALIZE_SAMPLE(without_velocyto)
    finalized_without_out = finalized_without.out
    finalized_without_summary = finalized_without.summary

    def all_finalized_out = finalized_with_out.mix(finalized_without_out)
    def all_finalized_summary = finalized_with_summary.mix(finalized_without_summary)

    // Publish results
    def published = PUBLISH_RESULTS(all_finalized_out)

    // Aggregate summaries
    def summary_bundle = all_finalized_summary.collect().filter { it && it.size() > 0 }
    AGGREGATE_SUMMARY(summary_bundle)

    return published.out
}
