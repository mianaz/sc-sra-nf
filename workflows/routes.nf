// Routing workflows for different data sources

include { EXTRACT_BAM_URLS; DOWNLOAD_BAMS; BAMTOFASTQ } from '../modules/convert'
include { DOWNLOAD_RUNS } from '../modules/download'

workflow PROCESS_BAM_ROUTE {
    take:
    format_info
    bin_ch

    main:
    def bam_urls = EXTRACT_BAM_URLS(format_info, bin_ch)
    def bams = DOWNLOAD_BAMS(bam_urls.out)
    def fastqs = BAMTOFASTQ(bams.out)

    emit:
    fastqs.out
}

workflow PROCESS_SRA_ROUTE {
    take:
    format_info
    bin_ch

    main:
    def downloads = DOWNLOAD_RUNS(format_info, bin_ch)

    emit:
    downloads.out
}

workflow PROCESS_MIXED_ROUTE {
    take:
    format_info
    bin_ch

    main:
    def downloads = DOWNLOAD_RUNS(format_info, bin_ch)

    emit:
    downloads.out
}
