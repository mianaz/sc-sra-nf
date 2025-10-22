// Helper workflows for sample grouping and organization

workflow MAKE_SAMPLE_GROUPS {
    take:
    map_ch

    main:
    def grouped_src = map_ch
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.SampleAccession as String, row.RunAccession as String ] }
        .groupTuple(by: 0)

    def sample_limit = params.sample_limit ? params.sample_limit.toString().toInteger() : null
    if( sample_limit ) {
        grouped_src = grouped_src.take(sample_limit)
    }

    grouped_ch = grouped_src.map { tup ->
        def sample = tup[0]
        def runs = tup[1]
        def of = file("${sample}.runs.tsv")
        of.text = runs.join('\n') + '\n'
        tuple(sample, of)
    }

    emit:
        grouped_ch
}
