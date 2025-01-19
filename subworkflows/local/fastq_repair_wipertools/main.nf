include { WIPERTOOLS_FASTQWIPER   } from '../../../modules/nf-core/wipertools/fastqwiper'
include { WIPERTOOLS_FASTQSCATTER } from '../../../modules/nf-core/wipertools/fastqscatter'
include { WIPERTOOLS_FASTQGATHER  } from '../../../modules/nf-core/wipertools/fastqgather'
include { WIPERTOOLS_REPORTGATHER  } from '../../../modules/nf-core/wipertools/reportgather'


workflow FASTQ_REPAIR_WIPERTOOLS {
    take:
    ch_fastq // channel: [ val(meta), [ .fastq ] ]

    main:
    // decouple fastq files [sample_id = id; id is the file name]
    ch_decoupled = ch_fastq.flatMap {
        meta, fastqList -> fastqList instanceof List
        ? fastqList.collect { fastq -> tuple("id": (fastq.name.endsWith(".gz") ? fastq.getBaseName(2): fastq.baseName), "sample_id":meta.id, "single_end": meta.single_end, fastq) }
        : [tuple("id": (fastqList.name.endsWith(".gz") ? fastqList.getBaseName(2): fastqList.baseName), "sample_id":meta.id, "single_end": meta.single_end, fastqList)]
    }

    // Split fastq files into chunks
    WIPERTOOLS_FASTQSCATTER(
        ch_decoupled,
        params.num_splits
    )

    /* decouple chunks
    id          = the chunk file name
    mate_id     = (previous) id, i.e., the original, not split, file name
    sample_id   = the id of the original meta
    single_end  = the single_end status of the original meta
    */
    ch_scattered_fastq = WIPERTOOLS_FASTQSCATTER.out.fastq_chunks.flatMap {
        meta, fastqList -> fastqList.collect {
            fastq -> tuple("id": (fastq.name.endsWith(".gz") ? fastq.getBaseName(2) : fastq.baseName),
            "mate_id":meta.id, "sample_id":meta.sample_id, "single_end": meta.single_end, fastq) } }

    // Wipe fastq files
    WIPERTOOLS_FASTQWIPER {
        ch_scattered_fastq
    }

    // group wiped chunks
    ch_cleaned_fastq = Channel.empty()
    WIPERTOOLS_FASTQWIPER.out.wiped_fastq
    | map { meta, fq -> [meta.subMap('mate_id', 'sample_id', 'single_end'), fq]}
    | groupTuple
    | map { meta, fq -> [['id':meta.mate_id+'_wiped', 'sample_id':meta.sample_id, 'single_end':meta.single_end], fq]}
    | set { ch_cleaned_fastq }

    // Gather fastq files
    WIPERTOOLS_FASTQGATHER(
        ch_cleaned_fastq
    )

    // group wiping reports
    ch_gathered_report = Channel.empty()
    WIPERTOOLS_FASTQWIPER.out.report
    | map { meta, report -> [meta.subMap('mate_id', 'sample_id', 'single_end'), report]}
    | groupTuple
    | map { meta, report -> [['id':meta.mate_id+'_wiped', 'sample_id':meta.sample_id, 'single_end':meta.single_end], report]}
    | set { ch_gathered_report }

    // Gather report files
    WIPERTOOLS_REPORTGATHER(
        ch_gathered_report
    )

    // gather versions
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(
        WIPERTOOLS_FASTQSCATTER.out.versions.first(),
        WIPERTOOLS_FASTQWIPER.out.versions.first(),
        WIPERTOOLS_FASTQGATHER.out.versions.first(),
        WIPERTOOLS_REPORTGATHER.out.versions.first()
    )

    emit:
    wiped_fastq = WIPERTOOLS_FASTQGATHER.out.gathered_fastq     // channel: [ val(meta), [ .fastq|.fastq.gz ] ]
    report      = WIPERTOOLS_REPORTGATHER.out.gathered_report   // channel: [ val(meta), [ .txt ] ]
    versions    = ch_versions                                   // channel: [ versions.yml ]
}
