/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                  } from '../modules/nf-core/fastqc/main'
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { GZRT                    } from '../modules/nf-core/gzrt/main'
include { BBMAP_REPAIR            } from '../modules/nf-core/bbmap/repair/main'
include { FASTQ_REPAIR_WIPERTOOLS } from '../subworkflows/local/fastq_repair_wipertools/main'
include { COLLECTRESULTS          } from '../modules/local/collectresults'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_fastqrepair_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FASTQREPAIR {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    // TODO: add integrity check in samplesheet (i.e., check that paired fastq files on each line have the same extensions: .gz or .fastq, .fq)

    main:
    ch_final = Channel.empty()      // channel: repaired fastq files
    ch_versions = Channel.empty()   // channel: versions of the software used in the pipeline

    // branch .gz and non gz files
    ch_fastq_ext = Channel.empty()
    ch_samplesheet
    | branch { _map, fq ->
        gz_files: fq.first().getExtension() == 'gz'
        non_gz_files: true }
    | set { ch_fastq_ext }

    //
    // Recover corrupted gz files
    //
    GZRT (ch_fastq_ext.gz_files)

    // Join recovered gz files with non-gz files
    ch_recovered_fastq = Channel.empty()
    GZRT.out.recovered
    | concat ( ch_fastq_ext.non_gz_files )
    | set { ch_recovered_fastq }

    //
    // Make fastq compliant and wipe bad characters
    //
    FASTQ_REPAIR_WIPERTOOLS (ch_recovered_fastq)

    //
    // Branch single- and paired-end reads for optional analyses
    ch_repaired_fastq = Channel.empty()
    FASTQ_REPAIR_WIPERTOOLS.out.wiped_fastq
    | branch {
        single_end: it[0].single_end == true
        paired_end: it[0].single_end == false }
    | set { ch_repaired_fastq }

    // Rename meta.id for single-end reads
    ch_repaired_fastq_single_end = Channel.empty()
    ch_repaired_fastq.single_end
    | map { meta, fq -> [meta.subMap('sample_id', 'single_end'), fq]}
    | map { meta, fq -> [['id':meta.sample_id + '_recovered_wiped', 'single_end':meta.single_end], fq]}
    | set { ch_repaired_fastq_single_end }

    // Group paired-reads by 'sample_id' and rename keys
    ch_repaired_fastq_paired_end = Channel.empty()
    ch_repaired_fastq.paired_end
    | map { meta, fq -> [meta.subMap('sample_id', 'single_end'), fq]}
    | groupTuple
    | map { meta, fq -> [['id':meta.sample_id + '_recovered_wiped', 'single_end':meta.single_end], fq]}
    | set { ch_repaired_fastq_paired_end }
    //

    //
    // Settle reads pairing (re-pair, optional)
    //
    if (!params.skip_bbmap_repair) {
        // Re-pair reads
        BBMAP_REPAIR (ch_repaired_fastq_paired_end, false)

        ch_repaired_fastq_paired_end_singleton = Channel.empty()
        BBMAP_REPAIR.out.repaired
        | concat ( BBMAP_REPAIR.out.singleton )
        | groupTuple
        | map { meta, fq -> [meta, fq.flatten()] }
        | set { ch_repaired_fastq_paired_end_singleton }

        ch_final = ch_repaired_fastq_paired_end_singleton.concat(ch_repaired_fastq_single_end)
        ch_versions = ch_versions.mix(BBMAP_REPAIR.out.versions.first())
    } else {
        ch_final = ch_repaired_fastq_paired_end.concat(ch_repaired_fastq_single_end)
    }

    //
    // local MODULE: COLLECTRESULTS
    //
    collected_fastq = ch_final.flatMap { meta, fastqList ->
        fastqList instanceof List ?
            fastqList.collect { fastq -> tuple(meta, fastq) } :
            [tuple(meta, fastqList)]
    }
    // collected_fastq.view()

    COLLECTRESULTS(
        collected_fastq
    )
    // COLLECTRESULTS.out.renamed_fastq.view()

    //
    // Assess QC of all fastq files (both single and paired end)
    //
    FASTQC ( ch_final )

    //
    // Collate and save software versions
    //
    ch_versions = ch_versions.mix(
        GZRT.out.versions.first(),
        FASTQ_REPAIR_WIPERTOOLS.out.versions,
        FASTQC.out.versions.first(),
        COLLECTRESULTS.out.versions.first()
    )

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf-core_fastqrepair_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})

    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions            = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
