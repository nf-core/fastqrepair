/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FQ_LINT                 } from '../modules/nf-core/fq/lint/main'
include { FASTQC                  } from '../modules/nf-core/fastqc/main'
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { GZRT                    } from '../modules/nf-core/gzrt/main'
include { BBMAP_REPAIR            } from '../modules/nf-core/bbmap/repair/main'
include { FASTQ_REPAIR_WIPERTOOLS } from '../subworkflows/local/fastq_repair_wipertools/main'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_fastqrepair_pipeline'
include { isFastqFileEmpty        } from '../subworkflows/local/utils_nfcore_fastqrepair_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FASTQREPAIR {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    ch_final = channel.empty()      // channel: repaired fastq files

    ch_samples = ch_samplesheet
    if (!params.skip_fq_lint) {
        FQ_LINT(ch_samplesheet)
        ch_versions = ch_versions.mix(FQ_LINT.out.versions.first())
        FQ_LINT.out.lint.map { meta, it -> [meta, file(it).text] }
            .filter {
                meta, text -> text.contains('ERROR') or text.contains('read 0 records') or !text.contains("fq-lint end")
                }
            .set { ch_lint_failed }

        ch_samplesheet
            .join(ch_lint_failed)
            .map {meta, it, text -> [meta, it]}
            .set{ ch_samples }
    }

    // branch .gz and non gz files
    ch_fastq_ext = channel.empty()
    ch_samplesheet
    | branch { _map, fq ->
        gz_files: fq.first().getExtension() == 'gz'
        non_gz_files: true }
    | set { ch_fastq_ext }

    //
    // Recover corrupted gz files
    //
    GZRT (ch_fastq_ext.gz_files)
    ch_versions = ch_versions.mix(GZRT.out.versions.first())

    // Join recovered gz files with non-gz files and filter empty files out
    ch_tobewiped_fastq = channel.empty()
    GZRT.out.recovered
    | concat ( ch_fastq_ext.non_gz_files )
    | filter { meta, fileList -> meta.single_end
        ? fileList instanceof List ? !isFastqFileEmpty(fileList[0]) : !isFastqFileEmpty(fileList)
        : !isFastqFileEmpty(fileList[0]) && !isFastqFileEmpty(fileList[1])}
    | set { ch_tobewiped_fastq }

    // If all input files are empty, then skip the rest of the pipeline
    ch_tobewiped_fastq
    | ifEmpty {
        log.warn "No non-empty FASTQ files found after GZRT. Skipping the rest of the pipeline!"
    }

    // If ch_tobewiped_fastq has a size < ch_samplesheet but > zero, then some files were empty and we need to log.warn them
    ch_samplesheet.map { meta, _fastq -> meta.id }
        .collect()
        .set { all_ids }

    // Extract meta.id from ch_subset and collect into a list
    ch_tobewiped_fastq.map { meta, _fastq -> meta.id }
        .collect()
        .set { subset_ids }

    all_ids.minus(subset_ids).subscribe{ id_list -> if (id_list.size > 0) log.warn "FASTQ files with the following meta.ids are empty: ${id_list.join(', ')}" }

    //
    // Make fastq compliant and wipe bad characters
    //
    ch_repaired_fastq = channel.empty()
    FASTQ_REPAIR_WIPERTOOLS (ch_tobewiped_fastq)
    ch_versions = ch_versions.mix(FASTQ_REPAIR_WIPERTOOLS.out.versions.first())

    FASTQ_REPAIR_WIPERTOOLS.out.wiped_fastq
    | map { meta, fq -> [meta.subMap('sample_id', 'single_end'), fq]}
    | map { meta, fq -> [['id':meta.sample_id, 'single_end':meta.single_end], fq]}
    | branch { item ->
        single_end: item[0].single_end == true
        paired_end: item[0].single_end == false }
    | set { ch_repaired_fastq }

    // Group paired-reads by 'sample_id' and rename keys
    ch_repaired_fastq_paired_end = channel.empty()
    ch_repaired_fastq.paired_end
    | groupTuple
    | set { ch_repaired_fastq_paired_end }

    // Settle reads pairing (re-pair, optional)
    //
    if (!params.skip_bbmap_repair) {
        // Re-pair reads
        BBMAP_REPAIR (ch_repaired_fastq_paired_end, false)
        ch_versions = ch_versions.mix(BBMAP_REPAIR.out.versions.first())

        ch_repaired_fastq_paired_end_singleton = channel.empty()
        BBMAP_REPAIR.out.repaired
        | concat ( BBMAP_REPAIR.out.singleton )
        | groupTuple
        | map { meta, fq -> [meta, fq.flatten()] }
        | set { ch_repaired_fastq_paired_end_singleton }

        ch_final = ch_repaired_fastq_paired_end_singleton.concat(ch_repaired_fastq.single_end)
    } else {
        ch_final = ch_repaired_fastq_paired_end.concat(ch_repaired_fastq.single_end)
    }

    //
    // Assess QC of all fastq files (both single and paired end)
    //
    FASTQC ( ch_final )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'fastqrepair_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'fastqrepair'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
