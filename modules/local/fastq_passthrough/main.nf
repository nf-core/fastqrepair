process FASTQ_PASSTHROUGH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'quay.io/nf-core/ubuntu:22.04' }"

    input:
    // staged in a sub-folder so that an input already named like the output cannot collide with it
    tuple val(meta), path(fastq, arity: '1..2', stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}*.fastq.gz"), emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    // Same names as the repaired output: <id>.fastq.gz (single-end), <id>_1.fastq.gz and <id>_2.fastq.gz (paired-end)
    def reads = fastq instanceof List ? fastq : [fastq]
    def names = reads.size() == 1 ? ["${meta.id}.fastq.gz"] : ["${meta.id}_1.fastq.gz", "${meta.id}_2.fastq.gz"]
    def cmds  = [reads, names].transpose().collect { read, name ->
        // gzip input is copied as it is, plain text input is compressed
        read.name.endsWith('.gz') ? "cp ${read} ${name}" : "gzip -c ${read} > ${name}"
    }.join('\n    ')
    """
    ${cmds}
    """

    stub:
    def reads = fastq instanceof List ? fastq : [fastq]
    def names = reads.size() == 1 ? ["${meta.id}.fastq.gz"] : ["${meta.id}_1.fastq.gz", "${meta.id}_2.fastq.gz"]
    """
    touch ${names.join(' ')}
    """
}
