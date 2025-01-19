process COLLECTRESULTS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(repaired_fastq)

    output:
    tuple val(meta), path("*.gz")       , emit: renamed_fastq
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # base_name=\$(basename "$repaired_fastq" .fastq.gz)
    # new_file="\$(echo "\$base_name" | sed -E 's/_recovered|_wiped|_repaired|_gather//g')_repaired.fastq.gz"
    # mv "$repaired_fastq" "\$new_file"

    # Extract the base name without path or extension
    base_name=\$(basename "$repaired_fastq" .fastq.gz)

    # Check if _1_ or _2_ is in the filename
    if [[ "\$base_name" == *_1_* ]]; then
        new_file="${prefix}_1.fastq.gz"
    elif [[ "\$base_name" == *_2_* ]]; then
        new_file="${prefix}_2.fastq.gz"
    else
        new_file="${prefix}.fastq.gz"
    fi

    # Rename the file
    mv "$repaired_fastq" "\$new_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collectresults: 1.0.0
    END_VERSIONS
    """
}
