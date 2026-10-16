process DORADO_DEMUX {
    tag "${meta.id}"
    label 'process_high'

    container "docker.io/nanoporetech/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(kit_name), path(sample_sheet)

    output:
    tuple val(meta), path("${meta.id}/*.bam"), emit: bam
    tuple val(meta), path("${meta.id}/*.bai"), emit: bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def kitName = "--kit_name ${kit_name}"
    def sampleSheet = sample_sheet ? "--sample_sheet ${sample_sheet}" : ""  // Optional

    """
    dorado demux \\
        -t ${task.cpus} \\
        ${kitName} \\
        ${sampleSheet} \\
        ${reads} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: "\$(dorado --version 2>&1 | head -n1)"
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}/${prefix}.bam
    touch ${prefix}/${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: "\$(dorado --version 2>&1 | head -n1)"
    END_VERSIONS
    """
}
