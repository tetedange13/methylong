process DORADO_DEMUX {
    tag "${meta.id}"
    label 'process_high'

    container "docker.io/nanoporetech/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f"

    input:
    tuple val(meta), path(reads), path(reference)
    tuple val(meta2), val(kit_name), path(sample_sheet)

    output:
    tuple val(meta), path("*.bam"), path(reference), emit: bam
    tuple val(meta), path("out_demux/sequencing_summary.txt")         , emit: summary, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def kitName = "--kit-name ${kit_name}"
    def sampleSheet = sample_sheet ? "--sample-sheet ${sample_sheet}" : ""  // Optional

    """
    dorado demux \\
        -t ${task.cpus} \\
        ${kitName} \\
        ${sampleSheet} \\
        ${reads} \\
        ${args} \\
        --output-dir ./out_demux
    # Get all 'barcodeNN' subdirs (WARN: With sample_sheet, they can be aliased)
    # FIXME: Should excluded unclassifed
    # Exit with error if only 'unclassifed' ???
    for file in \$(find ./out_demux -type f -name "*.bam");
    do
        echo \$file
        new_name=\$(basename \$(dirname \$file)).bam
        mv \$file ${meta.id}_\$new_name
    done

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
