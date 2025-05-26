process PROCESS_BEDS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.0--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(bed), path(roi)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path('*.bed'), emit: bed
    path  "versions.yml"          , emit: versions

    script:
    // Remove regions with no coverage from the callable regions BED file and intersect with an optional ROI file
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def unzip = bed.extension == "gz" ? "zcat" : "cat"
    def intersect = roi ? "| bedtools intersect -a ${roi} -b - ${args3} -g ${fai}" : ""
    """
    ${unzip} ${bed} \\
        | grep ${args} \\
        | bedtools merge ${args2} \\
        ${intersect} \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
