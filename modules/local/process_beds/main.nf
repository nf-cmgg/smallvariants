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
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), emit: versions_bedtools, topic: versions
    tuple val("${task.process}"), val('grep'), eval("grep --version |& sed -n 's/grep (GNU grep) *//p'"), emit: versions_grep, topic: versions
    tuple val("${task.process}"), val('cat'), eval(" cat --version |& sed -n 's/cat (GNU coreutils) *//p'"), emit: versions_cat, topic: versions
    tuple val("${task.process}"), val('zcat'), eval("zcat --version |& sed -n 's/zcat (gzip) *//p'"), emit: versions_zcat, topic: versions

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
        | bedtools sort -faidx ${fai} \\
        ${intersect} \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed
    """
}
