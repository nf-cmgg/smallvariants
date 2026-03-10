process PROCESS_BEDS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7df273d12f0c4d8539440b68876edf39b739cb78bb806418c5b5d057fe11bdbd/data' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"


    input:
    tuple val(meta), path(bed), path(roi)

    output:
    tuple val(meta), path('*.bed'), emit: bed
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), emit: versions_bedtools, topic: versions
    tuple val("${task.process}"), val('grep'), eval("grep --version |& sed -n 's/grep (GNU grep) *//p'"), emit: versions_grep, topic: versions
    tuple val("${task.process}"), val('cat'), eval("cat --version |& sed -n 's/cat (GNU coreutils) *//p'"), emit: versions_cat, topic: versions
    tuple val("${task.process}"), val('zcat'), eval("zcat --version |& sed -n 's/zcat (gzip) *//p'"), emit: versions_zcat, topic: versions
    tuple val("${task.process}"), val('sort'), eval("sort --version |& sed -n 's/sort (GNU coreutils) *//p'"), emit: versions_sort, topic: versions

    script:
    // Remove regions with no coverage from the callable regions BED file and intersect with an optional ROI file
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def unzip = bed.extension == "gz" ? "zcat" : "cat"
    def intersect = roi ? "| bedtools intersect -a ${roi} -b - ${args3}" : ""
    """
    ${unzip} ${bed} \\
        | grep ${args} \\
        | sort -k1,1 -k2,2n -S 2G --parallel=${task.cpus} \\
        | bedtools merge ${args2} \\
        ${intersect} \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed
    """
}
