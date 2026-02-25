process MERGE_BEDS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bed, stageAs: "?/*")
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path('*.bed'), emit: bed
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), emit: versions_bedtools, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Stream all .bed files without inflating on disk
    (
        for FILE in */*.bed; do
            cat "\$FILE"
        done
    ) \
    | awk '{print \$1"\\t"\$2"\\t"\$3 }' \
    | bedtools sort -faidx ${fai} -i - \
    | bedtools merge ${args} \
    > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
