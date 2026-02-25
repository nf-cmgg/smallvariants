process MERGE_BEDS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7df273d12f0c4d8539440b68876edf39b739cb78bb806418c5b5d057fe11bdbd/data' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"

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
    | sort -k1,1 -k2,2n -S 2G --parallel=${task.cpus} \
    | bedtools merge ${args} \
    > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
