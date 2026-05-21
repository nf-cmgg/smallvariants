process UPDIO {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    container "quay.io/cmgg/updio:1.0.0"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(cnv)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}"), emit: updio
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('updio'), val("1.0.0"), emit: versions_updio, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def common_cnv_file = cnv ? "--common_cnv_file $cnv" : "--common_cnv_file /usr/local/lib/updio/sample_data/common_dels_1percent_liftover.tsv"
    """
    UPDio \\
        --multisample_vcf $vcf \\
        --output_path $prefix \\
        $common_cnv_file \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix
    touch ${prefix}/${meta.child}.events_list
    touch ${prefix}/${meta.child}.log
    touch ${prefix}/${meta.child}.table
    touch ${prefix}/${meta.child}.upd
    """
}
