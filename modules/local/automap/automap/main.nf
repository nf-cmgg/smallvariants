process AUTOMAP_AUTOMAP {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    container "cmgg/automap:1.0.0"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(repeats)
    tuple val(meta3), path(panel)
    val(genome)

    output:
    tuple val(meta), path("${prefix}"), emit: automap
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('automap'), val("1.0.0"), emit: versions_automap, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def panel_file = panel ? "--panel $panel" : "--panel /usr/local/lib/automap/Resources/Biomodule_20220808_all_genes_hg38.txt"
    def hg_genome = genome ?: "hg38"

    """
    automap \\
        --vcf $vcf \\
        --genome $hg_genome \\
        --out $prefix/ \\
        --repeats $repeats \\
        $panel_file \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def panel_name = args.contains("--panelname") ? args.split("--panelname")[-1].trim().split(" ")[0] : ""
    prefix = task.ext.prefix ?: "${meta.id}"

    def create_outputs = meta.family_samples.tokenize(",").size() > 1 ? (1..meta.family_samples.tokenize(",").size()).collect { number ->
        def cmd_prefix = "touch ${prefix}/sample${number}"
        [
            "mkdir ${prefix}/sample${number}",
            "${cmd_prefix}/sample${number}.HomRegions.pdf",
            "${cmd_prefix}/sample${number}.HomRegions.${panel_name}.tsv",
            "${cmd_prefix}/sample${number}.HomRegions.tsv",
            "${cmd_prefix}/sample${number}.HomRegions.strict.${panel_name}.tsv"
        ].join(" && ")
    }.join(" && ") : [
        "touch ${prefix}/${meta.id}.HomRegions.pdf",
        "touch ${prefix}/${meta.id}.HomRegions.${panel_name}.tsv",
        "touch ${prefix}/${meta.id}.HomRegions.tsv",
        "touch ${prefix}/${meta.id}.HomRegions.strict.${panel_name}.tsv"
    ].join(" && ")

    """
    mkdir $prefix
    ${create_outputs}
    """
}
