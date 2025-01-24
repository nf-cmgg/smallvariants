process BCFTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.${extension}"), emit: vcf
    tuple val(meta), path("*.tbi")         , emit: tbi, optional: true
    tuple val(meta), path("*.csi")         , emit: csi, optional: true
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def last_args = args3 ?: args2 ?: args

    extension = last_args.contains("--output-type b") || last_args.contains("-Ob") ? "bcf.gz" :
                    last_args.contains("--output-type u") || last_args.contains("-Ou") ? "bcf" :
                    last_args.contains("--output-type z") || last_args.contains("-Oz") ? "vcf.gz" :
                    last_args.contains("--output-type v") || last_args.contains("-Ov") ? "vcf" :
                    "vcf"

    if ("$vcf" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    def filter_2 = args2 ? "| bcftools filter --threads ${task.cpus} ${args2}" : ""
    def filter_3 = args3 ? "| bcftools filter --threads ${task.cpus} ${args3}" : ""

    """
    bcftools filter \\
        --threads ${task.cpus} \\
        $args \\
        $vcf \\
        ${filter_2} \\
        ${filter_3} \\
        --output ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def last_args = args3 ?: args2 ?: args

    extension = last_args.contains("--output-type b") || last_args.contains("-Ob") ? "bcf.gz" :
                    last_args.contains("--output-type u") || last_args.contains("-Ou") ? "bcf" :
                    last_args.contains("--output-type z") || last_args.contains("-Oz") ? "vcf.gz" :
                    last_args.contains("--output-type v") || last_args.contains("-Ov") ? "vcf" :
                    "vcf"
    def index = last_args.contains("--write-index=tbi") || last_args.contains("-W=tbi") ? "tbi" :
                last_args.contains("--write-index=csi") || last_args.contains("-W=csi") ? "csi" :
                last_args.contains("--write-index") || last_args.contains("-W") ? "csi" :
                ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index}" : ""

    if ("$vcf" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
