process AUTOMAP_REPEATS {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    container "quay.io/cmgg/automap:1.0.0"

    input:
    tuple val(meta), val(genome)

    output:
    tuple val(meta), path("*.bed")  , emit: repeats
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('automap'), val("1.0.0"), emit: versions_automap, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Files are present in the container
    if (genome == "hg38" || genome == "GRCh38") {
        """
        zcat /usr/local/lib/automap/Resources/repeats_hg38.part*.bed.gz > ${prefix}.bed
        """
    } else if (genome == "hg19" || genome == "GRCh37") {
        """
        zcat /usr/local/lib/automap/Resources/repeats.part*.bed.gz > ${prefix}.bed
        """
    } else {
        error("No repeat regions can be found for genome '${genome}'. Available genomes are hg38/GRCh38 or hg19/GRCh37")
    }


    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (["hg38", "GRCh38", "hg19", "GRCh37"].contains(genome)) {
        """
        touch ${prefix}.bed
        """
    } else {
        error("No repeat regions can be found for genome '${genome}'. Available genomes are hg38/GRCh38 or hg19/GRCh37")
    }
}
