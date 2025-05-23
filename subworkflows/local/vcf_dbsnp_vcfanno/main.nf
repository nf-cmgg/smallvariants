include { VCFANNO            } from '../../../modules/nf-core/vcfanno/main'

workflow VCF_DBSNP_VCFANNO {
    take:
        ch_input             // channel: [mandatory] [ val(meta), path(vcf), path(tbi), ] => VCF files to be annotated
        ch_dbsnp             // channel: [optional]  [ val(meta), path(vcf) ] => the dbnsp vcf file
        ch_dbsnp_tbi         // channel: [optional]  [ val(meta), path(tbi) ] => the dbsnp vcf index file

    main:
    def ch_versions = Channel.empty()

    def ch_vcfanno_toml = ch_dbsnp.map { _meta, dbsnp -> [ get_vcfanno_config(dbsnp) ] }
        .collect()

    def ch_vcfanno_resources = ch_dbsnp.map { _meta, dbsnp -> dbsnp }
        .combine(ch_dbsnp_tbi.map { _meta, tbi -> tbi })
        .collect()

    VCFANNO(
        ch_input.map { meta, vcf, tbi -> [ meta, vcf, tbi, [] ] },
        ch_vcfanno_toml,
        [],
        ch_vcfanno_resources
    )
    ch_versions = ch_versions.mix(VCFANNO.out.versions.first())

    def ch_vcfs = VCFANNO.out.vcf
        .join(VCFANNO.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    emit:
    vcfs = ch_vcfs          // channel: [ val(meta), path(vcf), path(tbi) ]
    versions = ch_versions  // channel: [ path(versions.yml) ]

}

def get_vcfanno_config(vcf) {
    def old_toml = file("${projectDir}/assets/dbsnp.toml", checkIfExists: true)
    old_toml.copyTo("${workDir}/vcfanno/dbsnp.toml")
    def new_toml = file("${workDir}/vcfanno/dbsnp.toml")
    new_toml.text = old_toml.text.replace("DBSNP_FILE", vcf.getName())
    return new_toml
}
