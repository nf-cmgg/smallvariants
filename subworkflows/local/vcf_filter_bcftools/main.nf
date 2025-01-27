//
// Filter the VCFs
//

include { BCFTOOLS_FILTER } from '../../../modules/local/bcftools/filter/main'

workflow VCF_FILTER_BCFTOOLS {
    take:
        ch_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]

    main:

    def ch_versions = Channel.empty()

    BCFTOOLS_FILTER(
        ch_vcfs
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    def ch_filter_vcfs = BCFTOOLS_FILTER.out.vcf
        .join(BCFTOOLS_FILTER.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    emit:
    vcfs = ch_filter_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]

    versions = ch_versions // channel: [ versions.yml ]
}
