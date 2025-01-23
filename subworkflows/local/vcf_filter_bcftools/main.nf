//
// Filter the VCFs
//

include { BCFTOOLS_FILTER as FILTER_1 } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as FILTER_2 } from '../../../modules/nf-core/bcftools/filter/main'
include { TABIX_TABIX                 } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_FILTER_BCFTOOLS {
    take:
        ch_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]

    main:

    def ch_versions = Channel.empty()

    FILTER_1(
        ch_vcfs.map { meta, vcf, tbi=[] -> [ meta, vcf, tbi ]}
    )
    ch_versions = ch_versions.mix(FILTER_1.out.versions.first())

    FILTER_2(
        FILTER_1.out.vcf.join(FILTER_1.out.tbi, failOnDuplicate:true, failOnMismatch:true)
    )
    ch_versions = ch_versions.mix(FILTER_2.out.versions.first())

    def ch_filter_vcfs = FILTER_2.out.vcf.join(FILTER_2.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    emit:
    vcfs = ch_filter_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]

    versions = ch_versions // channel: [ versions.yml ]
}
