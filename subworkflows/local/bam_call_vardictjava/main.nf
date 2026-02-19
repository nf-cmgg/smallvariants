include { VARDICTJAVA                       } from '../../../modules/nf-core/vardictjava/main'

include { VCF_CONCAT_BCFTOOLS               } from '../vcf_concat_bcftools/main'
include { VCF_DBSNP_VCFANNO                 } from '../vcf_dbsnp_vcfanno/main'

workflow BAM_CALL_VARDICTJAVA {
    take:
        ch_input             // channel: [mandatory] [ val(meta), path(bam), path(bai), path(bed) ] => sample CRAM files and their indexes
        ch_fasta             // channel: [mandatory] [ val(meta), path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ val(meta), path(fai) ] => fasta reference index
        ch_dbsnp             // channel: [optional]  [ path(vcf) ] => the dbnsp vcf file
        ch_dbsnp_tbi         // channel: [optional]  [ path(tbi) ] => the dbsnp vcf index file

    main:
    VARDICTJAVA(
        ch_input.map { meta, bam, bai, bed ->
            def new_meta = meta + [caller:'vardict']
            [ new_meta, bam, bai, bed ]
        },
        ch_fasta,
        ch_fai
    )

    VCF_CONCAT_BCFTOOLS(
        VARDICTJAVA.out.vcf
    )

    def ch_annotated = channel.empty()
    if(!(ch_dbsnp instanceof List)) {
        VCF_DBSNP_VCFANNO(
            VCF_CONCAT_BCFTOOLS.out.vcfs,
            ch_dbsnp,
            ch_dbsnp_tbi
        )
        ch_annotated = VCF_DBSNP_VCFANNO.out.vcfs
    } else {
        ch_annotated = VCF_CONCAT_BCFTOOLS.out.vcfs
    }

    def ch_vcfs = ch_annotated
        .map { meta, vcf, tbi ->
            def new_meta = meta + [family_samples: meta.sample]
            [ new_meta, vcf, tbi ]
        }

    emit:
    vcfs = ch_vcfs          // channel: [ val(meta), path(vcf), path(tbi) ]
}
