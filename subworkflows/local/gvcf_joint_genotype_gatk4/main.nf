//
// GENOTYPE
//

include { MERGE_BEDS             } from '../../../modules/local/merge_beds'

include { GAWK                   } from '../../../modules/nf-core/gawk/main'
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { BEDTOOLS_SPLIT         } from '../../../modules/nf-core/bedtools/split/main'
include { BCFTOOLS_QUERY         } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_STATS         } from '../../../modules/nf-core/bcftools/stats/main'

include { INPUT_SPLIT_BEDTOOLS   } from '../input_split_bedtools/main'
include { VCF_CONCAT_BCFTOOLS    } from '../vcf_concat_bcftools/main'

workflow GVCF_JOINT_GENOTYPE_GATK4 {
    take:
        ch_gvcfs        // channel: [mandatory] [ val(meta), path(gvcf), path(tbi) ] => The GVCFs
        ch_fasta        // channel: [mandatory] [ path(fasta) ] => fasta reference
        ch_fai          // channel: [mandatory] [ path(fai) ] => fasta reference index
        ch_dict         // channel: [mandatory] [ path(dict) ] => sequence dictionary
        ch_dbsnp        // channel: [optional]  [ path(dbsnp) ] => The VCF containing the dbsnp variants
        ch_dbsnp_tbi    // channel: [optional]  [ path(dbsnp_tbi) ] => The index of the dbsnp VCF
        scatter_count   // integer: The amount of times each file should be scattered

    main:

    def ch_vcfs     = channel.empty()

    //
    // Get a BED file containing all contigs
    //

    GAWK(
        ch_fai,
        [],
        false
    )

    BCFTOOLS_QUERY(
        ch_gvcfs,
        [],
        [],
        []
    )

    def ch_merge_beds_input = BCFTOOLS_QUERY.out.output
        .map { meta, bed ->
            // Create the family meta
            def new_meta = meta.subMap(["family", "family_samples", "caller"]) + [id:meta.family]
            [ groupKey(new_meta, meta.family_samples.tokenize(",").size()), bed ]
        }
        .groupTuple()

    MERGE_BEDS(
        ch_merge_beds_input,
        ch_fai
    )
    def ch_beds = MERGE_BEDS.out.bed

    //
    // Split BED file into multiple BEDs specified by --scatter_count
    //

    INPUT_SPLIT_BEDTOOLS(
        MERGE_BEDS.out.bed.map { meta, bed ->
            // Multiply the scatter count by the family size to better scatter big families
            [meta, bed, (scatter_count * meta.family_samples.tokenize(",").size())]
        },
        ch_gvcfs.map { meta, gvcf, tbi ->
            // Create the family meta
            def new_meta = meta.subMap(["family", "family_samples", "caller"]) + [id:meta.family]
            [ groupKey(new_meta, new_meta.family_samples.tokenize(",").size()), gvcf, tbi ]
        }.groupTuple()
    )

    //
    // Create GenomicDBs for each family for each BED file
    //

    def ch_genomicsdbimport_input = INPUT_SPLIT_BEDTOOLS.out.split
        .map { meta, gvcfs, tbis, bed ->
            [ meta, gvcfs, tbis, bed, [], [] ]
        }

    GATK4_GENOMICSDBIMPORT(
        ch_genomicsdbimport_input,
        false,
        false,
        false
    )

    def ch_genotypegvcfs_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb
        .join(INPUT_SPLIT_BEDTOOLS.out.split, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, genomicsdb, _gvcf, _tbi, bed ->
            [ meta, genomicsdb, [], bed, [] ]
        }
    //
    // Genotype the genomicsDBs
    //

    GATK4_GENOTYPEGVCFS(
        ch_genotypegvcfs_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )

    def ch_gather_inputs = GATK4_GENOTYPEGVCFS.out.vcf
        .join(GATK4_GENOTYPEGVCFS.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    //
    // Combine the genotyped VCFs from each family back together
    //

    VCF_CONCAT_BCFTOOLS(
        ch_gather_inputs
    )
    ch_vcfs = VCF_CONCAT_BCFTOOLS.out.vcfs

    emit:
    vcfs = ch_vcfs                                      // [ val(meta), path(vcf), path(tbi) ]
    genomicsdb = GATK4_GENOMICSDBIMPORT.out.genomicsdb  // [ val(meta), path(genomicsdb) ]
    beds = ch_beds                                      // [ val(meta), path(bed) ]
}
