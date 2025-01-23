//
// ANNOTATION
//

include { VCFANNO                             } from '../../../modules/nf-core/vcfanno/main'

include { VCF_ANNOTATE_ENSEMBLVEP             } from '../../../subworkflows/local/vcf_annotate_ensemblvep/main'

workflow VCF_ANNOTATION {
    take:
        ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] => The post-processed VCFs
        ch_fasta                // channel: [mandatory] [ val(meta2), path(fasta) ] => fasta reference
        ch_vep_cache            // channel: [optional]  [ path(vep_cache) ] => The VEP cache to use
        ch_vep_extra_files      // channel: [optional]  [ path(file_1, file_2, file_3, ...) ] => All files necessary for using the desired plugins
        ch_vcfanno_config       // channel: [mandatory if params.vcfanno == true] [ path(toml_config_file) ] => The TOML config file for VCFanno
        ch_vcfanno_lua          // channel: [optional  if params.vcfanno == true] [ path(lua_file) ] => A VCFanno Lua file
        ch_vcfanno_resources    // channel: [mandatory if params.vcfanno == true] [ path(resource_dir) ] => The directory containing the reference files for VCFanno
        genome                  // string:  The genome needed for VEP
        species                 // string:  The species to use for VEP
        vep_cache_version       // integer: The version of the VEP cache to use
        vep_chunk_size          // integer: The size of each chunk to split VCFs into
        vcfanno                 // boolean: Use VCFanno for annotation

    main:

    def ch_annotated_vcfs   = Channel.empty()
    def ch_reports          = Channel.empty()
    def ch_versions         = Channel.empty()

    def ch_vep_input = ch_vcfs
        .map { meta, vcf, tbi ->
            [ meta, vcf, tbi, [] ]
        }

    //
    // Do the VEP annotation
    //

    VCF_ANNOTATE_ENSEMBLVEP(
        ch_vep_input,
        ch_fasta,
        genome == "hg38" ? "GRCh38" : genome,
        species,
        vep_cache_version,
        ch_vep_cache,
        ch_vep_extra_files,
        vep_chunk_size
    )

    ch_versions = ch_versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)
    ch_reports  = ch_reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vep_reports)

    //
    // Annotate the VCFs with VCFanno
    //

    if (vcfanno) {

        def ch_vcfanno_input = VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi
            .map { meta, vcf, tbi ->
                [ meta, vcf, tbi, [] ]
            }
            .dump(tag:'vcfanno_input', pretty:true)

        VCFANNO(
            ch_vcfanno_input,
            ch_vcfanno_config,
            ch_vcfanno_lua,
            ch_vcfanno_resources
        )
        ch_versions = ch_versions.mix(VCFANNO.out.versions.first())

        ch_annotated_vcfs = VCFANNO.out.vcf.join(VCFANNO.out.tbi, failOnDuplicate:true, failOnMismatch:true)
    }
    else {
        ch_annotated_vcfs = VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi
    }

    emit:
    annotated_vcfs  = ch_annotated_vcfs   // [ val(meta), path(vcf), path(tbi) ]
    reports         = ch_reports          // [ path(reports) ]
    versions        = ch_versions         // [ path(versions) ]
}
