//
// SAMPLE_PREPARATION
//

include { MERGE_BEDS as MERGE_ROI_PARAMS    } from '../../../modules/local/merge_beds'
include { MERGE_BEDS as MERGE_ROI_SAMPLE    } from '../../../modules/local/merge_beds'
include { PROCESS_BEDS                      } from '../../../modules/local/process_beds'

include { SAMTOOLS_MERGE                    } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                    } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_CONVERT                  } from '../../../modules/nf-core/samtools/convert/main'
include { MOSDEPTH                          } from '../../../modules/nf-core/mosdepth/main'

workflow CRAM_PREPARE_SAMTOOLS_BEDTOOLS {
    take:
        ch_crams             // channel: [mandatory] [ val(meta), path(cram), path(crai) ] => sample CRAM files and their optional indices
        ch_roi               // channel: [optional]  [ val(meta), path(roi) ] => ROI bed files for WES analysis
        ch_fasta             // channel: [mandatory] [ path(fasta) ] => fasta reference
        ch_fai               // channel: [mandatory] [ path(fai) ] => fasta reference index
        ch_default_roi       // channel: [optional]  [ path(roi) ] => bed containing regions of interest to be used as default
        output_bam           // boolean: Also output BAM files

    main:

    def ch_versions  = Channel.empty()
    def ch_reports   = Channel.empty()

    //
    // Merge the CRAM files if there are multiple per sample
    //

    def ch_cram_branch = ch_crams
        .map { meta, cram, crai ->
            [ groupKey(meta, meta.duplicate_count), cram, crai]
        }
        .groupTuple()
        .branch { meta, cram, crai ->
            multiple: cram.size() > 1
                return [meta.target, cram]
            single:   cram.size() == 1
                return [meta.target, cram[0], crai[0]]
        }

    SAMTOOLS_MERGE(
        ch_cram_branch.multiple,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    def ch_merged_crams = SAMTOOLS_MERGE.out.cram
        .join(SAMTOOLS_MERGE.out.crai, failOnDuplicate: true, failOnMismatch: true)

    //
    // Index the CRAM files which have no index
    //

    def ch_ready_crams_branch = ch_merged_crams
        .mix(ch_cram_branch.single)
        .branch { meta, cram, crai=[] ->
            not_indexed: crai == []
                return [ meta, cram ]
            indexed: crai != []
                return [ meta, cram, crai ]
        }

    SAMTOOLS_INDEX(
        ch_ready_crams_branch.not_indexed
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    def ch_ready_crams = ch_ready_crams_branch.not_indexed
        .join(SAMTOOLS_INDEX.out.crai, failOnDuplicate: true, failOnMismatch: true)
        .mix(ch_ready_crams_branch.indexed)

    //
    // Optionally convert the CRAM files to BAM
    //

    def ch_ready_bams = Channel.empty()
    if(output_bam) {
        SAMTOOLS_CONVERT(
            ch_ready_crams,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

        ch_ready_bams = SAMTOOLS_CONVERT.out.bam.join(SAMTOOLS_CONVERT.out.bai, failOnDuplicate:true, failOnMismatch:true)
    }

    //
    // Preprocess the ROI BED files => sort and merge overlapping regions
    //

    def ch_roi_branch = ch_roi
        .map { meta, roi ->
            [ groupKey(meta, meta.duplicate_count), roi ]
        }
        .groupTuple()
        .branch { meta, roi ->
            // Determine whether there is an ROI file given to the current sample
            // It's possible that a sample is given multiple times in the samplesheet, in which
            // case they have been merged earlier. This code checks if at least one entry of the same
            // sample contains an ROI file
            def output_roi = roi.findAll { entry -> entry != [] }
            found:      output_roi.size() > 0
                return [ meta.target, output_roi ]
            missing:    output_roi.size() == 0
                return [ meta.target, [] ]
        }

    MERGE_ROI_SAMPLE(
        ch_roi_branch.found,
        ch_fai
    )
    ch_versions = ch_versions.mix(MERGE_ROI_SAMPLE.out.versions.first())

    // Add the default ROI file to all samples without an ROI file
    // if an ROI BED file has been given through the --roi parameter
    def ch_missing_rois = Channel.empty()
    if (ch_default_roi) {
        MERGE_ROI_PARAMS(
            ch_default_roi.map { bed ->
                [[id:"default_roi"], bed]
            },
            ch_fai
        )
        ch_versions = ch_versions.mix(MERGE_ROI_PARAMS.out.versions)

        ch_missing_rois = ch_roi_branch.missing
            .combine(MERGE_ROI_PARAMS.out.bed.map { _meta, bed -> bed })
            .map { meta, _missing, default_roi ->
                [ meta, default_roi ]
            }
    } else {
        ch_missing_rois = ch_roi_branch.missing
    }

    def ch_ready_rois = ch_missing_rois.mix(MERGE_ROI_SAMPLE.out.bed)

    //
    // Create callable regions
    //

    def ch_mosdepth_input = ch_ready_crams
        .join(ch_ready_rois, failOnDuplicate:true, failOnMismatch:true)

    MOSDEPTH(
        ch_mosdepth_input,
        ch_fasta
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    def ch_mosdepth_reports = MOSDEPTH.out.summary_txt
        .mix(
            MOSDEPTH.out.global_txt,
            MOSDEPTH.out.regions_txt
        )
    ch_reports  = ch_reports.mix(
        ch_mosdepth_reports.map { _meta, report -> report }
    )
    def ch_perbase_beds = MOSDEPTH.out.per_base_bed
        .join(MOSDEPTH.out.per_base_csi, failOnMismatch: true, failOnDuplicate:true)

    def ch_beds_to_process = MOSDEPTH.out.quantized_bed
        .join(ch_ready_rois, failOnDuplicate:true, failOnMismatch:true)

    // Filter out the regions with no coverage
    PROCESS_BEDS(
        ch_beds_to_process,
        ch_fai
    )
    ch_versions = ch_versions.mix(PROCESS_BEDS.out.versions)

    def ch_ready_beds = PROCESS_BEDS.out.bed

    emit:
    ready_crams         = ch_ready_crams        // [ val(meta), path(cram), path(crai) ]
    merged_crams        = ch_merged_crams       // [ val(meta), path(cram), path(crai) ]
    ready_bams          = ch_ready_bams         // [ val(meta), path(bam), path(bai) ]
    ready_beds          = ch_ready_beds         // [ val(meta), path(bed) ]
    perbase_beds        = ch_perbase_beds       // [ val(meta), path(bed), path(csi) ]
    mosdepth_reports    = ch_mosdepth_reports   // [ val(meta), path(report) ]
    versions            = ch_versions           // [ path(versions) ]
    reports             = ch_reports            // [ path(reports) ]
}
