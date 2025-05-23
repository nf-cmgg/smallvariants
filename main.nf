#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/smallvariants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/smallvariants
----------------------------------------------------------------------------------------
*/

nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { getGenomeAttribute } from './subworkflows/local/utils_cmgg_smallvariants_pipeline'

// Take another look at this later!
params.fasta                = getGenomeAttribute('fasta', params.genomes, params.genome)
params.fai                  = getGenomeAttribute('fai', params.genomes, params.genome)
params.dict                 = getGenomeAttribute('dict', params.genomes, params.genome)
params.elfasta              = getGenomeAttribute('elfasta', params.genomes, params.genome)
params.strtablefile         = getGenomeAttribute('strtablefile', params.genomes, params.genome)
params.sdf                  = getGenomeAttribute('sdf', params.genomes, params.genome)
params.dbsnp                = getGenomeAttribute('dbsnp', params.genomes, params.genome)
params.dbsnp_tbi            = getGenomeAttribute('dbsnp_tbi', params.genomes, params.genome)
params.vep_cache            = getGenomeAttribute('vep_cache', params.genomes, params.genome)
params.dbnsfp               = getGenomeAttribute('dbnsfp', params.genomes, params.genome)
params.dbnsfp_tbi           = getGenomeAttribute('dbnsfp_tbi', params.genomes, params.genome)
params.spliceai_indel       = getGenomeAttribute('spliceai_indel', params.genomes, params.genome)
params.spliceai_indel_tbi   = getGenomeAttribute('spliceai_indel_tbi', params.genomes, params.genome)
params.spliceai_snv         = getGenomeAttribute('spliceai_snv', params.genomes, params.genome)
params.spliceai_snv_tbi     = getGenomeAttribute('spliceai_snv_tbi', params.genomes, params.genome)
params.mastermind           = getGenomeAttribute('mastermind', params.genomes, params.genome)
params.mastermind_tbi       = getGenomeAttribute('mastermind_tbi', params.genomes, params.genome)
params.eog                  = getGenomeAttribute('eog', params.genomes, params.genome)
params.eog_tbi              = getGenomeAttribute('eog_tbi', params.genomes, params.genome)
params.alphamissense        = getGenomeAttribute('alphamissense', params.genomes, params.genome)
params.alphamissense_tbi    = getGenomeAttribute('alphamissense_tbi', params.genomes, params.genome)
params.vcfanno_resources    = getGenomeAttribute('vcfanno_resources', params.genomes, params.genome)
params.vcfanno_config       = getGenomeAttribute('vcfanno_config', params.genomes, params.genome)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SMALLVARIANTS           } from './workflows/smallvariants'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_cmgg_smallvariants_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_cmgg_smallvariants_pipeline'
include { getWorkflowVersion      } from './subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // Check for dependencies between parameters
    //

    if(params.dbsnp_tbi && !params.dbsnp){
        error("Please specify the dbsnp VCF with --dbsnp VCF")
    }

    if (params.annotate) {
        // Check if a genome is given
        if (!params.genome) { error("A genome should be supplied for annotation (use --genome)") }

        // Check if the VEP versions were given
        if (!params.vep_version) { error("A VEP version should be supplied for annotation (use --vep_version)") }
        if (!params.vep_cache_version) { error("A VEP cache version should be supplied for annotation (use --vep_cache_version)") }

        // Check if a species is entered
        if (!params.species) { error("A species should be supplied for annotation (use --species)") }

        // Check if all vcfanno files are supplied when vcfanno should be used
        if (params.vcfanno && (!params.vcfanno_config || !params.vcfanno_resources)) {
            error("A TOML file and resource files should be supplied when using vcfanno (use --vcfanno_config and --vcfanno_resources)")
        }
    }

    def callers = params.callers.tokenize(",")
    callers.each { caller ->
        if(!(caller in GlobalVariables.availableCallers)) { error("\"${caller}\" is not a supported callers please use one or more of these instead: ${GlobalVariables.availableCallers}")}
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONFIG FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def multiqc_logo = params.multiqc_logo   ?: "$projectDir/assets/CMGG_logo.png"

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.input,
        params.ped,
        params.genomes,
        params.genome,
        params.watchdir
    )

    //
    // WORKFLOW: Run main workflow
    //

    SMALLVARIANTS (
        // Input channels
        PIPELINE_INITIALISATION.out.samplesheet,

        // File inputs
        params.fasta,
        params.fai,
        params.dict,
        params.elfasta,
        params.strtablefile,
        params.sdf,
        params.dbsnp,
        params.dbsnp_tbi,
        params.vep_cache,
        params.dbnsfp,
        params.dbnsfp_tbi,
        params.spliceai_indel,
        params.spliceai_indel_tbi,
        params.spliceai_snv,
        params.spliceai_snv_tbi,
        params.mastermind,
        params.mastermind_tbi,
        params.eog,
        params.eog_tbi,
        params.alphamissense,
        params.alphamissense_tbi,
        params.vcfanno_resources,
        params.vcfanno_config,
        params.multiqc_config,
        multiqc_logo,
        params.multiqc_methods_description,
        params.roi,
        params.somalier_sites,
        params.vcfanno_lua,
        params.updio_common_cnvs,
        params.automap_repeats,
        params.automap_panel,
        params.outdir,
        GlobalVariables.pedFiles,
        params.elsites,

        // Boolean inputs
        params.dragstr,
        params.annotate,
        params.vcfanno,
        params.only_call,
        params.only_merge,
        params.filter,
        params.normalize,
        params.add_ped,
        params.gemini,
        params.validate,
        params.updio,
        params.automap,
        params.vep_dbnsfp,
        params.vep_spliceai,
        params.vep_mastermind,
        params.vep_eog,
        params.vep_alphamissense,

        // Value inputs
        params.genome,
        params.species,
        params.vep_cache_version,
        params.vep_chunk_size,
        params.scatter_count,
        params.callers.tokenize(",")
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        SMALLVARIANTS.out.multiqc_report
    )

    publish:
    merged_crams        = SMALLVARIANTS.out.merged_crams
    mosdepth_reports    = SMALLVARIANTS.out.mosdepth_reports
    gvcfs               = SMALLVARIANTS.out.gvcfs.filter { _meta, gvcf, _tbi -> gvcf.startsWith(workflow.workDir) } // Filtering out input GVCFs from the output publishing fixes an issue in the current implementation of the workflow output definitions: https://github.com/nextflow-io/nextflow/issues/5480
    single_beds         = SMALLVARIANTS.out.single_beds
    perbase_beds        = SMALLVARIANTS.out.perbase_beds
    validation          = SMALLVARIANTS.out.validation
    gvcf_reports        = SMALLVARIANTS.out.gvcf_reports
    genomicsdb          = SMALLVARIANTS.out.genomicsdb
    vcfs                = SMALLVARIANTS.out.vcfs
    gemini              = SMALLVARIANTS.out.gemini
    peds                = SMALLVARIANTS.out.peds
    joint_beds          = SMALLVARIANTS.out.joint_beds
    final_reports       = SMALLVARIANTS.out.final_reports
    automap             = SMALLVARIANTS.out.automap
    updio               = SMALLVARIANTS.out.updio
    multiqc             = SMALLVARIANTS.out.multiqc_report
    multiqc_data        = SMALLVARIANTS.out.multiqc_data
}

output {
    merged_crams {
        path { meta, cram, crai ->
            cram >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.cram"
            crai >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.cram.crai"
        }
        enabled !params.skip_merged_cram_output
    }
    mosdepth_reports { path { meta, report ->
        report >> "${meta.family}/${meta.id}_${params.unique_out}/mosdepth/${report.name}"
    } }
    gvcfs { path { meta, gvcf, tbi ->
        gvcf >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.${meta.caller}.g.vcf.gz"
        tbi >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.${meta.caller}.g.vcf.gz.tbi"
    } }
    single_beds { path { meta, bed ->
        bed >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.bed"
    } }
    perbase_beds { path { meta, bed, csi ->
        bed >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.per-base.bed.gz"
        csi >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.per-base.bed.gz.csi"
    } }
    validation { path { meta, report ->
        report >> "${meta.family}/${meta.id}_${params.unique_out}/validation/${meta.caller}/${report.name}"
    } }
    gvcf_reports { path { meta, report ->
        report >> "${meta.family}/${meta.id}_${params.unique_out}/${meta.id}.${meta.caller}.bcftools_stats.txt"
    }}
    genomicsdb {
        enabled (params.output_genomicsdb || params.only_merge)
        path { meta, genomicsdb ->
            genomicsdb >> "${meta.family}/output_${params.unique_out}/${meta.id}_${meta.caller}_genomicsdb"
        }
    }
    vcfs { path { meta, vcf, tbi ->
        vcf >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.vcf.gz"
        tbi >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.vcf.gz.tbi"
    } }
    gemini { path { meta, db ->
        db >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.db"
    } }
    peds { path { meta, ped ->
        ped >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.ped"
    } }
    joint_beds { path { meta, bed ->
        bed >> "${meta.family}/output_${params.unique_out}/${meta.id}.${meta.caller}.bed"
    } }
    final_reports { path { meta, report ->
        report >> "${meta.family}/qc_${params.unique_out}/${report.name}"
    } }
    automap { path { meta, automap ->
        automap >> "${meta.family}/output_${params.unique_out}/automap/${meta.caller}"
    } }
    updio { path { meta, updio ->
        updio >> "${meta.family}/output_${params.unique_out}/updio/${meta.caller}"
    } }
    multiqc { path { report ->
        report[0] >> "${params.unique_out}/multiqc_report.html"
    } }
    multiqc_data { path { folder ->
        folder >> "${params.unique_out}/multiqc_data"
    } }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
