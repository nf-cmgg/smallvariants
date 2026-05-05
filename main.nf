#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/smallvariants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/smallvariants
----------------------------------------------------------------------------------------
*/

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
params.maxentscan           = getGenomeAttribute('maxentscan', params.genomes, params.genome)


params.unique_out = "v${workflow.manifest.version.replace('.', '_')}_${new java.util.Date().format( 'yyyy_MM_dd')}"

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

// Temporary type casting until Nextflow supports static typing (26.10+)

// Boolean parameters
params.igenomes_ignore = params.igenomes_ignore != null ? params.igenomes_ignore as Boolean : null
params.dragstr = params.dragstr != null ? params.dragstr as Boolean : null
params.validate = params.validate != null ? params.validate as Boolean : null
params.filter = params.filter != null ? params.filter as Boolean : null
params.annotate = params.annotate != null ? params.annotate as Boolean : null
params.add_ped = params.add_ped != null ? params.add_ped as Boolean : null
params.gemini = params.gemini != null ? params.gemini as Boolean : null
params.mosdepth_slow = params.mosdepth_slow != null ? params.mosdepth_slow as Boolean : null
params.only_call = params.only_call != null ? params.only_call as Boolean : null
params.only_merge = params.only_merge != null ? params.only_merge as Boolean : null
params.output_genomicsdb = params.output_genomicsdb != null ? params.output_genomicsdb as Boolean : null
params.normalize = params.normalize != null ? params.normalize as Boolean : null
params.only_pass = params.only_pass != null ? params.only_pass as Boolean : null
params.keep_alt_contigs = params.keep_alt_contigs != null ? params.keep_alt_contigs as Boolean : null
params.disable_updio = params.disable_updio != null ? params.disable_updio as Boolean : null
params.disable_automap = params.disable_automap != null ? params.disable_automap as Boolean : null
params.hc_phasing = params.hc_phasing != null ? params.hc_phasing as Boolean : null
params.disable_hc_dict_validation = params.disable_hc_dict_validation != null ? params.disable_hc_dict_validation as Boolean : null
params.skip_merged_cram_output = params.skip_merged_cram_output != null ? params.skip_merged_cram_output as Boolean : null
params.version = params.version != null ? params.version as Boolean : null
params.plaintext_email = params.plaintext_email != null ? params.plaintext_email as Boolean : null
params.monochrome_logs = params.monochrome_logs != null ? params.monochrome_logs as Boolean : null
params.validate_params = params.validate_params != null ? params.validate_params as Boolean : null
params.help_full = params.help_full != null ? params.help_full as Boolean : null
params.show_hidden = params.show_hidden != null ? params.show_hidden as Boolean : null
params.vep_merged = params.vep_merged != null ? params.vep_merged as Boolean : null
params.vep_dbnsfp = params.vep_dbnsfp != null ? params.vep_dbnsfp as Boolean : null
params.vep_spliceai = params.vep_spliceai != null ? params.vep_spliceai as Boolean : null
params.vep_spliceregion = params.vep_spliceregion != null ? params.vep_spliceregion as Boolean : null
params.vep_mastermind = params.vep_mastermind != null ? params.vep_mastermind as Boolean : null
params.vep_maxentscan = params.vep_maxentscan != null ? params.vep_maxentscan as Boolean : null
params.vep_eog = params.vep_eog != null ? params.vep_eog as Boolean : null
params.vep_alphamissense = params.vep_alphamissense != null ? params.vep_alphamissense as Boolean : null
params.vcfanno = params.vcfanno != null ? params.vcfanno as Boolean : null

// Integer parameters
params.scatter_count = params.scatter_count != null ? params.scatter_count as Integer : null
params.merge_distance = params.merge_distance != null ? params.merge_distance as Integer : null
params.min_callable_coverage = params.min_callable_coverage != null ? params.min_callable_coverage as Integer : null
params.vep_chunk_size = params.vep_chunk_size != null ? params.vep_chunk_size as Integer : null
params.vep_cache_version = params.vep_cache_version != null ? params.vep_cache_version as Integer : null

// Float parameters
params.vardict_min_af = params.vardict_min_af != null ? params.vardict_min_af as Float : null
params.vep_version = params.vep_version != null ? params.vep_version as Float : null

// Special parameters
params.help = params.help != null ? params.help.toString() == 'false' || params.help.toString() == 'true' ? params.help as Boolean : params.help : null

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

    def List<String> available_callers = ["haplotypecaller", "vardict", "elprep"]

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
        if(!(caller in available_callers)) { error("\"${caller}\" is not a supported callers please use one or more of these instead: ${available_callers.join(', ')}") }
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
        params.help,
        params.help_full,
        params.show_hidden,
        params.unique_out,
    )

    //
    // WORKFLOW: Run main workflow
    //

    SMALLVARIANTS (
        // Input channels
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.ped_files,

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
        params.maxentscan,
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
        params.elsites,
        params.msi_baseline,
        params.updio_regions,

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
        params.disable_updio,
        params.disable_automap,
        params.vep_dbnsfp,
        params.vep_spliceai,
        params.vep_mastermind,
        params.vep_eog,
        params.vep_alphamissense,
        params.vep_maxentscan,

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
        SMALLVARIANTS.out.multiqc_report.toList()
    )

    publish:
    merged_crams        = SMALLVARIANTS.out.merged_crams
    mosdepth_reports    = SMALLVARIANTS.out.mosdepth_reports
    gvcfs               = SMALLVARIANTS.out.gvcfs.filter { _meta, gvcf, _tbi -> gvcf.startsWith(workflow.workDir) } // Filtering out input GVCFs from the output publishing fixes an issue in the current implementation of the workflow output definitions: https://github.com/nextflow-io/nextflow/issues/5480
    msi                 = SMALLVARIANTS.out.msi
    single_beds         = SMALLVARIANTS.out.single_beds
    perbase_beds        = SMALLVARIANTS.out.perbase_beds
    validation          = SMALLVARIANTS.out.validation
    gvcf_reports        = SMALLVARIANTS.out.gvcf_reports
    genomicsdb          = SMALLVARIANTS.out.genomicsdb
    vcfs                = SMALLVARIANTS.out.vcfs.filter { _meta, vcf, _tbi -> vcf.startsWith(workflow.workDir) } // Filtering out input VCFs from the output publishing fixes an issue in the current implementation of the workflow output definitions: https://github.com/nextflow-io/nextflow/issues/5480
    gemini              = SMALLVARIANTS.out.gemini
    peds                = SMALLVARIANTS.out.peds
    joint_beds          = SMALLVARIANTS.out.joint_beds
    final_reports       = SMALLVARIANTS.out.final_reports
    automap             = SMALLVARIANTS.out.automap.map { meta, dir -> [ meta, files("${dir.toUriString()}/**") ] }.transpose(by:1)
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
    msi { path { meta, msi ->
        msi >> "${meta.family}/${meta.id}_${params.unique_out}/msi/${msi.name}"
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
    automap { path { meta, automap_file ->
        automap_file >> "${meta.family}/output_${params.unique_out}/automap/${meta.caller}/${automap_file.name}"
    } }
    updio { path { meta, updio ->
        updio >> "${meta.family}/output_${params.unique_out}/updio/${meta.caller}${params.updio_regions ? '.filtered' : ''}"
    } }
    multiqc { path { _meta, report ->
        report >> "${params.unique_out}/multiqc_report.html"
    } }
    multiqc_data { path { _meta, _folder ->
        "${params.unique_out}/"
    } }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
