//
// Subworkflow with functionality specific to the nf-cmgg/smallvariants pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { samplesheetToList         } from 'plugin/nf-schema'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { initializePed             } from 'plugin/nf-ped'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    pedFile           //  string: Path to the common PED file
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message
    unique_out        //  string: A unique name for the output folder to avoid overwriting previous runs

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )


    //
    // Validate parameters and generate parameter summary to stdout
    //

    def command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command
    )


    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    // Get the samplesheet list

    def List samplesheet_list = samplesheetToList(input, "assets/schema_input.json")

    // Pedigree handling
    def pedObject = initializePed()
    if(pedFile) {
        pedObject.importPed(file(pedFile))
    }
    samplesheet_list.each { _meta, _cram, _crai, _gvcf, _gtbi, _vcf, _tbi, _roi, ped, _truth_vcf, _truth_tbi, _truth_bed ->
        if (ped) {
            pedObject.importPed(ped, overwrite:true)
        }
    }
    pedObject.setEntries(pedObject.getEntries().collect { entry ->
        entry.setFamily(normalize_id(entry.getFamily()))
        entry.setIndividual(normalize_id(entry.getIndividual()))
        entry.setMother(normalize_id(entry.getMother()))
        entry.setFather(normalize_id(entry.getFather()))
        return entry
    })

    def List<String> errors = []

    def Map<String, List<String>> families = [:]
    def Map<String, Integer> sample_counts = [:]
    def Map<String, Map<String, Object>> sample_metas = [:]

    // Determine which files to watch for
    def List parsed_samplesheet_list = samplesheetToList(input, "assets/schema_input.json")
        // Do some calculations and manipulations here
        .collect { meta, cram, crai, gvcf, gtbi, vcf, tbi, roi, _ped, truth_vcf, truth_tbi, truth_bed ->
            // Replace dots with underscores in sample and family names to prevent breaking the multiqc report
            def new_meta = meta + [
                id:normalize_id(meta.id),
                sample:normalize_id(meta.sample),
                family:normalize_id(meta.family)
            ]

            // Pipeline logic
            if(!new_meta.family) {
                def Set<String> individualFamilies = pedObject.getFamiliesFromIndividual(new_meta.sample)
                if (individualFamilies.size() > 1) {
                    errors.add("Sample '${new_meta.sample}' is associated with multiple families (${individualFamilies.join(", ")}) in the PED files. Each sample can only belong to a single family.")
                } else {
                    new_meta.family = individualFamilies.size() == 1 ? individualFamilies.first() : new_meta.sample
                }
            }

            if (!families.containsKey(new_meta.family)) {
                families[new_meta.family] = [new_meta.sample]
            } else if(!families[new_meta.family].contains(new_meta.sample)) {
                families[new_meta.family].add(new_meta.sample)
            }


            if (!sample_counts.containsKey(new_meta.sample)) {
                sample_counts[new_meta.sample] = 1
                sample_metas[new_meta.sample] = new_meta
            } else {
                sample_counts[new_meta.sample] += 1
                if(sample_metas[new_meta.sample] != new_meta) {
                    def Map<String, Object> other_meta = sample_metas[new_meta.sample]
                    def List<String> diff_keys = []
                    other_meta.each { k,v ->
                        if (v != new_meta[k]) {
                            diff_keys.add(k)
                        }
                    }
                    errors.add("Found multiple entries for sample '${new_meta.sample}' in the samplesheet with differing meta values (`${diff_keys.join(' ')}`).")
                }
            }

            // Remove the PED file from output
            return [ new_meta, cram, crai, gvcf, gtbi, vcf, tbi, roi, truth_vcf, truth_tbi, truth_bed ]
        }

    // Stop the pipeline if extra validation errors have been detected
    if (errors.size() > 0) {
        error(errors.join("\n"))
    }

    def ch_samplesheet = channel.fromList(parsed_samplesheet_list)
        .map { meta, cram, crai, gvcf, gtbi, vcf, tbi, roi, truth_vcf, truth_tbi, truth_bed ->
            def new_meta = meta + [
                    family_samples:families[meta.family].sort(false).join(","),
                    duplicate_count:sample_counts[meta.id]
                ]
            return [ new_meta, cram, crai, gvcf, gtbi, vcf, tbi, roi, truth_vcf, truth_tbi, truth_bed ]
        }

    // Output the samplesheet
    file(input).copyTo("${outdir}/${unique_out}/samplesheet.${file(input).extension}")

    // Write PED files per family
    def pedFiles = channel.value(pedObject.getFamilies().collectEntries { family ->
        [ family, pedObject.writePed(families: [family]) ]
    })

    emit:
    samplesheet = ch_samplesheet
    ped_files = pedFiles
}

def normalize_id(String id) {
    if (!id) {
        return id
    }
    return id.replaceAll("\\.","_")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML

    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute, genomesMap, genome) {
    if (genomesMap && genome && genomesMap.containsKey(genome)) {
        if (genomesMap[ genome ].containsKey(attribute)) {
            return genomesMap[ genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "..."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
