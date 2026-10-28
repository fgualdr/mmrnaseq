//
// This file holds several functions specific to the workflow/rnaseq.nf in the nf-core/rnaseq pipeline
//

import nextflow.Nextflow
import groovy.json.JsonSlurper
import groovy.text.SimpleTemplateEngine

class WorkflowRnaseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        genomeExistsError(params, log)


        if (!params.fasta) {
            Nextflow.error("Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file.")
        }

        if (!params.gtf && !params.gff) {
            Nextflow.error("No GTF or GFF3 annotation specified! The pipeline requires at least one of these files.")
        }

        if (params.gtf) {
            if (params.gff) {
                gtfGffWarn(log)
            }
            if (params.genome == 'GRCh38' && params.gtf.contains('Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf')) {
                ncbiGenomeWarn(log)
            }
            if (params.gtf.contains('/UCSC/') && params.gtf.contains('Annotation/Genes/genes.gtf')) {
                ucscGenomeWarn(log)
            }
        }

        if (params.transcript_fasta) {
            transcriptsFastaWarn(log)
        }

        if (!params.skip_bbsplit && !params.bbsplit_index && !params.bbsplit_fasta_list) {
            Nextflow.error("Please provide either --bbsplit_fasta_list / --bbsplit_index to run BBSplit.")
        }

        if (params.remove_ribo_rna && !params.ribo_database_manifest) {
            Nextflow.error("Please provide --ribo_database_manifest to remove ribosomal RNA with SortMeRNA.")
        }


        if (params.with_umi && !params.skip_umi_extract) {
            if (!params.umitools_bc_pattern && !params.umitools_bc_pattern2) {
                Nextflow.error("UMI-tools requires a barcode pattern to extract barcodes from the reads.")
            }
        }

        if (!params.skip_trimming) {
            if (!valid_params['trimmers'].contains(params.trimmer)) {
                Nextflow.error("Invalid option: '${params.trimmer}'. Valid options for '--trimmer': ${valid_params['trimmers'].join(', ')}.")
            }
        }

        if (!params.skip_alignment) {
            if (!valid_params['aligners'].contains(params.aligner)) {
                Nextflow.error("Invalid option: '${params.aligner}'. Valid options for '--aligner': ${valid_params['aligners'].join(', ')}.")
            }
        } else {
            Nextflow.error("--skip_alignment is not supported by mmrnaseq because genomic multimapper resolution requires a STAR BAM.")
        }

        // Warn if --additional_fasta is combined with a pre-built STAR index.
        if (!params.skip_alignment && params.additional_fasta && params.star_index) {
            additionaFastaIndexWarn('star', log)
        }

        // Check which RSeQC modules we are running
        def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
        if ((valid_params['rseqc_modules'] + rseqc_modules).unique().size() != valid_params['rseqc_modules'].size()) {
            Nextflow.error("Invalid option: ${params.rseqc_modules}. Valid options for '--rseqc_modules': ${valid_params['rseqc_modules'].join(', ')}")
        }
    }

    //
    // Function to check whether biotype field exists in GTF file
    //
    public static Boolean biotypeInGtf(gtf_file, biotype, log) {
        def hits = 0
        gtf_file.eachLine { line ->
            def attributes = line.split('\t')[-1].split()
            if (attributes.contains(biotype)) {
                hits += 1
            }
        }
        if (hits) {
            return true
        } else {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" +
                "  Biotype QC will be skipped to circumvent the issue below:\n" +
                "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
                "  Amend '--featurecounts_group_type' to change this behaviour.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            return false
        }
    }

    //
    // Function to generate an error if contigs in genome fasta file > 512 Mbp
    //
    public static void checkMaxContigSize(fai_file, log) {
        def max_size = 512000000
        fai_file.eachLine { line ->
            def lspl  = line.split('\t')
            def chrom = lspl[0]
            def size  = lspl[1]
            if (size.toInteger() > max_size) {
                def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  Contig longer than ${max_size}bp found in reference genome!\n\n" +
                    "  ${chrom}: ${size}\n\n" +
                    "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n" +
                    "  Please see:\n" +
                    "  https://github.com/nf-core/rnaseq/issues/744\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.error(error_string)
            }
        }
    }

    //
    // Function that parses Salmon quant 'meta_info.json' output file to get inferred strandedness
    //
    public static String getSalmonInferredStrandedness(json_file) {
        def lib_type = new JsonSlurper().parseText(json_file.text).get('library_types')[0]
        def strandedness = 'reverse'
        if (lib_type) {
            if (lib_type in ['U', 'IU']) {
                strandedness = 'unstranded'
            } else if (lib_type in ['SF', 'ISF']) {
                strandedness = 'forward'
            } else if (lib_type in ['SR', 'ISR']) {
                strandedness = 'reverse'
            }
        }
        return strandedness
    }

    //
    // Function that parses TrimGalore log output file to get total number of reads after trimming
    //
    public static Integer getTrimGaloreReadsAfterFiltering(log_file) {
        def total_reads = 0
        def filtered_reads = 0
        log_file.eachLine { line ->
            def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
            def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
            if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
            if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
        }
        return total_reads - filtered_reads
    }

    //
    // Parse a two-column metric/value TSV emitted by the multimapper stages.
    //
    public static Map readMetricValueTsv(qc_file) {
        def metrics = [:]
        def line_no = 0
        qc_file.eachLine { line ->
            line_no += 1
            if (line_no > 1 && line.trim()) {
                def fields = line.split('\t', -1)
                if (fields.size() >= 2) {
                    metrics[fields[0]] = fields[1]
                }
            }
        }
        return metrics
    }

    //
    // STAR alignment availability. These are aligner-level categories only and are
    // deliberately kept separate from the final quantitative-usable fraction.
    // STAR reports paired-end input as read pairs/fragments, matching our fragment
    // semantics for the downstream multimapper contract.
    //
    public static Map getStarAlignmentMetrics(align_log) {
        def out = [
            input_fragments: 0L,
            unique_pct: 0.0,
            multiple_loci_pct: 0.0,
            too_many_loci_pct: 0.0
        ]
        align_log.eachLine { line ->
            def ni = line =~ /Number of input reads\s*\|\s*([0-9]+)/
            def up = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
            def mp = line =~ /% of reads mapped to multiple loci\s*\|\s*([\d\.]+)%/
            def tp = line =~ /% of reads mapped to too many loci\s*\|\s*([\d\.]+)%/
            if (ni) out.input_fragments = ni[0][1].toLong()
            if (up) out.unique_pct = up[0][1].toDouble()
            if (mp) out.multiple_loci_pct = mp[0][1].toDouble()
            if (tp) out.too_many_loci_pct = tp[0][1].toDouble()
        }
        if (out.input_fragments <= 0) {
            Nextflow.error("Could not parse 'Number of input reads' from STAR Log.final.out: ${align_log}")
        }
        return out
    }

    //
    // Quantitative mapping acceptance after multimapper handling.
    //
    // --skip_em false: usable = trusted unique + globally FDR-resolved multimappers
    // --skip_em true : usable = trusted unique only
    //
    // The denominator is STAR input fragments, not BAM records. This makes the
    // threshold compatible with PE fragment semantics and with the actual BAM used
    // for downstream counting.
    //
    public static ArrayList getQuantitativeFragmentUsability(params, align_log, contract_qc, realization_qc) {
        def star = getStarAlignmentMetrics(align_log)
        def contract = readMetricValueTsv(contract_qc)
        def realization = readMetricValueTsv(realization_qc)

        for (key in ['total_fragments','unique_fragments','multimapper_fragments','excluded_fragments']) {
            if (!contract.containsKey(key)) {
                Nextflow.error("Missing '${key}' in multimapper contract QC: ${contract_qc}")
            }
        }
        for (key in ['resolved_fdr_fragments','unresolved_fragments']) {
            if (!realization.containsKey(key)) {
                Nextflow.error("Missing '${key}' in BAM realization QC: ${realization_qc}")
            }
        }

        long input_fragments = star.input_fragments as long
        long contract_total  = contract.total_fragments.toDouble().round() as long
        long trusted_unique  = contract.unique_fragments.toDouble().round() as long
        long eligible_multi  = contract.multimapper_fragments.toDouble().round() as long
        long excluded        = contract.excluded_fragments.toDouble().round() as long
        long resolved        = realization.resolved_fdr_fragments.toDouble().round() as long
        long unresolved      = realization.unresolved_fragments.toDouble().round() as long

        if (contract_total > input_fragments) {
            Nextflow.error("Fragment accounting error for ${contract_qc}: mapped contract fragments (${contract_total}) exceed STAR input fragments (${input_fragments}).")
        }
        if (trusted_unique + eligible_multi + excluded != contract_total) {
            Nextflow.error("Fragment accounting error for ${contract_qc}: unique + multimapper + excluded != contract total.")
        }
        if (resolved + unresolved != eligible_multi) {
            Nextflow.error("Fragment accounting error for ${realization_qc}: resolved + unresolved != eligible multimappers.")
        }
        if (params.skip_em && resolved != 0) {
            Nextflow.error("--skip_em true must not resolve multimappers, but ${realization_qc} reports ${resolved} resolved fragments.")
        }

        long quantitative_usable = trusted_unique + resolved
        def pct = { long n -> input_fragments > 0 ? (100.0d * n / input_fragments) : 0.0d }
        double quantitative_usable_pct = pct(quantitative_usable)
        boolean pass = quantitative_usable_pct >= params.min_mapped_reads.toDouble()
        String mode = params.skip_em ? 'UNIQUE_ONLY' : 'BAYES_FDR'

        return [
            input_fragments,
            star.unique_pct,
            star.multiple_loci_pct,
            star.too_many_loci_pct,
            trusted_unique,
            eligible_multi,
            resolved,
            unresolved,
            excluded,
            quantitative_usable,
            pct(trusted_unique),
            pct(eligible_multi),
            pct(resolved),
            pct(unresolved),
            pct(excluded),
            quantitative_usable_pct,
            mode,
            pass
        ]
    }

    //
    // Function that parses and returns the predicted strandedness from the RSeQC infer_experiment.py output
    //
    public static ArrayList getInferexperimentStrandedness(inferexperiment_file, cutoff=30) {
        def sense        = 0
        def antisense    = 0
        def undetermined = 0
        inferexperiment_file.eachLine { line ->
            def undetermined_matcher = line =~ /Fraction of reads failed to determine:\s([\d\.]+)/
            def se_sense_matcher     = line =~ /Fraction of reads explained by "\++,--":\s([\d\.]+)/
            def se_antisense_matcher = line =~ /Fraction of reads explained by "\+-,-\+":\s([\d\.]+)/
            def pe_sense_matcher     = line =~ /Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)/
            def pe_antisense_matcher = line =~ /Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)/
            if (undetermined_matcher) undetermined = undetermined_matcher[0][1].toFloat() * 100
            if (se_sense_matcher)     sense        = se_sense_matcher[0][1].toFloat() * 100
            if (se_antisense_matcher) antisense    = se_antisense_matcher[0][1].toFloat() * 100
            if (pe_sense_matcher)     sense        = pe_sense_matcher[0][1].toFloat() * 100
            if (pe_antisense_matcher) antisense    = pe_antisense_matcher[0][1].toFloat() * 100
        }
        def strandedness = 'unstranded'
        if (sense >= 100-cutoff) {
            strandedness = 'forward'
        } else if (antisense >= 100-cutoff) {
            strandedness = 'reverse'
        }
        return [ strandedness, sense, antisense, undetermined ]
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Create MultiQC tsv custom content from a list of values
    //
    public static String multiqcTsvFromList(tsv_data, header) {
        def tsv_string = ""
        if (tsv_data.size() > 0) {
            tsv_string += "${header.join('\t')}\n"
            tsv_string += tsv_data.join('\n')
        }
        return tsv_string
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

    //
    // Print a warning if using GRCh38 assembly from igenomes.config
    //
    private static void ncbiGenomeWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if using a UCSC assembly from igenomes.config
    //
    private static void ucscGenomeWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if both GTF and GFF have been provided
    //
    private static void gtfGffWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Both '--gtf' and '--gff' parameters have been provided.\n" +
            "  Using GTF file as priority.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if using '--transcript_fasta'
    //
    private static void transcriptsFastaWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  '--transcript_fasta' parameter has been provided.\n" +
            "  Make sure transcript names in this file match those in the GFF/GTF file.\n\n" +
            "  Please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/753\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if --skip_alignment has been provided
    //
    private static void skipAlignmentWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  '--skip_alignment' parameter has been provided.\n" +
            "  Skipping alignment, genome-based quantification and all downstream QC processes.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if using '--additional_fasta' and '--<ALIGNER>_index'
    //
    private static void additionaFastaIndexWarn(index, log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  When using '--additional_fasta <FASTA_FILE>' the aligner index will not\n" +
            "  be re-built with the transgenes incorporated by default since you have \n" +
            "  already provided an index via '--${index}_index <INDEX>'.\n\n" +
            "  Set '--additional_fasta <FASTA_FILE> --${index}_index false --save_reference' to\n" +
            "  re-build the index with transgenes included and the index will be saved in\n" +
            "  'results/genome/index/${index}/' for re-use with '--${index}_index'.\n\n" +
            "  Ignore this warning if you know that the index already contains transgenes.\n\n" +
            "  Please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/556\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }
}
