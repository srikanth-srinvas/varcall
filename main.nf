#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/varcall
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/varcall
    Website: https://nf-co.re/varcall
    Slack  : https://nfcore.slack.com/channels/varcall
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC } from './modules/nf-core/fastqc/main'
include { CUTADAPT } from './modules/nf-core/cutadapt/main'
include { BWA_MEM } from './modules/nf-core/bwa/mem'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index'
include { FREEBAYES } from './modules/nf-core/freebayes/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { SNPEFF } from './modules/nf-core/snpeff/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_varcall_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_varcall_pipeline'
include { getGenomeAttribute } from './subworkflows/local/utils_nfcore_varcall_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARCALL {
    take:
    reads

    main:
    // Run FastQC
    FASTQC (reads)

    // Run Cutadapt
    CUTADAPT (FASTQC.out)

    // Run BWA
    BWA_MEM (CUTADAPT.out)

    // Run SAMtools Sort
    SAMTOOLS_SORT (BWA_MEM.out)

    // Run SAMtools Index
    SAMTOOLS_INDEX (SAMTOOLS_SORT.out)

    // Run FreeBayes
    FREEBAYES (SAMTOOLS_INDEX.out)

    // Run SnpEff
    SNPEFF (FREEBAYES.out)

    // Run MultiQC
    MULTIQC (FASTQC.out, CUTADAPT.out, BWA_MEM.out, SAMTOOLS_SORT.out, SAMTOOLS_INDEX.out, FREEBAYES.out, SNPEFF.out)

    emit:
    vcf = FREEBAYES.out.vcf
    annotated_vcf = SNPEFF.out.annotated_vcf
    multiqc_report = MULTIQC.out.report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    // SUBWORKFLOW: Run initialisation tasks
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    // WORKFLOW: Run main workflow
    VARCALL (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    // SUBWORKFLOW: Run completion tasks
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        VARCALL.out.vcf,
        VARCALL.out.annotated_vcf,
        VARCALL.out.multiqc_report
    )
}