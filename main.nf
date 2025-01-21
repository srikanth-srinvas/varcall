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
include { BWA_MEM } from './modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { FREEBAYES } from './modules/nf-core/freebayes/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { SNPEFF_SNPEFF } from './modules/nf-core/snpeff/snpeff/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_varcall_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_varcall_pipeline'
include { getGenomeAttribute } from './subworkflows/local/utils_nfcore_varcall_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta = getGenomeAttribute('fasta')
params.bwa_index = getGenomeAttribute('bwa_index')
params.sort_bam = true // Added parameter for BAM sorting

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARCALL {
    take:
    reads

    main:
    // Quality Control
    FASTQC(reads)

    // Adapter Trimming
    CUTADAPT(reads)

    // Alignment
    BWA_MEM(
        reads: CUTADAPT.out,
        index: params.bwa_index,
        fasta: params.fasta,
        sort_bam: params.sort_bam
    )

    // Continue the rest of the workflow
    SAMTOOLS_SORT(BWA_MEM.out.bam)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out)
    FREEBAYES(SAMTOOLS_INDEX.out)
    SNPEFF_SNPEFF(FREEBAYES.out)
    MULTIQC(
        FASTQC.out,
        CUTADAPT.out,
        BWA_MEM.out,
        SAMTOOLS_SORT.out,
        SAMTOOLS_INDEX.out,
        FREEBAYES.out,
        SNPEFF_SNPEFF.out
    )

    emit:
    vcf = FREEBAYES.out.vcf
    annotated_vcf = SNPEFF_SNPEFF.out.annotated_vcf
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