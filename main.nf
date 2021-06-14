#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC } from './modules/quality_control'
include { MULTIQC } from './modules/quality_control'
include { TRIMGALORE } from './modules/trimgalore'
//include { ALIGN } from './modules/bwa_mem'
include { ALIGN } from './modules/bowtie_align.nf'
include { FIXMATES } from './modules/fixmates.nf'
include { COORDINATE_SORT } from './modules/coordinate_sort.nf'
include { MARK_DUPLICATES } from './modules/mark_duplicates.nf'
include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
include { WRITE_TO_CRAM_FILE } from './modules/write_to_cram.nf'
include { MAKE_BQSR_TABLE } from './modules/make_bqsr_table.nf'
include { WRITE_BQSR_TABLE } from './modules/make_bqsr_table.nf'
include { ANALYSE_COVARIATES } from './modules/analyse_covariates.nf'
include { HAPLOTYPE_CALLER } from './modules/haplotype_caller.nf'
include { WRITE_VCF } from './modules/haplotype_caller.nf'
include { BGZIP_VCF } from './modules/rtg_tools.nf'
include { INDEX_GZIPPED_VCF } from './modules/rtg_tools.nf'
include { RTG_VCFEVAL } from './modules/rtg_tools.nf'


// sample id
params.id      = ''

//reads
params.reads = "$baseDir/data/*R{1,2}*.gz"

// reference genome prefix
	// --> there are many indecies/ files built from the ref_genome,
	// --- to cut down number of inputs needed in workflow and processes,
	// --- just supply prefix of ref_genome, and find relevant file in process
params.reference_genome_prefix = "$baseDir/data/GCA_000001405.15_GRCh38_full_analysis_set"

// reference genome
params.reference_genome_fasta = "$baseDir/data/GCA_000001405.15_GRCh38_full_analysis_set.fna"

// bowtie2_index
params.bowtie2_index_name = "$baseDir/data/bowtie2_index/GCA_000001405.15_GRCh38_full_analysis_set"
params.bowtie2_index_files = "$baseDir/data/bowtie2_index/"
//
params.known_sites_vcf = "$baseDir/data/Homo_sapiens_assembly38.dbsnp138.vcf"
params.known_sites_vcf_index = "$baseDir/data/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
//
params.reference_genome_dict = "$baseDir/data/GCA_000001405.15_GRCh38_full_analysis_set.dict"
params.reference_genome_fai_index = "$baseDir/data/GCA_000001405.15_GRCh38_full_analysis_set.fna.fai"
//
params.baseline_calls_vcf="$baseDir/data/NA12878.vcf.gz"
params.baseline_calls_vcf_index="$baseDir/data/NA12878.vcf.gz.tbi"
//
params.sdf="$baseDir/data/GCA_000001405.15_GRCh38_full_analysis_set.sdf/"

log.info """\
         C A N C E R   P I P E L I N E
         =============================
         index: ${params.ref}
         reads : ${params.reads}
         """
         .stripIndent()


workflow {
    main:
        // Always pair reads
        read_pairs_ch = channel.fromFilePairs(params.reads)
        FASTQC(read_pairs_ch)
        MULTIQC(FASTQC.out)
        TRIMGALORE(read_pairs_ch)
	ALIGN(TRIMGALORE.out.reads,params.reference_genome_fasta)
        FIXMATES(ALIGN.out)
        COORDINATE_SORT(FIXMATES.out)
        //SAMTOOLS_INDEX(COORDINATE_SORT.out)
	MARK_DUPLICATES(COORDINATE_SORT.out,params.reference_genome_fasta)
	MAKE_BQSR_TABLE(MARK_DUPLICATES.out,params.known_sites_vcf,params.known_sites_vcf_index,params.reference_genome_fasta,params.reference_genome_dict,params.reference_genome_fai_index)
	ANALYSE_COVARIATES(MAKE_BQSR_TABLE.out)
	SAMTOOLS_INDEX(MARK_DUPLICATES.out)
	HAPLOTYPE_CALLER(params.reference_genome_fasta,params.reference_genome_fai_index,params.reference_genome_dict,MARK_DUPLICATES.out,SAMTOOLS_INDEX.out)
	BGZIP_VCF(HAPLOTYPE_CALLER.out)
	INDEX_GZIPPED_VCF(BGZIP_VCF.out)
	RTG_VCFEVAL(params.baseline_calls_vcf,params.baseline_calls_vcf_index,BGZIP_VCF.out,INDEX_GZIPPED_VCF.out,params.sdf)
}

//workflow.onComplete {
//    println "workflow success: ${ workflow.success }"
//}

workflow.onError = {
    println "workflow success: ${ workflow.success }"
}
        
