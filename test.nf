#!/usr/bin/env nextflow
 
params.left_reads = "/alignment_data/trunc1.fq"
params.right_reads = "/alignment_data/trunc2.fq"
params.reference_genome_fasta = "/data/GCA_000001405.15_GRCh38_full_analysis_set.fna"
params.bowtie_index = "/data/bowtie2_index/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"
process ALIGN {
	
	input:
		path 'left_reads' from params.left_reads
		path 'right_reads' from params.right_reads
		val 'bowtie_index' from params.bowtie_index
	output:
		path 'aligned_bam' into aligned

	"""
	bowtie2 --rg-id 'definitelyarealid' --rg 'SM:samplemcsampleid\tLB:libraryname\tPL:ILLUMINA' -p 32 -q -x $bowtie_index -1 $left_reads -2 $right_reads > aligned_bam
	"""
}

process FIXMATES {

	input:
		path 'aligned_bam' from aligned
	output:
		path 'aligned_bam_f_mates' into fixed_mates

	"""
	samtools fixmate -m $aligned_bam aligned_bam_f_mates
	"""
}

process COORDINATE_SORT {

	input:
		path 'aligned_bam_f_mates' from fixed_mates
	output:
		path 'aligned_bam_f_mates_coord_sort' into coord_sorted

	"""
	samtools sort -T /data/sort123_ $aligned_bam_f_mates -o aligned_bam_f_mates_coord_sort
	"""
}

process MARK_DUPLICATES {

	input: 
		path 'aligned_bam_f_mates_coord_sort' from coord_sorted
		path 'reference_genome_fasta' from params.reference_genome_fasta
	output:
		path 'aligned_bam_f_mates_coord_sort_mrkd_dups.cram'
	"""
	samtools markdup --reference $reference_genome_fasta $aligned_bam_f_mates_coord_sort aligned_bam_f_mates_coord_sort_mrkd_dups.cram
	"""
}

