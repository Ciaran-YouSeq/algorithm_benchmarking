process HAPLOTYPE_CALLER {

        input:
                path reference_genome_fasta
                path reference_genome_fai_index
                path reference_genome_dict
                path aligned_bam_f_mates_coord_sort_mrkd_dups
		path 'aligned_bam_f_mates_coord_sort_mrkd_dups.bai'
        output:
                path 'haplotype_caller_output'

        """
        gatk HaplotypeCaller \
        --input $aligned_bam_f_mates_coord_sort_mrkd_dups \
        --reference $reference_genome_fasta \
        --output haplotype_caller_output
        """
}

process WRITE_VCF {

        input:
                path 'haplotype_caller_output'
        output:
                path 'haplotype_caller_output.vcf'
        """
        cat $haplotype_caller_output > haplotype_caller_output.vcf
        """
}
