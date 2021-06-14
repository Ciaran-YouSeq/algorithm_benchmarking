process MAKE_BQSR_TABLE {

        input:
                path aligned_bam_f_mates_coord_sort_mrkd_dups
                path known_sites_vcf
                path known_sites_vcf_index
                path reference_genome_fasta
                path reference_genome_dict
                path reference_genome_fai_index
        output:
                path 'BQSR_table'

        """
        gatk BaseRecalibrator \
        --input $aligned_bam_f_mates_coord_sort_mrkd_dups \
        --known-sites $known_sites_vcf \
        --reference $reference_genome_fasta \
        --output BQSR_table
        """
}

process WRITE_BQSR_TABLE {

        input:
                path 'BQSR_table'
        output:
                path 'BQSR.table'

        """
        cat $BQSR_table > BQSR.table
        """
}

