process SAMTOOLS_INDEX {

        input:
                path 'aligned_bam_f_mates_coord_sort_mrkd_dups'
        output:
                path 'aligned_bam_f_mates_coord_sort_mrkd_dups.bam.bai'

        """
        samtools index \
	$aligned_bam_f_mates_coord_sort_mrkd_dups \
	'aligned_bam_f_mates_coord_sort_mrkd_dups.bam.bai'

        """
}

