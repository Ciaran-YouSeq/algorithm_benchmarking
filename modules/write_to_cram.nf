process WRITE_TO_CRAM_FILE {

        input:
                path 'aligned_bam_f_mates_coord_sort_mrkd_dups' from marked_duplicates
        output:
                path 'aligned_bam_f_mates_coord_sort_mrkd_dups.cram'

        """
        cat $aligned_bam_f_mates_coord_sort_mrkd_dups > aligned_bam_f_mates_coord_sort_mrkd_dups.cram
        """
}

