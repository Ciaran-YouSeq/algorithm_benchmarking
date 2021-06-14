process MARK_DUPLICATES {

        input:
                path 'aligned_bam_f_mates_coord_sort'
                path 'reference_genome_fasta'
        output:
                path 'aligned_bam_f_mates_coord_sort_mrkd_dups'

        """
        samtools markdup \
        --reference $reference_genome_fasta \
        $aligned_bam_f_mates_coord_sort \
        aligned_bam_f_mates_coord_sort_mrkd_dups
        """
}
