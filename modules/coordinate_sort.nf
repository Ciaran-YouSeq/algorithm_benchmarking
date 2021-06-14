process COORDINATE_SORT {

        input:
                path 'aligned_bam_f_mates'
        output:
                path 'aligned_bam_f_mates_coord_sort'

        """
        samtools sort \
        -T /data/sort123_ \
        $aligned_bam_f_mates \
        -o aligned_bam_f_mates_coord_sort
        """
}
