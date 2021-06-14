process FIXMATES {

        input:
                path 'aligned_bam'
        output:
                path 'aligned_bam_f_mates'

        """
        samtools fixmate \
        -m $aligned_bam \
        aligned_bam_f_mates
        """
}

