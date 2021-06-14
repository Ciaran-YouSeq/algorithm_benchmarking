process ALIGN {

        input:
        	tuple val(pair_id), path(read_pairs_ch)

                val 'bowtie2_index'
        output:
                path 'pair_id'

	script:
        """
        bowtie2 \
        --rg-id 'definitelyarealid' \
        --rg 'SM:samplemcsampleid\tLB:libraryname\tPL:ILLUMINA' \
        -p 12 \
        -q \
        -x $bowtie2_index \
        -1  ${read_pairs_ch[0]}  \
        -2  ${read_pairs_ch[1]}  \
        > pair_id
        """
}
