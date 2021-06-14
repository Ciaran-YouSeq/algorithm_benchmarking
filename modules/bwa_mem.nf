
process ALIGN {
	
    input:
	val ref
	tuple val(pair_id), path(read_pairs_ch)

    output:
      file 'pair_id'

	script:
    """
    bwa-mem2 mem -v 1 -t 4 $ref ${read_pairs_ch[0]} ${read_pairs_ch[1]} > pair_id
    """
	
}
