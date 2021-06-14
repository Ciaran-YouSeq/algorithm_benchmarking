params.pub_dir = "$baseDir/data/trimgalore_logs"

process TRIMGALORE {
	
	publishDir params.pub_dir, mode: 'copy', overwrite: true

        input:
                tuple val(pair_id), path(read_pairs_ch)
        output:
                tuple val(pair_id), path("*.fq.gz")    , emit: reads
    		tuple val(pair_id), path("*report.txt"), emit: log

        """
        trim_galore \
	${read_pairs_ch[0]} \
	${read_pairs_ch[1]} \
	--fastqc \
	--paired \
	--stringency 3 \
	--gzip \
	--suppress_warn \
	> out 2>err
        """
}
