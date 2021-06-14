params.pub_dir = "$baseDir/data/qc_results"

process FASTQC {
 
	publishDir params.pub_dir, mode: 'copy', overwrite: true

 	input:
		tuple val(name), file(reads)

	output:
		file "*_fastqc.{zip,html}"

	script:
  	"""
  	fastqc $reads
 	"""
}

process MULTIQC {

       publishDir params.pub_dir, mode: 'copy', overwrite: true

    	input:
    		file ('fastqc/*')

    	output:
    		file "*multiqc_report.html"
    		file "*_data"

    	script:
    	"""
    	multiqc . -m fastqc
    	"""
}
