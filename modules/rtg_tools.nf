process BGZIP_VCF {

        input:
                path 'haplotype_caller_output'
        output:
                path  "${haplotype_caller_output}.gz"
        """
        rtg bgzip --stdout $haplotype_caller_output > "${haplotype_caller_output}.gz"
        """
}

process INDEX_GZIPPED_VCF {

        input:
                path 'haplotype_caller_output.gz'
        output:
                path 'haplotype_caller_output.gz.tbi'
        """
        rtg index --format=vcf 'haplotype_caller_output.gz'
        """
}

process RTG_VCFEVAL {

        input:
                path baseline_calls_vcf
                path baseline_calls_vcf_index
                path 'haplotype_caller_output.gz'
                path 'haplotype_caller_output.gz.tbi'
                path 'reference_genome_sdf'
        output:
                path 'rtg_output'
        """
        rtg vcfeval --baseline=$baseline_calls_vcf --calls='haplotype_caller_output.gz' --output=rtg_output --template=$reference_genome_sdf
        """
}
