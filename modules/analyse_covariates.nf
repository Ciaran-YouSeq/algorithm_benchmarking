process ANALYSE_COVARIATES {

        input:
                path 'BQSR_table'
        output:
                path 'analyze_covariates_plots.pdf'

        """
        gatk AnalyzeCovariates -bqsr $BQSR_table -plots analyze_covariates_plots.pdf
        """
}
