process normalize_counts {

    maxForks 8

    publishDir "${params.outdir}/04_normalized_counts", mode: 'copy'

    input:
    tuple val(sample), file(counts), file(counts_matrix), file(counts_summary)
    
    output:
    tuple val(sample), file("${sample}_rpkm.txt"), file("${sample}_tpm.txt")
   
    script:
    """
    Rscript $PWD/resources/rpkm_tpm_script.R ${counts} ${sample}
    """

}


