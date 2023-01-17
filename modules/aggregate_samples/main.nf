process aggregate_samples {
    
    publishDir "${params.outdir}/05_aggregate_samples", mode: 'copy'

    input:
    file(counts_matrix)
    file(rpkm)
    file(tpm)
    
    output:
    tuple file("rawCounts.csv"), file("rpkm.csv"), file("tpm.csv")

    script:
    """
    Rscript $PWD/resources/aggregate_samples/AggregateSample_rawCounts.R --folder ${counts_matrix}  --outputPath rawCounts.csv
    Rscript $PWD/resources/aggregate_samples/AggregateSample_rpkm.R      --folder ${rpkm}           --outputPath rpkm.csv
    Rscript $PWD/resources/aggregate_samples/AggregateSample_tpm.R       --folder ${tpm}            --outputPath tpm.csv
    """
}
