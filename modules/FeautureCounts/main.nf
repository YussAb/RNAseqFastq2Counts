
process FeautureCounts {
    maxForks 8

    publishDir "${params.outdir}/03_feauture_counts", mode: 'copy'

    input:
    tuple val(sample), file(ribo_in), file(ribo_ex), file(ribo_junk), file(ribo_summary)
    
    output:
    tuple val(sample), file("${sample}_counts"), file("${sample}_counts_matrix"), file("${sample}_counts.summary")
 
    script:
    """
   featureCounts -T 4 -t gene \
              -o ${sample}_counts  \
              -a ${ref_annotation} \
              -s ${params.strandness_featurecounts} \
              ${ribo_ex}


    cut -f1,7,8,9,10,11,12 ${sample}_counts >  ${sample}_counts_matrix
    """
}


