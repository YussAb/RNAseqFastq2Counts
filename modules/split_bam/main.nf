process split_bam {

    publishDir "${params.outdir}/02_splitBAM", mode:'copy'

    input:
    tuple val(sample), file(STAR_alignment) 
   
    output:
    tuple val(sample), file("${sample}.ribo.in.bam"), file("${sample}.ribo.ex.bam") , file("${sample}.ribo.junk.bam") , file("${sample}.ribo.ex.bam.summary") 

    script:
    """
    split_bam.py -i ${STAR_alignment}/${sample}*.bam \
        -r ${params.ribo_bed} \
        -o ${sample}.ribo \
        >  ${sample}.ribo.ex.bam.summary
    """
}


