process Qualimap {
    maxForks 4

    publishDir "${params.outdir}/QC/02_bam_qualimap", mode: 'copy'

    input:    
    tuple val(sample), file(STAR_alignment)
    
    output:
    tuple val(sample), file("${sample}_qm")

    script:
    """
    qualimap rnaseq \
    -outdir ${sample}_qm \
    -bam    ${STAR_alignment}/${sample}*.bam \
    -gtf    ${ref_annotation} \
    -p      ${params.standness_qualimap} \
    --java-mem-size=8G 
    """
}



