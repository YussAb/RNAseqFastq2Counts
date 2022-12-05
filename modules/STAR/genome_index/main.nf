process STAR_genome_index {
    publishDir "${params.outdir}/00_star_index", mode: 'copy'

    input:
    file(ref_genome)
    file(ref_annotation)

    output:
    file("STARgenome")

    script:
    //STAR Generate Index and create genome length file
    """
    mkdir STARgenome
    STAR --runMode genomeGenerate \
    --genomeDir STARgenome \
    --genomeFastaFiles ${ref_genome} \
    --sjdbGTFfile ${ref_annotation} \
    --runThreadN 16 \
    --outFileNamePrefix STARgenome 
    """
    }



