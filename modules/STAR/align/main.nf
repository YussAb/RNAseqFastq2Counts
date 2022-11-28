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


process STAR_aligment {
    maxForks 2
 
    tag "fastq: $sample"
    publishDir "${params.outdir}/01_star_mapped", mode: 'copy'

    input:
    tuple val(sample), file(fastq) 
    file STARgenome 

    output:
    tuple val(fastq.simpleName),file("STAR_${fastq.simpleName}")

    script:
    // STAR Mapper
    """
    STAR --genomeDir ${STARgenome} \
         --readFilesIn ${fastq} \
         --outFileNamePrefix ${fastq.simpleName} \
         --genomeLoad LoadAndRemove \
         --runThreadN 4 \
         --readFilesCommand zcat \
         --outSAMtype BAM Unsorted \
         --limitBAMsortRAM 10000000000 \
         --outSAMunmapped Within \
         --outFilterMultimapNmax 10 \
         --outFilterMultimapScoreRange 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.04

    #--genomeLoad NoSharedMemory
    samtools sort ${fastq.simpleName}Aligned.out.bam > ${fastq.simpleName}Aligned_sortedByCoord.out.bam
    samtools index ${fastq.simpleName}Aligned_sortedByCoord.out.bam
    rm ${fastq.simpleName}Aligned.out.bam 
    
    mkdir STAR_${fastq.simpleName}
    mv ${fastq.simpleName}Aligned* STAR_${fastq.simpleName}/.
    #mv ${fastq.simpleName}Signal* STAR_${fastq.simpleName}/.
    mv ${fastq.simpleName}SJ* STAR_${fastq.simpleName}/.
    mv ${fastq.simpleName}Log* STAR_${fastq.simpleName}/.
    """ 
}



