#!/usr/bin/env nextflow

/*
================================================================================
                yabili/rnaseq_fastq_to_counts
================================================================================
Started Apr 2022.
--------------------------------------------------------------------------------
yabili/rnaseq_fastq_to_counts:
  Generates tables with read counts for multiple samples 
  starting from fastq files
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/YussAb/...
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    The typical command for running the pipeline is as follows:
    nextflow .....nf --input  ...  -with-report report.html
    
    Mandatory arguments:
    --arg

   """.stripIndent()    
}

if (params.help) exit 0, helpMessage()

//params.input_tsv=''
//Channel
//    .fromPath(params.input_tsv)
//    .splitCsv(header:true, sep:'\t')
//    .map{ row-> tuple(row.tumor_id, row.normal_id, row.oncocode, row.setname, file(row.vcf)) }
//    .set { samples_ch }
//samples_ch.into {samples_ch1; samples_ch2; samples_ch3; samples_ch4 }



//REF:https://github.com/evanfloden/lncRNA-Annotation-nf/blob/master/lncRNA-Annotation.nf
/*
--------------------------------------------------------------------------------
*/
//input_list//

// Required Parameters
//params.reads = "/homenfs/yabili/giab/fastq/*_R{1,2}_*.fastq.gz"
params.reads = "/home/youssef/CANDIOLO_IRCCS/rnaseq/IsellaPipe/RNAseqFastq2Counts/Intro-to-rnaseq-hpc-O2/unix_lesson/raw_fastq/*.fq*"
params.genome = "/home/youssef/CANDIOLO_IRCCS/rnaseq/IsellaPipe/RNAseqFastq2Counts/Intro-to-rnaseq-hpc-O2/unix_lesson/reference_data/chr1.fa"
params.annotation="/home/youssef/CANDIOLO_IRCCS/rnaseq/IsellaPipe/RNAseqFastq2Counts/Intro-to-rnaseq-hpc-O2/unix_lesson/reference_data/chr1-hg19_genes.gtf"

//Input parameters validation//

ref_genome             = file(params.genome)
ref_annotation         = file(params.annotation) 

FEELnc_filter_options=params.feelnc_opts

//validate input files//

if( !ref_annotation.exists() ) exit 1, "Missing annotation file: ${ref_annotation}"
if( !ref_genome.exists() ) exit 1, "Missing genome directory: ${ref_genome}"
// Print some stuff here
println "reads: $params.reads"


//Create a channel for read files//
 
Channel.fromPath(params.reads, checkIfExists:true)
        .into{ reads_ch; reads_ch2 }

//Channel.fromFilePairs(params.reads, checkIfExists:true)
//       .set{ reads_pairs_ch }
      
//reads_pairs_ch.view()

process fastqc {    
    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: false
    
    input:
    file(fastq) from reads_ch
    
    output:
    file("${fastq.simpleName}_fastqc/*.zip") into fastqc_files

    script:
    """
    mkdir ${fastq.simpleName}_fastqc
    fastqc -o ${fastq.simpleName}_fastqc -t 4 ${fastq}
    """

}


process runMultiQC{
    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: false

    input:
    file('*') from fastqc_files.collect()

    output:
    file('multiqc_report.html')

    """
    multiqc .
    """
}

//option step add an if

process STAR_genome_index {
    publishDir "${params.outdir}/star_index", mode: 'copy'

    input:
    file ref_genome
    file ref_annotation

    output:
    file "STARgenome" into STARgenomeIndex_ch

    script:
    //STAR Generate Index and create genome length file
    """
    mkdir STARgenome
    STAR --runMode genomeGenerate \
    --genomeDir STARgenome \
    --genomeFastaFiles ${ref_genome} \
    --sjdbGTFfile ${ref_annotation} \
    --runThreadN 6 \
    --outFileNamePrefix STARgenome 
    """
}


process STAR_aligment {
    
    tag "reads: $fastq.simpleName"
    publishDir "${params.outdir}/star_mapped", mode: 'copy'

    input:
    file STARgenome from STARgenomeIndex_ch.first()
    file(fastq) from reads_ch2

    output:
    set val(fastq.simpleName), file("STAR_${fastq.simpleName}") into STARmappedReads_ch

    script:
    // STAR Mapper
    """
    STAR --genomeDir ${STARgenome} \
         --readFilesIn ${fastq} \
         --outFileNamePrefix ${fastq.simpleName} \
         --genomeLoad NoSharedMemory \
         --runThreadN 4 \
         --readFilesCommand zcat \
         --outSAMtype BAM Unsorted \
         --limitBAMsortRAM 10000000000 \
         --outSAMunmapped Within \
         --outFilterMultimapNmax 10 \
         --outFilterMultimapScoreRange 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.04
    
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

//STARmappedReads_ch.view()

process FeautureCounts {

    publishDir "${params.outdir}/feauture_counts", mode: 'copy'

    input:
    set val(name), file(STAR_alignment) from STARmappedReads_ch

    output:
    set val(name), file("${name}_counts"), file("${name}_counts_matrix") into featurecounts_ch
    
    script:
    """
    #mkdir featurecounts_${name}
    featureCounts -s 2 -T 4  \
              -o ${name}_counts  \
              -a ${ref_annotation} \
              ${STAR_alignment}/${name}*.bam

    cut -f1,7,8,9,10,11,12 ${name}_counts >  ${name}_counts_matrix
    """
}
