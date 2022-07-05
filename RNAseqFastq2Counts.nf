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
    nextflow RNAseqFastq2Counts.nf --fastq "path/to/fastq/fastq/*_R{1,2}_*.fastq.gz" \
                                   --genome path/to/ref.fa \
                                   --annotation path/to/annotation.gtf \
                                   --build_star_genome  \
                                   -with-report report.html
    
    Options:
    * --fastq        [help = path to fastq files to parse to counts ]
    * --genome       [help = path to reference genome.fa ]
    * --annotation   [help = path to annotation.gtf ]
    * --star_genome  [help = path to folder with STAR reference genome]

    If no STAR genome is has already built use this option
    * --build_star_genome  [help = optional to indicate to build STAR reference genome  ]
    
   """.stripIndent()    
}

if (params.help) exit 0, helpMessage()

//Future Development to provide input tsv
//params.input_tsv=''
//Channel
//    .fromPath(params.input_tsv)
//    .splitCsv(header:true, sep:'\t')
//    .map{ row-> tuple(row.tumor_id, row.normal_id, row.oncocode, row.setname, file(row.vcf)) }
//    .set { samples_ch }
//samples_ch.into {samples_ch1; samples_ch2; samples_ch3; samples_ch4 }


//Dev//Documentation//
//REF:https://github.com/evanfloden/lncRNA-Annotation-nf/blob/master/lncRNA-Annotation.nf
//REF:https://github.com/angelovangel/nxf-fastqc/blob/master/main.nf
//REF:Conditional process executions
//https://github.com/nextflow-io/patterns/blob/master/docs/conditional-process.adoc
//REF:Optional input
//https://github.com/nextflow-io/patterns/blob/master/docs/optional-input.adoc

/*
--------------------------------------------------------------------------------
*/
//debug_input_list//
//Required Parameters
//params.fastq = "/homenfs/yabili/giab/fastq/*_R{1,2}_*.fastq.gz"
//params.genome = "/home/youssef/CANDIOLO_IRCCS/rnaseq/yabiliPipe/01_RNAseqFastq2Counts/Intro-to-rnaseq-hpc-O2/unix_lesson/reference_data/chr1.fa"
//params.annotation="/home/youssef/CANDIOLO_IRCCS/rnaseq/yabiliPipe/01_RNAseqFastq2Counts/Intro-to-rnaseq-hpc-O2/unix_lesson/reference_data/chr1-hg19_genes.gtf"
//params.fastq = "/home/youssef/CANDIOLO_IRCCS/rnaseq/yabiliPipe/01_RNAseqFastq2Counts/Intro-to-rnaseq-hpc-O2/unix_lesson/raw_fastq/*.fq*"

//Set Flags//
params.build_star_genome = false

//Input parameters validation//
ref_genome             = file(params.genome)
ref_annotation         = file(params.annotation) 

//validate input files//
if( !ref_annotation.exists() ) exit 1, "Missing annotation file: ${ref_annotation}"
if( !ref_genome.exists() ) exit 1, "Missing genome directory: ${ref_genome}"

// Print some stuff here
println "reads: $params.fastq"

//Create a channel for read files// 
Channel.fromFilePairs(params.fastq, checkIfExists:true, size: -1) // default is 2, so set to -1 to allow any number of files
    .ifEmpty { error "Can not find any fastq matching ${fastq}" }
    .into{ fastq_ch; fastq_ch2; fastq_ch3 }

//fastq_ch3.view()


/*
--------------------------------------------------------------------------------
*/

process fastqc {    
    maxForks 4
    publishDir "${params.outdir}/01_fastqc", mode: 'copy', overwrite: false
    
    input:
    set val(sample), file(fastq) from fastq_ch
    
    output:
    file("${fastq.simpleName}_fastqc/*.zip") into fastqc_files

    script:
    """
    mkdir ${fastq.simpleName}_fastqc
    fastqc -o ${fastq.simpleName}_fastqc -t 6 ${fastq}
    """
}


process runMultiQC {
    publishDir "${params.outdir}/01_multiqc", mode: 'copy', overwrite: false

    input:
    file('*') from fastqc_files.collect()

    output:
    file('multiqc_report.html')

    """
    multiqc .
    """
}

if (params.build_star_genome) { 
    process STAR_genome_index {
    publishDir "${params.outdir}/00_star_index", mode: 'copy'

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
    --runThreadN 16 \
    --outFileNamePrefix STARgenome 
    """
    }
} else {
    params.star_genome=""
    Channel.fromPath(params.star_genome, checkIfExists:true)
        .ifEmpty ("STAR index need to be built")
        .set{ STARgenomeIndex_ch }
}
process STAR_aligment {
    maxForks 2
 
    tag "fastq: $sample"
    publishDir "${params.outdir}/02_star_mapped", mode: 'copy'

    input:
    file STARgenome from STARgenomeIndex_ch.first()
    set val(sample), file(fastq) from fastq_ch2

    output:
     set val(fastq.simpleName),file("STAR_${fastq.simpleName}") into STARmappedReads_ch, STARmappedReadsQualimap_ch

    script:
    // STAR Mapper
    """
    STAR --genomeDir ${STARgenome} \
         --readFilesIn ${fastq} \
         --outFileNamePrefix ${fastq.simpleName} \
         --genomeLoad LoadAndKeep \
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


process Qualimap {
    maxForks 4
    publishDir "${params.outdir}/03_qualimap", mode: 'copy'

    input:    
    set val(sample), file(STAR_alignment) from STARmappedReadsQualimap_ch
    
    output:
    set val(sample), file("${sample}_qm") into qualimap_ch

    script:
    """
    qualimap rnaseq \
    -outdir ${sample}_qm \
    -bam  ${STAR_alignment}/${sample}*.bam \
    -gtf ${ref_annotation} \
    --java-mem-size=8G 
    """

}


//STARmappedReads_ch.view()

process FeautureCounts {
    maxForks 8

    publishDir "${params.outdir}/04_feauture_counts", mode: 'copy'

    input:
    set val(sample), file(STAR_alignment) from STARmappedReads_ch

    output:
    set val(sample), file("${sample}_counts"), file("${sample}_counts_matrix") into featurecounts_ch
    
    script:
    """
    #mkdir featurecounts_${sample}
    featureCounts -s 2 -T 4  \
              -o ${sample}_counts  \
              -a ${ref_annotation} \
              ${STAR_alignment}/${sample}*.bam

    cut -f1,7,8,9,10,11,12 ${sample}_counts >  ${sample}_counts_matrix
    """
}

//Aggiungere una flag per aggregare i dati
//featurecounts_ch.view()


/*
if (params.aggregate_counts) { 
process AggregateCounts {
    }
}
*/

/*
process geneBodyCov {

    publishDir "${params.outdir}/05_gene_body_cov", mode: 'copy'


    input:
    set val(), file() from *_ch

    output:
    set val(), file("${sample}_counts"), file("${sample}_counts_matrix") into featurecounts_ch
 
    script:
    """
    geneBody_coverage.py -r hg19.HouseKeepingGenes.bed -i bam_path.txt | bam1,bam2  -o output_prefix
    """
}
*/
