/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Loaded from ./modules
//
include { fastqc            } from '../modules/fastqc/main.nf'
include { STAR_genome_index } from '../modules/STAR/align/main.nf'
include { STAR_aligment     } from '../modules/STAR/genome_index/main.nf'
include { Qualimap          } from '../modules/qualimap/main.nf'
include { split_bam         } from '../modules/split_bam/main.nf'
include { FeautureCounts    } from '../modules/FeautureCounts/main.nf'
include { normalize_counts  } from '../modules/normalize_counts/main.nf'


//
// WORKFLOW: Run main nf-core/rnaseq analysis pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {
    
    fastq_ch = Channel.fromFilePairs(params.fastq, checkIfExists:true, size: -1) // default is 2, so set to -1 to allow any number of files
                       .ifEmpty { error "Can not find any fastq matching ${fastq}" }
    fastqc(fastq_ch)

    if (params.build_star_genome){
        STAR_genome_index(ref_genome, ref_annotation)
        STARgenomeIndex_ch = STAR_genome_index.out
        STAR_aligment(fastq_ch, STARgenomeIndex_ch )
    }else{
        STARgenomeIndex_ch = Channel.fromPath(params.star_genome, checkIfExists:true).ifEmpty ("STAR index need to be built")
        STAR_aligment(fastq_ch, STARgenomeIndex_ch )}
    
    Qualimap(STAR_aligment.out)    
    split_bam(STAR_aligment.out)
    FeautureCounts(split_bam.out)
    normalize_counts(FeautureCounts.out)

}

