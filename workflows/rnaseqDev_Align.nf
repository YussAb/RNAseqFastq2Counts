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
include { STAR_aligment     } from '../modules/STAR/align/main.nf'
include { STAR_genome_index } from '../modules/STAR/genome_index/main.nf'
include { Qualimap          } from '../modules/qualimap/main.nf'
include { split_bam         } from '../modules/split_bam/main.nf'
include { FeautureCounts    } from '../modules/FeautureCounts/main.nf'
include { normalize_counts  } from '../modules/normalize_counts/main.nf'
include { aggregate_samples } from '../modules/aggregate_samples/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ_ALIGN_QUANT {

    fastq_ch = Channel.fromFilePairs(params.fastq, checkIfExists:true, size: -1) // default is 2, so set to -1 to allow any number of files
                       .ifEmpty { error "Can not find any fastq matching ${fastq}" }
    //fastq_ch.view()
    
 
    //DEFINE A SPECIFIC INPUT


    if (params.build_star_genome){
        STAR_genome_index(ref_genome, ref_annotation)
        STARgenomeIndex_ch = STAR_genome_index.out
        fastq_genome_ch = fastq_ch.combine(STARgenomeIndex_ch)
        STAR_aligment(fastq_genome_ch )
    }else{
        STARgenomeIndex_ch = Channel.fromPath(params.star_genome, checkIfExists:true).ifEmpty ("STAR index need to be built")
        fastq_genome_ch = fastq_ch.combine(STARgenomeIndex_ch)
        //fastq_genome_ch.view()
        STAR_aligment(fastq_genome_ch )}

    Qualimap(STAR_aligment.out)    
    split_bam(STAR_aligment.out)
    FeautureCounts(split_bam.out)
    normalize_counts(FeautureCounts.out)


    //DEV
    //CHECK EMITTED CH
    //FeautureCounts.out.view()
    //normalize_counts.out.view()

    //JOIN CHANNELS
    //aggregate_samples_ch = normalize_counts.out.join(FeautureCounts.out , by: [0] )

    //rpkm_ch = aggregate_samples_ch.map { it[1] }.collect()
    //tpm_ch = aggregate_samples_ch.map { it[2] }.collect()
    //rawCounts_ch = aggregate_samples_ch.map { it[4] }.collect()

    
    //rpkm_ch.view() 
    //tpm_ch.view()
    //rawCounts_ch.view()

    //aggregate_samples(rawCounts_ch,rpkm_ch,tpm_ch )
}
