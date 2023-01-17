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
include { fastqc as fastqc_post           } from '../modules/fastqc/main.nf'

//NF_CORE
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main.nf'
include { UMITOOLS_EXTRACT } from '../modules/nf-core/umitools/extract/main.nf' 
include { BBMAP_BBSPLIT } from '../modules/nf-core/bbsplit/main.nf'
include { SORTMERNA } from '../modules/nf-core/sortmerna/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ_QC {
    
    fastq_ch = Channel.fromFilePairs(params.fastq, checkIfExists:true, size: -1) // default is 2, so set to -1 to allow any number of files
                       .ifEmpty { error "Can not find any fastq matching ${fastq}" }
    fastq_ch.view()
    
    //FASTQC
    fastqc(fastq_ch)
    
    //TRIMGALORE
    //trimgalore_ch = Channel.of( [ id:'test' , single_end:true ])
    //trimgalore_fun =  fastq_ch.combine(trimgalore_ch)
    //trimgalore_fun.meta.id.view()
    //trimgalore_fun.view()

    //TRIMGALORE(trimgalore_fun)
    //TRIMGALORE.out[1].view()
    //TRIMGALORE.out[2].view()
    
    //UMITOOLS_EXTRACT(trimgalore_fun) 
 
    //BBMAP_BBSPLIT(trimgalore_fun)

    //BBMAP_BBSPLIT ( ch_filtered_reads, PREPARE_GENOME.out.bbsplit_index, [], [ [], [] ], false )

    //SORTMERNA(trimgalore_fun)

    //FASTQC_POST_CORRECTION 
    //fastqc_post(TRIMGALORE.out[0])

}

