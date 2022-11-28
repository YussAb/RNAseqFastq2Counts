process fastqc {     
     maxForks 4
 
     publishDir "${params.outdir}/QC/01_fastqc", mode: 'copy', overwrite: false
     
     input:
     tuple val(sample), file(fastq) 
     
     output:
     file("${fastq.simpleName}_fastqc/*.zip")
 
     script:
     """
     mkdir -p  ${fastq.simpleName}_fastqc
     fastqc -o ${fastq.simpleName}_fastqc -t 6 ${fastq}
     """
 }
