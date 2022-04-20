# RNAseq from fatsq to counts

Nextflow script to generate counts starting from fastq files.
Allow analysis of multiple samples at the same time

## RNAseqFastq2Counts.nf

Options:
* --fastq        [help = path to fastq files to parse to counts ]
* --genome       [help = path to reference genome.fa ]
* --annotation   [help = path to annotation.gtf ]
* --star_genome  [help = path to folder with STAR reference genome]

If no STAR genome is has already built use this option
* --build_star_genome  [help = optional to indicate to build STAR reference genome  ]


Example:

```bash
nextflow RNAseqFastq2Counts.nf --fastq "path/to/fastq"  --genome "path/to/ref.fa" --annotation "path/to/annotation.gtf" --star_genome "path/to/starGenome/dir" -with-report report.html
```

## AggregateSamples.R

Options:
* -f , --folder     [help = folder containing the file generated with feautureCounts *counts_matrix ]
* -o , --outputPath [help = output file path] 

Example:
```bash
Rscript AggregateSamples.R --folder results/04_feauture_counts  --outputPath test.csv
```  
 
