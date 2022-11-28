# RNAseq from fatsq to counts 

Nextflow script to generate counts starting from fastq files.
Allow analysis of multiple samples at the same time on different platforms.

This pipeline is implemented using **DSL2** from **Nextflow**.
In this repo you can find modules and workflows that you can use as wanted to create a new personalized workflow.

## Usage

Example:

```bash
nextflow main.nf --fastq "path/to/fastq"  --genome "path/to/ref.fa" --annotation "path/to/annotation.gtf" --star_genome "path/to/starGenome/dir" -with-report report.html
```

Options:
* --fastq        [help = path to fastq files to parse to counts ]
* --genome       [help = path to reference genome.fa ]
* --annotation   [help = path to annotation.gtf ]
* --star_genome  [help = path to folder with STAR reference genome]

If no STAR genome is has already built use this option
* --build_star_genome  [help = optional to indicate to build STAR reference genome  ]



