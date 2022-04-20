#DOCS:

##########################Inputs################################################
# Load libraries
suppressMessages(library(optparse))

# Parameters
option_list = list(optparse::make_option(c("-f", "--folder"), type = "character", default = F, help = "folder containing the file generated with feautureCounts *counts_matrix", metavar = ""),
                   optparse::make_option(c("-o", "--outputPath"), type = "character",  default = F , help = "output file path", metavar = ""))

opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)
folder = opt$folder ; 
outputPath = opt$outputPath

################################################################################
#DEBUG
#folder="~/CANDIOLO_IRCCS/rnaseq/yabiliPipe/01_RNAseqFastq2Counts/results/feauture_counts"
#outputPath="~/CANDIOLO_IRCCS/rnaseq/yabiliPipe/01_RNAseqFastq2Counts/output.tsv"


##########################STEP1#################################################
print("Reading Tables...")
#leggere una lista di csv
#https://stackoverflow.com/questions/11433432/how-to-import-multiple-csv-files-at-once
#folder = "~/yussab/filterGenesVCF/03_snpeff"
setwd(folder)
temp = list.files(pattern= "/*counts_matrix")
#myfiles = lapply(temp, read.csv)
#library(tools)
#https://stackoverflow.com/questions/39556188/how-to-get-named-list-when-reading-multiple-csv-files-from-given-folder
suppressMessages(library(data.table))
setDTthreads(32)
myfiles <- setNames(lapply(temp, data.table::fread), tools::file_path_sans_ext(basename(temp)))


##########################STEP2#################################################
#fare rbind della lista in un file unico
print("Merging Tables...")
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
myfiles <- myfiles[sapply(myfiles, nrow)>0]
df <- myfiles  %>% reduce(left_join, by = "Geneid")

print("Writing Table...")
write.table(df, file = outputPath , sep="\t", row.names = FALSE, quote = FALSE)
print("DONE!")
