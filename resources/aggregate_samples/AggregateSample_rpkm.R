#DOCS:

##########################Inputs################################################
# Load libraries
suppressMessages(library(optparse))

# Parameters
option_list = list(optparse::make_option(c("-f", "--folder"), type = "character", default = F, help = "folder containing the file rpkm", metavar = ""),
                   optparse::make_option(c("-o", "--outputPath"), type = "character",  default = F , help = "output file path", metavar = ""))

opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)
folder = opt$folder ; 
outputPath = opt$outputPath

################################################################################
#DEBUG
#folder="/data/yabili/RNAseqFastq2Counts/xpo1_3rd_run_seq/04_normalized_counts"
#outputPath="/data/yabili/RNAseqFastq2Counts/xpo1_3rd_run_seq/04_normalized_counts/xpo1_3rd_rpkm.tsv"


##########################STEP1#################################################
print("Reading Tables...")
#leggere una lista di csv
#https://stackoverflow.com/questions/11433432/how-to-import-multiple-csv-files-at-once
#folder = "~/yussab/filterGenesVCF/03_snpeff"
setwd(folder)
temp = list.files(pattern= "/*_rpkm.txt")
#myfiles = lapply(temp, read.csv)
#library(tools)
#https://stackoverflow.com/questions/39556188/how-to-get-named-list-when-reading-multiple-csv-files-from-given-folder
suppressMessages(library(data.table))
setDTthreads(32)
myfiles <- setNames(lapply(temp, data.table::fread), tools::file_path_sans_ext(basename(temp)))

#https://stackoverflow.com/questions/12664430/delete-a-column-in-a-data-frame-within-a-list
# Remove column by name in each element of the list
suppressMessages(library(purrr))
suppressMessages(library(dplyr))
myfiles <- map(myfiles, ~ (.x %>% select(-Chr, -Start , -End, -Strand, -Length)))


##########################STEP2#################################################
#fare rbind della lista in un file unico
print("Merging Tables...")


myfiles <- myfiles[sapply(myfiles, nrow)>0]
df <- myfiles  %>% reduce(left_join, by = "Geneid")

print("Writing Table...")
write.table(df, file = outputPath , sep="\t", row.names = FALSE, quote = FALSE)
print("DONE!")
