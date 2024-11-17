if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--reads_count","-i"),type = "character",
              default = F,
              help = "Please set the directory of reads count file."),
  make_option(c("--gene_length","-l"),type = "character",
              default = F,
              help = "Please set the directory of gene_length file"),
  make_option(c("--output","-o"),type = "character",
              default = F,
              help = "Please set the result name's prefix."))
opt_parser <- OptionParser(
  usage = "usage: comts rpkm [option]",
  option_list = option_list,
  add_help_option = TRUE,
  prog = NULL,
  description = "To calculate the RPKM abundance of geneset.")
opt <- parse_args(opt_parser)
reads_count <- opt$reads_count
gene_length <- opt$gene_length
output <- opt$output
if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)
}else{
  library(tidyr)}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}else{
  library(dplyr)}
l <- read.table(gene_length,header = T,sep = "\t")
r <- read.table(reads_count,header = T,sep = "\t",
                check.names = F)
colnames(r)[1] <- "GeneID"
colnames(l)[1] <- "GeneID"
colnames(l)[2] <- "GeneLen"
r <- left_join(r,l,by = "GeneID")
for(i in 2:(ncol(r)-1)){
    r[,i] <- as.numeric(r[,i])
    r[,i] <- r[,i]*10^9/(sum(r[,i])*r$GeneLen)}
r <- select(r,-GeneLen)
res <- paste0(output,".rpkm.txt")
write.table(r,res,sep = "\t",quote = F,row.names = F)
