if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--help","-h"),type = "character",default = F,
              help = "Show this help message and exit"))
opt_parser = OptionParser(
  usage = "comts geneset ribo     To calculate RPKM abandance of 14 universal single copy ribosomal genes (USCGs).\n       comts geneset res      To convert RPKM to AFG of single copy function gene through GeneSet.",
  add_help_option = TRUE,
  prog=NULL ,
  description = "This page is to show how to run the comts.")
opt <- parse_args(opt_parser)
