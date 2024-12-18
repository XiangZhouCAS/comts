if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--help","-h"),type = "character",default = F,
              help = "Show this help message and exit"))
opt_parser = OptionParser(
  usage = "usage: comts geneset	To calculate AFG and RPKM of GeneSet.\n       comts custom	To calculate single copy function enzyme AFG by customed database.",
  add_help_option = TRUE,
  prog=NULL ,
  description = "This page is to show how to run the program.")
opt <- parse_args(opt_parser)
