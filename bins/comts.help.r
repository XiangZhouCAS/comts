if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--help","-h"),type = "character",default = F,
              help = "Show this help message and exit"))
opt_parser = OptionParser(
  usage = "usage: comts geneset	To calculate community and RPKM abundance of GeneSet.\n       comts custom	To calculate single copy function enzyme gene community abandance by custom database.\n       comts stat	To generate a statistical table of differences and display the results graphically (Barplot and Heatmap).",
  add_help_option = TRUE,
  prog=NULL ,
  description = "This page is to show how to run the program.")
opt <- parse_args(opt_parser)
