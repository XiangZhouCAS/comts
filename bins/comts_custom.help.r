if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--help","-h"),type = "character",default = F,
              help = "Show this help message and exit"))
opt_parser <- OptionParser(
  usage = "usage: comts custom diy    To calculate single copy function enzyme gene community abandance by custom database.\n       comts custom ter    To calculate single copy terminal enzyme gene community abandance.\n       comts custom hyd    To calculate single copy Hydrogenase community abandance.",
  add_help_option = TRUE,
  prog=NULL ,
  description = "This page is to show how to run the program.")
opt <- parse_args(opt_parser)
