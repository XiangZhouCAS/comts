if (!require(optparse,quietly = TRUE)) {
  install.packages("optparse")
  library(optparse)
} else {
  library(optparse)}
option_list <- list(
  make_option(c("--input_reads","-i"),type = "character",default = F,
                                help = "please set the directory of reads"),
  make_option(c("--result","-o"),type = "character",default = F,
              help = "please set the result file name"),
  make_option(c("--threads","-t"),type = "numeric",default = 40,
              help = "Setting the threads of CPU,default is 40"),
  make_option(c("--singleM","-s"),type = "character",default = F,
              help = "please set the directory of singleM database"),
  make_option(c("--run_fastp","-f"),type = "character",default = "run",
              help = "If you have already filtered the reads (length >= 140), you can set `skip` to skip run fastp, default is run fastp"),
  make_option(c("--run_seqkit","-k"),type = "character",default = "run",
	      help = "If you have already counted the number of all reads with seqkit, you can set the directory of seqkit result to skip run seqkit, default is run seqkit"),
  make_option(c("--keep_samples","-e",type = "character",default = "keep"),
              help = "If you do not set 'keep', you can to delete these tmp results, default is keep."))
opt_parser = OptionParser(
  usage = "usage: comts rpkm ribo [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL ,
  description = "This Script is to calculate rpkm abandance through geneset.")
opt <- parse_args(opt_parser)
input_reads <- opt$input_reads
run_fastp <- opt$run_fastp
input_geneset <- opt$input_geneset
res <- opt$result
singleM <- opt$singleM
threads <- opt$threads
run_seqkit <- opt$run_seqkit
keep_samples <- opt$keep_samples
singleM_out <- paste0(res,".ribo.txt")
fastp_output <- paste0(res,".filtered.fq.gz")
seqkit_out <- paste0(res,".all.reads.txt")
seqkit <- sprintf("seqkit stat %s > %s",
                  fastp_output,seqkit_out)
seqkit2 <- sprintf("seqkit stat %s > %s",
                  input_reads,seqkit_out)
fastp <- sprintf("fastp -i %s -o %s --length_required 140 -w %d",
                 input_reads,fastp_output,threads)
diamond_singleM <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle length --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                           singleM, fastp_output, singleM_out, threads)
diamond_singleM2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle length --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                            singleM, input_reads, singleM_out, threads)
if(run_fastp != "skip"){
  system(fastp)}else{
    print("Not run fastp because you set the directory of filtered reads to skip it.")}
if(run_fastp != "skip"){
  system(diamond_singleM)}else{
    system(diamond_singleM2)
  }
if(run_fastp != "run"){
   if(run_seqkit != "run"){
   print("Not run seqkit because you set the directory of seqkit result to skip it")}else{system(seqkit2)}}else{
	   system(seqkit)}

if(!require(magrittr)){
  install.packages("magrittr")
  library(magrittr)
}else{
  library(magrittr)}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}else{
  library(dplyr)}

d2 <- read.table(singleM_out,sep = "\t")
colnames(d2)[c(1,2,3)] <- c("slen","sseq_id","aligned_length")
d2 <- filter(d2,aligned_length > 40)
d2 <- d2[,-3]
d2 <- d2%>%
  group_by_all()%>%
  count()
s_out <- paste0(res,".ribo.rpkm.txt")
if(run_seqkit == "run"){
  p <- read.table(seqkit_out,header = T)
  }else{p <- read.table(run_seqkit,header = T)}
d2$all_reads_nums <- p$num_seqs
d2$all_reads_nums <- gsub(",","",d2$all_reads_nums)
d2$all_reads_nums <- as.numeric(d2$all_reads_nums)
d2$RPKM <- d2$n*10^9/(d2$slen*d2$all_reads_nums*3)
d2$sample_name <- res
d2 <- d2[,c(6,2,1,3,4,5)]
d2$sseq_id <- gsub("-.*","",d2$sseq_id)
d2 <- d2%>%
  group_by(sseq_id)%>%
  mutate(T_RPKM = sum(RPKM))
d2 <- unique(d2[,c(1,2,7)])
write.table(d2,s_out,quote = F,sep = "\t",row.names = F)
system(sprintf("cat %s >> Ribo.rpkm.txt",s_out))
total <- read.table("Ribo.rpkm.txt",header = T,sep = "\t")
total <- filter(total,T_RPKM != "T_RPKM")
write.table(total,"Ribo.rpkm.txt",quote = F,sep = "\t",row.names = F)

system(sprintf("mkdir %s",res))
if(run_seqkit == "run"){
  system(sprintf("mv %s %s",seqkit_out,res))
  system(sprintf("mv %s %s",s_out,res))
  system(sprintf("mv %s %s",singleM_out,res))}else{
    system(sprintf("mv %s %s",s_out,res))
    system(sprintf("mv %s %s",singleM_out,res))}
if(run_fastp == "run"){
  system(sprintf("mv %s %s",fastp_output,res))}else{
    print("There is no fastp result.")}
wd1 <- paste0(getwd(),"/",res)
wd2 <- paste0(getwd())
if(keep_samples == "keep"){
  print("……")}else{if(file.exists("samples") == T){
    setwd(wd1)
    system(sprintf("mv %s ../samples",s_out))
    setwd(wd2)
    system(sprintf("rm %s -rf",res))}else{
      system(sprintf("mkdir samples"))
      setwd(wd1)
      system(sprintf("mv %s ../samples",s_out))
      setwd(wd2)
      system(sprintf("rm %s -rf",res))}}
