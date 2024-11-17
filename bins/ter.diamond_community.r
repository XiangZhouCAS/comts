if (!require(optparse,quietly = TRUE)) {
  install.packages("optparse")
  library(optparse)
} else {
  library(optparse)}
option_list <- list(
  make_option(c("--input_reads","-i"),type = "character",default = F,
              help = "please set the reads,only the forward reads if the data is PE"),
  make_option(c("--result","-o"),type = "character",default = F,
              help = "please set the result file name"),
  make_option(c("--threads","-t"),type = "numeric",default = 40,
              help = "Setting the threads of CPU,default is 40"),
  make_option(c("--diamond_db","-d"),type = "character",default = F,
              help = "please set the directory of diamond database"),
  make_option(c("--singleM","-s"),type = "character",default = F,
              help = "please set the directory of singleM database, but if you have already counted the RPKM of singleM marker genes, you can set the directory of result to skip it"),
  make_option(c("--run_fastp","-f"),type = "character",default = "run",
              help = "If you have already filtered the reads (length > 140), you can set `skip` to skip run fastp, default is run fastp"),
  make_option(c("--run_seqkit","-k"),type = "character",default = "run",
	      help = "If you have already counted the number of all reads with seqkit, you can set the directory of seqkit result to skip run seqkit, default is run seqkit"),
  make_option(c("--keep_samples","-e"),type = "character",default = "keeps",
	      help = "If you do not set 'keep', you can to delete these tmp results, default is keep."))
opt_parser = OptionParser(
  usage = "usage: comts custom ter [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL ,
  description = "This Script is to calculate terminal enzyme gene community abandance.")
opt <- parse_args(opt_parser)
input_reads <- opt$input_reads
outpath <- opt$result
threads <- opt$threads
fastp_output <- paste0(outpath,".filtered.fq.gz")
diamond_db <- opt$diamond_db
singleM <- opt$singleM
system(sprintf("file %s > tmp.txt",singleM))
system(sprintf("sed 's/.*://' -i tmp.txt"))
system(sprintf("sed 's/ ASCII //g' -i tmp.txt"))
tmp <- read.table("tmp.txt",sep = "\t")
run_fastp <- opt$run_fastp
run_seqkit <- opt$run_seqkit
keep_samples <- opt$keep_samples
diamond_out <- paste0(outpath,".hits.txt")
singleM_out <- paste0(outpath,".singleM.hits.txt")
seqkit_out <- paste0(outpath,".all.reads.txt")
fastp <- sprintf("fastp -i %s -o %s --length_required 140 -w %d",
                 input_reads,fastp_output,threads)
diamond <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle pident qcovhsp --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                   diamond_db, fastp_output, diamond_out, threads)
diamond_singleM <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle length --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                           singleM, fastp_output, singleM_out, threads)
seqkit <- sprintf("seqkit stat %s > %s",
                  fastp_output,seqkit_out)
diamond2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle pident qcovhsp --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
		   diamond_db, input_reads, diamond_out, threads)
diamond_singleM2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle length --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
			    singleM, input_reads, singleM_out, threads)
seqkit2 <- sprintf("seqkit stat %s > %s",
		   run_fastp,seqkit_out)
if(run_fastp != "skip"){
	system(fastp)}else{
		print("Not run the fastp because you set the directory of filtered reads to skip it.")}
if(run_fastp != "skip"){
	system(diamond)}else{system(diamond2)}
if(run_fastp != "skip"){
  if(tmp$V1 != "text"){
    system(diamond_singleM)
  }else{print("Not count the RPKM of singleM marker genes.")}
  }else{if(tmp$V1 != "text" ){
  system(diamond_singleM2)
    }else{print("Not count the RPKM of singleM marker genes.")}}
if(run_fastp != "skip"){
  if(run_seqkit == "run"){
    system(seqkit)}else{
      print("Not run the seqkit because you set the directory of seqkit result")}
  }else{if(run_seqkit == "run"){
  system(seqkit2)
}else{print("Not run the seqkit because you set the directory of seqkit result")}}

if (!require(magrittr)) {
  install.packages("magrittr")
  library(magrittr)
} else {
  library(magrittr)}
if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
} else {
  library(dplyr)}
if (!require(tidyr)) {
  install.packages("tidyr")
  library(tidyr)
} else {
  library(tidyr)}

if(run_seqkit == "run"){
  p <- read.table(seqkit_out,header = T)
  }else{p <- read.table(run_seqkit,header = T)}
d1 <- read.table(diamond_out,sep = "\t")
colnames(d1)[c(1:4)] <- c("slen","gene","pident","qcovhsp")
d1 <- d1%>%
  group_by_all()%>%
  count()
d1$all_reads_nums <- p$num_seqs
d1$all_reads_nums <- gsub(",","",d1$all_reads_nums)
d1$all_reads_nums <- as.numeric(d1$all_reads_nums)
d1$tmp_RPKM <- d1$n*10^9/(d1$slen*d1$all_reads_nums*3)
d1$sample_name <- basename(outpath)
d1$gene <- gsub(" .*","",d1$gene)
d12 <- filter(d1,gene %in% c("PsaA"))
d13 <- filter(d1,gene %in% c("HbsT"))
d14 <- filter(d1,gene %in% c("PsbA","IsoA","AtpA",
                             "YgfK","ARO"))
d15 <- filter(d1,gene %in% c("CoxL","MmoX","AmoA","NxrA","NuoF","RbcL"))
d11 <- filter(d1,!gene %in% c(d12$gene,d13$gene,d14$gene,d15$gene))%>%
  filter(pident > 50 & qcovhsp > 80)
d12 <- filter(d12,pident > 80 & qcovhsp > 80)
d13 <- filter(d13,pident > 75 & qcovhsp > 80)
d14 <- filter(d14,pident > 70 & qcovhsp > 80)
d15 <- filter(d15,pident > 60 & qcovhsp > 80)
d1 <- rbind(d11,d12,d13,d14,d15)
d1 <- d1[,c(8,2,7)]
d1 <- d1%>%
  group_by(gene)%>%
  mutate(RPKM = sum(tmp_RPKM))
d1 <- unique(d1[,-3])
if(tmp$V1 != "text"){
  d2 <- read.table(singleM_out,sep = "\t")
  }else{
  d2 <- read.table(singleM,sep = "\t")
  }
colnames(d2)[c(1,2,3)] <- c("slen","sseq_id","aligned_length")
d2 <- filter(d2,aligned_length > 40)
d2 <- d2[,-3]
d2 <- d2%>%
  group_by_all()%>%
  count()
d2$all_reads_nums <- p$num_seqs
d2$all_reads_nums <- gsub(",","",d2$all_reads_nums)
d2$all_reads_nums <- as.numeric(d2$all_reads_nums)
d2$RPKM <- d2$n*10^9/(d2$slen*d2$all_reads_nums*3)
d2$sample_name <- basename(outpath)
d2 <- d2[,c(6,2,1,3,4,5)]
d2$sseq_id <- gsub("-.*","",d2$sseq_id)
d2 <- d2%>%
  group_by(sseq_id)%>%
  mutate(T_RPKM = sum(RPKM))
d2 <- unique(d2[,c(1,2,7)])
s_out <- paste0(outpath,".singleM.hits.comb.txt")
write.table(d2,s_out,quote = FALSE,sep = "\t",row.names = F)
geo <- exp(mean(log(d2$T_RPKM)))
d1$community_abd <- d1$RPKM/geo*100
d1$community_abd <- ifelse(d1$community_abd > 100,100,d1$community_abd)
r_out <- paste0(outpath,".all.abd.txt")
write.table(d1,r_out,quote = FALSE,sep = "\t",row.names = F)
system(sprintf("mkdir %s",outpath))
if(run_fastp == "run"){
	system(sprintf("mv %s %s",fastp_output,outpath))
   }else{print("There is no need to move fastp result to outpath")}
system(sprintf("mv %s %s",diamond_out,outpath))
if(tmp$V1 != "text"){
  system(sprintf("mv %s %s",singleM_out,outpath))}else{
	  print("There is no need to move singleM result to outpath")}
if(run_seqkit == "run"){
  system(sprintf("mv %s %s",seqkit_out,outpath))}else{
    print("There is no need to move seqkit result to outpath")}
system(sprintf("cat %s >> all.abd.txt",r_out))
tmp <- read.table("all.abd.txt",sep = "\t",header = T)
tmp <- filter(tmp,sample_name != "sample_name")
write.table(tmp,"all.abd.txt",quote = F,sep = "\t",row.names = F)
system(sprintf("mv %s %s",s_out,outpath))
system(sprintf("mv %s %s",r_out,outpath))
if(run_fastp != "skip"){
	system("rm \"fastp.html\" \"fastp.json\" \"tmp.txt\"")}else{
		system("rm  \"tmp.txt\"")}
wd1 <- paste0(getwd(),"/",outpath)
wd2 <- paste0(getwd())
if(keep_samples == "keep"){
  print("……")}else{if(file.exists("samples") == T){
    setwd(wd1)
    system(sprintf("mv %s ../samples",r_out))
    setwd(wd2)
    system(sprintf("rm %s -rf",outpath))}else{
      system(sprintf("mkdir samples"))
      setwd(wd1)
      system(sprintf("mv %s ../samples",r_out))
      setwd(wd2)
      system(sprintf("rm %s -rf",outpath))}}
tmp <- read.table("all.abd.txt",header = T,sep = "\t")
rpkm <- pivot_wider(tmp,id_cols = gene,
              values_from = "RPKM",
              names_from = "sample_name")
rpkm[is.na(rpkm)] <- 0
write.table(rpkm,"rpkm.abd.txt",
                    sep = "\t",
                    quote = F,row.names = F)
com <- pivot_wider(tmp,id_cols = gene,
              values_from = "community_abd",
              names_from = "sample_name")
com[is.na(com)] <- 0
write.table(com,"community.abd.txt",
                    sep = "\t",
                    quote = F,row.names = F)
