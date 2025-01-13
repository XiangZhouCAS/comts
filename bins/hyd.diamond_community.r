if (!require(optparse,quietly = TRUE)) {
  install.packages("optparse")
  library(optparse)
} else {
  library(optparse)}
option_list <- list(
  make_option(c("--input_reads","-i"),type = "character",default = F,
              help = "Please set the reads to include only the forward reads if the data is paired-end (PE)"),
  make_option(c("--result","-o"),type = "character",default = F,
              help = "Please specify the name for the result file."),
  make_option(c("--threads","-t"),type = "numeric",default = 1,
              help = "Set the number of CPU threads, with the default value being 1."),
  make_option(c("--diamond_db","-d"),type = "character",default = F,
              help = "Please set the directory of diamond database (e.g.'hyddb.all.dmnd')."),
  make_option(c("--script","-c"),type = "character",default = F,
	      help = "Please set the directory of hyd_id-name.script"),
  make_option(c("--singleM","-s"),type = "character",default = F,
              help = "Please specify the directory for SingleM's universal single-copy ribosomal genes database (e.g., 'Ribo_14.dmnd'). If you have already calculated the RPKM for these genes, you may instead specify the directory for the results (e.g., 'sample_name.singleM.hits.txt') to skip this step."),
  make_option(c("--run_fastp","-f"),type = "character",default = "run",
              help = "If you have already filtered the reads (length > 140), you can set the value to 0 to skip running fastp. The default is to run fastp."),
  make_option(c("--run_seqkit","-k"),type = "character",default = "run",
	      help = "If you have already counted the total number of reads using seqkit, you can specify the directory of the seqkit results (e.g., 'sample_name.all.reads.txt') to skip running seqkit. By default, seqkit will be executed."),
  make_option(c("--keep_samples","-e"),type = "character",default = "keep",
              help = "By default, the temporary results will be deleted unless the 'keep' option is specified."))
opt_parser = OptionParser(
  usage = "usage: comts custom hyd [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL ,
  description = "This Script is to calculate Hydrogenase community abandance.")
opt <- parse_args(opt_parser)
input_reads <- opt$input_reads
outpath <- basename(opt$result)
threads <- opt$threads
fastp_output <- paste0(outpath,".filtered.fq.gz")
diamond_db <- opt$diamond_db
singleM <- opt$singleM
script <- opt$script
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
diamond_singleM <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                           singleM, fastp_output, singleM_out, threads)
seqkit <- sprintf("seqkit stat %s > %s",
                  fastp_output,seqkit_out)
diamond2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle pident qcovhsp --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
		   diamond_db, input_reads, diamond_out, threads)
diamond_singleM2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
			    singleM, input_reads, singleM_out, threads)
seqkit2 <- sprintf("seqkit stat %s > %s",
		   input_reads,seqkit_out)
if(run_fastp == "run"){
	print("fastp is running.")
	system(fastp)
        print("fastp is completed.")}else{
		print("Not run the fastp because you set the directory of filtered reads to skip it.")}
if(run_fastp == "run"){
	print("diamond is Running (functional genes).")
	system(diamond)
        print("diamond is completed (functional genes).")}else{
		print("diamond is Running (function genes).")
		system(diamond2)
		print("diamond is completed (function genes).")}
if(run_fastp == "run"){
  if(tmp$V1 != "text"){
    print("diamond is Running (USCGs).")
    system(diamond_singleM)
    print("diamond is completed (USCGs).")
  }else{print("Not count the RPKM of singleM marker genes.")}
  }else{if(tmp$V1 != "text" ){
  print("diamond is Running (USCGs).")
  system(diamond_singleM2)
  print("diamond is completed (USCGs).")
    }else{print("Not count the RPKM of singleM marker genes.")}}
if(run_fastp == "run"){
  if(run_seqkit == "run"){
    print("seqkit is Running.")
    system(seqkit)
    print("seqkit is completed.")}else{
      print("Not run the seqkit because you set the directory of seqkit result")}
  }else{if(run_seqkit == "run"){
  print("seqkit is Running.")
  system(seqkit2)
  print("seqkit is completed.")
}else{print("Not run the seqkit because you set the directory of seqkit result")}}
if (!require(magrittr)) {
  install.packages("magrittr")
  library(magrittr)
} else {
  library(magrittr)}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)}else{
  library(dplyr)}
if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)}else{
  library(tidyr)}
if(!require(data.table)){
  install.packages("data.table")
  library(data.table)}else{
  library(data.table)}
if(run_seqkit == "run"){
  p <- read.table(seqkit_out,header = T)
  }else{p <- read.table(run_seqkit,header = T)}
d1 <- fread(diamond_out,sep = "\t")
colnames(d1)[c(1:4)] <- c("slen","id","pident","qcovhsp")
d1 <- d1%>%
  group_by_all()%>%
  count()
d1$all_reads_nums <- p$num_seqs
d1$all_reads_nums <- gsub(",","",d1$all_reads_nums)
d1$all_reads_nums <- as.numeric(d1$all_reads_nums)
d1$tmp_RPKM <- d1$n*10^9/(d1$slen*d1$all_reads_nums*3)
d1$sample_name <- outpath
name <- read.table(script,sep = "\t",header = T)
d1$id <- gsub(" .*","",d1$id)
d1 <- left_join(d1,name,by = "id")
d12 <- as.data.frame(d1[grepl("NiFe4", d1$gene),])
d13 <- as.data.frame(d1[grepl("FeFe", d1$gene),])
d12 <- rbind(d12,d13)
d11 <- filter(d1,!id %in% d12$id)%>%
	as.data.frame()
d12 <- filter(d12,pident >= 60 & qcovhsp >= 80)
d11 <- filter(d11,pident >= 50 & qcovhsp >= 80)
d1 <- rbind(d11,d12)
d1 <- d1[,c(8,9,7)]
d1 <- d1%>%
  group_by(gene)%>%
  mutate(RPKM = sum(tmp_RPKM))
d1 <- unique(d1[,-3])
if(tmp$V1 != "text"){
  d2 <- read.table(singleM_out,sep = "\t")
  }else{
  d2 <- read.table(singleM,sep = "\t")
  }
colnames(d2)[c(1,2,3,4)] <- c("slen","sseq_id","qcovhsp","bitscore")
d2 <- filter(d2,qcovhsp >= 80 & bitscore >= 40)
d2 <- d2[,-c(3,4)]
d2 <- d2%>%
  group_by_all()%>%
  count()
d2$all_reads_nums <- p$num_seqs
d2$all_reads_nums <- gsub(",","",d2$all_reads_nums)
d2$all_reads_nums <- as.numeric(d2$all_reads_nums)
d2$RPKM <- d2$n*10^9/(d2$slen*d2$all_reads_nums*3)
d2$sample_name <- outpath
d2 <- d2[,c(6,2,1,3,4,5)]
d2$sseq_id <- gsub("-.*","",d2$sseq_id)
d2 <- d2%>%
  group_by(sseq_id)%>%
  mutate(T_RPKM = sum(RPKM))
d2 <- unique(d2[,c(1,2,7)])
s_out <- paste0(outpath,".singleM.hits.comb.txt")
write.table(d2,s_out,quote = FALSE,sep = "\t",row.names = F)
geo <- exp(mean(log(d2$T_RPKM)))
d1$AFG <- d1$RPKM/geo*100
d1$AFG <- ifelse(d1$AFG > 100,100,d1$AFG)
r_out <- paste0(outpath,".all.abd.txt")
write.table(d1,r_out,quote = FALSE,sep = "\t",row.names = F)
system(sprintf("mkdir %s",outpath))
if(run_fastp == "run"){
	system(sprintf("mv %s %s",fastp_output,outpath))}else{
		print("……")}
system(sprintf("mv %s %s",diamond_out,outpath))
if(tmp$V1 != "text"){
  system(sprintf("mv %s %s",singleM_out,outpath))}else{
	  print("There is no need to move singleM result to outpath")}
if(run_seqkit == "run"){
  system(sprintf("mv %s %s",seqkit_out,outpath))
  }else{
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
              values_from = "AFG",
              names_from = "sample_name")
com[is.na(com)] <- 0
com <- com[!duplicated(com[,c(1,2)]),]
write.table(com,"AFG.abd.txt",
                    sep = "\t",
                    quote = F,row.names = F)
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplt2)
}else{library(ggplot2)}
lv <- colnames(com)[-1]
lv <- as.factor(lv)
samplept <- com%>%
        pivot_longer(cols = !gene,values_to = "val",names_to = "sample")
samplept$sample <- factor(samplept$sample,levels = lv)
hp <- ggplot(samplept,aes(sample,gene))+
  geom_tile(aes(fill = log2(val+1)),color = "grey50")+
  scale_fill_gradientn(colours = c("white","#C12554","#31357F"),
                       name = expression("AFG [ Log"["2"]*"(%+1) ]"))+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15,color = "black"),
        legend.position = "bottom",
        legend.title = element_text(vjust = 0.75),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.ticks.y = element_line(color = "black"))
w <- length(unique(samplept$sample))
h <- length(unique(samplept$gene))
ggsave(plot = hp,"Heatmap_AFG.pdf",width = w,height = h/4)
print("All Completed!")
