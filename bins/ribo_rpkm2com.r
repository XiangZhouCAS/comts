if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--input_geneset","-i"),type = "character",default = F,
              help = "please set the directory of geneset RPKM abundance file"),
  make_option(c("--input_ribo_rpkm","-r"),type = "character",default = F,
	      help = "please set the directory of ribo_rpkm.txt"),
  make_option(c("--output","-o"),type = "character",default = F,
	      help = "Please set the directory of output"))
opt_parser = OptionParser(
  usage = "usage: comts rpkm com [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL ,
  description = "This Script is to convert rpkm to community abandance of function gene.")
opt <- parse_args(opt_parser)
geneset <- opt$input_geneset
rpkm <- opt$input_ribo_rpkm
output <- opt$output
all_data <- read.table(rpkm,header = T,sep = "\t")
geneset <- read.table(geneset,header = T,
                  sep = "\t",check.names = F)
if(!require(dplyr,quietly = TRUE)){
   install.packages("dplyr")
   library(dplyr)}else{
   library(dplyr)}
if(!require(tidyr,quietly = TRUE)){
   install.packages("tidyr")
   library(tidyr)}else{
   library(tidyr)}
all_data <- all_data%>%
  group_by(sample_name)%>%
  mutate(geomean = exp(mean(log(T_RPKM))))
all_data <- unique(all_data[,c(1,4)])
all_data <- all_data[match(colnames(geneset[,-1]),
		           all_data$sample_name),]
for(i in 1:nrow(all_data)){
  ribo_rpkm <- as.numeric(all_data[i,2])
  geneset[,i+1] <- geneset[,i+1]*100/ribo_rpkm}
df <- data.frame(lapply(geneset[,-1],
                        function(x)ifelse(x > 100, 100, x)),check.names = F)
res <- basename(output)
df$Entry <- geneset$Entry
df <- df[,c(ncol(df),1:ncol(df)-1)]
write.table(df,res,sep = "\t",quote = F,row.names = F)
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplt2)
}else{library(ggplot2)}
samplept <- df%>%
        pivot_longer(cols = !gene,values_to = "val",names_to = "sample")
lv <- colnames(samplept)[-1]
lv <- as.factor(lv)
samplept$sample <- factor(samplept$sample,levels = lv)
hp <- ggplot(samplept,aes(sample,gene))+
  geom_tile(aes(fill = val),color = "grey50")+
  scale_fill_gradientn(colours = c("white","#C12554","#31357F"),limits = c(0,100),
                       name = "Community (%)")+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15,color = "black"),
        legend.position = "bottom",
        legend.title = element_text(vjust = 0.75),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.ticks.y = element_line(color = "black"))
w <- length(unique(samplept$sample))
h <- length(unique(samplept$gene))
ggsave(plot = hp,"Heatmap_community_abd.pdf",width = w,height = h/4)
