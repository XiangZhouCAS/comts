if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--abundance_file","-i"),type = "character",default = F,
              help = "Please set the abundance file."),
  make_option(c("--metadata","-m"),type = "character",default = F,
              help = "Please set the metadata, col1 is sample name, col2 is group name."),
  make_option(c("--stat_method","-s"),type = "character",default = "kw",
              help = "Please set the stat method,default is Kruskal-Wallis test, Wilcoxon(do not set kw) is also allowed."),
  make_option(c("--output","-o"),type = "character",default = F,
              help = "Please set the output directory."),
  make_option(c("--padj","-p"),type = "character",default = "none",
              help = "You can choose the p-value adjust method,including FDR, BH and none, default is none."),
  make_option(c("--gene_list","-g"),type = "character",default = "all",
              help = "You can set the list of genes that you need to do differential analysis, and the default is all genes"))
opt_parser = OptionParser(
  usage = "usage: comts stat [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL ,
  description = "This Script is to calculate abundance difference.")
opt <- parse_args(opt_parser)
abundance_file <- opt$abundance_file
metadata <- opt$metadata
stat_method <- opt$stat_method
output <- opt$output
padj <- opt$padj
gene_list <- opt$gene_list
options(warn = -1)
if(!require("tidyr",quietly = T)){
  install.packages("tidyr")
  library(tidyr)
}else{library(tidyr)}
if(!require("dplyr",quietly = T)){
  install.packages("dplyr")
  library(dplyr)
}else{library(dplyr)}
if(!require("ggplot2",quietly = T)){
  install.packages("ggplot2")
  library(ggplot2)
}else{library(ggplot2)}
if(!require("ggnewscale",quietly = T)){
  install.packages("ggnewscale")
  library(ggnewscale)
}else{library(ggnewscale)}
if(!require("cowplot",quietly = T)){
  install.packages("cowplot")
  library(cowplot)
}else{library(cowplot)}
if(!require("ggpubr",quietly = T)){
  install.packages("ggpubr")
  library(ggpubr)
}else{library(ggpubr)}

ko <- read.table(abundance_file,header = T,sep = "\t")
colnames(ko)[1] <- "Entry"


if(gene_list != "all"){
  gene_list <- read.table(gene_list,header = T,sep = "\t")
  colnames(gene_list) <- "Entry"
}else{gene_list <- data.frame(Entry = ko$Entry)}
level1 <- gene_list$Entry

ko <- filter(ko,Entry %in% gene_list$Entry)
samplept <- ko

colnames(gene_list)[1] <- "Entry"
ko <- pivot_longer(ko,cols = !Entry,values_to = "val",
                   names_to = "sample")
ko$val <- ifelse(ko$val>100,100,ko$val)
sample <- read.table(metadata,header = T,sep = "\t")
colnames(sample)[c(1,2)] <- c("sample","group")
level2 <- sample$sample
ko <- left_join(ko,sample,by = "sample")

if(!require("rstatix",quietly = T)){
  install.packages("rstatix")
  library(rstatix)
}else{library(rstatix)}

if(padj != "none"){if(stat_method == "kw"){
  st <- ko%>%
    group_by(Entry)%>%
    kruskal_test(val~ group)%>%
    adjust_pvalue(method = padj)%>%
    add_significance()
}else{
  ko_out <- ko|>
    group_by(Entry)|>
    mutate(sum = sum(val))
  ko_out <- filter(ko_out,sum < 100*length(unique(ko$sample)))
  ko_out <- ko_out[,-5]
  st <- ko_out%>%
    group_by(Entry)%>%
    wilcox_test(val ~ group)%>%
    adjust_pvalue(method = padj)%>%
    add_significance()}}else{
      if(stat_method == "kw"){
        st <- ko%>%
          group_by(Entry)%>%
          kruskal_test(val ~ group)%>%
          add_significance()
      }else{
        ko_out <- ko|>
          group_by(Entry)|>
          mutate(sum = sum(val))
        ko_out <- filter(ko_out,sum < 100*length(unique(ko$sample)))
        ko_out <- ko_out[,-5]
        st <- ko_out%>%
          group_by(Entry)%>%
          wilcox_test(val ~ group)%>%
          add_significance()}}

ko <- ko%>%
  group_by(Entry,group)%>%
  mutate(values = mean(val),
         SEM = sd(val)/sqrt(length(val)))
coln1 <- unique(ko$group)[1]
coln2 <- unique(ko$group)[2]
ko$group <- ifelse(ko$group == coln1,"A","B")
ko <- unique(ko[,c(1,4:6)])
abd_long <- ko
matm <- ko[,c(1:3)]
mat_pre <- pivot_wider(ko,id_cols = Entry,
                       values_from = "values",
                       names_from = "group")
mat_pre$log2fc <- ifelse(mat_pre$B == 0,ifelse(mat_pre$A == 0,0,
                                               log2(mat_pre$A)),
                         log2(mat_pre$A/mat_pre$B))
group1 <- unique(ko$group)[1]
group2 <- unique(ko$group)[2]
ko1 <- filter(ko,group == group1)
ko2 <- filter(ko,group == group2)
colnames(ko1)[3] <- group1
ko1 <- ko1[,-c(2,4)]
colnames(ko2)[3] <- group2
ko2 <- ko2[,-2]

if(padj != "none"){if(stat_method == "kw"){
  ko <- left_join(ko1,ko2,by = "Entry")%>%
    left_join(st[,c(1,7,6,8,9)],by = "Entry")
}else{
  ko <- left_join(ko1,ko2,by = "Entry")%>%
    left_join(st[,c(1,8,9,10)],by = "Entry")
  ko$method <- "Wilcoxon"
  ko <- ko[,c(1:4,8,5:7)]}}else{
    if(stat_method == "kw"){
      ko <- left_join(ko1,ko2,by = "Entry")%>%
        left_join(st[,c(1,7,6,8)],by = "Entry")
    }else{
      ko <- left_join(ko1,ko2,by = "Entry")%>%
        left_join(st[,c(1,8,9)],by = "Entry")
      ko$method <- "Wilcoxon"
      ko <- ko[,c(1:4,7,5,6)]}}

ko_out <- ko
colnames(ko_out)[2] <- coln1
colnames(ko_out)[3] <- coln2
ko_out <- ko_out[match(level1,
                       ko_out$Entry),]
colnames(ko_out)[1] <- "Gene"
write.table(ko_out,output,sep = "\t",quote = F,row.names = F)

if(padj != "none"){
  (sig <- ko[,c(1:4,8)])}else{
    sig <- ko[,c(1:4,7)]}

if(padj == "none"){
  sig <- filter(sig,`p.signif` != "ns")
}else{sig <- filter(sig,`p.adj.signif`!= "ns")}

sig$x <- ifelse((sig$A+sig$SEM)>(sig$B+sig$SEM),
                sig$A+sig$SEM,sig$B+sig$SEM)
xmax <- max(sig$x)
sig$x2 <- ifelse(sig$x < 9,10,
                 ifelse(sig$x < 19,20,
                        ifelse(sig$x < 29,30,
                               ifelse(sig$x < 39,40,
                                      ifelse(sig$x < 49,50,
                                             ifelse(sig$x < 59,60,
                                                    ifelse(sig$x < 69,70,
                                                           ifelse(sig$x < 79,80,
                                                                  ifelse(sig$x < 89,90,
                                                                         ifelse(sig$x < 99,100,0))))))))))


sig <- sig[,c(1,7,5)]
abd_long <- left_join(abd_long,sig,by = "Entry")
abd_long$group <- ifelse(abd_long$group == "A",coln1,coln2)

abd_long$Entry <- factor(abd_long$Entry,levels = rev(gene_list$Entry))

p1h <- length(level1)
if(padj != "none"){p1 <- ggplot(abd_long,aes(values,Entry))+
  geom_errorbar(aes(xmin = values - SEM,
                    xmax = values + SEM,group = group),
                width = 0.2,
                position = position_dodge(0.9))+
  geom_bar(stat = "identity",aes(fill = group),
           position = position_dodge(0.9))+
  labs(x = "Community Abundance (%)",y = colnames(abundance_file)[1])+
  scale_fill_manual(values = c("#D8B646","#275AAC"))+
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0.01),
                     limits = c(0,100))+
  theme(legend.position = c(0.9,0.4),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.text.x = element_text(size = 15,color = "black"),
        axis.title = element_text(size = 15),
        legend.title = element_blank())+
  geom_text(label = abd_long$`p.adj.signif`,x = abd_long$x2,
            size = length(level1)*6/(3*p1h/5))}else{
              p1 <- ggplot(abd_long,aes(values,Entry))+
                geom_errorbar(aes(xmin = values - SEM,
                                  xmax = values + SEM,group = group),
                              width = 0.2,
                              position = position_dodge(0.9))+
                geom_bar(stat = "identity",aes(fill = group),
                         position = position_dodge(0.9))+
                labs(x = "Community Abundance (%)",y = colnames(abundance_file)[1])+
                scale_fill_manual(values = c("#D8B646","#275AAC"))+
                theme_classic()+
                scale_x_continuous(expand = c(0.01,0.01),
                                   limits = c(0,100))+
                theme(legend.position = c(0.9,0.4),
                      axis.text.y = element_text(size = 12,color = "black"),
                      axis.text.x = element_text(size = 15,color = "black"),
                      axis.title = element_text(size = 15),
                      legend.title = element_blank())+
                geom_text(label = abd_long$`p.signif`,x = abd_long$x2,
                          size = length(level1)*6/(3*p1h/5))}

ggsave(plot = p1,"gene_bar.pdf",width = 15,height = 3*p1h/5)

matfc <- mat_pre
matfc$group <- "A"
matfc <- left_join(matfc,sig[,c(1,3)],by = "Entry")
rangefc <- ifelse(abs(min(matfc$log2fc)) < abs(max(matfc$log2fc)),
                  -abs(max(matfc$log2fc)),-abs(min(matfc$log2fc)))
matfc$Entry <- substr(matfc$Entry,1,6)
matfc$Entry <- factor(matfc$Entry,levels = rev(gene_list$Entry))
if(padj != "none"){p2 <- ggplot(matfc)+
  geom_tile(aes(x = group,y = Entry,fill = log2fc),
            color = "grey",linewidth = 0.25)+
  scale_fill_gradient2(limits = c(rangefc,-rangefc),
                       low = "#275AAC",mid = "white",
                       high = "#D8B646")+
  theme(axis.title.y=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 15,color = "white"),
        legend.position = "bottom",
        legend.title = element_text(vjust = 0.75),
        axis.text.y = element_text(size = 12,color = "black"))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  xlab(NULL)+
  geom_text(aes(x = group,y = Entry,
                label = `p.adj.signif`),size = 8,vjust = 0.75)}else{
                  p2 <- ggplot(matfc)+
                    geom_tile(aes(x = group,y = Entry,fill = log2fc),
                              color = "grey",linewidth = 0.25)+
                    scale_fill_gradient2(limits = c(rangefc,-rangefc),
                                         low = "#275AAC",mid = "white",
                                         high = "#D8B646")+
                    theme(axis.title.y=element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.text.x = element_text(size = 15,color = "white"),
                          legend.position = "bottom",
                          legend.title = element_text(vjust = 0.75),
                          axis.text.y = element_text(size = 12,color = "black"))+
                    scale_x_discrete(expand = c(0,0))+
                    scale_y_discrete(expand = c(0,0))+
                    xlab(NULL)+
                    geom_text(aes(x = group,y = Entry,
                                  label = `p.signif`),size = 8,vjust = 0.75)}
matm1 <- filter(matm,group == "A")
matm2 <- filter(matm,group == "B")
triangle <- function(a, b, type = "up") {
  triangle_down <- function(a, b) {
    data.frame(x = c(0, 1, 1) + a, y = c(0, 0, 1) + b,
               group = paste0(a, "_", b), stringsAsFactors = FALSE)}
  triangle_up <- function(a, b) {
    data.frame(x = c(0, 0, 1) + a, y = c(0, 1, 1) + b,
               group = paste0(a, "_", b), stringsAsFactors = FALSE)
  }
  if (type == "up") {
    data <- do.call(rbind, lapply(1:b, function(i) {
      do.call(rbind, lapply(1:a, triangle_up, i))
    }))
  } else if (type == "down") {
    data <- do.call(rbind, lapply(1:b, function(i) {
      do.call(rbind, lapply(1:a, triangle_down, i))
    }))}
  return(data)}
updata <- triangle(1,length(matm1$Entry),"up")
downdata <- triangle(1,length(matm2$Entry),"down")
matm1 <- matm1[match(rev(level1),matm1$Entry),]
matm2 <- matm2[match(rev(level1),matm2$Entry),]
matm1$group <- paste0(1,"_",seq(1,nrow(matm1),1))
matm2$group <- paste0(1,"_",seq(1,nrow(matm2),1))
length <- nrow(matm1)
matm1$type <- "A"
matm2$type <- "B"
matm1 <- full_join(updata,matm1,by = "group")
matm2 <- full_join(downdata,matm2,by = "group")
matm1$values <- as.numeric(matm1$values)
matm2$values <- as.numeric(matm2$values)
ylabels <- matm1%>%filter(x==1)%>%
  distinct(Entry,.keep_all = T)%>%
  arrange(y)%>%
  pull(Entry)%>%
  as.character()
yellow  <- "#D8B646"
blue <- "#275AAC"
nake <- "white"
matm1 <- matm1[,c(4,1:3,5:6)]
matm2 <- matm2[,c(4,1:3,5:6)]
p3 <- ggplot(matm1,aes(x=x, y=y))+
  geom_polygon(data = matm1,aes(x=x, y=y,
                                group=group,fill=values),color="grey",
               linewidth = 0.25)+
  scale_fill_gradientn(colours = c(nake,yellow),name = "groupA",
                       limits = c(0,100)) +
  new_scale_fill()+
  geom_polygon(data = matm2,aes(x=x, y=y, group=group,
                                fill=values),color="grey",
               linewidth = 0.25) +
  scale_fill_gradientn(colours = c(nake,blue),name = "groupB",
                       limits = c(0,100))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits = c(1, NA), expand = c(0,0),
                     breaks = c(1:length)+0.5,labels=ylabels)+
  theme(axis.title.y=element_blank(),
        axis.text.x = element_text(size = 15,color = "white"),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(vjust = 0.75),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  xlab(NULL)
samplept <- pivot_longer(samplept,cols = !Entry,
                         values_to = "val",
                         names_to = "sample")
samplept$val <- ifelse(samplept$val>100,100,samplept$val)
samplept <- filter(samplept,Entry %in% level1)
samplept$Entry <- factor(samplept$Entry,levels = rev(level1))
samplept$sample <- factor(samplept$sample,levels = sample$sample)
p4 <- ggplot(samplept)+
  geom_tile(aes(sample,Entry,fill = val),color = "grey50")+
  scale_fill_gradientn(colours = c("white","#C12554","#31357F"),limits = c(0,100),
                       name = "Community (%)")+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 15,color = "black"),
        legend.position = "bottom",
        legend.title = element_text(vjust = 0.75))
leg2 <- get_legend(p2)
leg3 <- get_legend(p3)
leg4 <- get_legend(p4)
p2 <- p2+
  theme(legend.position = "none")
p3 <- p3+
  theme(legend.position = "none")
p4 <- p4+
  theme(legend.position = "none")
w1 <- (length(unique(samplept$sample))+4)*0.8
w2 <- (length(unique(samplept$sample))+4)*0.2
h1 <- length(unique(samplept$Entry))+4
h2 <- length(unique(samplept$sample))
p <- plot_grid(plot_grid(p2,p3,p4,ncol = 3,
                         rel_widths = c(w2*0.6,w2*0.4,w1)),
               plot_grid(leg2,leg3,leg4,ncol = 3),nrow = 2,
               rel_heights = c(9,1))
ggsave(plot = p,"Heatmap.pdf",height = 10*h1/50,
       width = h2*15/24)
if(file.exists("Rplots.pdf") == T){
  system(sprintf("rm -rf Rplots.pdf"))
}else{print("Done!")}
