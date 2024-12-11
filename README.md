# comts
## Comts: A pipeline for calculating single-copy genes’ community abundance in metagenome
## Introduction
Gene abundance in metagenome datasets is commonly represented by normalized values such as Reads Per Kilobase per Million Reads (RPKM), Fragments Per Kilobase per Million (FPKM), and Transcripts Per Million (TPM). However, the abundance of functional gene within the level of microbial community (AFG), representing proportion of the community that carries a specific metabolic function, remains underexplored and lacks a standardized methodology for estimation. In this study, we introduce Comts, a comprehensive framework for estimating AFG, and present a robust, user-friendly and efficient computational pipeline designed to calculate AFG from metagenomic sequencing data. The developed pipeline makes it accessible to researchers seeking to evaluate the metabolic capabilities of microbial communities, particularly with respect to single-copy functional genes.Gene abundance in metagenome datasets is commonly represented by normalized values such as Reads Per Kilobase per Million Reads (RPKM), Fragments Per Kilobase per Million (FPKM), and Transcripts Per Million (TPM). However, the abundance of functional gene within the level of microbial community (AFG), representing proportion of the community that carries a specific metabolic function, remains underexplored and lacks a standardized methodology for estimation. In this study, we introduce Comts, a comprehensive framework for estimating AFG, and present a robust, user-friendly and efficient computational pipeline designed to calculate AFG from metagenomic sequencing data. The developed pipeline makes it accessible to researchers seeking to evaluate the metabolic capabilities of microbial communities, particularly with respect to single-copy functional genes.
## The formula
AFG =  (RFG×100%)/MRUSCG

## Download and Installation
### The softwares listed below must have been installed before installation :robot:
> diamond  
> seqkit  
> fastp
### The R packages listed below must have been installed before installation  
> optparse
> dplyr
> tidyr
> data.table
> magrittr
> ggplot2
### Download throught `git clone` :
`git clone https://github.com/XiangZhouCAS/comts.git`
### Installation
1. `sh ./install.sh`
2. `source ~/.bashrc`
### DataBase
Ribo_14.dmnd  
hyddb_all.dmnd (Søndergaard, D., Pedersen, C. & Greening, C. HydDB: A web tool for hydrogenase classification and analysis. Sci Rep 6, 34212 (2016). https://doi.org/10.1038/srep34212)  
ter.dmnd.gz ([Hydrogen metabolism terminal enzyme's database providede by GreeningLab](https://github.com/GreeningLab/GreeningLab-database/blob/main/Original%20database%20(2020)))
## Usage
| Function | Description |
|-------|-------|
|`comts geneset`|To calculate community abandance of single copy function genes through GeneSet.|
|`comts custom`|To calculate community abandance of single copy function genes through custom database.|

### comts geneset
- `comts geneset ribo` To calculate RPKM abandance of 14 universal single copy ribosomal genes.  
```
comts geneset ribo -i sample1.1.fastq.gz -o sample1 -t 4 -s Ribo_14.dmnd
```
- `comts geneset res` To convert RPKM to community abandance of single copy function gene through GeneSet.  
```
comts geneset res -i geneset.rpkm.txt -r Ribo.rpkm.txt -o community.abd.txt
```
### comts custom
- `comts custom diy` To calculate single copy function enzyme gene community abandance by custom database.
```
comts custom diy -i sample1.1.fastq.gz -o sample1 -t 4 -d function_genes.dmnd -s Ribo_14.dmnd
```
- `comts custom ter` To calculate single copy terminal enzyme gene community abandance.  
```
comts custom ter -i sample1.1.fastq.gz -o sample1 -t 4 -d terminal_genes.dmnd -s Ribo_14.dmnd
```
- `comts custom hyd` To calculate single copy Hydrogenase community abandance.
```
comts custom hyd -i sample1.1.fastq.gz -o sample1 -t 4 -d hyddb.all.dmnd -s Ribo_14.dmnd -c hyd_id-name.script
```
