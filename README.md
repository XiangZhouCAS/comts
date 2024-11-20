# comts
## Comts: A simple pipeline for calculating single-copy genes’ community abundance in metagenome
## Introduction
Gene relative abundance algorithms such as RPKM (Reads Per Kilobase per Million mapped reads), FPKM (Fragments Per Kilobase of transcript per Million fragments mapped), and TPM (Transcripts Per Million) have been widely used in various scenarios, but in metagenome analysis, microbial communities are the object of analysis, the abundance of genes needs to correspond to the community to be more valuable for research, but the above algorithms cannot intuitively reflect to the percentage of community members of function genes in microbial community, i.e, community abundance. For this reason, we developed comts, a simple pipeline to calculate the community abundance of single-copy genes in metagenome.
## The formula
community abundance =  (Function Genes' RPKM×100%)/GeoMean( Total RPKM of universal single copy  Genes)

## Download and Installation
### The software listed below must have been installed before downloading :robot:
> diamond  
> seqkit  
> fastp
### Download throught `git clone` :
`git clone https://github.com/XiangZhouCAS/comts.git`
### Installation
1. `sh ./install.sh`
2. `source ~/.bashrc`
## Usage
| Function | Description |
|-------|-------|
|`comts geneset`|To calculate community abandance of single copy function genes through GeneSet.|
|`comts custom`|To calculate community abandance of single copy function genes through custom database.|

### comts geneset
- `comts geneset ribo` To calculate RPKM abandance of 14 universal single copy ribosomal genes.  
#### simple example
    comts geneset ribo -i sample1.1.fastq.gz -o sample1 -t 4 -s Ribo_14.dmnd
- `comts geneset res` To convert RPKM to community abandance of single copy function gene through GeneSet.  
#### simple example
    comts geneset res -i geneset.rpkm.txt -r Ribo.rpkm.txt -o community.abd.txt
### comts custom
- `comts custom diy` To calculate single copy function enzyme gene community abandance by custom database.
```comts custom diy -i sample1.1.fastq.gz -o sample1 -t 100 -d function_genes.dmnd -s Ribo_14.dmnd
- `comts custom ter` To calculate single copy terminal enzyme gene community abandance.  
    comts custom ter -i sample1.1.fastq.gz -o sample1 -t 100 -d terminal_genes.dmnd -s Ribo_14.dmnd
- `comts custom hyd` To calculate single copy Hydrogenase community abandance.
    comts custom hyd -i sample1.1.fastq.gz -o sample1 -t 100 -d hyddb.all.dmnd -s Ribo_14.dmnd -c hyd_id-name.script
