# 10/1 
library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)

#use pcadapt to read in uncompressed VCF data from class
vcf <- read.pcadapt("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf",
                    type = "vcf")

#compressed version of metadata 
#issue with reading in file so used the file in the class directory
#if wanted to use your own filted file need to provide path for vcf file that are uncompressed (no .gz extension) in the read,pcadapt() line and to the
# compressed version (ends in .gz) in the read.vcfR() line
# Be careful to only save the uncompressed version if you make one outside your github repo, because these
# are big and would cause github to not accept your push


#compressed version to subset metadata file 

vcfR <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf.gz")

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(vcfR@gt[,-1]),]


pcadapt.pca <- pcadapt(vcf,
                       K=2,
                       method = "componentwise",
                       min.maf = 0.01,# any loci less than 1% it will remove 
                       LD.clumping = list(size=500, thr=0.2)) #allows us to get function 2 parameters in size BP and thershold in correlation coeff.
#k is th enumebr of PC axes you want to test for selection 
#selection 'componentwise' returns a separate test for selection on each PC axis
# if not it will return a global selection test across all PC axes 
#min.maf is minor allele frequncey value to use when testing. SNPs below the value wont be counted
#ld.clumping removes loci in strong LD when intally fitting the PCs then tests fro slection using all loci 


#componentwise helps to get seperate test for each axis 
#. LD clumping helps to neutralize Null hypothesis so that it can get accruate data on each snp 
#allows us to fit inital PC axis and get test for each snp in data 
#threshold of .2 

summary(pcadapt.pca)

#Length Class  Mode   
#scores           1190  -none- numeric
#singular.values     2  -none- numeric
#loadings        31184  -none- numeric. how much each SNP is correlated with each axis. 2x the number of SNP for pc1 and pc2 
#zscores         31184  -none- numeric
#af              15592  -none- numeric. allele frequency 
#maf             15592  -none- numeric. allele frequency 
#chi2.stat       31184  -none- numeric. test to generate 
#stat            31184  -none- numeric
#gif                 2  -none- numeric
#pvalues         31184  -none- numeric. 
#pass             9987  -none- numeric. how many loci made it through the low frequncy filter 

#majority of variants are categorized as rare 

plot(pcadapt.pca, options = "scores",
     pop=meta2$region,
     i=1, j=2, K=2)
#generate plot after filtering out snps less than 1%
#loci in tails of dist are higher on graph 

#put data inton qqman format so that we cna look at manhattan plot organized by chromosome 

#look at the data in file 
View(head(vcfR@fix))

vcfR.fix <- as.data.frame(vcfR@fix[,1:2]) # all rows but first 2 columns 
#take first 2 columns for data 


chr.main <- unique(vcfR.fix$CHROM)[1:8]
#take first 8 chromsomes 

chrnum <- as.data.frame(cbind(chr.main,seq(1,8,1)))
#get all chromosomes with number associated with it 

pcadapt.pca$pvalues # look at pvalues 

Pval <- pcadapt.pca$pvalues

pcadapt.MHplot <- cbind(vcfR.fix, Pval)
pcadapt.MHplot <- left_join(chrnum, pcadapt.MHplot, join_by(chr.main==CHROM))

pcadapt.MHplot <- pcadapt.MHplot %>%
  mutate(SNP = paste0(chr.main,"_", POS))

pcadapt.MHplot$V2 = as.numeric(pcadapt.MHplot$V2)
pcadapt.MHplot$POS = as.numeric(pcadapt.MHplot$POS)
pcadapt.MHplot$pPC1 = as.numeric(pcadapt.MHplot[,4])
pcadapt.MHplot$pPC2 = as.numeric(pcadapt.MHplot[,5])

#eliminate NA so r doesnt crash 

pcadapt.MHplot <- pcadapt.MHplot %>% drop_na(pPC1)

#PC1
manhattan(pcadapt.MHplot,
          chr = "V2",
          bp = "POS",
          p = "pPC1",
          col = c("blue4", "orange"),
          logP=T,
          ylab = "-log10 p-value",
          genomewideline = F,
          main = "PCAadapt genome scane for selection (PC1)")


#PC2
manhattan(pcadapt.MHplot,
          chr = "V2",
          bp = "POS",
          p = "pPC2",
          col = c("blue4", "orange"),
          logP=T,
          ylab = "-log10 p-value",
          genomewideline = F,
          main = "PCAadapt genome scane for selection (PC2)")

#get the PC1 outliers that are less than .1%
#look where they are in theh refense genome and what genes are around them
View(pcadapt.MHplot %>%
       filter(pPC1<quantile(pcadapt.MHplot$pPC1, 0.001)) %>%
     select(chr.main,POS, pPC1))



