library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType = "cairo")

setwd("~/Projects/eco_genomics/population_genomics/")
vcf <-read.vcfR("outputs/vcf_final.filtered.vcf.gz")
#15380 variants 

#we need to thin the SNPs for LD before we run PCA and admixture 
#to satisfy the assumptions of independence among loci

# make new file to thin the distance to elimate any snp if its within 500 bp 
vcf.thin <- distance_thin(vcf, min.distance = 500)

#most snps are close or stacked close together bc lots of data was lost 

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)
#629 8

meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[,-1]) ,]

dim(meta2)
#593 8 

write.vcf(vcf.thin, "outputs/vcf.filtered.thinned.vcf.gz")
#hide uncompressed vcf file bc its too big for github so hibe outside repo 

#unzip file to view contents inside
#save to home directory bc file is too large 
system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")


geno <- vcf2geno(input.file = "/gpfs1/Home/clorourk/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.geno", scale = TRUE)






