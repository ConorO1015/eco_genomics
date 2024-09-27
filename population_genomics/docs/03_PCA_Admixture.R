library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

#helps with graphing 
options(bitmapType = "cairo")

setwd("~/Projects/eco_genomics/population_genomics/")
vcf <-read.vcfR("outputs/vcf_final.filtered.vcf.gz")
#15380 variants 

#we need to thin the SNPs for LD before we run PCA and admixture 
#to satisfy the assumptions of independence among loci

# make new file to thin the distance to elimate any snp if its within 500 bp 
vcf.thin <- distance_thin(vcf, min.distance = 500)
#3637 out of 15378 input snps were not located wihtin 500 base pairs of another SNP
#most snps are close or stacked close together bc lots of data was lost 

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)
#629 8

meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[,-1]) ,]

dim(meta2)
#593 8 

#have to uncompress vcf file so that the data within can be read. 

write.vcf(vcf.thin, "outputs/vcf.filtered.thinned.vcf.gz")
 

#unzip file to view contents inside
#save to home directory bc file is too large 
system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

#have to hide the file in our home dr bc its to big for r 
geno <- vcf2geno(input.file = "/users/c/l/clorourk/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.geno") #(scale = TRUE) didnt work 


#convert format so can be stored in ecogenomics 
CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale =TRUE) 

#used to reload outputs from previous command without running previous 
CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)
#593 indv
#3634 loci 

#pc1 is the frist and highest dot 
#most genetic info is going to be contained in pc1 than pc2 and cascading down in info
plot(CentPCA)

#plotting pca plot 
ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
       geom_point(alpha=1) + 
  labs(title = "Centaurea Genetic PCA", x = "PC2", y="PC3", color="Region", shape = "Continent")
 #xlim(-10,10)+
  #ylim(-10,10)
#
# add limits 
#shows that there is a lot of admixture btw each species in the middle of the plot or that there
#is very low divergance btw these species and have similar genetic makeup 
#CEU is the most divergent of the European species 
# acific northwest most non similar to all other samples and potentially have evolved 
#in allele frequinces that is much different than the other species 

#saving plot in demensions 
ggsave("figures/CentPCA_PC2vPC3.pdf", width = 6, height = 6, units = "in" )

# CentPCA$eigenvectors[1:5,1:2]. put into console. tells the vectors assoicated with each snp 
#the vector in seach dicetor of x and y or sn1 or sn2 is used to calculate eigenvalue for each snp
#the sum of all eigenvalue defines the whole data set 








