library(vcfR)
library(SNPfiltR)
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

#file of the pre filtered settings
vcf.filt <- read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.homework.vcf.gz")

#iorginal cuttoff was .75
#second cuttoff was .50 

#filter out by indv. missingness -- higher cutoff is more stringent 
vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap = meta2,
                                      cutoff=0.5)
#171 samples at .5 were cutout

#filtering data to get rid of monomrophic or multi-allelic sites 
vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5)

#assess clustering without MAC cutoff
#miss <- assess_missing_data_tsne(vcf.filt.indSNPMiss,
#popmap = meta,
#clustering = FALSE)

DP2 <- extract.gt(vcf.filt.indSNPMiss,
                  element = "DP",
                  as.numeric=T)

#thin data through LD 
vcf.filter.indSNPMiss.thin <- distance_thin(vcf.filt.indMiss,
                                            min.distance = 500)
#vcf.filter.indSNPMiss.thin

#vcf.filter.indSNPMiss.thin


#second file
library(tidyverse)
library(qqman)

X11.options(type="cairo")


#read in meta info on popualtion of orgin what region the file is in 
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(vcf.filter.indSNPMiss.thin@gt[,-1]),]



vcf.div <- genetic_diff(vcf.filter.indSNPMiss.thin,
                        pops=as.factor(meta2$region),
                        method="nei")

chr.main <- unique(vcf.div$CHROM)[1:8]


#number chrome
chrnun <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

#merge plots together of the numeric chrome and the full diversity results

vcf.div.MHplot <- left_join(chrnun, vcf.div, join_by(chr.main==CHROM))

#filtering gst values less than 0 and create new snp column that concatenates the 
#chrome ID with the bp postion 
vcf.div.MHplot <- vcf.div.MHplot %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main,"_", POS))


#have to make chromsome into values so that it can be plotted 
vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)

vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)

manhattan(vcf.div.MHplot,
          chr="V2" ,
          bp="POS",
          p="Gst",
          col=c("blue4", "orange3"),
          logp=F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999)
)


write.csv(vcf.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Gentic_Diff_byRegion75.csv",
          quote=F,
          row.names = F)

vcf.div.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, color = name)) + 
  geom_histogram(position = "identity", alpha = 0.5, bins=50) + 
  labs(title = "Genome Wide expected Heterozygosity (Hs)", fill = "Regions", 
       x = "Gene Diversity within Regions", y = "Counts of SNPs")

ggsave("Histogram_Genome_Diversity_byRegion.75.pdf",
       path = "~/Projects/eco_genomics/population_genomics/Homework/")





#summarise(avg=mean(value), SD=sd(value))


vcf.div.MHplot %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>%   
  group_by(name) %>%   
  summarise(avg=mean(value))


#Hs_table1 <- vcf.div.MHplot %>%   
 # as_tibble() %>%  
  #pivot_longer(c(4:9)) %>%  
  #group_by(name) %>%   
  #filter(value == 0) %>%   
  #summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n()) 

#view(Hs_table1)
#write.csv(Hs_table, "population_genomics/outputs/Hs_table_noZeros.csv",
                                                                         #quote=F,           
                                                                         #row.names=F
Hs_table2 <- vcf.div.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value==0) %>%
  summarise(avg=mean(value), StdDev_Hs=sd(value), N_Hs=n())

view(Hs_table2)


#filter(value=0 & value<0.25) %>%

#labels genome wide Hs 
# larger N has a larger avg snp 
#filter used get rid of bins in the table 
#first bin in plot is polymorphic within the entire data set bc they would not have survived the fiiltering process otherwise 

write.csv(Hs_table, "~/Projects/eco_genomics/population_genomics/Homework/Hs_Table_noZeros5.csv",
          quote = F,
          row.names = F)

#view(vcf.filter.indSNPMiss.thin)
#file 3 

library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

#helps with graphing 
options(bitmapType = "cairo")


#vcf <-read.vcfR("~/Projects/eco_genomics/population_genomics/Homework/vcf_final.filtered.5.vcf.gz")
#15380 variants 

#we need to thin the SNPs for LD before we run PCA and admixture 
#to satisfy the assumptions of independence among loci


#vcf.thin <- distance_thin(vcf.filter.indSNPMiss.thin, min.distance = 500)

#.5 led to 3937 out of 16741 input snps were not located within 500 bp

#read meta data 
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

#dim(meta)
#629 8

meta2 <- meta[meta$id %in% colnames(vcf.filter.indSNPMiss.thin@gt[,-1]) ,]

#dim(meta2)
#593 8 
setwd("~/Projects/eco_genomics/population_genomics")
vcf.thin <- vcf.filter.indSNPMiss.thin

write.vcf(vcf.thin,"outputs/vcf_final.filtered.homework.vcf.gz")
#view(meta2)


system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.homework.vcf.gz > ~/vcf_final.filtered.homework.vcf")

#convert file from uncompressed vcf to geno format 
geno <- vcf2geno(input.file = "/users/c/l/clorourk/vcf_final.filtered.homework.vcf",
                 output.file = "Homework/vcf_final.filtered.homework.geno") #(scale = TRUE) didnt work 


#"/users/c/l/clorourk/vcf_final.filtered.homework.vcf
CentPCA <- LEA::pca("Homework/vcf_final.filtered.homework.geno", scale =TRUE) 

#used to reload outputs from previous command without running previous 
CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)
#593 indv
#3634 loci 


#plot(CentPCA)

#plotting pca plot 
#ggplot(as.data.frame(CentPCA$projections),
#aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
#geom_point(alpha=1) + 
#labs(title = "Centaurea Genetic PCA", x = "PC2", y="PC3", color="Region", shape = "Continent")
#xlim(-10,10)+
#ylim(-10,10)





#saving plot in demensions 
#ggsave("figures/CentPCA_PC2vPC3.pdf", width = 6, height = 6, units = "in" )

# CentPCA$eigenvectors[1:5,1:2]. put into console. tells the vectors assoicated with each snp 
#the vector in seach dicetor of x and y or sn1 or sn2 is used to calculate eigenvalue for each snp

CentAdmix <- snmf("outputs/vcf_final.filtered.thinned.geno", 
                  K = 1:10, #computational cost as the model will run longer with greater K. fit of data will be better with greater K than level off but different for different data 
                  entropy = T, #cross validation of data 
                  repetitions = 3, #repeat 3 times to cross validate 
                  project = "new") #if youre adding to analysis later you could choose project = "continue" this would allow you to add more replicates of already completed runs 

#compare the evidence for differnet levels of K using the cross entropy from snmf and screenplot from PCA

par(mfrow=c(2,1)) #multi panel plot 
plot(CentAdmix, # plots the cross entropy score we can use for selecting models with K values that fit our data 
     col = "blue4",
     main = "SNMF"
) 

plot(CentPCA$eigenvalues[1:10],
     ylab = "Eigenvalues",
     xlab = "Number of PCs",
     col = "blue4",
     main = "PCA")

#take par setting off the two panel view 
dev.off()

#this is so we can change value of K so we dont have to rewtire code

myK = 5 


CE = cross.entropy(CentAdmix, K = myK)

best = which.min(CE)
#run is run 2 

myKQ = Q(CentAdmix, K=myK, run = best)
#get ancestory coefficent (Q score)



myKQmeta = cbind(myKQ, meta2)

#set colors 

my.colors = c("blue4","gold","tomato","lightblue","darkgreen")


myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region, pop, .by_group = TRUE) # sort by region or population in region 




#plot ancestor from 1 thorugh whatever K is selected 
barplot(as.matrix(t(myKQmeta[ , 1:myK])),
        border = NA, # no borders on bars
        space = 0, #no space btw bars
        col = my.colors[1:myK],
        xlab = "Geographic Region",
        ylab = "Ancestor Proportions",
        main = paste0("Ancestry Matrix K=",myK))


#add labels of region on x axis 
axis(1,
     at = 1:length(meta2$region),
     labels = myKQmeta$region,
     tick=F,
     cex.axis = 0.5,
     las = 3)
