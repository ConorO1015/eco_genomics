
install.packages("ape")
install.packages("vcfR")
library(vcfR)
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
list.files()

#looking at files
list.files("variants/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
#tells you the infromation in the file 
vcf
#629 samples 
#21 chrome
#18233 varients 

head(vcf)
#head of VCF 
# fixed section contains info of each snit. tells the alelles 
# gt genotype section has matrix of rows and columns 

# dp = depth 
# ad is breakdown of reads for each indv. 
# filtering 


#associate reference genome with filtered vcf

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")

#chromosome 1 as example 
chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)
# telling it that its vcf use ghe sequenced dna data and annotated

plot(chr1)
#example plot 

#export pdf to plots in project 
#
pdf(file="~/Projects/eco_genomics/population_genomics/figures/Chromoplot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
#x limit for base pair to base pair 

#mask poor quality variants 
chr1_masked <- masker(chr1,
                      min_QUAL = 50,
                      min_DP = 1000,
                      max_DP = 10000,
                      min_MQ = 30
                        )

plot(chr1_masked)
chromoqc(chr1_masked, xlim=c(1e1, 1.1e8))

#process the chromR object with proc.chromR
#default window size is 1000bp
chr1_proc <- proc.chromR(chr1_masked, win.size = 1e5)
plot(chr1_proc)
chromoqc(chr1_proc, xlim=c(1e1, 1.1e8))

#how to filter a vcf file for minDP
# look at depth from all indv
DP <- extract.gt(vcf, element = "DP" , as.numeric = T, 
                 convertNA = T)
#50% have 6 reads or fewer 
quantile(DP)

DP[DP==0] <- NA
#set the 0 depth genos to NA
# to get accurate date on genome 

quantile(DP, na.rm = T)
#removes na when true 

#dimensions of object 
#ensure loci are in rows; smaples in columns 
dim(DP)
#18233, 629

#first 5 rows and 10 indv. 
DP[1:5, 1:10]


# looking at the mean DP per indv.
avgDP = colMeans(DP, na.rm = T)
summary(avgDP)
hist(avgDP, breaks=50)
#mean is 24x range from 2.8-171x
#ideally want avg of 15-20X/ind


#visulaize the matrix of depth and missingness 

#looking at missingness through heatmap funciton 
heatmap.bp(DP[1:1000,], rlabels=F, clabels=F)
# if issue you can subset heat map by asking for differnt peramiters
dev.off()

#set indv genotypes with DP<X to 'NA'

#vcf@gt[,1][is.na(DP)==TRUE] <- NA

#check to see % missing data --> 
vcf


#start filtering by depth 
library(SNPfiltR)
#helps visualize and export filtering steps 

#for help with function ?package in console 

#first filter raw data by looking at hisotgram then assigning depth
# if not 3 reads then will assign NA to value 

hard_filter(vcf)

#look at mean depth per indv. 
vcf.filt <- hard_filter(vcf, depth=3)
# coudl assign other dp values depedning on goals such as 5 
#vcf.filt to view new values 

#look at allele balance
vcf.filt <- filter_allele_balance(vcf.filt,
                                  min.ratio = 0.15,
                                  max.ratio = 0.85
                                  )
vcf.filt <- max_depth(vcf.filt,
                      maxdepth = 60)

#usually 2x mean depth 




#across all indv. what is the avg depth of reads 
max_depth(vcf.filt)
#rule of thumb is twice the average as filter. our average is 30 so our max depth filter is 60 

#filter out genotypes with >60 reads/SNPs
max_depth(vcf.filt, maxdepth = 60)

#look for missing data 
meta <- read.csv("metadata/meta4vcf.csv", header=T)

#brakets seperate rows and columns 
meta2 <- meta[,c(1,4)]

#head() to view 


names(meta2) <- c("id","pop")

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)



#filter out by indv. missingness -- higher cutoff is more stringent 
vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap = meta2,
                                      cutoff=0.85)

#subset popmap to only include retained indv
#based on graph going to filter indv above 75% missing data 

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



#heatmap.bp(DP2[1:5000,],
           #rlabels=F, clabels=F)

#thin data through LD 
vcf.filter.indSNPMiss.thin <- distance_thin(vcf.filt.indMiss,
                                          min.distance = 500)


vcf.filter.indSNPMiss.thin

#write out files created during filtering steps 
#saving to directory 
write.vcf(vcf.filt.indSNPMiss, 
          "~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

mydiff <- vcfR::genetic_diff(vcf.filt.indMiss,
                             pops=as.factor(meta$region),
                             method="nei")
