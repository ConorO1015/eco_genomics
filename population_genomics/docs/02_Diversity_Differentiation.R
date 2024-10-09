#estimating diversity and gentic differentiation in the filtered centaurea data 

library(vcfR)
library(tidyverse)
library(qqman)

#helps solve plotting issue 
X11.options(type="cairo")

#read in our vcf file from repo outputs/directory 

vcf <-read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

#read in meta info on popualtion of orgin what region the file is in 
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) # vcf files has 595 samples 
dim(meta)

#need to subset the 'meta' dataframe so that it has teh sample ids that are present in the 
#filtered vcf.file
# %in% used to match 2 objects such as id and meta with the sample ids in the vcf file 
#vcf file stores vames as columns in the @gt slot minus the the first column whcih has INFO


meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]

# to call in code go to console and type in colnames(vcf@gt[,-1])[1] depending on what you want ot see 

#check dimesnions
dim(meta2)


#calculate diversity stats using the gentci_diff function in vcfR
# could use region or other grouping factor depending on what were intrested in 
#nei is name of guy who made stat 
vcf.div <- genetic_diff(vcf,
                        pops=as.factor(meta2$region),
                        method="nei")

str(vcf.div)
#used to get info on the characters to use later 
#tells CHromosome positon, Hs for each region 
#unique(vcf.div$CHROM) to look at each chromosome 
# numeric label set to each chrome
#first find how many chromosomes are present and take the first 8 values to 
#you can use number 1 to n if you know the number of chromosomes 
#1. finding out which unique chrome are present in vcf 
#2. take the 8 chromosomes and not the smaller fragment scafolds
chr.main <- unique(vcf.div$CHROM)[1:8]

#seq is sequence of things. number things 1 - 8 every interval
#number chrome
chrnun <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))


#merge plots together of the numeric chrome and the full diversity results

vcf.div.MHplot <- left_join(chrnun, vcf.div, join_by(chr.main==CHROM))

#filtering gst values less than 0 and create new snp column that concatenates the 
#chrome ID with the bp postion 
vcf.div.MHplot <- vcf.div.MHplot %>%
                filter(Gst>0) %>%
                mutate(SNP=paste0(chr.main,"_", POS))

# %>% pipe operator 

str(vcf.div.MHplot)



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

#value above is the quantile of the data so .5 is the average. this shows that there is very low divergance btw popualitons 
#this can mean that population is very similar and recently diverged 
#selection process can effect the snps taht are outliers 
#suggest sets a line at the top 99% of genome or the most differentiated parts of the genome 
#visualizing 
#dot = snp in the figrue and displayed as the inbreeding coeffiecnt 
#shows very low differention 
#line shows 99.9% quantile 

#steve notes 
#the suggestiveline is suggestive that those regions of the genome above the are experincing 
#high degreee of population differentation. 
#could be caused by low gene flow into that protion of the genome(genes contributing
# to reproductive isolation btw species) of genomic regions experincing section for
#local adaptations casuing regions ot diverge from eachother 

#9/24 code
view(vcf.div.MHplot)
write.csv(vcf.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Gentic_Diff_byRegion.csv",
          quote=F,
          row.names = F)

#looking at diversity wihtin a group of Hs 
names(vcf.div.MHplot)
#postion 4 - 9 are where Hs values are stored 
#prints values of Hs for each region 
#looks at the column names for each region 

#look at each subpopulaiton at the same time we use tidyverse that has each Hs value in row
#put them into 1 long column 

#tidy operation use %>% pipe one step into the other 

#histogram of the overall Hs values over all loci and overaly results 
#of each region to compare 


vcf.div.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, color = name)) + 
  geom_histogram(position = "identity", alpha = 0.5, bins=50) + 
  labs(title = "Genome Wide expected Heterozygosity (Hs)", fill = "Regions", 
       x = "Gene Diversity within Regions", y = "Counts of SNPs")

#ggsave("Histogram_Genome_Diversity_byRegion.5.pdf",
      # path = "~/Projects/eco_genomics/population_genomics/figures/")

#alpha is somehwhat transparent for overalying plots 
#distribution of each region 
# max hs value is .5 
#shows that one allele is almost compelety fixed in the popualiton as tehre is large bars near 0.0 and .1
#therefore low Hs value 
#large spike shows that there is only one allele present in the group by 10000 count 

#summerize the results in terms of avg, SD, and sample size. In tidy you can use the summerize,
# group_by, and filter commands to get the comparisons we want. 
#run commmand without assigning the outputs to a new variable and it will print results in console 
#or can assign to new variables 
#report in table SD sample size mean ect

Hs_table <- vcf.div.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0 & value<0.25) %>%
  summarise(avg=mean(value), StdDev_Hs=sd(value), N_Hs=n())
#labels genome wide Hs 
# larger N has a larger avg snp 
#filter used get rid of bins in the table 
#first bin in plot is polymorphic within the entire data set bc they would not have survived the fiiltering process otherwise 

write.csv(Hs_table, "~/Projects/eco_genomics/population_genomics/outputs/Hs_Table_noZeros5.csv",
          quote = F,
          row.names = F)

#steve notes 
#code was customized in class based on our discussion of looking at Hs values across all loci 
#or just ones that had non-zero values (added filtering step)
#You can use any analysis guided by what you want to learn and what you think 
#will be useful to look at such as what proportion of loci are Hs=0 in each group 

