


library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

#helps with graphing 
options(bitmapType = "cairo")

#set working directory 
setwd("~/Projects/eco_genomics/population_genomics/")
vcf <-read.vcfR("outputs/vcf_final.filtered.vcf.gz")
#15380 variants 

#we need to thin the SNPs for LD before we run PCA and admixture 
#to satisfy the assumptions of independence among loci

# make new file to thin the distance to elimate any snp if its within 500 bp 
#thin to keep only 1 SNP per 500 bp 

vcf.thin <- distance_thin(vcf, min.distance = 500)
#3637 out of 15378 input snps were not located wihtin 500 base pairs of another SNP
#most snps are close or stacked close together bc lots of data was lost 

#read meta data 
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)
#629 8

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]) ,]

dim(meta2)
#593 8 

#have to uncompress vcf file so that the data within can be read. 
#we have to save the thinned vcf to file so that it can be extracted and converted to a different format geno needed by the 
#LEA program

write.vcf(vcf.thin, "outputs/vcf.filtered.thinned.vcf.gz")
 

# We now need to uncompress this vcf file, but when we do it'll be too big for our
# repos (github caps individual files at something ridiculous like 25 Mb).
# So, we'll "hide" the uncompressed vcf file in our home directory (~/) but outside 
# of our repo...don't forget where it lives!  And if you make a mistake and save it within your
# repo, just be sure to move it outside your repo before you do a git commit > push.
system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

#convert file from uncompressed vcf to geno format 
geno <- vcf2geno(input.file = "/users/c/l/clorourk/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.geno") #(scale = TRUE) didnt work 



CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale =TRUE) 

#used to reload outputs from previous command without running previous 
CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)
#593 indv
#3634 loci 

#pc1 is the frist and highest dot 
#most genetic info is going to be contained in pc1 than pc2 and cascading down in info
plot(CentPCA)
#ggplot(data, aes(x = x_variable, y = y_variable)) +
  #geom_point(aes(colour = my_colour_vector, shape = my_shape_vector))
ggplot(data, aes(x = x_variable, y = y_variable, colour = group_variable, shape = group_variable)) +
  geom_point()


#plotting pca plot 
ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
       geom_point(alpha=1) + 
  labs(title = "Centaurea Genetic PCA", x = "PC1", y="PC2", color="Region", shape = "Continent")+
  xlim(-10,10)+
  ylim(-10,10)
#
# add limits 
#shows that there is a lot of admixture btw each species in the middle of the plot or that there
#is very low divergance btw these species and have similar genetic makeup 
#CEU is the most divergent of the European species 
# acific northwest most non similar to all other samples and potentially have evolved 
#in allele frequinces that is much different than the other species 

show(CentPCA)

#saving plot in demensions 
ggsave("figures/CentPCA_PC2vPC3.pdf", width = 6, height = 6, units = "in" )

# CentPCA$eigenvectors[1:5,1:2]. put into console. tells the vectors assoicated with each snp 
#the vector in seach dicetor of x and y or sn1 or sn2 is used to calculate eigenvalue for each snp
#the sum of all eigenvalue defines the whole data set


# 10/1 coding 
#run admixture analysis and create plots 
#for admixture going to use the LEA R package 
# fucntion inside LEA is called 'snmf'
# runs admixture like analysis bc the algorithm takes a long time to analyse data 
#snmf is quicker and is a comparable analysis 
#stocastisity of the solution that the algorithm generates so if you run it multiple
#times you could get a different outcome 

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
     ) #plots cross entropy score we can use for selecting k values that fit our data 
#plot of result of cross validation exp. 
#cross entropy is value of error 
#if data did not level off indicatign a better fit 
#k=1 to k=2 there is large area of crpss entropy showing that k=2 is better fit 
#k=2 to k=3 is a little jump but not as big as k=1 indicating a better fit 
# want to get a range of K values where elbow lies which is btw 4-6 in this data set 
#lower entropy score is better as it indicates how well the data fits together 
#look for range of values that bracket where the cross entropy begins to level out 

plot(CentPCA$eigenvalues[1:10],
     ylab = "Eigenvalues",
     xlab = "Number of PCs",
     col = "blue4",
     main = "PCA")

#take par setting off the two panel view 
dev.off()

#this is so we can change value of K so we dont have to rewtire code

myK = 5 

#calculate the cross-entropy (=model fit; lower values are better) for all reps
#then determine which rep has the lowest score so that we can use for plotting

CE = cross.entropy(CentAdmix, K = myK)
#from console 
#run 1 0.2609008
#run 2 0.2595418
#run 3 0.2636343

best = which.min(CE)
#best run is run 2 

myKQ = Q(CentAdmix, K=myK, run = best)
#get ancestory coefficent (Q score)


#bind columns together from the same rows
#NEEED TO HAVE ROWS IN SAME ORDER TO USE 
#cbind to metadata 

myKQmeta = cbind(myKQ, meta2)

#set colors 

my.colors = c("blue4","gold","tomato","lightblue","darkgreen")

#use tidyverse to sort matrix based on metadata to help with the clustering so that we can change how its clustered like by region  

#group by continent then in each continent group by region then by population 

#sort the dataset by features of intrest in the metadata before plotting 
#continent is the first group and is sorted by region and population within continent 

myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region, pop, .by_group = TRUE) # sort by region or population in region 




#make ancestory plot and save as a pdf to my figures
#pdf("figures/Admixture_K5.pdf", width = 10, height = 5)
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

#create multiple plots by varying myK values and running the code below it 

#PNW has much more different and distinct make up than the other regions 
#NE region is noisy but has very similar makeup over entire region 
#CEU has some separate area than the other regions 
#middle of plot is very admix which is the same as the PCA plot that has large clustering in the middle 
#overall not strong clustering 
#could suggest that there is a lot of gene flow that would prevent certain areas 
#from being distinct genetically 




