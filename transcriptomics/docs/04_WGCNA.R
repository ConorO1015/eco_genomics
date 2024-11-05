#02 WGCNA
#script fot anaylzing and visulazing gene corroleation 

library(DESeq2)
library(ggplot2)
library(WGCNA); options(stringAsFactors = FALSE);
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)


citation("C")

options(bitmapType = "cairo")

setwd("~/Projects/eco_genomics/transcriptomics/")

# STEP1 

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt")
tail(countsTable)
dim(countsTable) #119438     21

countTableRound <- round(countsTable)
tail(countTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)
#DevTemp FinalTemp
#N1C3      D18      BASE
#N2C2      D22      BASE
#sN1L3     D18       A28
#N1C4      D18      BASE
#sN2H2     D22       A33
#N1L2      D18       A28

traitData = read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", 
                       header = T, row.names = 1)

#filter the matrix to just BASE data becasye those are the dtaa for which we have traits measured 
filtered_count_matrix_BASEonly <- countsTable[, conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE", ]
rounded_filtered_counts_matrix <- round(filtered_count_matrix_BASEonly)

# STEP 2: detecting outliers 
#detect outlier genes 
gsg <- goodSamplesGenes(t(rounded_filtered_counts_matrix))
summary(gsg)

#Length Class  Mode   
#goodGenes   119438 -none- logical
#goodSamples      7 -none- logical
#allOK            1 -none- logical

table(gsg$goodGenes)
#FALSE  TRUE 
#37235 82203 

table(gsg$goodSamples)
#good genes determined by dispersion, counts to low the system flags 

# filter out bad genes 
data_WGCNA <- rounded_filtered_counts_matrix[gsg$goodGenes == TRUE,]
dim(data_WGCNA)
#82203     7

#use clustering with a tree dendrogram to ID outliers samples 
htree <- hclust(dist(t(data_WGCNA)), method = 'average')
plot(htree)


#PCA - outlier detection method 

pca <- prcomp(t(data_WGCNA))
pca_data <- pca$x

#transfrom into data frame
pca_data <- as.data.frame(pca_data)

pca.var <- pca$dev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

ggplot(pca_data, aes(PC1, PC2)) + 
  geom_point()+
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2: ', pca.var.percent[2]))

# STEP 3: normalization 

colData <- row.names(filtered_sample_metadata_BASEonly)
dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA,
                                    colData = filtered_sample_metadata_BASEonly, 
                                    design=~1) # no specified groups 

#filtering 
#sum of the counts for each gene has to be greater than 15 for at least 6 of the samples 
dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >=6,]
nrow(dds_WGCNA_75) #filtered down to 29559 transcripts 

dds_norm <- vst(dds_WGCNA_75) # perform variance stabilization 

# get and save normalized counts to use below 
norm.counts <- assay(dds_norm) %>% 
  t()

# STEP 4: Network construction 
#choose a set of soft threshold powers 
power <- c(c(1:10), seq(fro = 12, to = 50, by = 2))

#call the network topology analysis function (takes a couple minutes to run)
sft <- pickSoftThreshold(norm.counts, 
                         powerVector = power,
                         networkType = "signed", 
                         verbose = 5)

sft.data<- sft$fitIndices 

#plot to pick power 

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) + 
  geom_point()+
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") + 
  labs(x = "Power", y = "Scale free topography model fit, signed R^2") + 
  theme_classic() 



a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) + 
  geom_point()+
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") + 
  labs(x = "Power", y = "Mean Connectivity") + 
  theme_classic() 

grid.arrange(a1, a2, nrow = 2)
#how well do each power threshold affect to be the natural model of gene expression relatedness 
#shows that as you increase the power and strength of relationship you increase scale 
#have to consider the degree of conductivity 
#average number of connections a given gene has
# as you increase power you lower the connectivity 
#choose a threshold around 24-26 where we have cleared the .8 but still above the lowest connectivity 
#once you choose power you can then define modules and then you can test for correlation in the trait data 
#lower correlation the lower the connectivity btw genes 

# 10/29

soft_power <- 30 
temp_cor <- cor
cor <- WGCNA::cor #sets the temp_cor function to use WGCNA correlation function 

#covert data matrix to numeric 
norm.counts[] <- sapply(norm.counts, as.numeric)

#the command below creates the network and ID modules based on the parameters 
#that we chose 
bwnet26 <- blockwiseModules(norm.counts,
                            maxBlockSize = 30000,
                            TOMType = "signed",
                            power = soft_power,
                            mergeCutHeight = 0.25,
                            numericLabels = FALSE,
                            randomSeed = 1234,
                            verbose = 3 )
cor <- temp_cor #this resets the cor function to base R's cor function instead of using 
#WGCN's cor function 


# STEP 5: Explore Module Eigengenes 

module_eigengenes <- bwnet26$MEs

head(module_eigengenes)
dim(module_eigengenes)
#[1]  7 51

#get the number of genes for each module 
table(bwnet26$colors)

#plot the dendrogram and the module colors 
#produces a merged and unmerged format based on similarity cuttoff that we set with 0.25 above 
plotDendroAndColors(bwnet26$dendrograms[[1]], cbind(bwnet26$unmergedColors, bwnet26$colors),
                                                    c("unmerged", "merged"),
                                                    dendroLabels = FALSE,
                                                    addGuide = TRUE,
                                                    hang = 0.03,
                                                    guideHang = 0.05)


#save plot 
saveRDS(bwnet26, file = "outputs/bwnet26.rds")

#to load the bwnet file in later use:
#bwnet26 <- readRDS("bwnet26.rds")


# STEP 6: Correlation of modules with traits 
# Define the numbers of genes and samples 

nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

#Test for correlation between module eigengenes and trait data 

module.trait.corr <- cor(module_eigengenes, traitData, use = 'p')

#calculate values for each correlation 
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

#visulaize modul-trait association as a heatmap 

heatmap.data <- merge(module_eigengenes, traitData, by = 'row.names')
head(heatmap.data)

#adress error of row.names not being numeric 
heatmap.data <- heatmap.data %>%
  column_to_rownames(var = 'Row.names')

names(heatmap.data)

#make heat map of correlations 
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[42:44], #values may need to change based on 
             y = names(heatmap.data)[1:41], #number of eigengenes 
             col = c("blue2", "skyblue", "white", "pink", "red"))




                    
citation("goodSamplesGenes")
                    




