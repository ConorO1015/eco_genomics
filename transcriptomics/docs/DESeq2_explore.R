library(DESeq2)
#lib.loc = "/gpfs1/cl/pbiosw/R/4.4")
library(ggplot2)

setwd("~/Projects/eco_genomics/transcriptomics/")

#import counts matrix

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, 
                          row.names = 1)

tail(countsTable)



#get rid of extra zeros bc DESseq doesnt like decimals 
countTableRound <- round(countsTable)

tail(countTableRound) # look at data 


conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE,
                    stringsAsFactors = TRUE,
                    row.names = 1)
head(conds)


###########
#explore counts matrix 

colSums(countTableRound)
mean(colSums(countTableRound))
#18454529
#aim for 20m reads when designing study and around 18m when looking at data 

barplot(colSums(countTableRound),
        names.arg = colnames(countTableRound),
        cex.names = 0.5, 
        las = 2,
        ylim = c(0,30000000))

abline(h=mean(colSums(countTableRound)),
       col = "blue2",
       lwd=2)
#line is average 

#look at the average number of counts per gene 
rowSums(countTableRound)
mean(rowSums(countTableRound))
#3244.739

median(rowSums(countTableRound))
#64

apply(countTableRound,2,mean)
#N1C3     N2C2    sN1L3     N1C4    sN2H2     N1L2     N1H2     N2H1     N1H3     N1H4 
#138.9778 151.0551 251.5168 141.8671 138.1130 141.8121 184.8436 156.6249 150.6478 140.6399 
#sN1L1     N1C2     N2L4     N2C3    sN2C5     N2L2     N1L4     N2L3     N2H3     N2H4 
#152.8451 145.3552 156.2295 203.8649 102.6907 161.1628 175.2371 135.3535 143.7175 140.5130 
#N2C1 
#131.6714 
#gives a senses of variation in the sequencing effort across samples 


###################### start analysis in DESeq2 ###########################




#dds <- DESeqDataSetFromMatrix(countData = countTableRound,
                             # colData = conds,
                              #design = ~DevTemp + FinalTemp)

dds <- DESeqDataSetFromMatrix(countData = countTableRound,
                              colData = conds,
                              design = ~DevTemp + FinalTemp)
#set demsnions 
dim(dds)
#119438     21

dds <- dds[rowSums(counts(dds) >= 10) >=15]

nrow(dds)
 # 35527 = number of trasncripts with more than 10 reads in more than or equal to 15 samples 

#run the DESeq model to test for global differential gene expression 
dds <- DESeq(dds)

#list the results youve generated 
resultsNames(dds)
#"Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"
 

#visualize our global gene expression patterns using PCA
#first we need to transfrom the data for plotting using variance stablization 

vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup=c("DevTemp", "FinalTemp"),
                   returnData = TRUE)

percentVar <- round(100*attr(pcaData,"precentVar"))

finalTempColor <- c("BASE" = "grey", "A28" = "hotpink", "A33" = "red")
shapesChoose <- c("D18" = 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, shape = DevTemp)) + 
  geom_point(size = 5) + 
  scale_shape_manual(values = shapesChoose)+
  scale_color_manual(values = finalTempColor)+
  labs(x = paste0('PC1: ', percentVar[1], '%'),
       y = paste0('PC2: ', percentVar[2], '%')) + 
  theme_bw(base_size = 16)
  
  
p
  








