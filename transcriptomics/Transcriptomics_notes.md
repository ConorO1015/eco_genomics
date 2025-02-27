---
editor_options: 
  markdown: 
    wrap: 72
---

# Transcriptomics

#10/8 experimental questions regarding copepod development and placidity
at different environments (physiological mechanisms of development
plasticity) Experimental questions - does the temp that they experience
growing up affect there UTL? - how does gene expression response differ
between 28c and 33c and does this differ w baseline? - what genes are
deferentially expressed at developmental temperature 22 compared to DT
18?

Factors 1. developmental temperature - 18 - 22 2 levels

2.  final temperature

-   baseline
-   28 (A28) -33 (A33) 3 levels

### coding

# 10/10

Dev temp is 18 and 22 from the data the final temps are listed on top of
matrix with the number of transcripts/SNPS in the matrix - some genes
are very highly expressed - many genes that have very low expression
this leads to an inverse exponential graph over dispersion is where the
variance is greater than the mean in the data - this is why gene
expression uses negative binomial distribution - could also use a poison
distribution but doesn't account for some variation in the data

vcf file

1\. look at matric

2\. improt into DESeq

3\. use linear model to relate our variables such as gene expression to
our explanatory variables (treatment conditions) - gene expression can
be explained by dev temp and final temp (2 factorsDT and 3 factors in
FT) size factors

## DESeq 2 explore code

File was read in and the data was rounded so that it could be analyzed
and plotted. The average number of reads was 18454529 which is an
accurate amount of data to have to be able to analyze the gene
expression in Copepods. A bar graph was made to analyze the number of
reads of each sample. Next the average number of counts per gene was
3244.739 and a median of 64. Demensions of the matrix were set and the
data the data was filtered using the 'nrow' function where counts are at
least 10 in 15 or more samples. This led to 35527 counts in the data for
analysis. A PCA plot was generated to asses the gene expression pattern
occurring between each group. The plot shows that there is large
variation in gene expression based on final temperature and
developmental temperature. The Base line expression and the 28 degree
final temperature were plotted very close to each other PC1 and PC2
space. This could indicate that there is similar gene expression
occurring during these environmental conditions. However, between these
two groups and the Final temperature of 33 degrees there is a lot of
variation in gene expression. Overall, developmental temperature seems
to play less of a factor in gene expression than the final temperatures
as the different developmental temperatures are mostly plotted next near
each other regardless of final temperature. However, for the baseline
expression and final temperature of 33 degrees there is some variation
in gene expression between development temperature.

# 10/15
## Trancsript 02

Fist we compared the each gene with an adjusted p-value. It was ordered from lowest to highest p-value. Next used the 'head' function to look at the the genes with the lowest p-value(most significant) and the log2fold change that compares the developmental temperature of 22 v 18. It also shows directionality with negative or positive value in comparing gene expression. Next we look at the most gene with the most significant p-value. We used plotCounts function to look at the counts of each data point comparing developmental temperature and final temperature. We then created a plot to compare the counts with the developmental temperature. The plot showed that there is a clear difference between the developmental and final temperature separated by a count of roughly 1000 as D18 was centered around 1800 counts and D22 was centered around 800 counts. Then we made a MA plot that showed the log fold change (y -axis) and the mean normalized counts(x-axis). Showed that there is a lot of up regulation in gene expression in D22 compared to D18. Pots that are to the right of the plot show genes that are constantly expressed. Next a volcano plot was generated to determine which genes were being up and down regulated. Points in the right quadrant showed that there is D22 is being up regulated and D18 is being down regulated. Finally, a heat map was generated to assess the gene expression between the developmental temperature and final temperature. Showed that there is much less regulation in the D18 group as corroborated by the previous plots. For specific genes this trend in up and down regulation is flipped in which D22 is down regulated and D18 is up regulated. 


# 10/17 
In this set of code we are going to be comparing the gene expression of different levels of developmental temperature to the different levels of final temperature. First we manipulated the data to group each development temperature with the final temperature. This created 5 different groups to compare. We compared D18 and D22 to the Baseline temp, D18 and D22 to A28, and D18 and D22 to A33. Three separate portions of code were generated for each of these comparisons to determine the up regulation and down regulation of gene expression under each condition and to generated an MA plot. 
Next we determined the number of genes expressed for each developmental/final temperature comparison as well as the overlap of genes that are expressed in multiple different groups. This was in order to generate a ven diagram of the gene expression. There were 23 genes that were expressed among all treatment groups. 

#10/22
overall question: 
- look at each group 
- BASE (d18 v d22)
- a28(d18 v d22)
- a33(d18 v d22)
- want to create ven diagram to compare all three 
- need 7 different numbers. one for each section of the ven diagram 


After a ven diagram was generated, a scatter plot was made to compare specific treatment groups. We compared the D18 with the Baseline temperature and A28 and D22 with the Baseline temp and A33. This was also done to compare A33 with both dev temps. The first scatter plot comparing the response to 28c v Baseline varying by DevTemp. It showed that there was up regulation in D22 (r272 in red) under A28 conditions compared to baseline conditions as it is clustered around the 2.5 with some samples reading as high as 7.5. There is very little down regulation in D22(blue). There is some up regulation (pink) in D18 in response to 28c compared to base but not as large as that of D22. There is some down regulation in D18 compared to A28(turquoise) however, it is spread over a large area which may demonstrate that the baseline treatment has little effect on down regulation. 

The second scatter plot is comparing baseline vs 33c at developmental temperatures shows more values than the first scatter plot. There is a great number of genes that are unregulated (1127 in red) but also a large number of genes down regulated (240 turquoise) when comparing the change
Comparing both scatter plots there are many more genes differential expressed at 33 degrees than 28 degrees. 


# 10/24 
# WGCNA 
- if low p value there will be large test stat 
- if negative log fold change there will be negetive test stat 
grid arrange 
 - allows you to take objects already saved and put them in order 

In this script we are using WGCNA to visualize and gene correlation. First the files were read in an d working directory set. We created a matrix of just BASE reads as this is the control of the experiment and all samples will be compared to this matrix. The second step was to filter and detect any outliers. We used the 'goodSamplesGenes' function to assess the quality of gene expression counts. It looks in the matrix and returns objects that are missing values, low variability or dispersion. R returns true (good) false (bad) genes that do or do not need to be filtered out based on criteria. Based on this the bad genes were filtered out leaving us with 82,203 genes in 7 samples. Next outliners were detected using the 'htree' function (tree dendrogram) that assesses clustering. We also used PCA that showed how much variation between samples existed. Helps to visualize any outliers. Step three was to normalize the data. Genes were filtered based on read depth of counts of at least 15 in 6 or more samples. 'vst' was used for variance stabilization and the data was read back into a matrix. Step 4 was to construct the network for WGCNA analysis. A soft power threshold was also set as more stringent power can be applied later. After this was ran, 2 plots were generated that compared the scale free topography model v power and the mean connectivity(average connectivity between genes) vs power. Regarding the first plot, you want to select a power that meets or just exceeds the R^2 value(high R^2). Generally, as the power increases the connectivity decreases as stronger filtering reduces the number of connections. You want to choose a value that retains sufficient connectivity but not to much as to filter out too many connections. A soft power of 26 was choosen. 

 
#10/29
We continued with 






