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

Fist we compared the each gene with an adjusted p-value. It was ordered from lowest to highest p-value. Next used the 'head' function to look at the the genes with the lowest p-value(most signficant) and the log2fold change that compares the developmental temperature of 22 v 18. It also shows directionality with negative or positive value in comparing gene expression. Next we look at the most gene with the most significant p-value. We used plotCounts function to look at the counts of each data point comparing developmental temperature and final temperature. We then created a plot to compare the counts with the developmental temperature. The plot showed that there is a clear difference between the developmental and final temperature separated by a count of roughly 1000 as D18 was centered around 1800 counts and D22 was centered around 800 counts. Then we made a MA plot that showed the log fold change (y -axis) and the mean normalized counts(x-axis). Showed that there is a lot of up regulation in gene expression in D22 compared to D18. Pots that are to the right of the plot show genes that are constantly expressed. Next a volcano plot was generated to determine which genes were being up and down regulated. Points in the right quadrant showed that there is D22 is being up regulated and D18 is being down regulated. Finally, a heat map was generated to assess the gene expression between the developmental temperature and final temperature. Showed that there is much less regulation in the D18 group as corroborated by the previous plots. For specific genes this trend in up and down regulation is flipped in which D22 is down regulated and D18 is up regulated. 






Looked at first result name but could choose other names to look at
differential expression between groups generated 2 plots. 
