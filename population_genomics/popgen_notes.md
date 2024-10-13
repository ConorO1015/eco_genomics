# coding and data notes for pop genomics module

## Author: Conor O'Rourke

### 09/10/24 - Intro to Centaurea GBS data and working with VCF files

We'll be analyzing the GBS data from 3 regions (EU, NE, PNW) starting today with variant call format files (VCF's)

zcat file.gz \| head to get token go to settings the at the bottom is developer setting. cd is a change in directory ls - l is a long list can be abbreviated at ll #zcat to view a zipper file head is used to view the first 10 lines of the file pipe is used to send output of one fxn to another command line/terminal #letters at the end of the sequence is the specific barcode for that indv.

looked at the fastq file using zcat

# 09/12/24

## viewing VCF files and talking about filtering

higher number of contigs the more assembled the genome is to acess things on the vac use spack load

use tview to view insdide doc o

q = quit cd .. or - to go back from whatever file you were in

# 9/17

./. is NA when looking at head of file

# filtering strats

depth, low frequency alleles, and missingness at indv. level and SNP. level depth: low dp = limited accuracy of genotype call. if high DP assembly error or gene paraolgy which is when you get reads from duplicated genes mapping to same position. If SNP level is 10% then only 10% of indv. have that SNP. Low frequncy allele is an allele that is less than 1% of population. Key is to find a sweet spot btw all 3 filtering peramiters df[rows, columns] by leaving a postion blank it will take all parameters

Went through how to filter data and the different ways to cotinually name new filerting graphs. this can be done by looking at the raw data and determining at which point based on q score, missingness, depth, and allele frequencies that data should be cut out ot generate the most accurate data possible.

looking at vcfr file \## Header Section

The first eight columns of the fixed region are titled CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO. This is per variant information which is fixed, or the same, over all samples. The first two columns indicate the location of the variant by chromosome and position within that chromosome. Here, the ID field has not been used, so it consists of missing data (NA). The REF and ALT columns indicate the reference and alternate allelic states for a diploid sample. When multiple alternate allelic states are present they are delimited with commas. The QUAL column attempts to summarize the quality of each variant over all samples. The FILTER field is not used here but could contain information on whether a variant has passed some form of quality assessment.

## Genotype section

The gt (genotype) region contains information about each variant for each sample. The values for each variant and each sample are colon delimited. Multiple types of data for each genotype may be stored in this manner. The format of the data is specified by the FORMAT column (column nine). Here we see that we have information for GT, AD, DP, GQ and PL. The definition of these acronyms can be referenced by querying the the meta region, as demonstrated previously. Every variant does not necessarily have the same information (e.g., SNPs and indels may be handled differently), so the rows are best treated independently. Different variant callers may include different information in this region.

<https://grunwaldlab.github.io/Population_Genetics_in_R/reading_vcf.html#>:\~:text=The%20fix%20region%20contains%20information,the%20same%2C%20over%20all%20samples.

# 09/19

### `%in%` used to match bojects

Made Manhattan plot after filtering data so that each chromosome has correct syntax for plot

Made a Manhattan plot showing FST. First the vcf file was manipulated so that it can be used for data analysis. data was grouped by region but could be grouped in other ways. Next Chromes were changed to be numbered from 1-8 and fst values associated with each position on set chrome. Data was manipulated so that pure numbers could be read by code and plotted to show fst. In addition, any 0 fst values were filtered out of the data.

Analysis of Fst of the groups using ggplot showed that there is very little differentiation within the data with the exception of chromosome 8 as it has many indv(dots) that are higher in variability compared to all other chromosomes. This could be do to changes in this populations environment causing selection pressures. When line of quadrilles is changed to .5 or the average of the data the line is almost at 0 futhur showing that there is very little divergence btw populations

#09/24

Continuing with analysis of the data from 9/19. Diversity within each group was assessed using the Hs values. Hs. values from the data were manipulated so that it could be assigned to each group and plotted appropriate. Tidyverse was used to group each Hs to its group using '%\>%'. A histogram of each group was generated with each group overlayed on top to compare populations heterozgosity. Showed that one allele is almost completely fixed in the populations as there is an extreme skew to the left of the graph which is the homozygous dominant alleles. Could be some error within the data as there is an extreme value of 10,000 snps in the Southern European population. Finally a table was generated of each subpopulation average Hs, Standard deviation of Hs and the sample size.

# PCA and admixture

# 09/26

admixture analysis (structure) has a genetic model behind it which PCA does not. The model is one of Hardy-W equilibrium. Good way to predict classes of genotypes when random mating is at play. Can also predict Hs if you know allele frequencies. Tries to group mixture into a certain number of groups into K. K has to be chosen 1. choose value of K. typical value 1-10 2. Assign individuals to one of the K groups. can be random or with prior knowledge 3. Calculate allele frequencies in each group 4. Calculate 2pq based on these frequencies and compare to the observed frequencies of Het 5. Does process over again to increase the prob of matching the observed heterozygosity ind. could be a member of multiple groups that leads to Q Q is the fractional ancestor of each indv in a group. Made in a matrix in K columns and row indv. some indv might be 0.5 0.5(completely separate) or 1 0(Homozygous for one allele) Plot is created based on this Then a cross validation model is generated to assess how good the model is working If K 'works' or is close to true value of K there will be a dip in the CV model as there is less error in predicting the genotype. Assumes that population structure is the only thing that is at work affecting heterozygosity

### coding

Files were input and read. The data was filtered using the min.distance function to get rid of any SNPs that were to close together to remove bias in the data. A file that was very large was opened converted from vcf to geno and moved to our home directory because if was to large for R. After a PCA graph was generated using the info contained in the thinned geno file using the pcaProject function. The eigenvalues were on the y and the index was on the x. The graph showed each PCA that will be input into a full PCA graph. The highest dot has the highest genetic value at PCA. the second highest has the second highest genetic info and so forth. The graph is steep and becomes almost linear because the lower dots have less accurate info regarding genetic makeup. Finally, a PCA plot was generated for PC1 v PC2 and PC2 v PC3. Showed that the pacific North West has the most divergent population as it is clustered closely together with little overlap to the other populations. Overall, the other populations in the study had lots of overlap on the graph showing that there is much less divergence between these populations and that they are very similar. There is slight divergence of the central European population but there is more overlap with other populations than the Pacific Northwest

#10/1 continued admixture/PCA First we used the 'snmf' function to run an admixture like analysis to compare each population. An admixture algorithm would take to long for our purposes. The analysis was run 3 times but can be ran as many times as needed to generate most accurate scores. Next we generated a plot that compared number of ancestreal popualitons and the cross entropy scores. We also generate a plot that compared the Eigenvalues and the number of PCs and compared both to interment he best number of K. The K value was picked based on the 'elbow' of both graphs as this showed the most accurate number for K based on the data. K value of 5 was used to generate a bar plot that was separated by continent and then further by region and population within each continent. The plot corroborated what the PCA plot generated in which the Pacific Northeast had the most genetic divergence from the other populations tested. It also showed that there was a degree of divergence in the Central European population.

#10/4 selection First pcaadapt function was ran with a k value of two and a countwise method to apply PCA directly to the SNP data and returns a separate test for selection on each PC axis instead of the entire data. A manhattan plot was generated with p-values were then plotted on the y axis and the SNP on the x axis to determine how statistically significant each SNP is.There were some SNPs with very significant p-values over 20 and 30. This could indicate that there is some selection going on at particular locations between 0-5000 and 10000 - 15000 which could indicate association of traits or under selection pressures.Next PC1 and PC2 plots were generated for each chromosome to assess any selection in the data on a chromosomal level. Almost all chromosomes have SNPS that are above a p-value of 20. The most significant peaks are in chromosome 1,7, and 8. These are locations that are of interests in the genome to determine if selection is at play. If the location in the genome is known the SNP can potentially be linked to a specific function. Since this is PC1 the data could also indicate population structure. PC2 has much less variation in p-values as the highest with the highest p-value being at 15. The most significant finding is that from chromosome 1 and 8 as both have very large p-values in both PC1 and PC2. This indicates that in these chromosomes specifically there is linked section occurring. It could also mean that there is local adaption occurring with a specific regional population.
