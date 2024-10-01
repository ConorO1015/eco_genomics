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

Analysis of Fst of the groups using ggplot showed that there is very little differentiation within the data with the exception of chromosome 8 as it has many indv(dots) that are higher in variablity compared to all other chromosomes. This could be do to changes in this populations enviornment causing selection pressures. When line of quartiles is changed to .5 or the average of the data the line is almost at 0 futhur showing that there is very little divergence btw populations

#09/24

Continuing with analysis of the data from 9/19. Diversity within each group was assessed using the Hs values. Hs. values from the data were manipulated so that it could be assigned to each group and plotted appropriate. Tidyverse was used to group each Hs to its group using '%\>%'. A histogram of each group was generated with each group overlayed on top to compare populations heterozgosity. Showed that one allele is almost completely fixed in the populations as there is an extreme skew to the left of the graph which is the homozygous dominant alleles. Could be some error within the data as there is an extreme value of 10,000 snps in the Southern European population. Finally a table was generated of each subpopulation average Hs, Standard deviation of Hs and the sample size.

# PCA and admixture

# 09/26

admixture analysis (structure) has a genetic model behind it which PCA does not. The model is one of Hardy-W equilibrium. Good way to predict classes of genotypes when random mating is at play. Can also predict Hs if you know allele frequencies. Tries to group mixture into a certain number of groups into K. K has to be chosen 1. choose value of K. typical value 1-10 2. Assign individuals to one of the K groups. can be random or with prior knowledge 3. Calculate allele frequencies in each group 4. Calculate 2pq based on these frequencies and compare to the observed frequencies of Het 5. Does process over again to increase the prob of matching the observed heterozygosity ind. could be a member of multiple groups that leads to Q Q is the fractional ancestor of each indv in a group. Made in a matrix in K columns and row indv. some indv might be 0.5 0.5(completely separate) or 1 0(Homozygous for one allele) Plot is created based on this Then a cross validation model is generated to assess how good the model is working If K 'works' or is close to true value of K there will be a dip in the CV model as there is less error in predicting the genotype. Assumes that population structure is the only thing that is at work affecting heterozygosity

### coding 

Files were input and read. The data was filtered using the min.distance function to get rid of any SNPs that were to close together to remove bias in the data. A file that was very large was opened converted and moved to our home directory because if was to large for R. After a PCA graph was generated using the info contained in the thinned geno file using the pcaProject function. The eigenvalues were on the y and the index was on the x. The graph showed each PCA that will be input into a full PCA graph. The highest dot has the highest genetic value at PCA. the second highest has the second highest genetic info and so forth. The graph is steep and becomes almost linear because the lower dots have less accurate info regarding genetic makeup. Finally, a PCA plot was generated for PC1 v PC2 and PC2 v PC3. Showed that the pacific North West has the most divergent population as it is clustered closely together with little overlap to the other populations. Overall, the other populations in the study had lots of overlap on the graph showing that there is much less divergence between these populations and that they are very similar. There is slight divergence of the central European population but there is more overlap with other populations than the Pacific Northwest


# 10/1
### agenda 
admixture analysis 
- calulations
- plotting 

