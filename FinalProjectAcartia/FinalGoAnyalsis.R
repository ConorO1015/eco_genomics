getlibrary(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)
library(pheatmap)
library(DOSE)
library(AnnotationDbi)
library(clusterProfiler)
library(tibble)
library(GOfuncR)
library(topGO)
library(GO.db)
options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/GroupProjects/Acartia/")

################################################
#
#normalize hudsonica data 
#
###############################################
ahuddata <- read.delim("ahud_samples_R.txt")
ahuddata_subset <- subset(ahuddata,treatment=="AM" & generation=="F0")
ahud_subset <- subset(ahuddata, subset = (treatment %in% c("AM","OW")& generation == "F0"))


countsTable <- read.table("/gpfs1/cl/pbio3990/GroupProjects/Acartia/salmon.isoform.counts.matrix.filteredAssembly",
                          header = TRUE, row.names = 1)

##########################
#
#subset the tonsa data files
#
#########################

tonsadata <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                        header = TRUE,
                        stringsAsFactors = TRUE,
                        row.names = 1)
head(tonsadata)
treatment <- tonsadata[which(tonsadata$FinalTemp =="BASE"), ]
print(treatment)

#keeping these colums because they are the correct conditions for this study 
columns_to_keep <- grep("^AA_F0|HA_F0", colnames(countsTable), value=TRUE)
subset.counts_hudsonica <- countsTable[, columns_to_keep]

#round the subsetted data becase DeSeq does not like decimals
subset.counts_hudsonica <- round(subset.counts_hudsonica)

head(subset.counts_hudsonica)


################################################
#
#subsetting the tonsa counts matrix
#
###############################################
countsTable2 <- read.table("/gpfs1/cl/pbio3990/GroupProjects/Acartia/tonsa_counts.txt",
                           header = TRUE, row.names = 1)
subset.counts_tonsa <- countsTable2[, c("N1C3","N2C2","N1C4","N1C2","N2C3","N2C1")]


#round the subsetted data becase DeSeq does not like decimals
subset.counts_tonsaround<- round(subset.counts_tonsa)

head(subset.counts_tonsaround)


##################################################
#
# filtering BLAST query genes for best evalue 
#
##################################################

library(dplyr)

blast_data <- read.table("acartia_blast_max1", header = FALSE, sep = "\t", 
                         stringsAsFactors = FALSE)

colnames(blast_data) <- c("query", "subject", "identity", "length", "mismatch", 
                          "gapopen", "q.start", "q.end", "s.start", "s.end", 
                          "evalue", "bitscore")
best_blast_matches <- blast_data %>%
  group_by(query)%>%
  filter(evalue == min(evalue))%>%
  ungroup()

best_blast_matches <- best_blast_matches %>%
  group_by(query) %>%
  slice(1) %>%
  ungroup()

print(best_blast_matches)


############################################################
#
#call genes from counts matrix and BLAST 
#merge the filtered BLAST file with the counts matricies 
#
###########################################################

merged_ahud1 <- merge(best_blast_matches, subset.counts_hudsonica, by.x = "query", by.y = "row.names", all = TRUE)

merged_tonsa <- merge(best_blast_matches, subset.counts_tonsaround, by.x = "subject", by.y = "row.names", all = TRUE)

merged_all <- merge (merged_ahud1, merged_tonsa, by.x = "query", by.y = "query")

############################################################
#
# subset merged table
#
#############################################################

subset.merged_all <- merged_all[, c("subject.x","query","AA_F0_Rep1_","AA_F0_Rep2_","AA_F0_Rep3_","HA_F0_Rep1_","HA_F0_Rep2_","HA_F0_Rep3_","N2C2","N1C4","N1C2","N2C3","N2C1","N1C3")]

subset.merged_all$combined <- paste(subset.merged_all[,1], subset.merged_all[,2], sep = "_")

#write.csv(subset.merged_all, "/gpfs1/cl/pbio3990/GroupProjects/Acartia/subset.merged_all", row.names = FALSE)

#make column 1 the row names 
rownames(subset.merged_all) <- subset.merged_all$combined
#merge the names of the genes so they are all unique and then make them the row names 

subset.merged_all <- subset.merged_all[, -c(1,2, ncol(subset.merged_all))]

################################
#
#Normalize subsetted merged file
#
##############################
head(subset.merged_all)
dim(subset.merged_all)
# 5669     12
# subset.merged_all[acartia_counts < 0] <- 0
# any(acartia_counts < 0)
# which(acartia_counts < 0, arr.ind = TRUE)

acartia_conds <- read.delim("/gpfs1/cl/pbio3990/GroupProjects/Acartia/acartia_conditions.txt")
head(acartia_conds)

#create dds object 
dds <- DESeqDataSetFromMatrix(countData = subset.merged_all, colData = acartia_conds,
                              design = ~ Treatment + Species + Treatment:Species)
dds <- DESeq(dds)

#DESeq results, use these names to run contrasts 
resultsNames(dds)
#[1] "Intercept"                 "Treatment_OW_vs_AM"        "Species_Atonsa_vs_Ahud"   
#[4] "TreatmentOW.SpeciesAtonsa"


# Filter data to include only common genes
#common_genes <- intersect(tonsa_go$geneID, normalized_counts$subject)
#tonsa_go <- tonsa_go %>% filter(geneID %in% common_genes)
#norm_counts <- normalized_counts %>% filter(subject %in% common_genes)


############### GO by Species ##################
library(stringr)


#call specific contrasts and filter by signifcantly differentially expressed genes 
res_species <- results(dds, name = "Species_Atonsa_vs_Ahud", alpha = 0.05)
sig_genes <- res_species[!is.na(res_species$padj) & res_species$padj < 0.05, ]

sig_genes <- sig_genes[order(sig_genes$padj),] #sort by p-val
view(sig_genes)

sig_gened <- as.data.frame(sig_genes)

sig_gened <- rownames_to_column(sig_gened, var = "BothNames")
view(sig_gened)

sig_genes1 <- sig_gened %>% 
  mutate(ExtractGene = str_extract(sig_gened[[1]], "(?<=::)[^:]+(?=::)"))

view(sig_genes1)
#summary(sig_genes)

common_genes <- intersect(tonsa_go$geneID, sig_genes1$ExtractGene) 
tonsa_go <- tonsa_go %>% filter(geneID %in% common_genes)
sig_genes1 <- sig_genes1 %>% filter(ExtractGene %in% common_genes)

sig_genes1 <- sig_genes1 %>%
  dplyr::select(ExtractGene, log2FoldChange, padj) %>%
  mutate(Group = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Up_Tonsa",
    log2FoldChange < 0 & padj < 0.05 ~ "Down_Tonsa",
    TRUE ~ "nonsig"
  ))


#load in GO file from tonsa
tonsa_go <- read.delim("Gene_GO_noNA.tsv", header = TRUE) 
#gene2go <- split(tonsa_go$GO, tonsa_go$geneID)
#view(tonsa_go)

gene2GO <- tonsa_go$GO %>% strsplit(";") %>% setNames(tonsa_go$geneID)
geneList_UPelevated <- setNames(as.integer(sig_genes1$ExtractGene == "Up_Tonsa"), sig_genes1$ExtractGene)
geneList_UPambient <- setNames(as.integer(sig_genes1$ExtractGene == "Up_Hudsonica"), sig_genes1$ExtractGene)



# Run GO enrichment analysis (Biological Process as default)
#go_results <- enrichGO(
#  gene          = "ExtractGene",                # Vector of gene names
#  OrgDb         = gene2go,             # Annotation database (e.g., human)
#  keyType       = "SYMBOL",                 # Use gene symbols; change as needed
#  ont           = "BP",                     # Ontology: BP (Biological Process), MF, CC
#  pAdjustMethod = "BH",                     # Adjust for multiple testing
#  pvalueCutoff  = 0.05,                     # p-value threshold
#  qvalueCutoff  = 0.2                      # FDR threshold
#)



########## GO #############


# TopGO analysis function
#dont know if i need this 
run_topgo <- function(geneList, description) {
  GOdata <- new("topGOdata",
                description = description,
                ontology = "BP",
                allGenes = geneList,
                geneSel = function(x) x == 1,
                nodeSize = 10,
                annot = annFUN.gene2GO,
                gene2GO = gene2GO)
  result <- runTest(GOdata, algorithm = "parentChild", statistic = "fisher")
  res_table <- GenTable(GOdata, classicFisher = result, topNodes = length(usedGO(GOdata))) %>%
    mutate(classicFisher = as.numeric(classicFisher),
           Annotated = as.numeric(Annotated),
           Significant = as.numeric(Significant)) %>%
    filter(classicFisher < 0.05)
  return(res_table)
}

#everything works up to here :(

A_up_results <- run_topgo(geneList_UPambient, "Ambient_up") %>%
  mutate(group = "Ambient_up")
E_up_results <- run_topgo(geneList_UPelevated, "Elevated_up") %>%
  mutate(group = "Elevated_up")

# Filter broad terms
filter_terms <- function(df) {
  df %>% filter(Annotated >= 15 & Annotated < 500)
}

setwd(path.expand("~"))
















######### Scarp code ##################
#Run GO 
ego <- enricher(sig_genes1, TERM2GENE = gene2go)

library(topGO) #need to have teacher load
all_genes <- unique(tonsa_go$geneID)
gene_vector <- factor(as.integer(all_genes %in% sig_genes1))
names(gene_vector) <- all_genes

go_data <- new("topGodata",
               ontology = "BP", 
               allGenes = gene_vector,
               annot = annFUN.gene2GO,
               gene2GO = gene2go)


result_fish <- runTest(go_data, algorithm = "classic", statistic = "fisher" )

res_table <- GenTable(go_data, 
                      classicFisher = result_fish,
                      orderBy = "classicFisher",
                      topNodes = length(usedGo(go_data))) %>% 
  mutate(classicFisher = as.numeric(classicFisher),
         Annotated = as.numeric(Annotated),
         Significant = as.numeric(Significant)) %>% 
  filte(classicFisher <0.05)


############ plot ############
pval <- score(result_fish)
log_pval <- -log10(pval)
top_terms <- sort(log_pval, decreasing = TRUE)[1:10] # top 10 GO terms 


barplot(top_terms,
        las = 2,
        main = "Top Go Terms for Species Contrast",
        ylab = "-log10(p-values)",
        names.arg = names(top_terms))








######### GO by treatment ##################
res_treatment <- results(dds, name = "Treatment_OW_vs_AM", alpha = 0.05)
res_treatmentP <- res_treatment[order(res_treatment$padj),]