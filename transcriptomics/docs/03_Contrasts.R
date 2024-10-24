# 10/17 

library(eulerr)

#set up grouaps within DESeq object 

dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)
dim(dds)
# 35527    21

resultsNames(dds)
#[1] "Intercept"               "group_D18A33_vs_D18A28" 
#[3] "group_D18BASE_vs_D18A28" "group_D22A28_vs_D18A28" 
#[5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28

#create files for each contrast into one object into our environment  
#this is so we can compare gene expression among groups 

#taking the design and putting them into a new factor called group 
#so the matrix has all the possible group combinations so that you can compare
#each group to another

# 1. compare baseline gene expression between developmental treatment groups 
res_D18_BASE_D22_BASE <- results(dds, contrast = c("group", "D18BASE", "D22BASE"), alpha = 0.05)
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),]
head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE)
#up regulated 556 1.7%
#down regulated 1378 4.2% 

#make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj < 0.05,])

plotMA(res_D18_BASE_D22_BASE, ylim = c(-4,4))

# 2. compare gene expression between developmental treatment groups 
res_D18_A28_D22_A28 <- results(dds, contrast = c("group", "D18A28", "D22A28"), alpha = 0.05)
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),]
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),]
head(res_D18_A28_D22_A28)
summary(res_D18_A28_D22_A28)
# up reg 146 .49%
#down reg 150 .5%

#make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj < 0.05,])

plotMA(res_D18_A28_D22_A28, ylim = c(-4,4))

# 3. compare A33 gene expression between developmental treatment groups 
res_D18_A33_D22_A33 <- results(dds, contrast = c("group", "D18A33", "D22A33"), alpha = 0.05)
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),]
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),]
head(res_D18_A33_D22_A33)
summary(res_D18_A33_D22_A33)
#up reg 47, .22%
#down reg 31 .15%

#make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj < 0.05,])

plotMA(res_D18_A33_D22_A33, ylim = c(-4,4))

length(degs_D18_BASE_D22_BASE) #1935
length(degs_D18_A28_D22_A28) #296
length(degs_D18_A33_D22_A33) #78

#look at the overlap in which genes are differentially expressed in multiple contrasts 

length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A28_D22_A28)) #107
length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A33_D22_A33)) #44
length(intersect(degs_D18_A33_D22_A33,degs_D18_A28_D22_A28)) #29


nested_intersection <- intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28)
length(intersect(degs_D18_A33_D22_A33, nested_intersection)) #23
#23 genes that should be in the middle of the diagram 

#calculate the number of unique genes in each portion of the Euler plot 
#A portion is Base 
#B portion is A28
#c portion is A33
#D portion is Base and A28 intersection 
#E portion is A28 and A33 intersection 
#F portion is Base and A33 intersect
#G portion is all three 



1935 - 107 - 44 + 23 #number of genes in A portion of diagram (1807)
# 1807 genes diff expressed uniquely at baseline btw 18 v 22 

296 - 107 - 29 + 23 # 183 uniquely expressed when exposed to 28 

78 - 44 - 29 + 23 #28 genes uniquely expressed when exposed to 33 

#D 
107 - 23 #84 genes unique to BASE and A28 

#F
44 - 23 #21 genes unique to Base and A33 

#E
29 - 23 #6 genes unique to A28 and A33 

myEuler <- euler(c("BASE" = 1807, "A28" = 183, "A33" = 28, "BASE&A28" = 84, 
                   "BASE&A33" = 21, "A28&A33" = 6, "BASE&A28&A33" = 23))


#plot 
plot(myEuler, lty = 1:3, quantities= TRUE)


`##################################

#make a scatter plot of responses to A28 when copepods develop at 18 vs 22
#

#contrast D18_A28vsBASE
res_D18_BASEvsA28 <- as.data.frame(results(dds,contrast = c("group", "D18BASE", "D18A28"), alpha = 0.05))

#contrast D22_A28vsBASE
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast = c("group", "D22BASE", "D22A28"), alpha = 0.05))

#merge dataframes 
res_df28 <- merge(res_D18_BASEvsA28, res_D22_BASEvsA28, by = "row.names", suffixes = c(".18", ".22"))

#change to row names
rownames(res_df28) <- res_df28$Row.names
res_df28 <- res_df28[,-1]

library(dplyr)
library(tidyr)

#define color mapping logic with mutate function 

res_df28 <- res_df28 %>% 
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise2", #greater than 0 it is unregulated compared to 28 relative to base 
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta1",
    padj.22 < 0.05 & stat.22 < 0 ~ "blue2",
    padj.22 < 0.05 & stat.22 > 0 ~ "red1"
    ))

#count the number of points per fill color 
color_counts <- res_df28 %>%
  group_by(fill) %>% 
  summarise(count = n())

label_position28 <- data_frame(
  fill = c("blue2", "magenta1", "red1", "turquoise2"),
  x_pos = c(1, 5, 0, -7.5),
  y_pos = c(-5, 0, 9, 3)
)

label_data28 <- merge(color_counts, label_position28, by = "fill")




plot28 <- ggplot(res_df28, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) + 
  geom_point(alpha = 0.8) + 
  scale_color_identity() + 
  geom_text(data = label_data28, aes(x = x_pos, y = y_pos, label = count, color = fill),
            size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") + 
  xlim(-10, 10) + ylim(-10, 10) + 
  labs(x = "Log2FoldChange 28 vs. BASE at 18", 
       y = "Log2FoldChange 28 v BASE at 22", 
       title = "How does response to 28C vary by DevTemp?") + 
  theme_minimal()

#plot shows
#many genes unregulated at 22 degrees relative to baseline 
#when developing at 18 very small response (magenta)
#down reg when exposed to 22 
#



#make a scatter plot of responses to A33 when copepods develop at 18 vs 22


#contrast D18_A33vsBASE
res_D18_BASEvsA33 <- as.data.frame(results(dds,contrast = c("group", "D18BASE", "D18A33"), alpha = 0.05))

#contrast D22_A33vsBASE
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast = c("group", "D22BASE", "D22A33"), alpha = 0.05))

#merge dataframes 
res_df33 <- merge(res_D18_BASEvsA33, res_D22_BASEvsA33, by = "row.names", suffixes = c(".18", ".22"))

#change to row names
rownames(res_df33) <- res_df33$Row.names
res_df33 <- res_df33[,-1]


#define color mapping logic with mutate function 

res_df33 <- res_df33 %>% 
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise2", #greater than 0 it is unregulated compared to 28 relative to base 
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta1",
    padj.22 < 0.05 & stat.22 < 0 ~ "blue2",
    padj.22 < 0.05 & stat.22 > 0 ~ "red1"
  ))
#count the number of points per fill color 
color_counts <- res_df33 %>%
  group_by(fill) %>% 
  summarise(count = n())

label_position33 <- data_frame(
  fill = c("blue2", "magenta1", "red1", "turquoise2"),
  x_pos = c(1, 5, 0, -7.5),
  y_pos = c(-5, 0, 9, 3)
)

label_data33 <- merge(color_counts, label_position33, by = "fill")




plot33 <- ggplot(res_df33, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) + 
  geom_point(alpha = 0.8) + 
  scale_color_identity() + 
  geom_text(data = label_data33, aes(x = x_pos, y = y_pos, label = count, color = fill),
            size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") + 
  xlim(-10, 10) + ylim(-10, 10) + 
  labs(x = "Log2FoldChange 33 vs. BASE at 18", 
       y = "Log2FoldChange 33 vs. BASE at 22", 
       title = "How does response to 33C vary by DevTemp?") + 
  theme_minimal()

#put the two plots together in a two panel plot 

library(gridExtra)

combined_plot <- grid.arrange(plot28, plot33, ncol = 2)

ggsave("~/Projects/eco_genomics/transcriptomics/figures/combined_plot.png", combined_plot, width = 12, height = 6) 

# many more genes differential epressed at 33 








