#10/15 

#load in day 1 transcriptomic script and run DESeq object 
library(pheatmap)
resultsNames(dds)
options(bitmapType = "cairo")
#pull out the results for Devleopmental temp 22 v 218 
res_D22vsD18 <- results(dds, name="DevTemp_D22_vs_D18", alpha = 0.05)

#order by significance 

res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]

head(res_D22vsD18)
#shows the measure of differences in gene expression 
#compares 22 v 18 
#what is the gene expression in 22 relative to 18 
#log2fold change does that in doubling (log 2 of 1 = 2 times much expression)
#shows directionality with signs. 
#large log2fold change is a greater the effect size 
#validate what two genes are getting compared expression wise.

#look at counds of a specific top gene that we're intrested in to validate that the model is working 
d <- plotCounts(dds, gene = "TRINITY_DN140854_c0_g5_i2", int =(c("DevTemp", "FinalTemp")),
                returnData = TRUE)
d

p <- ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp)) +
  theme_minimal() + theme(text=element_text(size=20), panel.grid.major=element_line(colour = "grey"))

        
        
p <- p + geom_point(position = position_jitter(w=.2, h=0),size=3)
p  


#plot shows that there is a differnce in count btw the dev temp of 22 and 18.
#18 group was around 1800 counts and the 22 group was around 800 groups 


#making manhattan plot 
plotMA(res_D22vsD18, ylim=c(-4,4))
#lots of up regulation in development in 22 degrees compared to 18 
#MA plot is a logfold chand v the aveage gene expression 
#genes very highly expressed are all the way to right and do not sure much differential gene expression 
#these could show replication/metablims etc
#genes that are on all the time 

#volcano plot 
#convert our deseq results object to a data frame to plot 

res_df <- as.data.frame(res_D22vsD18)

#add column to df to label if gene is significantly differential expressed or not 

res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >1, "Significant", 
                             "Not Significant")

#plot 
ggplot(res_df, aes(x = log2FoldChange,y = -log10(padj),color = Significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("slateblue", "tomato")) +
  labs(x = "Log2 Fold Change", y = "log10 Adjusted P-value", title = "Volcano Plot") +
  theme_minimal() +
  theme(legend.position = "top") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") + 
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "orange")

# shows more significance in down regulation 
#specific gene in red quadrant gene is up regulated 
#some genes in the other quadrant are down regulated 
#upregulated in 22 and down reg in 18 in the right quadrant 
  

#heat map 

vsd <-vst(dds, blind = FALSE)

#to much data to anaylze so going to take the top genes in the results tab 
topgenes <- head(rownames(res_D22vsD18), 20) 
mat <- assay(vsd)[topgenes,]
df <- as.data.frame(colData(dds)[,c("DevTemp", "FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=FALSE, cluster_cos=T, cluster_rows=T)

#red is highly expressed 
#blue is low expression 
#block of blue in the middle many samples corresponding to d18 
#higher expression corresponds to d22









