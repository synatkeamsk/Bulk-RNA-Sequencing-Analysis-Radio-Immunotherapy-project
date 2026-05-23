
#' Differential gene expression
library(DESeq2)
library(tidyverse)
library("readxl")

#' Write the csv file for now 
write.csv(counts, file = "Raw_Cound_data.csv") ## this is raw count a.k.a unnormalized

#' read the raw_count back in 
Raw.Count<- read.csv("Raw.Count.csv", stringsAsFactors = TRUE, row.names=1)
Raw.Count<- Raw.Count %>% mutate_if(is.numeric, round, digits=0)
Raw.Count<- na.omit(Raw.Count)

#' Read Synat Pheno 
Phenotype<- read.csv("Pheno_data.csv", stringsAsFactors = TRUE)
Phenotype<- Phenotype[,c("Mouse","Group", "Treatment", "Time")]
phenotype

#' another pheno with exce file 
meta<- read_excel("Pheno_sk.xlsx")
meta<- meta[, c("Mouse", "Group", "Treatment")]
all(colnames(Raw.Count) == meta$Mouse)
meta

#' Make sure the phenotypic data and count matrix are matched 
all(colnames(Raw.Count) == Phenotype$Mouse)

#' Create DESeq2 object 
#' Subset the phenotype data frame with just two groups of 3 hrs

# Select sample of 3 hrs 
Raw.Count_3hrs<- Raw.Count %>% select(M267, M273, M276, M306, M271, M285, M288, M301)

#'=====================================================================================================
                                      # 3 hours post-radiotherapy
#' ====================================================================================================

#' select relevant metatdata at 3hr
Pheno_3hrs<- Phenotype %>% 
  filter(Group %in% c("0Gy_3h", "RT_3h"))

#' Deseq2 object
dds <- DESeqDataSetFromMatrix(countData = Raw.Count_3hrs,
                              colData = Pheno_3hrs,
                              design = ~ Group)


# Prefiltering: not necessary for DESeq, but good to reduce memory size and increase speed of transformation 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)  ## DESeq2 did not put this one as a workflow =
## spare code 
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,] # https://support.bioconductor.org/p/125071/ 

# Set the reference group, which is untreated in our case. if not DESeq2 will pick the alphabetical order

#' Option 1 

#' 3 hrs 
dds$Group <- factor(dds$Group, levels = c("0Gy_3h","RT_3h"))

# DESeq
dds <- DESeq(dds)

# Examine result
resultsNames(dds)
R_3hrs <- results(dds, alpha = 0.0001)
summary(R_3hrs)

## Save differentially expressed genes at 3 hrs
# Upregulated genes 
up_gene_3hrs<- R_3hrs[which(R_3hrs$log2FoldChange > 1.3 & R_3hrs$padj < 0.0001),]
up<- as.data.frame(up_gene_3hrs)
write.csv(up, file = "up_gene_3hrs.csv")

# Down regulated genes 
down_gene_3hrs<- R_3hrs[which(R_3hrs$log2FoldChange < 0 & R_3hrs$padj < .05),]
upgene<- as.data.frame(up_gene)
write.csv(as.data.frame(down_gene_3hrs), file="down_3hr.csv")

#'=====================================================================================================
                                      # 1 day post-radiotherapy
#' ====================================================================================================
# Subset the phenotype data frame with just two groups of 1 day

# Select sample of 1 day 
Raw.Count_1d<- Raw.Count %>% select(M287, M289, M294, M305, M269, M277, M290, M293)
meta
# select relevant metatdata at day 1
Pheno_1d<- meta %>% 
  filter(Group %in% c("0Gy:1 day", "2Gy×5:1 day"))
Pheno_1d
## DEseq2 object 
library(tidyverse)
library(DESeq2)
dds_1d <- DESeqDataSetFromMatrix(countData = Raw.Count_1d,
                              colData = Pheno_1d,
                              design = ~ Group)


## Prefiltering: not necessary for DESeq, but good to reduce memory size and increase speed of transformation 
keep <- rowSums(counts(dds_1d)) >= 10
dds_1d <- dds_1d[keep,]
dds_1d <- estimateSizeFactors(dds_1d)  ## DESeq2 did not put this one as a workflow 
dds_1d$Group <- factor(dds_1d$Group, levels = c("0Gy:1 day","2Gy×5:1 day"))

#dds <- dds[, dds$Group %in% c("0Gy_1d","RT_1d")]

# DESeq
dds_1d <- DESeq(dds_1d)

# Examine result
resultsNames(dds_1d)
R_1day <- results(dds_1d, alpha = 0.05)
summary(R_1day)

#' Save upregulated genes @ 1day
#' Upregulated genes 
up_gene_1day<- R_1day[which(R_1day$log2FoldChange > 0 & R_1day$padj < .05),]
upgene<- as.data.frame(up_gene)
write.csv(as.data.frame(up_gene_1day), file="up_1day.csv")

#' Down regulated genes 
down_gene_1day<- R_1day[which(R_1day$log2FoldChange < 0 & R_1day$padj < .05),]
upgene<- as.data.frame(up_gene)
write.csv(as.data.frame(down_gene_1day), file="down_1day.csv")

#' Volcano plot @ 1 day 
Volcano.1day <- R_1day[order(R_1day$padj),]
Volcano.1day<- na.omit(Volcano.1day)

# create data.frame
results.1day<- as.data.frame(mutate(as.data.frame(Volcano.1day), 
                                 significant=ifelse(Volcano.1day$padj<0.05, "Up", "Down")), 
                          row.names=rownames(Volcano.1day))
results.1day<- na.omit(results.1day)
is.na(results.a)  ## check missing value 

#' ggplot plus ggrepel
options(ggrepel.max.overlaps = Inf)
library(ggrepel)
results.1day$significant<- factor(results.1day$significant, levels = c('Up', 'Down'))

#' plot
volcano<- ggplot(results.1day, aes(log2FoldChange, -log10(pvalue))) + 
  theme_bw() +
  geom_point(aes(col = significant)) +
  scale_color_manual(values = c("red", "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold"), 
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size=10), 
        legend.title = element_text(size = 11, face = "bold"), 
        legend.text = element_text(size = 9)) + 
  geom_text_repel(data=results.1day[1:30,], 
                  aes(label=rownames(results.1day[1:30,]))) +
  ggtitle("Radiotherapy vs. sham-radiotherapy__ [1 day]") 

volcano
#'=====================================================================================================
                                      # 4 days post-radiotherapy
#'=====================================================================================================

#' Select sample of 4 days 
Raw.Count_4d<- Raw.Count %>% select(M265, M275, M283, M303, M266, M279, M280, M296)

#' select relevant metatdata at 4 days
Pheno_4d<- Phenotype %>% 
  filter(Group %in% c("0Gy_4d", "RT_4d"))

# DEseq2 object 
dds_4d <- DESeqDataSetFromMatrix(countData = Raw.Count_4d,
                                 colData = Pheno_4d,
                                 design = ~ Group)


#' Prefiltering: not necessary for DESeq, but good to reduce memory size and increase speed of transformation 
keep <- rowSums(counts(dds_4d)) >= 10
dds_4d <- dds_4d[keep,]
dds_4d <- estimateSizeFactors(dds_4d)  ## DESeq2 did not put this one as a workflow 
dds_4d$Group <- factor(dds_4d$Group, levels = c("0Gy_4d", "RT_4d" )) # did not work
dds <- dds[, dds$Group %in% c("0Gy_4d","RT_4d")]

# DESeq
dds_4d <- DESeq(dds_4d)

# Examine result
resultsNames(dds_4d)
R_4day <- results(dds_4d, alpha = 0.05)
summary(R_4day)

#' save differentially expressed genes at 4 days

#' Upregulated genes 
up_gene_4day<- R_4day[which(R_4day$log2FoldChange > 0 & R_4day$padj < .05),]
#upgene<- as.data.frame(up_gene)
write.csv(as.data.frame(up_gene_4day), file="up_4day.csv")

#' Down regulated genes 
down_gene_4day<- R_4day[which(R_4day$log2FoldChange < 0 & R_4day$padj < .05),]
#upgene<- as.data.frame(up_gene)
write.csv(as.data.frame(down_gene_4day), file="down_4day.csv")

#'=====================================================================================================
                                      # 8 days post-radiotherapy
#'=====================================================================================================

# Select sample of 8 days 
Raw.Count_8d<- Raw.Count %>% select(M268, M278, M295, M309, M270, M286, M292, M298)

# select relevant metatdata at 4 days
Pheno_8d<- Phenotype %>% 
  filter(Group %in% c("0Gy_8d", "RT_8d"))
dds_8d <- DESeqDataSetFromMatrix(countData = Raw.Count_8d,
                              colData = Pheno_8d,
                              design = ~ Group)


#' Prefiltering: not necessary for DESeq, but good to reduce memory size and increase speed of transformation 

keep <- rowSums(counts(dds_8d)) >= 10
dds_8d <- dds_8d[keep,]
dds_8d <- estimateSizeFactors(dds_8d)  ## DESeq2 did not put this one as a workflow 

#' set reference group
dds_8d$Group <- factor(dds_8d$Group, levels = c("0Gy_8d", "RT_8d" )) # did not work
dds <- dds[, dds$Group %in% c("0Gy_8d","RT_8d")]

#' DESeq
dds_8d <- DESeq(dds_8d)

#' Examine result
resultsNames(dds_8d)
R_8day <- results(dds_8d, alpha = 0.05)
summary(R_8day)

#' Save differentially expressed genes at 8 days

#' Upregulated genes 
up_gene_8day<- R_8day[which(R_8day$log2FoldChange > 0 & R_8day$padj < .05),]
#upgene<- as.data.frame(up_gene)
write.csv(as.data.frame(up_gene_8day), file="up_8day.csv")

#' Down regulated genes 
down_gene_8day<- R_8day[which(R_8day$log2FoldChange < 0 & R_8day$padj < .05),]
#upgene<- as.data.frame(up_gene)
write.csv(as.data.frame(down_gene_8day), file="down_8day.csv")

#' whole comparison so fit new model with interaction between Treatment and Time
dds_all <- DESeqDataSetFromMatrix(countData = Raw.Count, 
                              colData = Phenotype,
                              design = ~ Treatment + Time + Treatment:Time)

Phenotype
#' relevel 

#' set reference group
dds_all$Treatment <- relevel(dds_all$Treatment, ref = "Sham-radiotherapy")

#' DESeq
dds_all <- DESeq(dds_all)

#' Examine result
resultsNames(dds_all)
all_time_interaction <- results(dds_all, alpha = 0.05)
summary(all_time_interaction)

#' volcano plot of model with interaction 
Volcano.inter <- all_time_interaction[order(all_time_interaction$padj),]
Volcano.inter<- na.omit(Volcano.inter)


#' create data.frame
results.inter<- as.data.frame(mutate(as.data.frame(Volcano.inter), 
                                    significant=ifelse(Volcano.inter$padj<0.05, "Up", "Down")), 
                             row.names=rownames(Volcano.inter))
results.inter<- na.omit(results.inter)
is.na(results.a)  ## check missing value 

#' ggplot plus ggrepel
options(ggrepel.max.overlaps = Inf)
library(ggrepel)
results.inter$significant<- factor(results.inter$significant, levels = c('Up', 'Down'))

#' plot
ggplot(results.inter, aes(log2FoldChange, -log10(pvalue))) + 
  theme_classic() +
  geom_point(aes(col = significant)) +
  scale_color_manual(values = c("red", "black")) +
  ggtitle("Radiotherapy vs. sham-radiotherapy") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size=10), 
        legend.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 15)) + 
  geom_text_repel(data=results.inter[1:20,], 
                  aes(label=rownames(results.inter[1:20,]))) + 
  xlim(c(-1, 1))


#' Upregulated genes 
up_gene_all<- all_time_interaction[which(all_time_interaction$log2FoldChange > 0 & all_time_interaction$padj < .05),]
#upgene<- as.data.frame(up_gene)
write.csv(as.data.frame(up_gene_all), file="up_all_interaction.csv")

#' Down regulated genes 
down_gene_all<- all_time_interaction[which(all_time_interaction$log2FoldChange < 0 & all_time_interaction$padj < .05),]
#' upgene<- as.data.frame(up_gene)
write.csv(as.data.frame(down_gene_8day), file="down_all_interaction.csv")

#' ===================================================================================================================
                                         # fit additive model 
#' ===================================================================================================================

#' Additive model 
dds_additive <- DESeqDataSetFromMatrix(countData = Raw.Count, 
                                  colData = Phenotype,
                                  design = ~ Treatment + Time)


#' set reference group
dds_additive$Treatment <- relevel(dds_additive$Treatment, ref = "Sham-radiotherapy")

#' DESeq
dds_additive <- DESeq(dds_additive)

#' Examine result
resultsNames(dds_additive)
all_time_additive <- results(dds_additive, alpha = 0.05)
summary(all_time_additive)

#' Save differentially expressed genes in csv files 

#' Filter genes<0.05
library(tidyverse)
additive.model <- subset(all_time_additive, padj < 0.05)

#' csv file 
write.csv(as.data.frame(additive.model), file="DEG_additive_model.csv")

#' Hierachical clustring based on differentially expressed genes 
#' Data visualization such as clustering and dimensionality reduction 

#' Data transformation for the downstream visualization 

#' Extracting transformed values  
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds_all, blind=FALSE)
head(assay(vsd), 3)

#' Effects of transformations on the variance 
#' this gives log2(n + 1) ====
ntd <- normTransform(dds_all)
ntd<- normTransform(dds)

#' load vsn library 
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#' Data quality assessment by sample clustering and visualization
#' Heatmap of the count matrix
library("pheatmap")

#' my classic code from previous session 
#' Convert Phenotype to factor  
#' read metadata

library(tidyverse)
library(ggsignif)
library(ggpubr)
require(readxl)
library(ggthemes)

#' Read data 
metadata<- read_excel("metadata.xlsx")
metadata

#' load pheatmap library  
library(pheatmap)

# Extract the matrix of transformed counts
ntd.matrix <- assay(ntd)

# Compute the correlation values between samples
ntd.corr <- cor(ntd.matrix) 

# match row name 
rownames(metadata) <- colnames(ntd.corr)  

# Euclidean/complete 

# default distance and linkage 

pheatmap(ntd.corr, annotation= select(Phenotype, Group), show_rownames = TRUE, 
         cex=1, border_color=TRUE, annotation_col = select(Phenotype, Group), 
         clustering_distance_cols = "maximum", clustering_distance_rows = "maximum",
         clustering_method = "ward.D2", 
         cluster_cols = TRUE, cluster_rows = FALSE)

#' vsd  
# Extract the matrix of transformed counts
vsd.matrix <- assay(vsd)

#' Compute the correlation values between samples
vsd.corr <- cor(vsd.matrix) 

#' match row name 
rownames(Phenotype) <- colnames(vsd.corr)  

#' Euclidean/complete =================
# default distance and linkage =============================

pheatmap(vsd.corr, annotation= select(Phenotype, Group), show_rownames = TRUE, 
         cex=1, border_color=TRUE, annotation_col = select(Phenotype, Group), 
         clustering_distance_cols = "maximum", clustering_distance_rows = "maximum",
         clustering_method = "ward.D2", 
         cluster_cols = TRUE, cluster_rows = FALSE)

#' rld  ===================
# Extract the matrix of transformed counts
rld.matrix <- assay(rld)

# Compute the correlation values between samples
rld.corr <- cor(rld.matrix) 

# match row name 
rownames(Phenotype) <- colnames(rld.corr)  

#' Euclidean/complete =================
# default distance and linkage =============================
pheatmap(rld.corr, annotation= select(Phenotype, Group), show_rownames = TRUE, 
         cex=1, border_color=TRUE, annotation_col = select(Phenotype, Group), 
         clustering_distance_cols = "maximum", clustering_distance_rows = "maximum",
         clustering_method = "ward.D2", 
         cluster_cols = TRUE, cluster_rows = FALSE)


#' PCA plot 

#' Treat vs. not treated

plotPCA(vsd, intgroup=c("condition", "Group"))
plotPCA(vsd, intgroup=c("Treatment")) + 
  theme_bw() + 
  labs(title = "Principle component analysis")

#' Different groups
plotPCA(vsd, intgroup=c("Group")) + 
  theme_bw() + 
  labs(title = "Principle component analysis")

# Advanced plot ====================================
pcaData <- plotPCA(vsd, intgroup=c("Group", "Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$Group<- factor(pcaData$Group, levels = c('0Gy_3h', '0Gy_1d', '0Gy_4d', '0Gy_8d', 
                                                     'RT_3h', 'RT_1d', 'RT_4d', 'RT_8d'))
ggplot(pcaData, aes(PC1, PC2, color=Group, shape=Treatment)) +
  geom_point(size=5) + theme_minimal() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_point(size=4) +
labs(title = "Principle Component Analysis") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size=12), axis.title = element_text(face = "bold", size=13),
        legend.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 12))


# link= https://stackoverflow.com/questions/34869768/how-to-make-a-pca-plots-as-i-posted-here 


# Save high quality plot 
ggsave("PCA.tiff", plot = pca, height=10, width = 6, units = "in", dpi = 300)

## sample distance heatmap 

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Mouse, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


#' annotation
Phenotype
annotation.row <- data.frame(Treatment = factor(vsd@colData$Treatment,
                                                levels = c("Radiotherapy", "Sham-radiotherapy")), 
                             Group = factor(vsd@colData$Group,
                                            levels = c("0Gy (0.1 day)","0Gy (1 day)","0Gy (4 days)",
                                                       "0Gy (8 days)", "2Gy × 5 (0.1 day)", "2Gy × 5 (1 day)", 
                                                       "2Gy × 5 (4 days)", "2Gy × 5 (8 days)")))

meta

#' mapping
rownames(annotation.row)<- rownames(sampleDistMatrix)

#' plot the pheatmap ====================
library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation_row = annotation.row, 
         border_color = TRUE, cluster_cols = F)

#' Code chunk for pheat-map from ASI worshop ================
sampleDists.vsd.blind <- dist(t(assay(vsd)))
sampleDistMatrix.vsd.blind <- as.matrix(sampleDists.vsd.blind)

#' Annotate the samples for the heatmap
annotation.row <- data.frame(Treatment = factor(vsd@colData$Treatment,
                                                levels = c("Radiotherapy", "Sham-radiotherapy")), 
                             Group = factor(vsd@colData$Group,
                                            levels = c("0Gy","2Gy × 1","2Gy × 3",
                                                       "2Gy × 5", "6Gy × 1", "6Gy × 2")))

#' mapping the data =========================================
rownames(annotation.row) <- rownames(sampleDistMatrix.vsd.blind)
head(annotation.row)

#' plot pheatmap =================================
# Create a colour palette for the heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))
pheatmap(sampleDistMatrix.vsd.blind, clustering_distance_rows = sampleDists.vsd.blind,
         clustering_distance_cols = sampleDists.vsd.blind, show_colnames = FALSE,
         show_rownames = TRUE, annotation_row = annotation.row, legend = FALSE, 
         clustering_method = "ward.D2")


#' Another heatmap 
library("pheatmap")
select <- order(rowMeans(counts(dds_all,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_all)[,c("Group","Treatment")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)



#' Volcano plot 
Volcano.a <- all_time_additive[order(all_time_additive$padj),]
Volcano.a<- na.omit(Volcano.a)


#' create data.frame
results.a<- as.data.frame(mutate(as.data.frame(Volcano.a), 
                                 sig=ifelse(Volcano.a$padj<0.05, "padj<0.05", "Not Sig")), 
                          row.names=rownames(Volcano.a))
results.a<- na.omit(results.a)
is.na(results.a)  ## check missing value 

#' ggplot plus ggrepel
options(ggrepel.max.overlaps = Inf)
library(ggrepel)
ggplot(results.a, aes(log2FoldChange, -log10(pvalue))) + 
  theme_classic() +
  geom_point(aes(col = sig)) +
  scale_color_manual(values = c("black", "red")) +
  ggtitle("Radiotherapy vs. sham-radiotherapy") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size=10), 
        legend.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 15)) + 
  geom_text_repel(data=results.a[1:60,], 
  aes(label=rownames(results.a[1:60,]))) 


#' exporting normalized count 
# Does normalized count impact by design formula: https://support.bioconductor.org/p/106548/
# sample distance heatmap ==========================================================

#' 1). read the raw_count back in 

Raw.Count<- read.csv("Raw.Count.csv", stringsAsFactors = TRUE, row.names=1)
Raw.Count<- Raw.Count %>% mutate_if(is.numeric, round, digits=0)
Raw.Count

# another pheno with exce file 
meta<- read_excel("Pheno_sk.xlsx")
meta<- meta[, c("Mouse", "Group", "Treatment")]
all(colnames(Raw.Count) == meta$Mouse)

#' model for heatmap 
meta$Group<- as.factor(meta$Group)
meta$Mouse<- as.factor(meta$Mouse)
meta$Treatment<- as.factor(meta$Treatment)

#' model object 

dds <- DESeqDataSetFromMatrix(countData = Raw.Count,
                              colData = meta,
                              design = ~ Group)


#' filter gene with at least 10 counts
keep <- rowSums(counts(dds)) >= 10
dds<- dds[keep,]
dds<- estimateSizeFactors(dds)  ## DESeq2 did not put this one as a workflow 

dds_1d$Group <- factor(dds_1d$Group, levels = c("0Gy_1d","RT_1d"))
#' dds <- dds[, dds$Group %in% c("0Gy_1d","RT_1d")]
#' DESeq
dds <- DESeq(dds)

#' Estimate size factor 
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts<- as.data.frame(normalized_counts)

# Write csv files 
write.csv(normalized_counts, file = "normalized_count.csv")

#' Write csv file == good tip to generate dataframe from vst transformation 
vsd <- vst(dds)
count<- assay(vsd)
write.csv(count, file = "transform_count.csv")

#' PCA 
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

#' plot 
plotPCA(vsd, intgroup=c("Group"))

#' another plot type 
pcaData <- plotPCA(vsd, intgroup=c("Group", "Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca<- ggplot(pcaData, aes(PC1, PC2, color=Group, shape=Treatment)) +
  geom_point(size=4) +
  theme_bw()+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title = "Principal Component Analysis") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        axis.text = element_text(face = "bold", size=12),
        axis.title = element_text(face = "bold", size=12),
        legend.title = element_text(size = 9, face = "bold"), 
        legend.text = element_text(size = 6, face = "bold")) 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pca
ggsave("PCA.tiff", plot = pca, height=4, width =7, units = "in", dpi = 300)

#' Heatmap 
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Mouse, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)

#' annotation ============================ 
annotation.row <- data.frame(Treatment = factor(vsd@colData$Treatment,
                                                levels = c("Radiotherapy", "Sham-radiotherapy")), 
                             Group = factor(vsd@colData$Group,
                                            levels = c("0Gy:0.1 day","0Gy:1 day","0Gy:4 days",
                                                       "0Gy:8 days", "2Gy×5:0.1 day", "2Gy×5:1 day", 
                                                       "2Gy×5:4 days", "2Gy×5:8 days")))


# mapping =======================
rownames(annotation.row)<- rownames(sampleDistMatrix)

# plot the pheatmap ====================
library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation_row = annotation.row, 
         border_color = TRUE, cluster_cols = F)