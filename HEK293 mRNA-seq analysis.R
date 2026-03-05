########################################################################################################
# Code for performing Differential Gene Expression Analysis on Exp. 13-1 (HEK SR64 vs SR71 mRNA-seq)
# row means >=10 filter before dds and splicing
# Christina Ramelow, PhD
# Date: 08/27/2025
########################################################################################################
# Exp. 13-1 Comparisons
#Current comparisons (C) and hypotheses: 

#1.HEK global vs SR64 pulldown
#Hypothesis: Little DEGs between groups

#2.HEK global vs SR71 pulldown
#Hypothesis: Very little DEGs between groups.

#3.SR64 pulldown vs SR71 pulldown
#Hypothesis: More cytosolic stuff in SR64 and more nuclear stuff in SR71


##################################################################################################################
##################################################################################################################
## Part 1: Data processing and clean up
# steps:
# 1. Filter out low counts
  # a. Venns: keep a gene if row mean count >= 10/sample group
  # b. PCAs, correlations, diffex: keep a gene if row mean count >= 10 across all samples
# 2. DESeq2 normalization
  # a. median of ratios approach (one file with all samples and 1 file w/o neg.pulldowns)
##################################################################################################################
#load Libraries
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)
library(ggplot2)
library(ggrepel)

stringsAsFactors = FALSE

rootdir <- setwd("/Users/christina/Desktop/Exp. 13-1 HEK SR64 vs SR71 RNA-seq analysis")

# load feature_counts.txt file
raw_counts <- read.table("Exp. 13-1_feature_counts.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
dim(raw_counts)
# [1] 50116    19

# remove the first 5 columns that are primarily categorical
raw_counts_2 <- raw_counts[,c(6:19)]
dim(raw_counts_2)
# [1] 50116    14

# check if there are NAs in the counts df
sum(is.na(raw_counts_2))
#[1] 0

# load RNA traits file
traits <- read.csv("Exp. 13-1 traits_RNA.csv",  header = TRUE, row.names = 1)
dim(traits)
# [1] 14  3

# rename the column name of raw_counts_2 df to match traits df
# Extract row names from traits
new_col_names <- rownames(traits)

# Assign row names of traits as column names of raw_counts_3
names(raw_counts_2) <- new_col_names

# check that the row names of the traits df and the column names of the counts df are the same
raw_counts_clean <- raw_counts_2[,match(rownames(traits), colnames(raw_counts_2))]

#sanity check 1 colnames in counts file = rownames in traits file
all(colnames(raw_counts_clean) %in% rownames(traits))
#TRUE

#sanity check 2 colnames in counts file = rownames in traits file
all(colnames(raw_counts_clean) == rownames(traits))
#TRUE

# write .csv file for raw_counts_clean for Exp. 13-1
write.csv(raw_counts_clean, file=paste0("0.Exp. 13-1_RNA-seq_raw_counts_clean_unfiltered-",dim(raw_counts_clean)[1],"x",dim(raw_counts_clean)[2],"_CR_08272025.csv"))

################################################################################
################################################################################
## Step 1a: Filter out low counts using a row mean >/= 10 by group for Venn diagrams
# Exp. 13-1
# 1. subset the traits and counts dfs by sample group to create list of genes for Venn diagrams
# 2. calculate the row mean
# 3. filter for genes with a row mean >= 10
# 4. write a .csv file for each sample group
################################################################################

# row mean >= 10 HEK.SR71.global
counts_HEK.SR71.global <- raw_counts_clean[,c(1:3)]
dim(counts_HEK.SR71.global)
# [1] 50116     3

traits_HEK.SR71.global  <- traits[c(1:3),]
traits_HEK.SR71.global
#             Sample.ID..      Group Replicate
# 23238R-11-01         CR1 HEK.global         1
# 23238R-11-02         CR2 HEK.global         2
# 23238R-11-03         CR3 HEK.global         3

rownames(traits_HEK.SR71.global)
# [1] "23238R-11-01" "23238R-11-02" "23238R-11-03"

colnames(counts_HEK.SR71.global)
# [1] "23238R-11-01" "23238R-11-02" "23238R-11-03"

# calculate row mean
HEK.SR71.global_row_mean <- rowMeans(counts_HEK.SR71.global)

# subset rows where row means >= 10
HEK.SR71.global_filtered_rows <- counts_HEK.SR71.global[HEK.SR71.global_row_mean >= 10, ]

dim(HEK.SR71.global_filtered_rows)
# [1] 18017     3

# write .csv file 
write.csv(HEK.SR71.global_filtered_rows, file=paste0("1a. Exp. 13-1_RNA-seq_HEK.SR71.global_Venn filter-",dim(HEK.SR71.global_filtered_rows)[1],"x",dim(HEK.SR71.global_filtered_rows)[2],"_CR_08272025.csv"))

################################################################################
# row mean >= 10 HEK.SR64.pulldn
counts_HEK.SR64.pulldn <- raw_counts_clean[,c(9:11)]
dim(counts_HEK.SR64.pulldn)
# [1] 50116     3
traits_HEK.SR64.pulldn <- traits[c(9:11),]
traits_HEK.SR64.pulldn
#              Sample.ID..           Group Replicate
# 23238R-11-09         CR9 HEK.SR64.pulldn         1
# 23238R-11-10        CR10 HEK.SR64.pulldn         2
# 23238R-11-11        CR11 HEK.SR64.pulldn         3

rownames(traits_HEK.SR64.pulldn)
# [1] "23238R-11-09" "23238R-11-10" "23238R-11-11"

colnames(counts_HEK.SR64.pulldn)
# [1] "23238R-11-09" "23238R-11-10" "23238R-11-11"

# calculate row mean
HEK.SR64.pulldn_row_mean <- rowMeans(counts_HEK.SR64.pulldn)

# subset rows where row means >= 10
HEK.SR64.pulldn_filtered_rows <- counts_HEK.SR64.pulldn[HEK.SR64.pulldn_row_mean >= 10, ]

dim(HEK.SR64.pulldn_filtered_rows)
# [1] 18094     3

# write .csv file 
write.csv(HEK.SR71.global_filtered_rows, file=paste0("1a. Exp. 13-1_RNA-seq_HEK.SR64.pulldn_Venn filter-",dim(HEK.SR64.pulldn_filtered_rows)[1],"x",dim(HEK.SR64.pulldn_filtered_rows)[2],"_CR_08272025.csv"))

################################################################################
# row mean >= 10 HEK.SR71.pulldn
counts_HEK.SR71.pulldn <- raw_counts_clean[,c(12:14)]
dim(counts_HEK.SR71.pulldn)
# [1] 50116     3
traits_HEK.SR71.pulldn <- traits[c(12:14),]
traits_HEK.SR71.pulldn
#             Sample.ID..           Group Replicate
# 23238R-11-12        CR12 HEK.SR71.pulldn         1
# 23238R-11-13        CR16 HEK.SR71.pulldn         2
# 23238R-11-14        CR17 HEK.SR71.pulldn         3

rownames(traits_HEK.SR71.pulldn)
# [1] "23238R-11-12" "23238R-11-13" "23238R-11-14"

colnames(counts_HEK.SR71.pulldn)
# [1] "23238R-11-12" "23238R-11-13" "23238R-11-14"

# calculate row mean
HEK.SR71.pulldn_row_mean <- rowMeans(counts_HEK.SR71.pulldn)

# subset rows where row means >= 10
HEK.SR71.pulldn_filtered_rows <- counts_HEK.SR71.pulldn[HEK.SR71.pulldn_row_mean >= 10, ]

dim(HEK.SR71.pulldn_filtered_rows)
# [1] 17210     3

# write .csv file 
write.csv(HEK.SR71.pulldn_filtered_rows, file=paste0("1a. Exp. 13-1_RNA-seq_HEK.SR71.pulldn_Venn filter-",dim(HEK.SR71.pulldn_filtered_rows)[1],"x",dim(HEK.SR71.pulldn_filtered_rows)[2],"_CR_08272025.csv"))

################################################################################
################################################################################
# clean up filtered row dfs to only include gene lists

# HEK.SR71.global_filtered_rows
HEK.SR71.global_genes <- as.matrix(rownames(HEK.SR71.global_filtered_rows))
dim(HEK.SR71.global_genes)
# [1] 18017     1

# HEK.SR64.pulldn_filtered_rows
HEK.SR64.pulldn_genes <- as.matrix(rownames(HEK.SR64.pulldn_filtered_rows))
dim(HEK.SR64.pulldn_genes)
# [1] 18094     1

# HEK.SR71.pulldn_filtered_rows
HEK.SR71.pulldn_genes <- as.matrix(rownames(HEK.SR71.pulldn_filtered_rows))
dim(HEK.SR71.pulldn_genes)
# [1] 17210     1

################################################################################
# Create Venn diagrams for each comparison
list_z <- NULL

# #import fonts to use Arial in output pdf
# https://r-coder.com/custom-fonts-r/
library(extrafont)
font_import()
fonts()
loadfonts(device = "pdf", quiet = TRUE)

# load library to create Venn's
# GUI: https://www.biovenn.nl/
# R package: https://cran.r-project.org/web/packages/BioVenn/index.html
library(BioVenn)

## Create Venn diagrams and output as a pdf ##
pdf(file = "1a. Exp. 13-1_RNA-seq_Venns_CR_08272025.pdf", height = 11, width = 8.5, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# HEK.SR71.global vs HEK.SR64.pulldn
biovenn1 <- draw.venn(HEK.SR71.global_genes, HEK.SR64.pulldn_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="HEK.SR71.global vs HEK.SR64.pulldn RNA", 
                      xtitle = "global", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "HEK.SR64.pulldn",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "black",
                      y_c = "magenta",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x global: 18017"
# [1] "y global: 18094"
# [1] "z global: 0"
# [1] "x only: 844"
# [1] "y only: 921"
# [1] "z only: 0"
# [1] "x-y global overlap: 17173"

# HEK.SR71.global vs HEK.SR71.pulldn
biovenn2 <- draw.venn(HEK.SR71.global_genes, HEK.SR71.pulldn_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="HEK.SR71.global vs HEK.SR71.pulldn RNA", 
                      xtitle = "global", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "HEK.SR71.pulldn",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "black",
                      y_c = "magenta2",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)
# [1] "x global: 18017"
# [1] "y global: 17210"
# [1] "z global: 0"
# [1] "x only: 1440"
# [1] "y only: 633"
# [1] "z only: 0"
# [1] "x-y global overlap: 16577"


# HEK.SR64.pulldn vs HEK.SR71.pulldn
biovenn3 <- draw.venn(HEK.SR64.pulldn_genes, HEK.SR71.pulldn_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="HEK SR64 vs SR71 pulldown RNA", 
                      xtitle = "HEK.SR64.pulldn", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "HEK.SR71.pulldn",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "lightgrey",
                      y_c = "magenta3",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x global: 18094"
# [1] "y global: 17210"
# [1] "z global: 0"
# [1] "x only: 1136"
# [1] "y only: 252"
# [1] "z only: 0"
# [1] "x-y global overlap: 16958"


dev.off()

################################################################################
################################################################################
## Step 1b: Create row mean filter >=10 across all samples for downstream diffex analysis

# calculate row mean
counts_clean_row_mean <- rowMeans(raw_counts_clean)

# subset rows where row means >= 10
counts_clean_filtered_rows <- raw_counts_clean[counts_clean_row_mean  >= 10, ]

dim(counts_clean_filtered_rows)
# 1] 17998    14

# write .csv file 
write.csv(counts_clean_filtered_rows, file=paste0("1b. Exp. 13-1_RNA-seq_clean_counts_diffex filter-",dim(counts_clean_filtered_rows)[1],"x",dim(counts_clean_filtered_rows)[2],"_CR_08272025.csv"))

################################################################################
################################################################################
## Step 2: DESeq2 normalization of counts based on median of ratios

# load library
library(DESeq2)

# create DESeq data set matrix (dds) for all 14 samples
dds_all<- DESeqDataSetFromMatrix(countData = counts_clean_filtered_rows,
                                     colData = traits,
                                     design = ~Group)
dim(dds_all)
# [1] 17998    14

# set the factor level
dds_all$Group <- relevel(dds_all$Group, ref = "HEK.global")

# run DESeq
dds_all<- DESeq(dds_all)

View(counts(dds_all))

dds <- dds_all

dds <- estimateSizeFactors(dds)

sizeFactors(dds)
# 23238R-11-01 23238R-11-02 23238R-11-03 23238R-11-04 23238R-11-05 23238R-11-06 23238R-11-07 23238R-11-08 23238R-11-09 
# 1.0890903    1.1754968    1.2605920    1.0208653    0.8150499    1.5939941    1.1929057    1.4791863    0.8989332 
# 23238R-11-10 23238R-11-11 23238R-11-12 23238R-11-13 23238R-11-14 
# 0.8644851    0.9698238    0.7939392    0.7066748    0.7060415 

normalized_counts_1 <- counts(dds, normalized=TRUE)

write.csv(normalized_counts_1, file=paste0("2a. Exp. 13-1_RNA-seq_clean_filtered_DESeq2 norm counts_diffex filter-",dim(normalized_counts_1)[1],"x",dim(normalized_counts_1)[2],"_CR_08272025.csv"))

##################################################################################
# norm counts without negative pulldowns

# create DESeq data set matrix dds for all 20 samples
dds_nonegs <- DESeqDataSetFromMatrix(countData = counts_clean_filtered_rows_no_negs,
                                     colData = traits_no_negs,
                                     design = ~Group)
dim(dds_nonegs)
# [1] 17998    11

# set the factor level
dds_nonegs$Group <- relevel(dds_nonegs$Group, ref = "HEK.global")

# run DESeq
dds_nonegs <- DESeq(dds_nonegs)

dds_nonegs <- estimateSizeFactors(dds_nonegs)

sizeFactors(dds_nonegs)
# 23238R-11-01 23238R-11-02 23238R-11-03 23238R-11-04 23238R-11-05 23238R-11-09 23238R-11-10 23238R-11-11 23238R-11-12 
# 1.1864666    1.2826704    1.3739823    1.1172610    0.8977354    0.9814663    0.9416285    1.0572049    0.8707841 
# 23238R-11-13 23238R-11-14 
# 0.7768515    0.7763584 

normalized_counts_nonegs <- counts(dds_nonegs, normalized=TRUE)

write.csv(normalized_counts_nonegs, file=paste0("2a. Exp. 13-1_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-",dim(normalized_counts_1)[1],"x",dim(normalized_counts_1)[2],"_CR_08272025.csv"))

################################################################################

## end of data processing and clean up ##

################################################################################
################################################################################
## Part II: Data analysis and visualization
# Exp. 13-1 mRNA-seq analysis

# steps:
# 1. PCAs
  # 1a. all groups
  # 1b. SR64 vs SR71 pulldowns
# 2. Correlations
  # a. HEK.SR71.global vs HEK.SR64.pulldn
  # b. HEK.SR71.global vs HEK.SR71.pulldn
# 3. Diffex volcanoes & GSEA
  # a. HEK.SR71.global vs HEK.SR64.pulldn
  # b. HEK.SR71.global vs HEK.SR71.pulldn
################################################################################
## Step 1: Generate PCA plots
library(DESeq2)
library(ggplot2)

# 1a: all groups
# create DESeq data set matrix dds for all 20 samples
dds_all<- DESeqDataSetFromMatrix(countData = counts_clean_filtered_rows,
                                     colData = traits,
                                     design = ~Group)
dim(dds_all)
# [1] 17998    24

# set the factor level
dds_all$Group <- relevel(dds_all$Group, ref = "HEK.global")

# run DESeq
dds_all<- DESeq(dds_all)

# create pdf output file
pdf(file = "1a-1. Exp. 13-1_RNA-seq_all groups_DESeq2_PCA_CR_08272025.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

# perform DESeq2-basedPCA on rlog transformed data of all 20 samples
rldr <- rlog(dds_all, blind = TRUE)
plotPCA(rldr, intgroup = "Group") # plots based on normalized counts

pcaData <- plotPCA(rldr, intgroup = "Group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca1 <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("DESeq2 PCA plot for log2(mean mRNA counts)") +
  theme(panel.background = element_rect(fill = "transparent", color = "black"), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5)) +  # Center the title
  labs(color = "Transcriptome groups") +
  coord_fixed()

dev.off()

library(limma)

pdf(file = "1a-2. Exp. 13-1_RNA-seq_all groups_Limma_MDS_CR_08272025.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- traits$Group

pch_values <- ifelse(Grouping_vector == "HEK.global", 16,
                     ifelse(Grouping_vector == "HEK.SR64.global", 16,
                            ifelse(Grouping_vector == "HEK.SR71.global", 16,
                                   ifelse(Grouping_vector == "HEK.neg.pulldn", 16, 
                                          ifelse(Grouping_vector == "HEK.SR64.pulldn", 16, 
                                                 ifelse(Grouping_vector == "HEK.SR71.pulldn", 16, NA))))))


pt_colors <- ifelse(Grouping_vector == "HEK.global", "coral",
                    ifelse(Grouping_vector == "HEK.SR64.global", "gold3",
                           ifelse(Grouping_vector == "HEK.SR71.global", "green3",
                                  ifelse(Grouping_vector == "HEK.neg.pulldn", "darkturquoise", 
                                       ifelse(Grouping_vector == "HEK.SR64.pulldn", "dodgerblue", 
                                         ifelse(Grouping_vector == "HEK.SR71.pulldn", "magenta", NA))))))
  
legend_groups <- c("HEK.global", "HEK.SR64.global", "HEK.SR71.global", "HEK.neg.pulldn", "HEK.SR64.pulldn", "HEK.SR71.pulldn")

plotMDS_allgroups <- plotMDS(log2(normalized_counts_1), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for log2(mean mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("coral", "gold3", "green3", "darkturquoise", "dodgerblue", "magenta"), 
       pch = 16, title = "Transcriptome groups",cex=1.4)

dev.off()


################################################################################
# 1b: all groups w/o neg.pulldowns

counts_clean_filtered_rows_no_negs <- counts_clean_filtered_rows[,c(1:5, 9:14)]
dim(counts_clean_filtered_rows_no_negs)
# [1] 17998    11

colnames(counts_clean_filtered_rows_no_negs)
# [1] "23238R-11-01" "23238R-11-02" "23238R-11-03" "23238R-11-04" "23238R-11-05" "23238R-11-09" "23238R-11-10"
# [8] "23238R-11-11" "23238R-11-12" "23238R-11-13" "23238R-11-14"

traits_no_negs <- traits[c(1:5, 9:14),]
traits_no_negs
#              Sample.ID..           Group Replicate
# 23238R-11-01         CR1      HEK.global         1
# 23238R-11-02         CR2      HEK.global         2
# 23238R-11-03         CR3      HEK.global         3
# 23238R-11-04         CR4 HEK.SR64.global         1
# 23238R-11-05         CR5 HEK.SR71.global         1
# 23238R-11-09         CR9 HEK.SR64.pulldn         1
# 23238R-11-10        CR10 HEK.SR64.pulldn         2
# 23238R-11-11        CR11 HEK.SR64.pulldn         3
# 23238R-11-12        CR12 HEK.SR71.pulldn         1
# 23238R-11-13        CR16 HEK.SR71.pulldn         2
# 23238R-11-14        CR17 HEK.SR71.pulldn         3

pdf(file = "1b-1. Exp. 13-1_RNA-seq_all groups_nonegs_DESeq2_PCA_CR_08272025.pdf", width = 8, height = 8.5)
# perform DESeq2-basedPCA on rlog transformed data of 18 samples w/o neg pulldowns

rldr <- rlog(dds_nonegs, blind = TRUE)
plotPCA(rldr, intgroup = "Group") # plots based on normalized counts

pcaData <- plotPCA(rldr, intgroup = "Group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca2 <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("DESeq2 PCA plot for log2(mean mRNA counts)") +
  theme(panel.background = element_rect(fill = "transparent", color = "black"), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5)) +  # Center the title
  labs(color = "Transcriptome groups") +
  coord_fixed()

dev.off()

# perform limma:plotMDS-based PCA of 22 samples w/o neg pulldowns
library(limma)

pdf(file = "1b-2. Exp. 13-1_RNA-seq_all groups_nonegs_Limma_MDS_CR_08272025.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- traits$Group[c(1:5, 9:14)]

pch_values <- ifelse(Grouping_vector == "HEK.global", 16,
                     ifelse(Grouping_vector == "HEK.SR64.global", 16,
                            ifelse(Grouping_vector == "HEK.SR71.global", 16,
                                          ifelse(Grouping_vector == "HEK.SR64.pulldn", 16, 
                                                 ifelse(Grouping_vector == "HEK.SR71.pulldn", 16, NA)))))



pt_colors <- ifelse(Grouping_vector == "HEK.global", "coral",
                    ifelse(Grouping_vector == "HEK.SR64.global", "gold3",
                           ifelse(Grouping_vector == "HEK.SR71.global", "green3",
                                         ifelse(Grouping_vector == "HEK.SR64.pulldn", "dodgerblue", 
                                                ifelse(Grouping_vector == "HEK.SR71.pulldn", "magenta", NA)))))

legend_groups <- c("HEK.global", "HEK.SR64.global", "HEK.SR71.global", "HEK.SR64.pulldn", "HEK.SR71.pulldn")

plotMDS_nonegs <- plotMDS(log2(normalized_counts_nonegs), #top = 500, 
                             labels = NULL, pch = pch_values, col = pt_colors, 
                             cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                             xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for log2(mean mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("coral", "gold3", "green3", "darkturquoise", "dodgerblue", "magenta"), 
       pch = 16, title = "Transcriptome groups",cex=1.4)

dev.off()

################################################################################
# 1c. SR64 vs SR71 pulldowns
# subset DESeq2 norm counts to only include pulldowns

#Updated clean_counts
clean_counts_pos_pulldns <- counts_clean_filtered_rows[,c(9:14)]

#Updated traits
pos_pulldns_traits <- traits[c(9:14),]
pos_pulldns_traits
# Sample.ID..           Group Replicate
# 23238R-11-09         CR9 HEK.SR64.pulldn         1
# 23238R-11-10        CR10 HEK.SR64.pulldn         2
# 23238R-11-11        CR11 HEK.SR64.pulldn         3
# 23238R-11-12        CR12 HEK.SR71.pulldn         1
# 23238R-11-13        CR16 HEK.SR71.pulldn         2
# 23238R-11-14        CR17 HEK.SR71.pulldn         3


colnames(clean_counts_pos_pulldns)
# [1] "23238R-11-09" "23238R-11-10" "23238R-11-11" "23238R-11-12" "23238R-11-13" "23238R-11-14"

#create DESeq data set (dds) matrix for the 3 Aldh1l1 and 3 Camk2a pos.pulldowns
pos_pulls_dds <- DESeqDataSetFromMatrix(countData = clean_counts_pos_pulldns, #filtered, unnorm counts is the input here
                                              colData = pos_pulldns_traits,
                                              design = ~Group)
dim(pos_pulls_dds)
# [1] 17998     6

# set the factor level
pos_pulls_dds$Group <- relevel(pos_pulls_dds$Group, ref = "HEK.SR71.pulldn")

#Run DESeq
pos_pulls_dds <- DESeq(pos_pulls_dds)

# Create PCA plots
pdf(file = "1c-1. Exp. 13-1_RNA-seq_SR64 and SR71 pull groups_DESeq2_PCA_CR_08272025.pdf", width = 8, height = 8.5)
# perform DESeq2-basedPCA on rlog transformed data of two positive pulldown groups

rldr <- rlog(pos_pulls_dds, blind = TRUE)
plotPCA(rldr, intgroup = "Group") # plots based on normalized counts

pcaData <- plotPCA(rldr, intgroup = "Group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca3 <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(panel.background = element_rect(fill = "transparent", color = "black"), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent")) +
  labs(color = "Transcriptome groups") +
  coord_fixed()

dev.off()

# perform limma:plotMDS-based PCA of 18 samples w/o neg pulldowns
pdf(file = "1c-2. Exp. 13-1_RNA-seq_SR64 and SR71 pulldn groups_Limma_MDS_CR_08272025.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector_pulls <- pos_pulldns_traits$Group

pos_pull_norm_counts <- normalized_counts_1[,c(9:14)]

pch_values <- ifelse(Grouping_vector_pulls== "HEK.SR64.pulldn", 16,
                     ifelse(Grouping_vector_pulls== "HEK.SR71.pulldn", 16, NA))


pt_colors <- ifelse(Grouping_vector_pulls== "HEK.SR64.pulldn", "dodgerblue",
                    ifelse(Grouping_vector_pulls== "HEK.SR71.pulldn", "magenta", NA))

legend_groups <- c("HEK.SR64.pulldn", "HEK.SR71.pulldn")

plotMDS_pos_pulldns <- plotMDS(log2(pos_pull_norm_counts), #top = 500, 
                          labels = NULL, pch = pch_values, col = pt_colors, 
                          cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                          xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)
mtext(side=3, text="MDS Plot for log2(mean mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("dodgerblue", "magenta"), 
       pch = 16, title = "Transcriptome groups", cex = 1.4)

dev.off()

################################################################################
################################################################################
# Step 2. Correlations of mean mRNA counts
  # a. HEK.global vs HEK.SR64.pulldn
  # b. HEK.global vs HEK.SR71.pulldn

# check for finite values
any(is.finite(normalized_counts_1))
# [1] TRUE zeroes are considered finite

# CR: Because there are zeroes in the DESeq2 normalized counts data matrix, 
# the log2 transformation changes these values to infinite values (-ifn).
# A correlation coefficient cannot be calculated with infinite values, so
# I am going to add a "pseudocount" of + 1 to all the counts values to
# avoid taking a log of 0.

norm_counts_pseudo <- normalized_counts_1 + 1
any(is.finite(norm_counts_pseudo))

norm_counts_log <- log2(norm_counts_pseudo)
any(is.infinite(norm_counts_log))
# [1] FALSE

# load library to generate plot
library(ggplot2)

################################################################################
# 2a. HEK.global vs HEK.SR64.pulldn

pdf(file = "2a-e. Exp. 13-1_RNA-seq_counts correlations_CR_08272025.pdf", width = 6, height = 6, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# subset into groups of interest and take the row mean
HEK.global <- rowMeans(norm_counts_log[,c(1:3)])
HEK.neg.pulldn <- rowMeans(norm_counts_log[,c(6:8)])
HEK.SR64.pulldn <- rowMeans(norm_counts_log[,c(9:11)])
HEK.SR71.pulldn <- rowMeans(norm_counts_log[,c(12:14)])

# generate a linear correlation
correlation_2a <- cor(HEK.global, HEK.SR64.pulldn)
correlation_2a
# 1] 0.9511533

correlation_2a_test <- cor.test(HEK.global, HEK.SR64.pulldn, method = "pearson")
correlation_2a_test
# 	Pearson's product-moment correlation

# data:  HEK.global and HEK.SR64.pulldn
# t = 413.31, df = 17996, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9497413 0.9525266
# sample estimates:
#  cor 
# 0.9511533

# convert data to a data frame
df_2a <- data.frame(HEK.global = HEK.global, HEK.SR64.pulldn = HEK.SR64.pulldn)

# create correlation plot
correlation_plot_2a <- ggplot(df_2a, aes(x = HEK.global, y = HEK.SR64.pulldn)) +
  geom_point(color = "black", fill = "magenta3", pch = 21, alpha = 0.5, size = 3) +
  labs(
    x = "HEK.global log2(mean mRNA counts)",
    y = "HEK.SR64.pulldn log2(mean mRNA counts)",
    caption = paste("Number of genes:", nrow(df_2a))
  ) +
  ggtitle("2a. HEK.global vs HEK.SR64.pulldn") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0)  # left align caption
  ) +
  annotate(
    "text",
    x = max(df_2a$HEK.global),
    y = min(df_2a$HEK.SR64.pulldn) + 1,
    label = paste("Correlation coefficient:", round(correlation_2a, 2)),
    hjust = 0, vjust = 1, color = "black"
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2a)

################################################################################
# generate a linear correlation
correlation_2b <- cor(HEK.global, HEK.SR71.pulldn)
correlation_2b
# [1] 0.9379205

correlation_2b_test <- cor.test(HEK.global, HEK.SR71.pulldn, method = "pearson")
correlation_2b_test
# 	Pearson's product-moment correlation

# data:  HEK.global and HEK.SR71.pulldn
# t = 362.75, df = 17996, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9361385 0.9396544
# sample estimates:
#   cor 
# 0.9379205 

# convert data to a data frame
df_2b <- data.frame(HEK.global = HEK.global, HEK.SR71.pulldn = HEK.SR71.pulldn)

# create correlation plot
correlation_plot_2b <- ggplot(df_2b, aes(x = HEK.global, y = HEK.SR71.pulldn)) +
  geom_point(color = "black", fill = "magenta4", pch = 21, alpha = 0.5, size = 3) +
  labs(
    x = "HEK.global log2(mean mRNA counts)",
    y = "HEK.SR71.pulldn log2(mean mRNA counts)",
    caption = paste("Number of genes:", nrow(df_2b))
  ) +
  ggtitle("2b. HEK.global vs HEK.SR71.pulldn") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0)  # left align caption
  ) +
  annotate(
    "text",
    x = max(df_2b$HEK.global),
    y = min(df_2b$HEK.SR71.pulldn) + 1,
    label = paste("Correlation coefficient:", round(correlation_2b, 2)),
    hjust = 0, vjust = 1, color = "black"
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2b)

################################################################################
# generate a linear correlation
correlation_2c <- cor(HEK.SR64.pulldn, HEK.SR71.pulldn)
correlation_2c
# [1] 0.9887139

correlation_2c_test <- cor.test(HEK.SR64.pulldn, HEK.SR71.pulldn, method = "pearson")
correlation_2c_test
# 	Pearson's product-moment correlation

# data:  HEK.SR64.pulldn and HEK.SR71.pulldn
# t = 885.32, df = 17996, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9883812 0.9890371
# 3 sample estimates:
#  cor 
# 0.9887139 

# convert data to a data frame
df_2c <- data.frame(HEK.SR64.pulldn = HEK.SR64.pulldn, HEK.SR71.pulldn = HEK.SR71.pulldn)

# create correlation plot
correlation_plot_2c <- ggplot(df_2c, aes(x = HEK.SR64.pulldn, y = HEK.SR71.pulldn)) +
  geom_point(color = "black", fill = "magenta", pch = 21, alpha = 0.5, size = 3) +
  labs(
    x = "HEK.SR64.pulldn log2(mean mRNA counts)",
    y = "HEK.SR71.pulldn log2(mean mRNA counts)",
    caption = paste("Number of genes:", nrow(df_2c))
  ) +
  ggtitle("2c. HEK.SR64.pulldn vs HEK.SR71.pulldn") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0)  # left align caption
  ) +
  annotate(
    "text",
    x = max(df_2c$HEK.SR64.pulldn),
    y = min(df_2c$HEK.SR71.pulldn) + 1,
    label = paste("Correlation coefficient:", round(correlation_2c, 2)),
    hjust = 0, vjust = 1, color = "black"
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2c)

dev.off()

################################################################################
################################################################################
# Step 3 Part I: Diffex volcanoes

# 3. Diffex volcanoes & GSEA
  # C1 HEK.global vs HEK.SR64.pulldn
  # C2 HEK.global vs HEK.SR71.pulldn
  # C3 HEK.SR64.pulldn vs HEK.SR71.pulldn

#load Libraries
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)
library(ggplot2)
library(ggrepel)

#use base pdf function to extract plots into a single pdf
pdf(file = "C1_HEK.SR64.pulldn vs HEK.global RNA_rowmeans_17998_Diffex_08272025.pdf", height = 8, width = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

final_clean_counts <- counts_clean_filtered_rows

#Updated clean_counts for comparison 1
C1_clean_counts <- final_clean_counts[,c(1:3, 9:11)]

#Updated traits for comparison 1
C1_traits <- traits[c(1:3, 9:11),]
C1_traits
#              Sample.ID..           Group Replicate
# 23238R-11-01         CR1      HEK.global         1
# 23238R-11-02         CR2      HEK.global         2
# 23238R-11-03         CR3      HEK.global         3
# 23238R-11-09         CR9 HEK.SR64.pulldn         1
# 23238R-11-10        CR10 HEK.SR64.pulldn         2
# 23238R-11-11        CR11 HEK.SR64.pulldn         3

#Deseq2 data matrix is needed for downstream analysis
dds1 <- DESeqDataSetFromMatrix(countData = C1_clean_counts,
                               colData = C1_traits,
                               design = ~Group)   
dim(dds1)
# [1] 17998     6

dds1

# set the factor level
dds1$Group <- relevel(dds1$Group, ref = "HEK.global")

#Run DESeq 
dds1 <- DESeq(dds1)

#diff p adj value 
res05 <- results(dds1, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 17998 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5754, 32%
# LFC < 0 (down)     : 6003, 33%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 6)

dim(res05)
#[1] 17998     6

#contrasts
resultsNames(dds1)
# [1] "Intercept"                           "Group_HEK.SR64.pulldn_vs_HEK.global"

#plot MA
#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res05, main = "Comparison 1: HEK.SR64.pulldn vs HEK.global")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds1, blind=TRUE)
de1 <- as.data.frame(res05)
de1$padj[is.na(de1$padj)] <- 1

write.csv(de1, file = "C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv")

var <- order(de1$padj, decreasing=F)
#topAdjpGenes <- which(var <= 20)
#allows you to hold your enrs ids and sample names in a df 
#Most sig genes! Top 50
matrix <- assay(rld)[ var[1:50], ]
matrix <- matrix - rowMeans(matrix)

matrixx <- matrix
#thresh_val <- min(abs(max(matrix)),abs(min(matrix)))
thresh_val <- 0.5
matrixx[matrixx < -thresh_val] <- -thresh_val
matrixx[matrixx > thresh_val] <- thresh_val

annotation_data <- as.data.frame(colData(rld)[c("Group")])
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998") 

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998_Diffex_Up_08272025.csv")
write.csv(down, file = "C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998_Diffex_Down_08272025.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998_Diffex_Up_08272025.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998_Diffex_Down_08272025.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 3749    6

#volcano plot of diffex with padj df
de1 <- read.csv("C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv")

de1$diffexpressed <- "NS"
de1$diffexpressed[de1$log2FoldChange > 1.0 & de1$padj < 0.05] <- "UP"
de1$diffexpressed[de1$log2FoldChange < -1.0 & de1$padj < 0.05] <- "DOWN"

de1$delabel <- NA
de1$delabel[de1$diffexpressed != "NS"] <- de1$X[de1$diffexpressed != "NS"]

ggplot(data=de1, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(-5, 300)) +
  scale_x_continuous(name = "log2(difference) 
HEK.SR64.pulldn - HEK.global",
limits = c(-5, 5))

library(ggbreak)

ggplot(data=de1, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle("C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_x_continuous(name = "log2(difference)\nHEK.SR64.pulldn - HEK.global",
                     limits = c(-5, 5)) +
  scale_y_continuous(name = "-log10(padj)") +
  scale_y_break(c(300, 1000))   # axis break between 50 and 5000

summary(-log10(de1$padj))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000012  0.729581  2.908815       Inf 11.469127       Inf 

quantile(-log10(de1$padj), probs = c(0.95, 0.99, 0.995, 1), na.rm = TRUE)
#       95%       99%     99.5%      100% 
# 54.05938 119.66406 151.56848       Inf 

#Perform PCA on rlog transformed data
rldr <- rlog(dds1, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

######################################################################################################################################
######################################################################################################################################
#Comparison 2 - C2 signifies the below mentioned comparison!

# C2 HEK.global vs HEK.SR71.pulldn

#use base pdf function to extract plots into a single pdf
pdf(file = "C2_HEK.SR71.pulldn vs HEK.global RNA_rowmeans_17998_Diffex_08272025.pdf", height = 8, width = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

#Updated clean_counts for comparison 1
C2_clean_counts <- final_clean_counts[,c(1:3, 12:14)]

#Updated traits for comparison 1
C2_traits <- traits[c(1:3, 12:14),]
C2_traits
#              Sample.ID..           Group Replicate
# 23238R-11-01         CR1      HEK.global         1
# 23238R-11-02         CR2      HEK.global         2
# 23238R-11-03         CR3      HEK.global         3
# 23238R-11-12        CR12 HEK.SR71.pulldn         1
# 23238R-11-13        CR13 HEK.SR71.pulldn         2
# 23238R-11-14        CR14 HEK.SR71.pulldn         3

#Deseq2 data matrix is needed for downstream analysis
dds2 <- DESeqDataSetFromMatrix(countData = C2_clean_counts,
                               colData = C2_traits,
                               design = ~Group)   
dim(dds2)
# [1] 17998     6

dds2

# set the factor level
dds2$Group <- relevel(dds2$Group, ref = "HEK.global")

#Run DESeq 
dds2 <- DESeq(dds2)

#diff p adj value 
res05 <- results(dds2, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 17998 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5695, 32%
# LFC < 0 (down)     : 5700, 32%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 4)

dim(res05)
#[1] 17998     6

#contrasts
resultsNames(dds2)
# [1] "Intercept"                           "Group_HEK.SR71.pulldn_vs_HEK.global"

#plot MA
#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res05, main = "Comparison 2: HEK.SR71.pulldn vs HEK.global")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds2, blind=TRUE)
de2 <- as.data.frame(res05)
de2$padj[is.na(de2$padj)] <- 1

write.csv(de2, file = "C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv")

var <- order(de2$padj, decreasing=F)
#topAdjpGenes <- which(var <= 20)
#allows you to hold your enrs ids and sample names in a df 
#Most sig genes! Top 50
matrix <- assay(rld)[ var[1:50], ]
matrix <- matrix - rowMeans(matrix)

matrixx <- matrix
#thresh_val <- min(abs(max(matrix)),abs(min(matrix)))
thresh_val <- 0.5
matrixx[matrixx < -thresh_val] <- -thresh_val
matrixx[matrixx > thresh_val] <- thresh_val

annotation_data <- as.data.frame(colData(rld)[c("Group")])
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998") 

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998_Diffex_Up_08272025.csv")
write.csv(down, file = "C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998_Diffex_Down_08272025.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998_Diffex_Up_08272025.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998_Diffex_Down_08272025.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 4557    6

#volcano plot of diffex with padj df
de2 <- read.csv("C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv")

de2$diffexpressed <- "NS"
de2$diffexpressed[de2$log2FoldChange > 1.0 & de2$padj < 0.05] <- "UP"
de2$diffexpressed[de2$log2FoldChange < -1.0 & de2$padj < 0.05] <- "DOWN"

de2$delabel <- NA
de2$delabel[de2$diffexpressed != "NS"] <- de2$X[de2$diffexpressed != "NS"]

ggplot(data=de2, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(-5, 300)) +
  scale_x_continuous(name = "log2(difference) 
HEK.SR71.pulldn - HEK.global",
                     limits = c(-5, 5))

library(ggbreak)

ggplot(data=de2, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle("C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_x_continuous(name = "log2(difference)\nHEK.SR71.pulldn - HEK.global",
                     limits = c(-5, 5)) +
  scale_y_continuous(name = "-log10(padj)") +
  scale_y_break(c(300, 1000))   # axis break between 50 and 5000

summary(-log10(de2$padj))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000004 0.651862 2.500473      Inf 8.265511      Inf 

quantile(-log10(de2$padj), probs = c(0.95, 0.99, 0.995, 1), na.rm = TRUE)
#       95%       99%     99.5%      100% 
# 34.42094 72.24957 92.09444      Inf 

#Perform PCA on rlog transformed data
rldr <- rlog(dds2, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()
######################################################################################################################################
######################################################################################################################################
#Comparison 3 - C3 signifies the below mentioned comparison!
# C3 HEK.SR64.pulldn vs HEK.SR71.pulldn
#Hypothesis: 

#use base pdf function to extract plots into a single pdf
pdf(file = "C3_HEK.SR64.pulldn vs HEK.SR71.pulldn RNA_rowmeans_17998_Diffex_08272025.pdf", height = 8, width = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

#Updated clean_counts for comparison 1
C3_clean_counts <- final_clean_counts[,c(9:14)]

#Updated traits for comparison 1
C3_traits <- traits[c(9:14),]
C3_traits
#              Sample.ID..           Group Replicate
# 23238R-11-09         CR9 HEK.SR64.pulldn         1
# 23238R-11-10        CR10 HEK.SR64.pulldn         2
# 23238R-11-11        CR11 HEK.SR64.pulldn         3
# 23238R-11-12        CR12 HEK.SR71.pulldn         1
# 23238R-11-13        CR13 HEK.SR71.pulldn         2
# 23238R-11-14        CR14 HEK.SR71.pulldn         3

#Deseq2 data matrix is needed for downstream analysis
dds3 <- DESeqDataSetFromMatrix(countData = C3_clean_counts,
                               colData = C3_traits,
                               design = ~Group)   
dim(dds3)
# [1] 17998     6

dds3

# set the factor level
dds3$Group <- relevel(dds3$Group, ref = "HEK.SR71.pulldn")

#Run DESeq 
dds3 <- DESeq(dds3)

#diff p adj value 
res05 <- results(dds3, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 17998 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1258, 7%
# LFC < 0 (down)     : 2056, 11%
# outliers [1]       : 0, 0%
# low counts [2]     : 2443, 14%
# (mean count < 19)

dim(res05)
#[1] 17998     6

#contrasts
resultsNames(dds3)
# [1] "Intercept"                           "Group_HEK.SR64.pulldn_vs_HEK.SR71.pulldn"

#plot MA
#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res05, main = "Comparison 3: HEK.SR64.pulldn vs HEK.SR71.pulldn")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds3, blind=TRUE)
de3 <- as.data.frame(res05)
de3$padj[is.na(de3$padj)] <- 1

write.csv(de3, file = "C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998_Diffex_08272025.csv")

var <- order(de3$padj, decreasing=F)
#topAdjpGenes <- which(var <= 20)
#allows you to hold your enrs ids and sample names in a df 
#Most sig genes! Top 50
matrix <- assay(rld)[ var[1:50], ]
matrix <- matrix - rowMeans(matrix)

matrixx <- matrix
#thresh_val <- min(abs(max(matrix)),abs(min(matrix)))
thresh_val <- 0.5
matrixx[matrixx < -thresh_val] <- -thresh_val
matrixx[matrixx > thresh_val] <- thresh_val

annotation_data <- as.data.frame(colData(rld)[c("Group")])
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998") 

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998_Diffex_Up_08272025.csv")
write.csv(down, file = "C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998_Diffex_Down_08272025.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998_Diffex_Up_08272025.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998_Diffex_Down_08272025.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 161   6

#volcano plot of diffex with padj df
de3 <- read.csv("C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998_Diffex_08272025.csv")

de3$diffexpressed <- "NS"
de3$diffexpressed[de3$log2FoldChange > 1.0 & de3$padj < 0.05] <- "UP"
de3$diffexpressed[de3$log2FoldChange < -1.0 & de3$padj < 0.05] <- "DOWN"

de3$delabel <- NA
de3$delabel[de3$diffexpressed != "NS"] <- de3$X[de3$diffexpressed != "NS"]

ggplot(data=de3, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(-5, 80)) +
  scale_x_continuous(name = "log2(difference) 
HEK.SR64.pulldn - HEK.SR71.pulldn",
                     limits = c(-5, 5))

#Perform PCA on rlog transformed data 
rldr <- rlog(dds3, blind = TRUE)
plotPCA(rldr, intgroup = "Group")


dev.off()

######################################################################################################################################
######################################################################################################################################
dev.off()

################################################################################
################################################################################
# Step 3 Part II: Diffex GSEA

# C1 HEK.global vs HEK.SR64.pulldn
# C2 HEK.global vs HEK.SR71.pulldn
# C3 HEK.SR64.pulldn vs HEK.SR71.pulldn

################################################################################
# Gene set erichment analysis with Eric Dammer's GOparallel pipeline
# https://github.com/edammer/GOparallel
# performs hypergeometric overlap with Fisher Exact Test for enrichment p<0.05 and 5 minimum genes per ontology
# currently outputs files that can be used as GO-Elite input, not read by this script...
# outputs Z score barplots for each input list (module, or DiffEx Protein list) across 6 ontology types in the larger GMT files available from the Bader Lab website.


# Parameters for GOparallel GSEA on volcano-defined DEG lists  
source("GOparallel-FET.R")

options(stringsAsFactors=FALSE)
library(NMF)

## Parameters saved to variables in memory.
filePath <- getwd() 
#Folder that (may) contain the input file specified above, and which will contain the outFilename project Folder.
modulesInMemory=FALSE
ANOVAgroups=TRUE  #change to FALSE if using inputFile
#if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
#inputFile will be ignored
maxBarsPerOntology=6

panelDimensions=c(2,2)    #dimensions of the individual parblots within a page of the main barplot PDF output
pageDimensions=c(8.5,11)  #main barplot PDF output page dimensions, in inches
parallelThreads=7 # need at least 2 threads to run
removeRedundantGOterms=TRUE
#if true, the 3 GO ontology types are collapased into a minimal set of less redundant terms using the below OBO file
cocluster=TRUE
#If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)

# "Clean up" to run before GOPar
rm(ANOVAout)

###############################################################################

# Relevant diffex comparisons:
  # C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv
  # C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv
  # C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998_Diffex_08272025.csv

###############################################################################
# C1 HEK.SR64.pulldn vs HEK.global

# organize DESeq2 diffex output file to mirror LFQ-MS ANOVAout file
C1_ANOVA<-read.csv("C1_HEK.SR64.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv", header=TRUE, row.names=1)
C1_ANOVAout<-C1_ANOVA[,c(4,1,5,2)]
colnames(C1_ANOVAout)[c(3,4)]<-c("HEK.SR64.pulldn-HEK.global","diff HEK.SR64.pulldn.HEK.global")
write.csv(C1_ANOVAout, file = "C1_RNA_ANOVAout_HEK.SR64.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv")

# GOpar
flip=c(0)
outFilename <- "C1_RNA_HEK.SR64.pulldn_vs_HEK.global_pval_DEGs-GO"
GMTdatabaseFile="Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt"
GOparallel(C1_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

`###############################################################################
# C2 HEK.global vs HEK.SR71.pulldn

# organize DESeq2 diffex output file to mirror LFQ-MS ANOVAout file
C2_ANOVA<-read.csv("C2_HEK.SR71.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv", header=TRUE, row.names=1)
C2_ANOVAout<-C2_ANOVA[,c(4,1,5,2)]
colnames(C2_ANOVAout)[c(3,4)]<-c("HEK.SR71.pulldn-HEK.global","diff HEK.SR71.pulldn.HEK.global")
write.csv(C2_ANOVAout, file = "C2_RNA_ANOVAout_HEK.SR71.pulldn vs HEK.global_rowmeans_17998_Diffex_08272025.csv")

# GOpar
flip=c(0)
outFilename <- "C2_RNA_HEK.SR71.pulldn_vs_HEK.global_pval_DEGs-GO"
GMTdatabaseFile="Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt"
GOparallel(C2_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.


###############################################################################
# C3 HEK.SR64.pulldn vs HEK.SR71.pulldn

# organize DESeq2 diffex output file to mirror LFQ-MS ANOVAout file
C3_ANOVA<-read.csv("C3_HEK.SR64.pulldn vs HEK.SR71.pulldn_rowmeans_17998_Diffex_08272025.csv", header=TRUE, row.names=1)
C3_ANOVAout<-C3_ANOVA[,c(4,1,5,2)]
colnames(C3_ANOVAout)[c(3,4)]<-c("HEK.64.pulldn-HEK.71.pulldn","diff HEK.64.pulldn.HEK.71.pulldn")
write.csv(C3_ANOVAout, file = "C3_RNA_ANOVAout_HEK.64.pulldn vs HEK.71.pulldn_rowmeans_17998_Diffex_08272025.csv")

# GOpar
flip=c(0)
outFilename <- "C3_RNA_HEK.64.pulldn_vs_HEK.71.pulldn_pval_DEGs-GO"
GMTdatabaseFile="Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt"
GOparallel(C3_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

###############################################################################

# end of data analysis and visualization
