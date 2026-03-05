########################################################################################################
# Code for performing Differential Gene Expression Analysis on Exp. 1-4 (BV2)
# row means >=10 filter before dds and splicing
# Christina Ramelow, MS
# Date: 02/03/2024
########################################################################################################

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

setwd("/Users/christina/Desktop/TurboID dual-omics/Exp. 1-4/3. RNA-seq")

# load feature_counts.txt file
raw_counts <- read.table("Exp. 1-4, 9-2 and 10-1_feature_counts.txt", sep = "\t", header = TRUE, row.names = 1)
dim(raw_counts)
# [1] 35976    55

# remove the first 5 columns that are primarily categorical
raw_counts_2 <- raw_counts[,c(6:55)]
dim(raw_counts_2)
# [1] 35976    50

# remove the smRNA-seq samples that were accidentally included
raw_counts_3 <- raw_counts_2[,c(1:20, 27:50)]
dim(raw_counts_3)
# [1] 35976    44
# 21047FL-134-01-22 - 21047FL-134-01-27 should be removed. CR: checked

# check if there are NAs in the counts df
sum(is.na(raw_counts_3))
#[1] 0

# load RNA traits file
traits <- read.csv("Exp. 1-4, 9-2 and & 10-1 combined traits_RNA.csv",  header = TRUE, row.names = 1)
dim(traits)
# [1] 44  2

# rename the column name of raw_counts_3 df to match traits df
# Extract row names from traits
new_col_names <- rownames(traits)

# Assign row names of traits as column names of raw_counts_3
names(raw_counts_3) <- new_col_names

# check that the row names of the traits df and the column names of the counts df are the same
raw_counts_clean <- raw_counts_3[,match(rownames(traits), colnames(raw_counts_3))]

#sanity check 1 colnames in counts file = rownames in traits file
all(colnames(raw_counts_clean) %in% rownames(traits))
#TRUE

#sanity check 2 colnames in counts file = rownames in traits file
all(colnames(raw_counts_clean) == rownames(traits))
#TRUE

# splice the raw_counts_clean df to separate Exp. 1-4, 9-2 and 10-1
# Exp. 1-4 spliced
raw_counts_clean_1_4 <- raw_counts_clean[,c(1:20)]
dim(raw_counts_clean_1_4)
# [1] 35976    20

colnames(raw_counts_clean_1_4)
# should be "21047FL-134-01-01" - "21047FL-134-01-21" CR: checked

# write .csv file for raw_counts_clean for Exp. 1-4
write.csv(raw_counts_clean_1_4, file=paste0("0.Exp. 1-4_RNA-seq_raw_counts_clean_unfiltered-",dim(raw_counts_clean_1_4)[1],"x",dim(raw_counts_clean_1_4)[2],"_CR_02032024.csv"))

# Exp. 9-2 and 10-1 splced
raw_counts_clean_9_2_10_1 <- raw_counts_clean[,c(21:44)]
dim(raw_counts_clean_9_2_10_1)
# [1] 35976    24

colnames(raw_counts_clean_9_2_10_1)
# should be "21047FL-134-01-28" - "22082R-12-03"  CR: checked

# write .csv file for raw_counts_clean for Exp. 9-2 and 10-1
write.csv(raw_counts_clean_9_2_10_1, file=paste0("0.Exp. 9-2 and 10-1_RNA-seq_raw_counts_clean_unfiltered-",dim(raw_counts_clean_9_2_10_1)[1],"x",dim(raw_counts_clean_9_2_10_1)[2],"_CR_02032024.csv"))

################################################################################
################################################################################
## Step 1a: Filter out low counts using a row mean >/= 10 by group for Venn diagrams
# Exp. 1-4
# 1. subset the traits and counts dfs by sample group to create list of genes for Venn diagrams
# 2. calculate the row mean
# 3. filter for genes with a row mean >= 10
# 4. write a .csv file for each sample group
################################################################################
# subset the traits file to only include Exp. 1-4
traits_1_4 <- traits[c(1:20),]
dim(traits_1_4)
# [1] 20  2
################################################################################
# row mean >= 10 BV2T.global
counts_BV2T.global <- raw_counts_clean_1_4[,c(7:9)]
dim(counts_BV2T.global)
# [1] 35976     3
traits_BV2T.global <- traits_1_4[c(7:9),]
rownames(traits_BV2T.global)
# [1] "21047FL-134-01-07" "21047FL-134-01-08" "21047FL-134-01-09"
colnames(counts_BV2T.global)
# [1] "21047FL-134-01-07" "21047FL-134-01-08" "21047FL-134-01-09"

# calculate row mean
BV2T.global_row_mean <- rowMeans(counts_BV2T.global)

# subset rows where row means >= 10
BV2T.global_filtered_rows <- counts_BV2T.global[BV2T.global_row_mean >= 10, ]

dim(BV2T.global_filtered_rows)
# [1] 12371     3

# write .csv file 
write.csv(BV2T.global_filtered_rows, file=paste0("1a.Exp. 1-4_RNA-seq_BV2T.global_Venn filter-",dim(BV2T.global_filtered_rows)[1],"x",dim(BV2T.global_filtered_rows)[2],"_CR_02032024.csv"))

################################################################################
# row mean >= 10 BV2T.LPS.global
counts_BV2T.LPS.global <- raw_counts_clean_1_4[,c(10:12)]
dim(counts_BV2T.LPS.global)
# [1] 35976     3
traits_BV2T.LPS.global <- traits_1_4[c(10:12),]
rownames(traits_BV2T.LPS.global)
# [1] "21047FL-134-01-10" "21047FL-134-01-11" "21047FL-134-01-12"
colnames(counts_BV2T.LPS.global)
# [1] "21047FL-134-01-10" "21047FL-134-01-11" "21047FL-134-01-12"

# calculate row mean
BV2T.LPS.global_row_mean <- rowMeans(counts_BV2T.LPS.global)

# subset rows where row means >= 10
BV2T.LPS.global_filtered_rows <- counts_BV2T.LPS.global[BV2T.LPS.global_row_mean >= 10, ]

dim(BV2T.LPS.global_filtered_rows)
# [1] 12183     3

# write .csv file 
write.csv(BV2T.global_filtered_rows, file=paste0("1a.Exp. 1-4_RNA-seq_BV2T.LPS.global_Venn filter-",dim(BV2T.LPS.global_filtered_rows)[1],"x",dim(BV2T.LPS.global_filtered_rows)[2],"_CR_02032024.csv"))

################################################################################
# row mean >= 10 BV2T.pos.pulldown
counts_BV2T.pos.pulldown <- raw_counts_clean_1_4[,c(15:17)]
dim(counts_BV2T.pos.pulldown)
# [1] 35976     3
traits_BV2T.pos.pulldown <- traits_1_4[c(15:17),]
rownames(traits_BV2T.pos.pulldown)
# [1] "21047FL-134-01-16" "21047FL-134-01-17" "21047FL-134-01-18"
colnames(counts_BV2T.pos.pulldown)
# [1] "21047FL-134-01-16" "21047FL-134-01-17" "21047FL-134-01-18"

# calculate row mean
BV2T.pos.pulldown_row_mean <- rowMeans(counts_BV2T.pos.pulldown)

# subset rows where row means >= 10
BV2T.pos.pulldown_filtered_rows <- counts_BV2T.pos.pulldown[BV2T.pos.pulldown_row_mean >= 10, ]

dim(BV2T.pos.pulldown_filtered_rows)
# [1] 11964     3

# write .csv file 
write.csv(BV2T.pos.pulldown_filtered_rows, file=paste0("1a.Exp. 1-4_RNA-seq_BV2T.pos.pulldown_Venn filter-",dim(BV2T.pos.pulldown_filtered_rows)[1],"x",dim(BV2T.pos.pulldown_filtered_rows)[2],"_CR_02032024.csv"))
################################################################################
# row mean >= 10 BV2T.LPS.pos.pulldown
counts_BV2T.LPS.pos.pulldown <- raw_counts_clean_1_4[,c(18:20)]
dim(counts_BV2T.LPS.pos.pulldown)
# [1] 35976     3
traits_BV2T.LPS.pos.pulldown<- traits_1_4[c(18:20),]
rownames(traits_BV2T.LPS.pos.pulldown)
# 1] "21047FL-134-01-19" "21047FL-134-01-20" "21047FL-134-01-21"
colnames(counts_BV2T.LPS.pos.pulldown)
# 1] "21047FL-134-01-19" "21047FL-134-01-20" "21047FL-134-01-21"

# calculate row mean
BV2T.LPS.pos.pulldown_row_mean <- rowMeans(counts_BV2T.LPS.pos.pulldown)

# subset rows where row means >= 10
BV2T.LPS.pos.pulldown_filtered_rows <- counts_BV2T.LPS.pos.pulldown[BV2T.LPS.pos.pulldown_row_mean >= 10, ]

dim(BV2T.LPS.pos.pulldown_filtered_rows)
# [1] 12346     3

# write .csv file 
write.csv(BV2T.LPS.pos.pulldown_filtered_rows, file=paste0("1a.Exp. 1-4_RNA-seq_BV2T.LPS.pos.pulldown_Venn filter-",dim(BV2T.LPS.pos.pulldown_filtered_rows)[1],"x",dim(BV2T.LPS.pos.pulldown_filtered_rows)[2],"_CR_02032024.csv"))

################################################################################
################################################################################
# clean up filtered row dfs to only include gene lists

# BV2T.global_filtered_rows
BV2T.global_genes <- as.matrix(rownames(BV2T.global_filtered_rows))
dim(BV2T.global_genes)
# [1] 12371     1

# BV2T.LPS.global_filtered_rows
BV2T.LPS.global_genes <- as.matrix(rownames(BV2T.LPS.global_filtered_rows))
dim(BV2T.LPS.global_genes)
# [1] 12183     1

# BV2T.pos.pulldown_filtered_rows
BV2T.pos.pulldown_genes <- as.matrix(rownames(BV2T.pos.pulldown_filtered_rows))
dim(BV2T.pos.pulldown_genes)
# [1] 11964     1

# BV2T.LPS.pos.pulldown_filtered_rows
BV2T.LPS.pos.pulldown_genes <- as.matrix(rownames(BV2T.LPS.pos.pulldown_filtered_rows))
dim(BV2T.LPS.pos.pulldown_genes)
# [1] 12346     1

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
pdf(file = "1a. Exp. 1-4__RNA-seq_Venns_CR_02042024.pdf", height = 11, width = 8.5, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# BV2T.global vs BV2T.LPS.global
biovenn1 <- draw.venn(BV2T.global_genes, BV2T.LPS.global_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="BV2T global vs BV2T LPS global RNA", 
                      xtitle = "global", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "LPS global",
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
                      y_c = "darkgreen",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x global: 12371"
# [1] "y global: 12183"
# [1] "z global: 0"
# [1] "x only: 584"
# [1] "y only: 396"
# [1] "z only: 0"
# [1] "x-y global overlap: 11787"
# [1] "x-z global overlap: 0"
# [1] "y-z global overlap: 0"
# [1] "x-y only overlap: 11787"

# BV2T.global vs BV2T.pos.pulldown
biovenn2 <- draw.venn(BV2T.global_genes, BV2T.pos.pulldown_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="BV2T global vs BV2T pulldown RNA", 
                      xtitle = "global", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "pulldown",
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
                      y_c = "lightgrey",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)
# [1] "x global: 12371"
# [1] "y global: 11964"
# [1] "z global: 0"
# [1] "x only: 649"
# [1] "y only: 242"
# [1] "z only: 0"
# [1] "x-y global overlap: 11722"
# [1] "x-z global overlap: 0"
# [1] "y-z global overlap: 0"
# [1] "x-y only overlap: 11722"


# BV2T.pos.pulldown vs BV2T.LPS pulldown
biovenn3 <- draw.venn(BV2T.pos.pulldown_genes, BV2T.LPS.pos.pulldown_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="BV2T vs BV2T LPS pulldown RNA", 
                      xtitle = "pulldown", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "LPS pulldown",
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
                      y_c = "lightgreen",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x global: 11964"
# [1] "y global: 12346"
# [1] "z global: 0"
# [1] "x only: 295"
# [1] "y only: 677"
# [1] "z only: 0"
# [1] "x-y global overlap: 11669"
# [1] "x-z global overlap: 0"
# [1] "y-z global overlap: 0"
# [1] "x-y only overlap: 11669"


# BV2T.LPS.global vs BV2T.LPS.pos.pulldown
biovenn4 <- draw.venn(BV2T.LPS.global_genes, BV2T.LPS.pos.pulldown_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="BV2T LPS global vs BV2T LPS pulldown RNA", 
                      xtitle = "LPS global", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "LPS pulldown",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "darkgreen",
                      y_c = "lightgreen",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x global: 12183"
# [1] "y global: 12346"
# [1] "z global: 0"
# [1] "x only: 236"
# [1] "y only: 399"
# [1] "z only: 0"
# [1] "x-y global overlap: 11947"
# [1] "x-z global overlap: 0"
# [1] "y-z global overlap: 0"
# [1] "x-y only overlap: 11947"

dev.off()

################################################################################
################################################################################
## Step 1b: Create row mean filter >=10 across all samples for downstream diffex analysis

# calculate row mean
counts_clean_row_mean <- rowMeans(raw_counts_clean_1_4)

# subset rows where row means >= 10
counts_clean_filtered_rows_1_4 <- raw_counts_clean_1_4[counts_clean_row_mean  >= 10, ]

dim(counts_clean_filtered_rows_1_4)
# [1] 12742    20

# write .csv file 
write.csv(counts_clean_filtered_rows_1_4, file=paste0("1b.Exp. 1-4_RNA-seq_clean_counts_diffex filter-",dim(counts_clean_filtered_rows_1_4)[1],"x",dim(counts_clean_filtered_rows_1_4)[2],"_CR_02032024.csv"))

################################################################################
################################################################################
## Step 2: DESeq2 normalization of counts based on median of ratios

# load library
library(DESeq2)

# create DESeq data set matrix (dds) for all 20 samples
dds_all1_4 <- DESeqDataSetFromMatrix(countData = counts_clean_filtered_rows_1_4,
                                     colData = traits_1_4,
                                     design = ~Group)
dim(dds_all1_4)
# [1] 12742    20

# set the factor level
dds_all1_4$Group <- relevel(dds_all1_4$Group, ref = "BV2.global")

# run DESeq
dds_all1_4<- DESeq(dds_all1_4)

View(counts(dds_all1_4))

dds <- dds_all1_4

dds <- estimateSizeFactors(dds)

sizeFactors(dds)
# 21047FL-134-01-01 21047FL-134-01-02 21047FL-134-01-03 21047FL-134-01-04 21047FL-134-01-05 21047FL-134-01-06 
# 1.7276857         1.2474688         1.3044619         1.1194946         1.9957233         1.3747657 
# 21047FL-134-01-07 21047FL-134-01-08 21047FL-134-01-09 21047FL-134-01-10 21047FL-134-01-11 21047FL-134-01-12 
# 1.2796322         0.8957683         0.8664624         0.9231729         0.7790115         0.6659619 
# 21047FL-134-01-13 21047FL-134-01-14 21047FL-134-01-16 21047FL-134-01-17 21047FL-134-01-18 21047FL-134-01-19 
# 0.9239079         0.7433480         0.8352943         0.7899724         0.7779676         1.0221982 
# 21047FL-134-01-20 21047FL-134-01-21 
# 1.0045302         0.8827356 

normalized_counts_1 <- counts(dds, normalized=TRUE)

write.csv(normalized_counts_1, file=paste0("2a.Exp. 1-4_RNA-seq_clean_filtered_DESeq2 norm counts_diffex filter-",dim(normalized_counts_1)[1],"x",dim(normalized_counts_1)[2],"_CR_02032024.csv"))

################################################################################
# updated clean_counts with neg.pulldowns removed
all1_4_no_negs_clean_counts <- counts_clean_filtered_rows_1_4[,c(1:12, 15:20)]

# udated traits with neg.pulldowns removed
all1_4_no_negs_traits <- traits[c(1:12, 15:20),]

# create DESeq data set matrix dds for all 18 samples
dds_all1_4_no_negs <- DESeqDataSetFromMatrix(countData = all1_4_no_negs_clean_counts,
                                             colData = all1_4_no_negs_traits,
                                             design = ~Group)
dim(dds_all1_4_no_negs)
# [1] 12742    18

# set the factor level
dds_all1_4_no_negs$Group <- relevel(dds_all1_4_no_negs$Group, ref = "BV2.global")

# run DESeq
dds_all1_4_no_negs<- DESeq(dds_all1_4_no_negs)

View(counts(dds_all1_4_no_negs))

dds <- dds_all1_4_no_negs

dds <- estimateSizeFactors(dds)

sizeFactors(dds)
# 21047FL-134-01-01 21047FL-134-01-02 21047FL-134-01-03 21047FL-134-01-04 21047FL-134-01-05 21047FL-134-01-06 
# 1.6876309         1.2211140         1.2732130         1.0953255         1.9517175         1.3441513 
# 21047FL-134-01-07 21047FL-134-01-08 21047FL-134-01-09 21047FL-134-01-10 21047FL-134-01-11 21047FL-134-01-12 
# 1.2481733         0.8759875         0.8472266         0.9007616         0.7610231         0.6489440 
# 21047FL-134-01-16 21047FL-134-01-17 21047FL-134-01-18 21047FL-134-01-19 21047FL-134-01-20 21047FL-134-01-21 
# 0.8176543         0.7708591         0.7594104         1.0010039         0.9801711         0.8630847 

normalized_counts_2 <- counts(dds, normalized=TRUE)

write.csv(normalized_counts_2, file=paste0("2a.Exp. 1-4_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-",dim(normalized_counts_2)[1],"x",dim(normalized_counts_2)[2],"_CR_02032024.csv"))

## end of data processing and clean up ##

################################################################################
################################################################################
## Part II: Data analysis and visualization
# Exp. 1-4 mRNA-seq analysis

# steps:
# 1. PCAs
    # 1a. all groups
    # 1b. all groups w/o neg.pulldowns
    # 1c. BV2T.global +/- LPS (n=3/group) vs BV2T.pos.pulldown (n=3/group)
# 2. Correlations
    # a. BV2T.global vs BV2T.LPS.global
    # b. BV2T.global vs BV2T.pos.pulldown
    # c. BV2T.global.LPS vs BV2T.LPS.pos.pulldown
    # d. BV2T.LPS.pos.pulldown vs BV2T.pos.pulldown
# 3. Diffex volcanoes & GSEA
    # a. BV2.LPS.global vs BV2.global
    # b. BV2T.global vs BV2.global
    # c. BV2T.LPS.global vs BV2.LPS.global
    # d. BV2T.LPS.global vs BV2T.global
    # e. BV2T.pos.pulldown vs BV2T.global
    # f. BV2T.LPS.pos.pulldown vs BV2T.pos.pulldown
    # g. BV2T.LPS.pos.pulldown vs BV2.LPS.neg.pulldown
    # h. BV2T.LPS.pos.pulldown vs BV2T.LPS.global

################################################################################
## Step 1: Generate PCA plots
library(DESeq2)
library(ggplot2)

# 1a: all groups
# create DESeq data set matrix dds for all 20 samples
dds_all1_4 <- DESeqDataSetFromMatrix(countData = counts_clean_filtered_rows_1_4,
                                             colData = traits_1_4,
                                             design = ~Group)
dim(dds_all1_4)
# [1] 12742    20

# set the factor level
dds_all1_4$Group <- relevel(dds_all1_4$Group, ref = "BV2.global")

# run DESeq
dds_all1_4<- DESeq(dds_all1_4)

# create pdf output file
pdf(file = "1a-1. Exp. 1-4_RNA-seq_all groups_DESeq2_PCA_CR_02032024.pdf", width = 8, height = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

# perform DESeq2-basedPCA on rlog transformed data of all 20 samples
rldr <- rlog(dds_all1_4, blind = TRUE)
plotPCA(rldr, intgroup = "Group") # plots based on normalized counts

pcaData <- plotPCA(rldr, intgroup = "Group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca1 <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("DESeq2 PCA plot for log2(mRNA counts)") +
  theme(panel.background = element_rect(fill = "transparent", color = "black"), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5)) +  # Center the title
  labs(color = "Transcriptome groups") +
  coord_fixed()

dev.off()

pdf(file = "1a-2. Exp. 1-4_RNA-seq_all groups_Limma_PCA_CR_02032024.pdf", width = 8, height = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- read.csv("Exp. 1-4_Grouping vector_with negs.csv")

pch_values <- ifelse(Grouping_vector == "BV2.global", 16,
                     ifelse(Grouping_vector == "BV2.LPS.global", 16,
                            ifelse(Grouping_vector == "BV2T.global", 16,
                                   ifelse(Grouping_vector == "BV2T.LPS.global", 16, 
                                          ifelse(Grouping_vector == "BV2.LPS.neg.pulldown", 16, 
                                          ifelse(Grouping_vector == "BV2T.pos.pulldown", 16, 
                                                 ifelse(Grouping_vector == "BV2T.LPS.pos.pulldown", 16, NA)))))))


pt_colors <- ifelse(Grouping_vector == "BV2.global", "orange",
                    ifelse(Grouping_vector == "BV2.LPS.global", "purple",
                           ifelse(Grouping_vector == "BV2T.global", "limegreen",
                                  ifelse(Grouping_vector == "BV2T.LPS.global", "darkturquoise", 
                                         ifelse(Grouping_vector == "BV2.LPS.neg.pulldown", "black", 
                                            ifelse(Grouping_vector == "BV2T.pos.pulldown", "dodgerblue", 
                                                ifelse(Grouping_vector == "BV2T.LPS.pos.pulldown",  "magenta", NA)))))))

legend_groups <- c("BV2.global", "BV2.LPS.global", "BV2T.global", "BV2T.LPS.global", "BV2.LPS.neg.pulldown", "BV2T.pos.pulldown","BV2T.LPS.pos.pulldown")

plotMDS_1_4_allgroups_nonegs <- plotMDS(log2(normalized_counts_2), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for log2(mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("orange", "purple", "limegreen", "darkturquoise", "black", "dodgerblue", "magenta"), 
       pch = 16, title = "Transcriptome groups",cex=1.4)

dev.off()


################################################################################
# 1b: all groups w/o neg.pulldowns
pdf(file = "1b-1. Exp. 1-4_RNA-seq_all groups_nonegs_DESeq2_PCA_CR_02032024.pdf", width = 8, height = 10.5)
# perform DESeq2-basedPCA on rlog transformed data of 18 samples w/o neg pulldowns

rldr <- rlog(dds_all1_4_no_negs, blind = TRUE)
plotPCA(rldr, intgroup = "Group") # plots based on normalized counts

pcaData <- plotPCA(rldr, intgroup = "Group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca2 <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("DESeq2 PCA plot for log2(mRNA counts)") +
  theme(panel.background = element_rect(fill = "transparent", color = "black"), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5)) +  # Center the title
  labs(color = "Transcriptome groups") +
  coord_fixed()

dev.off()

# perform limma:plotMDS-based PCA of 18 samples w/o neg pulldowns
library(limma)

pdf(file = "1b-2. Exp. 1-4_RNA-seq_all groups_nonegs_Limma_PCA_CR_02032024.pdf", width = 8, height = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- read.csv("Exp. 1-4_Grouping vector.csv", header = TRUE)

pch_values <- ifelse(Grouping_vector == "BV2.global", 16,
                     ifelse(Grouping_vector == "BV2.LPS.global", 16,
                            ifelse(Grouping_vector == "BV2T.global", 16,
                                   ifelse(Grouping_vector == "BV2T.LPS.global", 16, 
                                          ifelse(Grouping_vector == "BV2T.pos.pulldown", 16, 
                                                 ifelse(Grouping_vector == "BV2T.LPS.pos.pulldown", 16, NA))))))


pt_colors <- ifelse(Grouping_vector == "BV2.global", "orange",
                    ifelse(Grouping_vector == "BV2.LPS.global", "purple",
                           ifelse(Grouping_vector == "BV2T.global", "limegreen",
                                  ifelse(Grouping_vector == "BV2T.LPS.global", "darkturquoise", 
                                         ifelse(Grouping_vector == "BV2T.pos.pulldown", "dodgerblue", 
                                                ifelse(Grouping_vector == "BV2T.LPS.pos.pulldown",  "magenta", NA))))))

legend_groups <- c("BV2.global", "BV2.LPS.global", "BV2T.global", "BV2T.LPS.global", "BV2T.pos.pulldown","BV2T.LPS.pos.pulldown")

plotMDS_1_4_allgroups_nonegs <- plotMDS(log2(normalized_counts_2), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for log2(mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("orange", "purple", "limegreen", "darkturquoise", "dodgerblue", "magenta"), 
       pch = 16, title = "Transcriptome groups",cex=1.4)

dev.off()

################################################################################
# 1c BV2T.global +/- LPS vs BV2T.pos.pulldown +/- LPS
# subset DESeq2 norm counts to only include BVT2 samples

#Updated clean_counts with neg.pulldowns removed
BV2T_clean_counts <- counts_clean_filtered_rows_1_4[,c(7:12, 15:20)]

#Updated traits with neg.pulldowns removed
BV2T_traits <- traits[c(7:12, 15:20),]

#create DESeq data set (dds) matrix for all 12 BV2T samples
BV2T_dds <- DESeqDataSetFromMatrix(countData = BV2T_clean_counts, #filtered, unnorm counts is the input here
                                             colData = BV2T_traits,
                                             design = ~Group)
dim(BV2T_dds)
# [1] 12742    12 

# set the factor level
BV2T_dds$Group <- relevel(BV2T_dds$Group, ref = "BV2T.global")

#Run DESeq
BV2T_dds<- DESeq(BV2T_dds)


# Create PCA plots
pdf(file = "1c-1. Exp. 1-4_RNA-seq_BV2T groups_DESeq2_PCA_CR_02032024.pdf", width = 8, height = 10.5)
# perform DESeq2-basedPCA on rlog transformed data of 18 samples w/o neg pulldowns

rldr <- rlog(BV2T_dds, blind = TRUE)
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
pdf(file = "1c-2. Exp. 1-4_RNA-seq_BV2T groups_Limma_PCA_CR_02032024.pdf", width = 8, height = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector_BV2T <- Grouping_vector[c(7:18),]
BV2T_norm_counts <- normalized_counts_2[,c(7:18)]

pch_values <- ifelse(Grouping_vector_BV2T == "BV2T.global", 16,
                     ifelse(Grouping_vector_BV2T == "BV2T.LPS.global", 16,
                            ifelse(Grouping_vector_BV2T == "BV2T.pos.pulldown", 16,
                                   ifelse(Grouping_vector_BV2T == "BV2T.LPS.pos.pulldown", 16, NA))))


pt_colors <- ifelse(Grouping_vector_BV2T == "BV2T.global", "limegreen",
                      ifelse(Grouping_vector_BV2T == "BV2T.LPS.global", "darkturquoise", 
                              ifelse(Grouping_vector_BV2T == "BV2T.pos.pulldown", "dodgerblue", 
                                      ifelse(Grouping_vector_BV2T == "BV2T.LPS.pos.pulldown",  "magenta", NA))))

legend_groups <- c("BV2T.global", "BV2T.LPS.global", "BV2T.pos.pulldown","BV2T.LPS.pos.pulldown")

plotMDS_BV2T <- plotMDS(log2(BV2T_norm_counts), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)
mtext(side=3, text="MDS Plot for log2(mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("limegreen", "darkturquoise", "dodgerblue", "magenta"), 
       pch = 16, title = "Transcriptome groups", cex = 1.4)

dev.off()

################################################################################
################################################################################
# Step 2. Correlations of mRNA counts
  # a. BV2T.global vs BV2T.LPS.global
  # b. BV2T.global vs BV2T.pos.pulldown
  # c. BV2T.LPS.global vs BV2T.LPS.pos.pulldown
  # d. BV2T.LPS.pos.pulldown vs BV2T.pos.pulldown

# check for finite values
any(is.finite(normalized_counts_2))
# [1] TRUE zeroes are considered finite

# CR: Because there are zeroes in the DESeq2 normalized counts data matrix, 
# the log2 transformation changes these values to infinite values (-ifn).
# A correlation coefficient cannot be calculated with infinite values, so
# I am going to add a "pseudocount" of + 1 to all the counts values to
# avoid taking a log of 0.

norm_counts_pseudo <- normalized_counts_2 + 1
any(is.finite(norm_counts_pseudo))

norm_counts_log <- log2(norm_counts_pseudo)
any(is.finite(norm_counts_log))

# load library to generate plot
library(ggplot2)

################################################################################
# 2a. BV2T.global vs BV2T.LPS.global

pdf(file = "2a-d. Exp. 1-4__RNA-seq_counts correlations_CR_02062024.pdf", width = 6, height = 6, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# subset into groups of interest and take the row mean
BV2T.global <- rowMeans(norm_counts_log[,c(7:9)])
BV2T.LPS.global <- rowMeans(norm_counts_log[,c(10:12)])

# generate a linear correlation
correlation_2a <- cor(BV2T.global, BV2T.LPS.global)

# convert data to a data frame
df_2a <- data.frame(BV2T.global = BV2T.global, BV2T.LPS.global = BV2T.LPS.global)

# create correlation plot
correlation_plot_2a <- ggplot(df_2a, aes(x = BV2T.global, y = BV2T.LPS.global)) +
  geom_point(color = "black", fill = "darkgreen", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "BV2T global log2(mRNA counts)",
       y = "BV2T LPS global log2(mRNA counts)") +
  ggtitle("Linear correlation for mRNA counts") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(BV2T.global), y = min(BV2T.LPS.global), 
           label = paste("Correlation coefficient:", round(correlation_2a, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2d$BV2T.pos.pulldown), y = min(df_2d$BV2T.LPS.pos.pulldown) - 2, 
           label = paste("Number of genes:", nrow(df_2a)), 
           hjust = 0, vjust = 1, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2a)

################################################################################
# 2b.  BV2T.global vs BV2T.pos.pulldown
# subset into groups of interest and take the row mean
BV2T.global <- rowMeans(norm_counts_log[,c(7:9)])
BV2T.pos.pulldown <- rowMeans(norm_counts_log[,c(13:15)])

# convert data to a data frame
df_2b <- data.frame(BV2T.global = BV2T.global, BV2T.pos.pulldown = BV2T.pos.pulldown)

# calculate correlation coefficient
correlation_2b <- cor(BV2T.global, BV2T.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others

correlation_plot_2b <- ggplot(df_2b, aes(x = BV2T.global, y = BV2T.pos.pulldown)) +
  geom_point(color = "black", fill = "darkgrey", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "BV2T global log2(mRNA counts)",
       y = "BV2T pulldown log2(mRNA counts)") +
  ggtitle("Linear correlation for mRNA counts") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(BV2T.global), y = min(BV2T.pos.pulldown), 
           label = paste("Correlation coefficient:", round(correlation_2b, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2d$BV2T.pos.pulldown), y = min(df_2d$BV2T.LPS.pos.pulldown) - 2, 
           label = paste("Number of genes:", nrow(df_2a)), # # of genes are the same for each df
           hjust = 0, vjust = 1, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2b)


################################################################################
# 2c. BV2T.LPS.global vs BV2T.LPS.pos.pulldown
# subset into groups of interest and take the row mean
BV2T.LPS.global <- rowMeans(norm_counts_log[,c(10:12)])
BV2T.LPS.pos.pulldown <- rowMeans(norm_counts_log[,c(16:18)])

# convert data to a data frame
df_2c <- data.frame(BV2T.LPS.global = BV2T.LPS.global, BV2T.LPS.pos.pulldown = BV2T.LPS.pos.pulldown)

# calculate correlation coefficient
correlation_2c <- cor(BV2T.LPS.global, BV2T.LPS.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others

correlation_plot_2c <- ggplot(df_2c, aes(x = BV2T.LPS.global, y = BV2T.LPS.pos.pulldown)) +
  geom_point(color = "black", fill = "lightgreen", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "BV2T LPS global log2(mRNA counts)",
       y = "BV2T LPS pulldown log2(mRNA counts)") +
  ggtitle("Linear correlation for mRNA counts") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(BV2T.LPS.global), y = min(BV2T.LPS.pos.pulldown), 
           label = paste("Correlation coefficient:", round(correlation_2b, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2d$BV2T.pos.pulldown), y = min(df_2d$BV2T.LPS.pos.pulldown) - 2, 
           label = paste("Number of genes:", nrow(df_2a)), # # of genes are the same for each df
           hjust = 1, vjust = 0, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2c)


################################################################################
# 2d. BV2T.LPS.pos.pulldown vs BV2T.pos.pulldown
# subset into groups of interest and take the row mean
BV2T.pos.pulldown <- rowMeans(norm_counts_log[,c(13:15)])
BV2T.LPS.pos.pulldown <- rowMeans(norm_counts_log[,c(16:18)])

# convert data to a data frame
df_2d <- data.frame(BV2T.pos.pulldown = BV2T.pos.pulldown, BV2T.LPS.pos.pulldown  = BV2T.LPS.pos.pulldown )

# calculate correlation coefficient
correlation_2d <- cor(BV2T.pos.pulldown, BV2T.LPS.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others

correlation_plot_2d <- ggplot(df_2d, aes(x = BV2T.pos.pulldown, y = BV2T.LPS.pos.pulldown)) +
  geom_point(color = "black", fill = "grey", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "BV2T pulldown log2(mRNA counts)",
       y = "BV2T LPS pulldown log2(mRNA counts)") +
  ggtitle("Linear correlation for mRNA counts") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(BV2T.pos.pulldown), y = min(BV2T.LPS.pos.pulldown), 
           label = paste("Correlation coefficient:", round(correlation_2d, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2d$BV2T.pos.pulldown), y = min(df_2d$BV2T.LPS.pos.pulldown) - 2, 
           label = paste("Number of genes:", nrow(df_2a)), # # of genes are the same for each df
           hjust = 1, vjust = 0, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2d)

dev.off()


################################################################################
################################################################################
# Step 3 Part I: Diffex volcanoes

# C1. BV2.LPS.global vs BV2.global
# C2. BV2T.global vs BV2.global
# C3. BV2T.LPS.global vs BV2.LPS.global
# C4. BV2T.LPS.global vs BV2T.global
# C5. BV2T.pos.pulldown vs BV2T.global
# C6. BV2T.LPS.pos.pulldown vs BV2T.pos.pulldown
# C7. BV2T.LPS.pos.pulldown vs BV2.LPS.neg.pulldown
# C8. BV2T.LPS.pos.pulldown vs BV2T.LPS.global

#load Libraries
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)
library(ggplot2)
library(ggrepel)

#use base pdf function to extract plots into a single pdf
pdf(file = "C1_RNA_BV2 LPS global vs BV2 global_rowmeans_12742_Diffex_02082024.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

#Updated clean_counts for comparison 1
final_clean_counts <- all1_4_no_negs_clean_counts
C1_clean_counts <- final_clean_counts[,c(1:6)]

#Updated traits for comparison 1
traits <- traits_1_4[c(1:12, 15:20),]
C1_traits <- traits[c(1:6),]

#Deseq2 data matrix is needed for downstream analysis
dds1 <- DESeqDataSetFromMatrix(countData = C1_clean_counts,
                               colData = C1_traits,
                               design = ~Group)   
dim(dds1)
# [1] 12742     6

dds1

# set the factor level
dds1$Group <- relevel(dds1$Group, ref = "BV2.global")

#Run DESeq 
dds1 <- DESeq(dds1)

#diff p adj value 
res05 <- results(dds1, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 12742 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2444, 19%
# LFC < 0 (down)     : 2414, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 2)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds1)
# [1] "Intercept"                        "Group_BV2.LPS.global_vs_BV2.global"

#plot MA
#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(dds1, main = "Comparison 1 RNA: BV2 LPS global vs BV2 global")

#use the log transform on the data set
rld <- rlog(dds1, blind=TRUE)
de1 <- as.data.frame(res05)
de1$padj[is.na(de1$padj)] <- 1

write.csv(de1, file = "C1_RNA_BV2 LPS global vs BV2 global_rowmeans_12742_Diffex_02082024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=FALSE, main = "Top 50 padj C1_RNA_BV2 LPS global vs BV2 global_rowmeans_12742_Diffex_02082024") 

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C1_RNA_BV2 LPS global vs BV2 global_rowmeans_12742_Diffex_Up_02082024.csv")
write.csv(down, file = "C1_RNA_BV2 LPS global vs BV2 global_rowmeans_12742_Diffex_Down_02082024.csv")

up_genes <- read.csv("C1_RNA_BV2 LPS global vs BV2 global_rowmeans_12742_Diffex_Up_02082024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C1_RNA_BV2 LPS global vs BV2 global_rowmeans_12742_Diffex_Down_02082024.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 1322   6

#volcano plot of diffex with padj df
de1 <- read.csv("C1_RNA_BV2 LPS global vs BV2 global_rowmeans_12742_Diffex_02082024.csv")

de1$diffexpressed <- "NS"
de1$diffexpressed[de1$log2FoldChange > 1.0 & de1$padj < 0.05] <- "UP"
de1$diffexpressed[de1$log2FoldChange < -1.0 & de1$padj < 0.05] <- "DOWN"

de1$delabel <- NA
de1$delabel[de1$diffexpressed != "NS"] <- de1$X[de1$diffexpressed != "NS"]

ggplot(data=de1, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C1_RNA_BV2 LPS global vs BV2 global_rowmeans_12742") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)", limits = c(-100, 500), breaks = seq(0, 500, by = 100)) +
  scale_x_continuous(name = "log2(difference) BV2 LPS global - BV2 global", limits = c(-10, 10), breaks = seq(-15, 15, by = 5))

# CR: there are genes that are super diffex on the volcano (e.g., Ccl5, Nfkbia, Gstm1, etc.) that have a p/padj value of zero.
# I am not sure how these are are showing up on the volcano since the -log10(0) = infinite value.
# Is it possible that the de df shows a zero value because the actually value is crazy small and close to zero?
  # Looking online, the consensus I have found is that DESeq2 does not report very small p-values as exactly 0 in the output data frame. 
  # Instead, it represents them as values close to 0 that are machine-representable.

#Perform PCA on rlog transformed data
rldr <- rlog(dds1, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

######################################################################################################################################
######################################################################################################################################
#Comparison 2 - C2 signifies the below mentioned comparison!
#2.BV2T global vs BV2 global

#use base pdf function to extract plots into a single pdf
pdf(file = "C2_RNA_BV2T global vs BV2 global_rowmeans_12742_Diffex_02082024.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

#Updated clean_counts for comparison 2
C2_clean_counts <- final_clean_counts[,c(1:3, 7:9)]

#Updated traits for comparison 2
C2_traits <- traits[c(1:3, 7:9),]

#Deseq2 data matrix
dds2 <- DESeqDataSetFromMatrix(countData = C2_clean_counts,
                               colData = C2_traits,
                               design = ~Group)   

dim(dds2)
# [1] 12742     6

dds2

# set the factor level
dds2$Group <- relevel(dds2$Group, ref = "BV2.global")

#Run DESeq
dds2 <- DESeq(dds2)

#diff p adj value 
res05 <- results(dds2, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 12735 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 180, 1.4%
# LFC < 0 (down)     : 161, 1.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 1976, 16%
# (mean count < 43)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds2)
# [1] "Intercept"                     "Group_BV2T.global_vs_BV2.global"

#plot MA
plotMA(res05, main = "Comparison 2 RNA: BV2T global vs BV2 global")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds2, blind=TRUE)
de2 <- as.data.frame(res05)
de2$padj[is.na(de2$padj)] <- 1

write.csv(de2, file = "C2_RNA_BV2T global vs BV2 global_rowmeans_12742_Diffex_02082024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C2_RNA_BV2T global vs BV2 global_rowmeans_12742_Diffex_02082024") 

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C2_RNA_BV2T global vs BV2 global_rowmeans_12742_Diffex_Up_02082024.csv")
write.csv(down, file = "C2_RNA_BV2T global vs BV2 global_rowmeans_12742_Diffex_Down_02082024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C2_RNA_BV2T global vs BV2 global_rowmeans_12742_Diffex_Up_02082024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C2_RNA_BV2T global vs BV2 global_rowmeans_12742_Diffex_Down_02082024.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 42  6

#volcano plot of diffex with padj df
de2 <- read.csv("C2_RNA_BV2T global vs BV2 global_rowmeans_12742_Diffex_02082024.csv")

de2$diffexpressed <- "NS"
de2$diffexpressed[de2$log2FoldChange > 1.0 & de2$padj < 0.05] <- "UP"
de2$diffexpressed[de2$log2FoldChange < -1.0 & de2$padj < 0.05] <- "DOWN"

de2$delabel <- NA
de2$delabel[de2$diffexpressed != "NS"] <- de2$X[de2$diffexpressed != "NS"]

ggplot(data=de2, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C2_RNA_BV2T global vs BV2 global_rowmeans") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
    scale_y_continuous(name = "-log10(padj)", limits = c(-10, 300), breaks = seq(0, 300, by = 100)) +
    scale_x_continuous(name = "log2(difference) BV2T global - BV2 global", limits = c(-10, 10), breaks = seq(-10, 10, by = 5))  
  
  
#Perform PCA on rlog transformed data 
rldr <- rlog(dds2, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

######################################################################################################################################
######################################################################################################################################
#Comparison 3 - C3 signifies the below mentioned comparison!
#3.BV2T LPS global vs BV2 LPS global

#use base pdf function to extract plots into a single pdf
pdf(file = "C3_RNA_BV2T LPS global vs BV2 LPS global_rowmeans_12742_Diffex_02082024.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)


#Updated clean_counts for comparison 3
C3_clean_counts <- final_clean_counts[,c(4:6, 10:12)]

#Updated traits for comparison 3
C3_traits <- traits[c(4:6, 10:12),]

#Deseq2 data matrix
dds3 <- DESeqDataSetFromMatrix(countData = C3_clean_counts,
                               colData = C3_traits,
                               design = ~Group)   

dim(dds3)
# [1] 12742     6

dds3

# set the factor level
dds3$Group <- relevel(dds3$Group, ref = "BV2.LPS.global")

#Run DESeq
dds3 <- DESeq(dds3)

#diff p adj value 
res05 <- results(dds3, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 12742 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1091, 8.6%
# LFC < 0 (down)     : 1072, 8.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 742, 5.8%
# (mean count < 16)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds3)
# [1] "Intercept"                             "Group_BV2T.LPS.global_vs_BV2.LPS.global"

#plot MA
plotMA(res05, main = "Comparison 3 RNA: BV2T LPS global vs BV2 LPS global")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds3, blind=TRUE)
de3 <- as.data.frame(res05)
de3$padj[is.na(de3$padj)] <- 1

write.csv(de3, file = "C3_RNA_BV2T LPS global vs BV2 LPS global_rowmeans_12742_Diffex_02082024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C3_RNA_BV2T LPS global vs BV2 LPS global_rowmeans_12742_Diffex_02082024")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C3_RNA_BV2T LPS global vs BV2 LPS global_rowmeans_12742_Diffex_Up_02082024.csv")
write.csv(down, file = "C3_RNA_BV2T LPS global vs BV2 LPS global_rowmeans_12742_Diffex_Down_02082024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C3_RNA_BV2T LPS global vs BV2 LPS global_rowmeans_12742_Diffex_Up_02082024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C3_RNA_BV2T LPS global vs BV2 LPS global_rowmeans_12742_Diffex_Down_02082024.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 135   6

#volcano plot of diffex with padj df
de3 <- read.csv("C3_RNA_BV2T LPS global vs BV2 LPS global_rowmeans_12742_Diffex_02082024.csv")

de3$diffexpressed <- "NS"
de3$diffexpressed[de3$log2FoldChange > 1.0 & de3$padj < 0.05] <- "UP"
de3$diffexpressed[de3$log2FoldChange < -1.0 & de3$padj < 0.05] <- "DOWN"

de3$delabel <- NA
de3$delabel[de3$diffexpressed != "NS"] <- de3$X[de3$diffexpressed != "NS"]

ggplot(data=de3, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C3_RNA_BV2T LPS global vs BV2 LPS global_rowmeans_12742") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)", limits = c(-10, 200), breaks = seq(0, 200, by = 100)) +
  scale_x_continuous(name = "log2(difference) BV2T LPS global - BV2 LPS global", limits = c(-10, 10), breaks = seq(-10, 10, by = 5))  


#Perform PCA on rlog transformed data 
rldr <- rlog(dds3, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

######################################################################################################################################
######################################################################################################################################
#Comparison 4 - C4 signifies the below mentioned comparison!
#4.BV2T LPS global vs BV2T global

#use base pdf function to extract plots into a single pdf
pdf(file = "C4_RNA_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_02082024.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)


#Updated clean_counts for comparison 4
C4_clean_counts <- final_clean_counts[,c(7:9, 10:12)]

#Updated traits for comparison 4
C4_traits <- traits[c(7:9, 10:12),]

#Deseq2 data matrix
dds4 <- DESeqDataSetFromMatrix(countData = C4_clean_counts,
                               colData = C4_traits,
                               design = ~Group)   

dim(dds4)
# [1] 12742     6

dds4

# set the factor level
dds4$Group <- relevel(dds4$Group, ref = "BV2T.global")

#Run DESeq
dds4 <- DESeq(dds4)

#diff p adj value 
res05 <- results(dds4, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 12742 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2339, 18%
# LFC < 0 (down)     : 2220, 17%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds4)
# [1] "Intercept"                          "Group_BV2T.LPS.global_vs_BV2T.global"

#plot MA
plotMA(res05, main = "Comparison 4 RNA: BV2T LPS global vs BV2T global")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds4, blind=TRUE)
de4 <- as.data.frame(res05)
de4$padj[is.na(de4$padj)] <- 1

write.csv(de4, file = "C4_RNA_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_02082024.csv")

var <- order(de4$padj, decreasing=F)
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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C4_RNA_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_02082024")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C4_RNA_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_Up_02082024.csv")
write.csv(down, file = "C4_RNA_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_Down_02082024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C4_RNA_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_Up_02082024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C4_RNA_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_Down_02082024.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 1242    6

#volcano plot of diffex with padj df
de4 <- read.csv("C4_RNA_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_02082024.csv")

de4$diffexpressed <- "NS"
de4$diffexpressed[de4$log2FoldChange > 1.0 & de4$padj < 0.05] <- "UP"
de4$diffexpressed[de4$log2FoldChange < -1.0 & de4$padj < 0.05] <- "DOWN"

de4$delabel <- NA
de4$delabel[de4$diffexpressed != "NS"] <- de4$X[de4$diffexpressed != "NS"]

ggplot(data=de4, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C4_RNA_BV2T LPS global vs BV2T global_rowmeans_127425") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)", limits = c(-10, 500), breaks = seq(0, 500, by = 100)) +
  scale_x_continuous(name = "log2(difference) BV2T LPS global - BV2T global", limits = c(-10, 10), breaks = seq(-10, 10, by = 5))  

#Perform PCA on rlog transformed data 
rldr <- rlog(dds4, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

#################################################################################
#################################################################################
#Comparison 5 - C5 signifies the below mentioned comparison!
#5. BVT2 pulldown vs BV2T global

#use base pdf function to extract plots into a single pdf
pdf(file = "C5_BVT2 pulldown vs BV2T global_rowmeans_12742_Diffex_02082024.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

#Updated clean_counts for comparison 5
C5_clean_counts <- final_clean_counts[,c(7:9, 13:15)]

#Updated traits for comparison 5
C5_traits <- traits[c(7:9, 13:15),]

#Deseq2 data matrix
dds5 <- DESeqDataSetFromMatrix(countData = C5_clean_counts,
                               colData = C5_traits,
                               design = ~Group)   
dim(dds5)
# [1] 12742     6

dds5

# set the factor level
dds5$Group <- relevel(dds5$Group, ref = "BV2T.global")

#Run DESeq
dds5 <- DESeq(dds5)

#diff p adj value 
res05 <- results(dds5, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 12730 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 590, 4.6%
# LFC < 0 (down)     : 730, 5.7%
# outliers [1]       : 0, 0%
# low counts [2]     : 3702, 29%
# (mean count < 108)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds5)
# [1] "Intercept"                             "Group_BV2T.pos.pulldown_vs_BV2T.global"

#plot MA
plotMA(res05, main = "Comparison 5 RNA: BV2T pulldown vs BV2T global")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds5, blind=TRUE)
de5<- as.data.frame(res05)
de5$padj[is.na(de5$padj)] <- 1

write.csv(de5, file = "C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_02082024.csv")

var <- order(de5$padj, decreasing=F)
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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_02082024.csv")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_Up_02082024.csv")
write.csv(down, file = "C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_Down_02082024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_Up_02082024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_Down_02082024.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 59  6

#volcano plot of diffex with padj df
de5 <- read.csv("C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_02082024.csv")

de5$diffexpressed <- "NS"
de5$diffexpressed[de5$log2FoldChange > 1.0 & de5$padj < 0.05] <- "UP"
de5$diffexpressed[de5$log2FoldChange < -1.0 & de5$padj < 0.05] <- "DOWN"

de5$delabel <- NA
de5$delabel[de5$diffexpressed != "NS"] <- de5$X[de5$diffexpressed != "NS"]

ggplot(data=de5, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)", limits = c(-10, 300), breaks = seq(0, 300, by = 100)) +
  scale_x_continuous(name = "log2(difference) BV2T pulldown - BV2T global", limits = c(-10, 10), breaks = seq(-10, 10, by = 5))  

#Perform PCA on rlog transformed data 
rldr <- rlog(dds5, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

#################################################################################
#################################################################################
#Comparison 6 - C6 signifies the below mentioned comparison!
#6. BVT2 LPS pulldown vs BV2T pulldown

#use base pdf function to extract plots into a single pdf
pdf(file = "C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_02082024.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

#Updated clean_counts for comparison 6
C6_clean_counts <- final_clean_counts[,c(13:18)]

#Updated traits for comparison 6
C6_traits <- traits[c(13:18),]

#Deseq2 data matrix
dds6 <- DESeqDataSetFromMatrix(countData = C6_clean_counts,
                               colData = C6_traits,
                               design = ~Group)   
dim(dds6)
# [1] 12742     6

dds6

# set the factor level
dds6$Group <- relevel(dds6$Group, ref = "BV2T.pos.pulldown")

#Run DESeq
dds6 <- DESeq(dds6)

#diff p adj value 
res05 <- results(dds6, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 12741 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2014, 16%
# LFC < 0 (down)     : 1941, 15%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds6)
# [1] "Intercept"                                        "Group_BV2T.LPS.pos.pulldown_vs_BV2T.pos.pulldown"

#plot MA
plotMA(res05, main = "Comparison 6 RNA: BV2T LPS pulldown vs BV2T pulldown")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds6, blind=TRUE)
de6 <- as.data.frame(res05)
de6$padj[is.na(de6$padj)] <- 1

write.csv(de6, file = "C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_02082024.csv")

var <- order(de6$padj, decreasing=F)
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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_02082024")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_Up_02082024.csv")
write.csv(down, file = "C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_Down_02082024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_Up_02082024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_Down_02082024.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 1088    6

#volcano plot of diffex with padj df
de6 <- read.csv("C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_02082024.csv")

de6$diffexpressed <- "NS"
de6$diffexpressed[de6$log2FoldChange > 1.0 & de6$padj < 0.05] <- "UP"
de6$diffexpressed[de6$log2FoldChange < -1.0 & de6$padj < 0.05] <- "DOWN"

de6$delabel <- NA
de6$delabel[de6$diffexpressed != "NS"] <- de6$X[de6$diffexpressed != "NS"]

ggplot(data=de6, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)", limits = c(-10, 500), breaks = seq(0, 500, by = 100)) +
  scale_x_continuous(name = "log2(difference) BV2T LPS pulldown - BV2T pulldown", limits = c(-10, 10), breaks = seq(-10, 10, by = 5))  

#Perform PCA on rlog transformed data 
rldr <- rlog(dds6, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

#################################################################################
#################################################################################
#Comparison 7 - C7 signifies the below mentioned comparison!
#7.BV2T LPS pulldown vs BV2 LPS pulldown

# CR: skipped for now since I don't have a row mean filtered norm counts matrix that includes the neg pulldowns

#use base pdf function to extract plots into a single pdf
pdf(file = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_02082024.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)


#Updated clean_counts for comparison 7
C7_clean_counts <- final_clean_counts[,c(13:14, 18:20)]

#Updated traits for comparison 7
C7_traits <- traits_1_4[c(13:14, 18:20),]

#Deseq2 data matrix
dds7 <- DESeqDataSetFromMatrix(countData = C7_clean_counts,
                               colData = C7_traits,
                               design = ~Group)   
dim(dds7)
# [1] 12742     5

dds7

# set the factor level
dds7$Group <- relevel(dds7$Group, ref = "BV2.LPS.neg.pulldown")

#Run DESeq
dds7 <- DESeq(dds7)

#diff p adj value 
res05 <- results(dds7, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 14364 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2275, 16%
# LFC < 0 (down)     : 2378, 17%
# outliers [1]       : 0, 0%
# low counts [2]     : 1381, 9.6%
# (mean count < 4)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds7)
# [1] "Intercept"                                           "Group_BV2T.LPS.pos.pulldown_vs_BV2.LPS.neg.pulldown"

#plot MA
plotMA(res05, main = "Comparison 7 RNA: BV2T LPS pulldown vs BV2 LPS pulldown")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds7, blind=TRUE)
de7 <- as.data.frame(res05)
de7$padj[is.na(de7$padj)] <- 1

write.csv(de7, file = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_02082024.csv")

var <- order(de7$padj, decreasing=F)
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

#select the 'contrast' you want
annotation_data <- as.data.frame(colData(rld)[c("Group")])
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_02082024")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_Up_02082024.csv")
write.csv(down, file = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_Down_02082024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_Up_02082024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_Down_02082024.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 1753    6

#volcano plot of diffex with padj df
de7 <- read.csv("C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_02082024.csv")

de7$diffexpressed <- "NS"
de7$diffexpressed[de7$log2FoldChange > 1.0 & de7$padj < 0.05] <- "UP"
de7$diffexpressed[de7$log2FoldChange < -1.0 & de7$padj < 0.05] <- "DOWN"

de7$delabel <- NA
de7$delabel[de7$diffexpressed != "NS"] <- de7$X[de7$diffexpressed != "NS"]

ggplot(data=de7, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)", limits = c(-10, 200), breaks = seq(0, 200, by = 50)) +
  scale_x_continuous(name = "log2(difference) BV2T LPS pulldown - BV2 LPS pulldown", limits = c(-10, 10), breaks = seq(-10, 10, by = 5))  

#Perform PCA on rlog transformed data 
rldr <- rlog(dds7, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()


## 03/07/2025
#Comparison 7 - C7 signifies the below mentioned comparison!
#7.BV2T LPS pulldown vs BV2 LPS pulldown

# CR: skipped for now since I don't have a row mean filtered norm counts matrix that includes the neg pulldowns

#use base pdf function to extract plots into a single pdf
pdf(file = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_03072025.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)


#Updated clean_counts for comparison 7
C7_clean_counts <- counts_clean_filtered_rows_1_4[,c(13:14, 18:20)]

#Updated traits for comparison 7
C7_traits <- traits_1_4[c(13:14, 18:20),]

#Deseq2 data matrix
dds7 <- DESeqDataSetFromMatrix(countData = C7_clean_counts,
                               colData = C7_traits,
                               design = ~Group)   
dim(dds7)
# [1] 12742     5

dds7

# set the factor level
dds7$Group <- relevel(dds7$Group, ref = "BV2.LPS.neg.pulldown")

#Run DESeq
dds7 <- DESeq(dds7)

#diff p adj value 
res05 <- results(dds7, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 12742 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2301, 18%
# LFC < 0 (down)     : 2408, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds7)
# [1] "Intercept"                                           "Group_BV2T.LPS.pos.pulldown_vs_BV2.LPS.neg.pulldown"

#plot MA
plotMA(res05, main = "Comparison 7 RNA: BV2T LPS pulldown vs BV2 LPS pulldown")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds7, blind=TRUE)
de7 <- as.data.frame(res05)
de7$padj[is.na(de7$padj)] <- 1

write.csv(de7, file = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_03072025.csv")

var <- order(de7$padj, decreasing=F)
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

#select the 'contrast' you want
annotation_data <- as.data.frame(colData(rld)[c("Group")])
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_03072025")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_Up_03072025.csv")
write.csv(down, file = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_Down_03072025.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_Up_03072025.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_Down_03072025.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 1772    6

#volcano plot of diffex with padj df
de7 <- read.csv("C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_03072025.csv")

de7$diffexpressed <- "NS"
de7$diffexpressed[de7$log2FoldChange > 1.0 & de7$padj < 0.05] <- "UP"
de7$diffexpressed[de7$log2FoldChange < -1.0 & de7$padj < 0.05] <- "DOWN"

de7$delabel <- NA
de7$delabel[de7$diffexpressed != "NS"] <- de7$X[de7$diffexpressed != "NS"]

ggplot(data=de7, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)", limits = c(-10, 200), breaks = seq(0, 200, by = 50)) +
  scale_x_continuous(name = "log2(difference) BV2T LPS pulldown - BV2 LPS pulldown", limits = c(-10, 10), breaks = seq(-10, 10, by = 5))  

#Perform PCA on rlog transformed data 
rldr <- rlog(dds7, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()



#################################################################################
#################################################################################
#Comparison 8 - C8 signifies the below mentioned comparison!
#8. BV2T LPS pulldown vs BV2T LPS global

#use base pdf function to extract plots into a single pdf
pdf(file = "C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_02082024.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

#Updated clean_counts for comparison 8
C8_clean_counts <- final_clean_counts[,c(10:12, 16:18)]

#Updated traits for comparison 8
C8_traits <- traits[c(10:12, 16:18),]

#Deseq2 data matrix
dds8 <- DESeqDataSetFromMatrix(countData = C8_clean_counts,
                               colData = C8_traits,
                               design = ~Group)   
dim(dds8)
# [1] 12742     6

dds8

# set the factor level
dds8$Group <- relevel(dds8$Group, ref = "BV2T.LPS.global")

#Run DESeq
dds8 <- DESeq(dds8)

#diff p adj value 
res05 <- results(dds8, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 12742 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1085, 8.5%
# LFC < 0 (down)     : 1196, 9.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 2965, 23%
# (mean count < 71)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds8)
# [1] "Intercept"                                     "Group_BV2T.LPS.pos.pulldown_vs_BV2T.LPS.global"

#plot MA
plotMA(res05, main = "Comparison 8 RNA: BV2T LPS pulldown vs BV2T LPS global")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds8, blind=TRUE)
de8<- as.data.frame(res05)
de8$padj[is.na(de8$padj)] <- 1

write.csv(de8, file = "C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_02082024.csv")

var <- order(de8$padj, decreasing=F)
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

#select the 'contrast' you want
annotation_data <- as.data.frame(colData(rld)[c("Group")])
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_02082024")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_Up_02082024.csv")
write.csv(down, file = "C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_Down_02082024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_Up_02082024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_Down_02082024.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 44  6

#volcano plot of diffex with padj df
de8 <- read.csv("C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_02082024.csv")

de8$diffexpressed <- "NS"
de8$diffexpressed[de8$log2FoldChange > 1.0 & de8$padj < 0.05] <- "UP"
de8$diffexpressed[de8$log2FoldChange < -1.0 & de8$padj < 0.05] <- "DOWN"

de8$delabel <- NA
de8$delabel[de8$diffexpressed != "NS"] <- de8$X[de8$diffexpressed != "NS"]

ggplot(data=de8, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)", limits = c(-10, 200), breaks = seq(0, 200, by = 100)) +
  scale_x_continuous(name = "log2(difference) BV2T LPS pulldown - BV2T LPS global", limits = c(-10, 10), breaks = seq(-10, 10, by = 5))  

#Perform PCA on rlog transformed data 
rldr <- rlog(dds8, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

################################################################################
#Comparison 9 - C9 signifies the below mentioned comparison!
#9. BV2 LPS global vs BV2 LPS pulldown

#use base pdf function to extract plots into a single pdf
pdf(file = "C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_09122025.pdf", height = 8, width = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

#Updated clean_counts for comparison 7
C9_clean_counts <- counts_clean_filtered_rows_1_4[,c(4:6, 13:14)]

#Updated traits for comparison 7
C9_traits <- traits_1_4[c(4:6, 13:14),]

library(DESeq2)

#Deseq2 data matrix
dds9 <- DESeqDataSetFromMatrix(countData = C9_clean_counts,
                               colData = C9_traits,
                               design = ~Group)   
dim(dds9)
# [1] 12742     5

dds9

# set the factor level
dds9$Group <- relevel(dds9$Group, ref = "BV2.LPS.neg.pulldown")

#Run DESeq
dds9 <- DESeq(dds9)

#diff p adj value 
res05 <- results(dds9, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 12742 with nonzero global read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2602, 20%
# LFC < 0 (down)     : 2516, 20%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

dim(res05)
# [1] 12742     6

#contrasts
resultsNames(dds9)
# [1] "Intercept"                                           "Group_BV2.LPS.global_vs_BV2.LPS.neg.pulldown"

#plot MA
plotMA(res05, main = "Comparison 9 RNA: BV2 LPS global vs BV2 LPS pulldown")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds9, blind=TRUE)
de9 <- as.data.frame(res05)
de9$padj[is.na(de9$padj)] <- 1

write.csv(de9, file = "C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_09122025.csv")

var <- order(de9$padj, decreasing=F)
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

#select the 'contrast' you want
annotation_data <- as.data.frame(colData(rld)[c("Group")])
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_09122025")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_Up_09122025.csv")
write.csv(down, file = "C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_Down_09122025.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_Up_09122025.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_Down_09122025.csv", header =  TRUE, row.names = 1)

#global df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 2217    6

#volcano plot of diffex with padj df
de9 <- read.csv("C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_09122025.csv")

de9$diffexpressed <- "NS"
de9$diffexpressed[de9$log2FoldChange > 1.0 & de9$padj < 0.05] <- "UP"
de9$diffexpressed[de9$log2FoldChange < -1.0 & de9$padj < 0.05] <- "DOWN"

de9$delabel <- NA
de9$delabel[de9$diffexpressed != "NS"] <- de9$X[de9$diffexpressed != "NS"]

ggplot(data=de9, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)", limits = c(-10, 200), breaks = seq(0, 200, by = 50)) +
  scale_x_continuous(name = "log2(difference) BV2 LPS global - BV2 LPS pulldown", limits = c(-10, 10), breaks = seq(-10, 10, by = 5))  

#Perform PCA on rlog transformed data 
rldr <- rlog(dds9, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()



################################################################################
################################################################################
# Step 3 Part II: Diffex GSEA

# C4. BV2T.LPS.global vs BV2T.global
# C5. BV2T.pos.pulldown vs BV2T.global
# C6. BV2T.LPS.pos.pulldown vs BV2T.pos.pulldown
# C8. BV2T.LPS.pos.pulldown vs BV2T.LPS.global

################################################################################
# Gene set erichment analysis with Eric Dammer's GOparallel pipeline
# https://github.com/edammer/GOparallel
  # performs hypergeometric overlap with Fisher Exact Test for enrichment p<0.05 and 5 minimum genes per ontology
  # currently outputs files that can be used as GO-Elite input, not read by this script...
  # outputs Z score barplots for each input list (module, or DiffEx Protein list) across 6 ontology types in the larger GMT files available from the Bader Lab website.


# Parameters for GOparallel GSEA on volcano-defined DEG lists  
source("GOparallel-FET.R")

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
cocluster=FALSE
#If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)

# "Clean up" to run before GOPar
rm(ANOVAout)

###############################################################################
# Comparison 4: BV2T LPS global vs BV2T global
# organize DESeq2 diffex output file to mirror LFQ-MS ANOVAout file
C4_ANOVA<-read.csv("C4_RNA_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_02082024.csv", header=TRUE, row.names=1)
C4_ANOVAout<-C4_ANOVA[,c(4,1,5,2)]
colnames(C4_ANOVAout)[c(3,4)]<-c("BV2T.LPS.global-BV2T.global","diff BV2T.LPS.global.BV2T.global")
write.csv(C4_ANOVAout, file = "C4_RNA_ANOVAout_BV2T LPS global vs BV2T global_rowmeans_12742_Diffex_02082024.csv")


# GOpar
flip=c(0)
outFilename <- "C4_RNA_BV2T.LPS.global_vs_BV2T.global_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C4_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

# Comparison 5: BV2T.pos.pulldown vs BV2T.global
#C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_02082024.csv
C5_ANOVA<-read.csv("C5_RNA_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_02082024.csv", header=TRUE, row.names=1)
C5_ANOVAout<-C5_ANOVA[,c(4,1,5,2)]
colnames(C5_ANOVAout)[c(3,4)]<-c("BV2T.pulldown-BV2T.global","diff BV2T.pulldown.BV2T.global")
write.csv(C5_ANOVAout, file = "C5_RNA_ANOVAout_BV2T pulldown vs BV2T global_rowmeans_12742_Diffex_02082024.csv" )


# GOpar
flip=c(0)
outFilename <- "C5_RNA_BV2T.pulldown_vs_BV2T.global_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C5_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

# Comparison 6: BV2T LPS pulldown bs BV2T pulldown
#C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_02082024.csv
C6_ANOVA<-read.csv("C6_RNA_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_02082024.csv", header=TRUE, row.names=1)
C6_ANOVAout<-C6_ANOVA[,c(4,1,5,2)]
colnames(C6_ANOVAout)[c(3,4)]<-c("BV2T.LPS.pulldown-BV2T.pulldown","diff BV2T.LPS.pulldown.BV2T.pulldown")
write.csv(C6_ANOVAout, file = "C6_RNA_ANOVAout_BVT2 LPS pulldown vs BV2T pulldown_rowmeans_12742_Diffex_02082024.csv")

# GOpar
flip=c(0)
outFilename <- "C6_RNA_BV2T.LPS.pulldown_vs_BV2T.pulldown_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C6_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

# Comparison 7: BV2T.LPS.pos.pulldown vs BV2.LPS.neg.pulldown
C7_ANOVA<-read.csv("C7_RNA_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_03072025.csv", header=TRUE, row.names=1)
C7_ANOVAout<-C7_ANOVA[,c(4,1,5,2)]
colnames(C7_ANOVAout)[c(3,4)]<-c("BV2T.LPS.pulldown-BV2.LPS.neg.pulldown","diff BV2T.LPS.pulldown.BV2.LPS.neg.pulldown")
write.csv(C7_ANOVAout, file = "C7_RNA_ANOVAout_BV2T LPS pulldown vs BV2 LPS pulldown_rowmeans_12742_Diffex_03072025.csv" )

# GOpar
flip=c(0)
outFilename <- "C7_RNA_BV2T.LPS.pulldown_vs_BV2.LPS.neg.pulldown_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C7_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.


# Comparison 8: BV2T.LPS.pos.pulldown vs BV2T.LPS.global
#C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_02082024.csv
C8_ANOVA<-read.csv("C8_RNA_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_02082024.csv", header=TRUE, row.names=1)
C8_ANOVAout<-C8_ANOVA[,c(4,1,5,2)]
colnames(C8_ANOVAout)[c(3,4)]<-c("BV2T.LPS.pulldown-BV2T.LPS.global","diff BV2T.LPS.pulldown.BV2T.LPS.global")
write.csv(C8_ANOVAout, file = "C8_RNA_ANOVAout_BV2T LPS pulldown vs BV2T LPS global_rowmeans_12742_Diffex_02082024.csv" )

# GOpar
flip=c(0)
outFilename <- "C8_RNA_BV2T.LPS.pos.pulldown_vs_BV2T.LPS.global_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C8_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

# Comparison 9: BV2.LPS.global vs BV2.LPS.neg.pulldowb
# C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_09122025.csv
C9_ANOVA<-read.csv("C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_09122025.csv", header=TRUE, row.names=1)
C9_ANOVAout<-C9_ANOVA[,c(4,1,5,2)]
colnames(C9_ANOVAout)[c(3,4)]<-c("BV2.LPS.global-BV2.LPS.neg.pulldown","diff BV2.LPS.global.BV2.LPS.neg.pulldown")
write.csv(C9_ANOVAout, file = "C9_RNA_BV2 LPS global vs BV2 LPS pulldown_rowmeans_12742_Diffex_02082024.csv" )

# GOpar
flip=c(0)
outFilename <- "C9_RNA_BV2.LPS.global_vs_BV2.LPS.neg.pulldown_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C9_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

#################################################################################

# end of data analysis and visualization


