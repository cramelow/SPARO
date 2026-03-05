########################################################################################################
# Code for performing Differential Gene Expression Analysis on Exp. 8-1 (Aldh1l1-RiboTag) and Exp. 9-2 (Aldh1l1-TurboID)
# row means >=10 filter before dds and splicing
# Christina Ramelow, MS
# Date: 06/18/2024
##################################################################################################################
##################################################################################################################
## Part 1: Data processing and clean up
# steps:
# 1. Filter out low counts
#   a. Venns: keep a gene if row mean count >= 10/sample group
#   b. PCAs, correlations, diffex: keep a gene if row mean count >= 10 across all samples
# 2. DESeq2 normalization
#   a. median of ratios approach (one file with all samples and 1 file w/o neg.pulldowns)
##################################################################################################################
#load Libraries
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)
library(ggplot2)
library(ggrepel)

stringsAsFactors = FALSE

rootdir <- setwd("/Users/christina/Desktop/Exp. 8-1 mRNA-seq analysis")

# load Exp. 9-2 feature_counts.txt file
raw_counts_Aldh_Turbo <- read.table("Exp. 1-4, 9-2 and 10-1_feature_counts.txt", sep = "\t", header = TRUE, row.names = 1)
dim(raw_counts_Aldh_Turbo)
# [1] 35976    55

# remove the first 5 columns that are primarily categorical
raw_counts_TurboID_2 <- raw_counts_Aldh_Turbo[,c(6:55)]
dim(raw_counts_TurboID_2)
# [1] 35976    50

# remove the smRNA-seq samples that were accidentally included
raw_counts_TurboID_3 <- raw_counts_TurboID_2[,c(1:20, 27:50)]
dim(raw_counts_TurboID_3)
# [1] 35976    44
# 21047FL-134-01-22 - 21047FL-134-01-27 should be removed. CR: checked

# check if there are NAs in the counts df
sum(is.na(raw_counts_TurboID_3))
#[1] 0

# create a geneID column and move it to the front of the df
raw_counts_TurboID_3$geneID <- rownames(raw_counts_TurboID_3)
raw_counts_TurboID_3 <- raw_counts_TurboID_3[, c("geneID", setdiff(names(raw_counts_TurboID_3), "geneID"))]

# load Exp. 8-1 feature_counts.txt file
raw_counts_Aldh_RiboTag <- read.table("Exp. 8-1 feature_counts.txt", sep = "\t", header = TRUE, row.names = 1)
dim(raw_counts_Aldh_RiboTag)
# [1] 35976    15

# remove the first 5 columns that are primarily categorical
raw_counts_RiboTag_2 <- raw_counts_Aldh_RiboTag[,c(6:15)]
dim(raw_counts_RiboTag_2)
# [1] 35976    10

# check if there are NAs in the counts df
sum(is.na(raw_counts_RiboTag_2))
#[1] 0

# create a geneID column and move it to the front of the df
raw_counts_RiboTag_2$geneID <- rownames(raw_counts_RiboTag_2)
raw_counts_RiboTag_2 <- raw_counts_RiboTag_2[, c("geneID", setdiff(names(raw_counts_RiboTag_2), "geneID"))]

# combine the two featurecounts.txt files by shared geneIDs
TurboID_RiboTag_merged_raw_counts <- merge(raw_counts_TurboID_3, raw_counts_RiboTag_2, by = "geneID")
dim(TurboID_RiboTag_merged_raw_counts)
# [1] 35976    55

# set the geneID column back as row names
rownames(TurboID_RiboTag_merged_raw_counts) <- TurboID_RiboTag_merged_raw_counts$geneID

# remove the geneID column
TurboID_RiboTag_merged_raw_counts <- TurboID_RiboTag_merged_raw_counts[, -which(names(TurboID_RiboTag_merged_raw_counts) == "geneID")]

dim(TurboID_RiboTag_merged_raw_counts)
# [1] 35976    54

# load Exp. 1-4, 9-2, 10-1 and 8-1 RNA traits file
traits <- read.csv("Exp. 1-4, 9-2, 10-1 and 8-1 combined traits_RNA.csv",  header = TRUE, row.names = 1)
dim(traits)
# [1] 54  2

# rename the column name of the merged_raw_counts df to match traits df
# Extract row names from traits
new_col_names <- rownames(traits)

# Assign row names of traits as column names
names(TurboID_RiboTag_merged_raw_counts) <- new_col_names

# check that the row names of the traits df and the column names of the counts df are the same
raw_counts_clean <- TurboID_RiboTag_merged_raw_counts[,match(rownames(traits), colnames(TurboID_RiboTag_merged_raw_counts))]

#sanity check 1 colnames in counts file = rownames in traits file
all(colnames(raw_counts_clean) %in% rownames(traits))
#TRUE

#sanity check 2 colnames in counts file = rownames in traits file
all(colnames(raw_counts_clean) == rownames(traits))
#TRUE

# Exp. 9-2 and 8-1 splced
raw_counts_clean_9_2_8_1 <- raw_counts_clean[,c(21:41, 45:54)]
dim(raw_counts_clean_9_2_8_1)
# [1] 35976    31

colnames(raw_counts_clean_9_2_8_1)
# should be "21047FL-134-01-28" - "21047FL-134-01-48" and  "23238R-04-01" - "23238R-04-10" CR: checked

# write .csv file for raw_counts_clean for Exp. 9-2 and 8-1
write.csv(raw_counts_clean_9_2_8_1, file=paste0("0.Exp. 9-2 and 8-1_RNA-seq_raw_counts_clean_unfiltered-",dim(raw_counts_clean_9_2_8_1)[1],"x",dim(raw_counts_clean_9_2_8_1)[2],"_CR_06182024.csv"))

################################################################################
################################################################################
## Step 1a: Filter out low counts using a row mean >/= 10 by group for Venn diagrams
# Exp. 9-2 and 8-1
# 1. subset the traits and counts dfs by sample group to create list of genes for Venn diagrams
# 2. calculate the row mean
# 3. filter for genes with a row mean >= 10
# 4. write a .csv file for each sample group
################################################################################
# subset the traits file to only include Exp. 9-2 and 8-1
traits_9_2_8_1 <- traits[c(21:41, 45:54),]
dim(traits_9_2_8_1)
# [1] 31  2
write.csv(traits_9_2_8_1, file = "Exp. 9-2 and 8-1 RNA traits.csv")
################################################################################
# row mean >= 10 Aldh1l1.bulk.cortex
counts_Aldh1l1.bulk.cortex <- raw_counts_clean_9_2_8_1[,c(4:6)]
dim(counts_Aldh1l1.bulk.cortex)
# [1] 35976     3

traits_Aldh1l1.bulk.cortex  <- traits_9_2_8_1[c(4:6),]
traits_Aldh1l1.bulk.cortex
#                 Sample.ID..         Group
# 21047FL-134-01-31   4-RNA 30L Aldh1l1.bulk.cortex
# 21047FL-134-01-32   5-RNA 31B Aldh1l1.bulk.cortex
# 21047FL-134-01-33   6-RNA 24L Aldh1l1.bulk.cortex

rownames(traits_Aldh1l1.bulk.cortex)
# [1] "21047FL-134-01-31" "21047FL-134-01-32" "21047FL-134-01-33"

colnames(counts_Aldh1l1.bulk.cortex)
# [1] "21047FL-134-01-31" "21047FL-134-01-32" "21047FL-134-01-33"

# calculate row mean
Aldh1l1.bulk.cortex_row_mean <- rowMeans(counts_Aldh1l1.bulk.cortex)

# subset rows where row means >= 10
Aldh1l1.bulk.cortex_filtered_rows <- counts_Aldh1l1.bulk.cortex[Aldh1l1.bulk.cortex_row_mean >= 10, ]

dim(Aldh1l1.bulk.cortex_filtered_rows)
# [1] 16737     3

# write .csv file 
write.csv(Aldh1l1.bulk.cortex_filtered_rows, file=paste0("1a.Exp. 9-2 and 8-1_RNA-seq_Aldh1l1.bulk.cortex_Venn filter-",dim(Aldh1l1.bulk.cortex_filtered_rows)[1],"x",dim(Aldh1l1.bulk.cortex_filtered_rows)[2],"_CR_06182024.csv"))

################################################################################
# row mean >= 10 RiboTag bulk.cortex
counts_RiboTag.bulk.cortex <- raw_counts_clean_9_2_8_1[,c(22:23)]
dim(counts_RiboTag.bulk.cortex)
# [1] 35976     2

traits_RiboTag.bulk.cortex <- traits_9_2_8_1[c(22:23),]
traits_RiboTag.bulk.cortex
#          Sample.ID..         Group
# 23238R-04-01         CR1 RiboTag.bulk.cortex
# 23238R-04-02         CR2 RiboTag.bulk.cortex

rownames(traits_RiboTag.bulk.cortex)
# [1] "23238R-04-01" "23238R-04-02"

colnames(counts_RiboTag.bulk.cortex)
# [1] "23238R-04-01" "23238R-04-02"

# calculate row mean
RiboTag.bulk.cortex_row_mean <- rowMeans(counts_RiboTag.bulk.cortex)

# subset rows where row means >= 10
RiboTag.bulk.cortex_filtered_rows <- counts_RiboTag.bulk.cortex[RiboTag.bulk.cortex_row_mean >= 10, ]

dim(RiboTag.bulk.cortex_filtered_rows)
# [1] 17584     2

# write .csv file 
write.csv(RiboTag.bulk.cortex_filtered_rows, file=paste0("1a.Exp. 9-2 and 8-1_RNA-seq_RiboTag.bulk.cortex_Venn filter-",dim(RiboTag.bulk.cortex_filtered_rows)[1],"x",dim(RiboTag.bulk.cortex_filtered_rows)[2],"_CR_06182024.csv"))

###############################################################################
# row mean >= 10 Aldh1l1.RiboTag.bulk.cortex
counts_Aldh1l1.RiboTag.bulk.cortex <- raw_counts_clean_9_2_8_1[,c(24:26)]
dim(counts_Aldh1l1.RiboTag.bulk.cortex)
# [1] 35976     3
traits_Aldh1l1.RiboTag.bulk.cortex <- traits_9_2_8_1[c(24:26),]
traits_Aldh1l1.RiboTag.bulk.cortex
#              Sample.ID..                 Group
# 23238R-04-03         CR3 Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-04         CR4 Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-05         CR5 Aldh1l1.RiboTag.bulk.cortex

rownames(traits_Aldh1l1.RiboTag.bulk.cortex)
# [1] "23238R-04-03" "23238R-04-04" "23238R-04-05"

colnames(counts_Aldh1l1.RiboTag.bulk.cortex)
# [1] "23238R-04-03" "23238R-04-04" "23238R-04-05"

# calculate row mean
Aldh1l1.RiboTag.bulk.cortex_row_mean <- rowMeans(counts_Aldh1l1.RiboTag.bulk.cortex)

# subset rows where row means >= 10
Aldh1l1.RiboTag.bulk.cortex_filtered_rows <- counts_Aldh1l1.RiboTag.bulk.cortex[Aldh1l1.RiboTag.bulk.cortex_row_mean >= 10, ]

dim(Aldh1l1.RiboTag.bulk.cortex_filtered_rows)
# [1] 17583     3

# write .csv file 
write.csv(Aldh1l1.RiboTag.bulk.cortex_filtered_rows, file=paste0("1a.Exp. 9-2 and 8-1_RNA-seq_Aldh1l1.RiboTag.bulk.cortex_Venn filter-",dim(Aldh1l1.RiboTag.bulk.cortex_filtered_rows)[1],"x",dim(Aldh1l1.RiboTag.bulk.cortex_filtered_rows)[2],"_CR_06182024.csv"))

################################################################################
# row mean >= 10 Aldh1l1.pos.pulldown
counts_Aldh1l1.pos.pulldown <- raw_counts_clean_9_2_8_1[,c(14:16)]
dim(counts_Aldh1l1.pos.pulldown)
# [1] 35976     3
traits_Aldh1l1.pos.pulldown <- traits_9_2_8_1[c(14:16),]
traits_Aldh1l1.pos.pulldown
#                   Sample.ID..                Group
# 21047FL-134-01-41  15-RNA 30L Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L Aldh1l1.pos.pulldown

rownames(traits_Aldh1l1.pos.pulldown)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43"

colnames(counts_Aldh1l1.pos.pulldown)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43"

# calculate row mean
Aldh1l1.pos.pulldown_row_mean <- rowMeans(counts_Aldh1l1.pos.pulldown)

# subset rows where row means >= 10
Aldh1l1.pos.pulldown_filtered_rows <- counts_Aldh1l1.pos.pulldown[Aldh1l1.pos.pulldown_row_mean >= 10, ]

dim(Aldh1l1.pos.pulldown_filtered_rows)
# [1] 16632     3

# write .csv file 
write.csv(Aldh1l1.pos.pulldown_filtered_rows, file=paste0("1a.Exp. 9-2 and 8-1_RNA-seq_Aldh1l1.pos.pulldown_Venn filter-",dim(Aldh1l1.pos.pulldown_filtered_rows)[1],"x",dim(Aldh1l1.pos.pulldown_filtered_rows)[2],"_CR_06182024.csv"))

################################################################################
# row mean >= 10 Aldh1l1.RiboTag.pos.pulldown
counts_Aldh1l1.RiboTag.pos.pulldown<- raw_counts_clean_9_2_8_1[,c(29:31)]
dim(counts_Aldh1l1.RiboTag.pos.pulldown)
# [1] 35976     3

traits_Aldh1l1.RiboTag.pos.pulldown<- traits_9_2_8_1[c(29:31),]
traits_Aldh1l1.RiboTag.pos.pulldown
#              Sample.ID..                        Group
# 23238R-04-08         CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09         CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10        CR10 Aldh1l1.RiboTag.pos.pulldown

rownames(traits_Aldh1l1.RiboTag.pos.pulldown)
# [1] "23238R-04-08" "23238R-04-09" "23238R-04-10"

colnames(counts_Aldh1l1.RiboTag.pos.pulldown)
# [1] "23238R-04-08" "23238R-04-09" "23238R-04-10"

# calculate row mean
Aldh1l1.RiboTag.pos.pulldown_row_mean <- rowMeans(counts_Aldh1l1.RiboTag.pos.pulldown)

# subset rows where row means >= 10
Aldh1l1.RiboTag.pos.pulldown_filtered_rows <- counts_Aldh1l1.RiboTag.pos.pulldown[Aldh1l1.RiboTag.pos.pulldown_row_mean >= 10, ]

dim(Aldh1l1.RiboTag.pos.pulldown_filtered_rows)
# [1] 16737     3

# write .csv file 
write.csv(Aldh1l1.RiboTag.pos.pulldown_filtered_rows, file=paste0("1a.Exp. 9-2 and 8-1_RNA-seq_Aldh1l1.RiboTag.pos.pulldown_Venn filter-",dim(Aldh1l1.RiboTag.pos.pulldown_filtered_rows)[1],"x",dim(Aldh1l1.RiboTag.pos.pulldown_filtered_rows)[2],"_CR_06182024.csv"))


################################################################################
################################################################################
# clean up filtered row dfs to only include gene lists

# Aldh1l1.bulk.cortex_filtered_rows
Aldh1l1.bulk.cortex_genes <- as.matrix(rownames(Aldh1l1.bulk.cortex_filtered_rows))
dim(Aldh1l1.bulk.cortex_genes)
# [1] 16737     1

# Aldh1l1.RiboTag.bulk.cortex_filtered_rows
Aldh1l1.RiboTag.bulk.cortex_genes <- as.matrix(rownames(Aldh1l1.RiboTag.bulk.cortex_filtered_rows))
dim(Aldh1l1.RiboTag.bulk.cortex_genes)
# [1] 17583     1

# RiboTag.bulk.cortex_filtered_rows
RiboTag.bulk.cortex_genes <- as.matrix(rownames(RiboTag.bulk.cortex_filtered_rows))
dim(RiboTag.bulk.cortex_genes)
# [1] 17584     1

# Aldh1l1.pos.pulldown_filtered_rows
Aldh1l1.pos.pulldown_genes <- as.matrix(rownames(Aldh1l1.pos.pulldown_filtered_rows))
dim(Aldh1l1.pos.pulldown_genes)
# [1] 16632     1

# Aldh1l1.RiboTag.pos.pulldown_filtered_rows
Aldh1l1.RiboTag.pos.pulldown_genes <- as.matrix(rownames(Aldh1l1.RiboTag.pos.pulldown_filtered_rows))
dim(Aldh1l1.RiboTag.pos.pulldown_genes)
# [1] 16737     1

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
pdf(file = "1a. Exp. 9-2 and 8-1_RNA-seq_Venns_CR_06182024.pdf", height = 11, width = 8.5, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# Aldh1l1.bulk.cortex vs Aldh1l1.RiboTag.bulk.cortex
biovenn1 <- draw.venn(Aldh1l1.bulk.cortex_genes, Aldh1l1.RiboTag.bulk.cortex_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Aldh1l1.TurboID vs Aldh1l1.RiboTag bulk.cortex RNA", 
                      xtitle = "Aldh1l1.bulk.cortex", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "Aldh1l1.RiboTag.bulk.cortex",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "purple",
                      y_c = "red",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE)
# [1] "x bulk.cortex: 16737"
# [1] "y bulk.cortex: 17583"
# [1] "z bulk.cortex: 0"
# [1] "x only: 245"
# [1] "y only: 1091"
# [1] "z only: 0"
# [1] "x-y bulk.cortex overlap: 16492"

# Aldh1l1.RiboTag.bulk.cortex vs Aldh1l1.RiboTag.pos.pulldown
biovenn2 <- draw.venn(Aldh1l1.RiboTag.bulk.cortex_genes, Aldh1l1.RiboTag.pos.pulldown_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Aldh1l1.RiboTag.bulk.cortex vs Aldh1l1.RiboTag.pos.pulldown RNA", 
                      xtitle = "Aldh1l1.RiboTag.bulk.cortex", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "Aldh1l1.RiboTag.pos.pulldown",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "grey",
                      y_c = "red",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE)
# [1] "x bulk.cortex: 17583"
# [1] "y bulk.cortex: 16737"
# [1] "z bulk.cortex: 0"
# [1] "x only: 1319"
# [1] "y only: 473"
# [1] "z only: 0"
# [1] "x-y bulk.cortex overlap: 16264"


# Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown
biovenn3 <- draw.venn(Aldh1l1.pos.pulldown_genes, Aldh1l1.RiboTag.pos.pulldown_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Aldh1l1.TurboID vs Aldh1l1.RiboTag pulldown RNA", 
                      xtitle = "Aldh1l1.TurboID", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "Aldh1l1.RiboTag",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "purple4",
                      y_c = "red4",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE)

# [1] "x bulk.cortex: 16632"
# [1] "y bulk.cortex: 16737"
# [1] "z bulk.cortex: 0"
# [1] "x only: 666"
# [1] "y only: 771"
# [1] "z only: 0"
# [1] "x-y bulk.cortex overlap: 15966"
# [1] "x-z bulk.cortex overlap: 0"
# [1] "y-z bulk.cortex overlap: 0"
# [1] "x-y only overlap: 15966"


dev.off()

################################################################################
################################################################################
## Step 1b: Create row mean filter >=10 across all samples for downstream diffex analysis

# calculate row mean
counts_clean_row_mean <- rowMeans(raw_counts_clean_9_2_8_1)

# subset rows where row means >= 10
counts_clean_filtered_rows_9_2_8_1 <- raw_counts_clean_9_2_8_1[counts_clean_row_mean  >= 10, ]

dim(counts_clean_filtered_rows_9_2_8_1)
# [1] 17054    31

# write .csv file 
write.csv(counts_clean_filtered_rows_9_2_8_1, file=paste0("1b.Exp. 9-2 and 8-1_RNA-seq_clean_counts_diffex filter-",dim(counts_clean_filtered_rows_9_2_8_1)[1],"x",dim(counts_clean_filtered_rows_9_2_8_1)[2],"_CR_06182024.csv"))

################################################################################
################################################################################
## Step 2: DESeq2 normalization of counts based on median of ratios

# load library
library(DESeq2)

# create DESeq data set matrix (dds) for all 31 samples
dds_all9_2_8_1 <- DESeqDataSetFromMatrix(countData = counts_clean_filtered_rows_9_2_8_1,
                                     colData = traits_9_2_8_1,
                                     design = ~Group)
dim(dds_all9_2_8_1)
# [1] 17054    31

# set the factor level
dds_all9_2_8_1$Group <- relevel(dds_all9_2_8_1$Group, ref = "Aldh1l1.bulk.cortex")

# run DESeq
dds_all9_2_8_1<- DESeq(dds_all9_2_8_1)

View(counts(dds_all9_2_8_1))

dds <- dds_all9_2_8_1

dds <- estimateSizeFactors(dds)

sizeFactors(dds)
# 21047FL-134-01-28 21047FL-134-01-29 21047FL-134-01-30 21047FL-134-01-31 21047FL-134-01-32 21047FL-134-01-33 
# 0.8893720         0.8049612         0.7749720         0.8572228         1.0256973         1.1482553 
# 21047FL-134-01-34 21047FL-134-01-35 21047FL-134-01-36 21047FL-134-01-37 21047FL-134-01-38 21047FL-134-01-39 
# 0.9989146         0.6364181         0.7855086         0.9811747         0.9203750         1.2266314 
# 21047FL-134-01-40 21047FL-134-01-41 21047FL-134-01-42 21047FL-134-01-43 21047FL-134-01-44 21047FL-134-01-45 
# 1.5524386         1.0284815         1.0962545         0.8934085         1.0540818         0.7051570 
# 21047FL-134-01-46 21047FL-134-01-47 21047FL-134-01-48      23238R-04-01      23238R-04-02      23238R-04-03 
# 0.8453669         0.8415983         0.8265650         1.1347803         1.2736475         1.3087446 
# 23238R-04-04      23238R-04-05      23238R-04-06      23238R-04-07      23238R-04-08      23238R-04-09 
# 1.2948424         1.0700477         1.2484744         1.2103432         0.8805125         1.1037337 
# 23238R-04-10 
# 1.2092379 

normalized_counts_1 <- counts(dds, normalized=TRUE)

write.csv(normalized_counts_1, file=paste0("2a.Exp. 9-2 and 8-1_RNA-seq_clean_filtered_DESeq2 norm counts_diffex filter-",dim(normalized_counts_1)[1],"x",dim(normalized_counts_1)[2],"_CR_06182024.csv"))

################################################################################
# updated clean_counts with neg.pulldowns removed
all9_2_8_1_no_negs_clean_counts <- counts_clean_filtered_rows_9_2_8_1[,c(1:11, 14:26, 29:31)]

# udated traits with neg.pulldowns removed
all9_2_8_1_no_negs_traits <- traits_9_2_8_1[c(1:11, 14:26, 29:31),]
all9_2_8_1_no_negs_traits

# create DESeq data set matrix dds for all 22 samples
dds_all9_2_8_1_no_negs <- DESeqDataSetFromMatrix(countData = all9_2_8_1_no_negs_clean_counts,
                                             colData = all9_2_8_1_no_negs_traits,
                                             design = ~Group)
dim(dds_all9_2_8_1_no_negs)
# [1] 17054    27

# set the factor level
dds_all9_2_8_1_no_negs$Group <- relevel(dds_all9_2_8_1_no_negs$Group, ref = "Aldh1l1.bulk.cortex")

# run DESeq
dds_all9_2_8_1_no_negs<- DESeq(dds_all9_2_8_1_no_negs)

View(counts(dds_all9_2_8_1_no_negs))

dds <- dds_all9_2_8_1_no_negs

dds <- estimateSizeFactors(dds)

sizeFactors(dds)
# 21047FL-134-01-28 21047FL-134-01-29 21047FL-134-01-30 21047FL-134-01-31 21047FL-134-01-32 21047FL-134-01-33 
# 0.9213783         0.8349613         0.8023031         0.8899217         1.0636470         1.1919804 
# 21047FL-134-01-34 21047FL-134-01-35 21047FL-134-01-36 21047FL-134-01-37 21047FL-134-01-38 21047FL-134-01-41 
# 1.0358285         0.6600755         0.8153691         1.0210895         0.9576344         1.0654161 
# 21047FL-134-01-42 21047FL-134-01-43 21047FL-134-01-44 21047FL-134-01-45 21047FL-134-01-46 21047FL-134-01-47 
# 1.1313599         0.9255927         1.0888734         0.7342807         0.8767168         0.8755041 
# 21047FL-134-01-48      23238R-04-01      23238R-04-02      23238R-04-03      23238R-04-04      23238R-04-05 
# 0.8588437         1.1870658         1.3348204         1.3673661         1.3577919         1.1199568 
# 23238R-04-08      23238R-04-09      23238R-04-10 
# 0.9121930         1.1422150         1.2556182 

normalized_counts_2 <- counts(dds, normalized=TRUE)

write.csv(normalized_counts_2, file=paste0("2a.Exp. 9-2 and 8-1_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-",dim(normalized_counts_2)[1],"x",dim(normalized_counts_2)[2],"_CR_06182024.csv"))

## end of data processing and clean up ##

################################################################################
################################################################################
## Part II: Data analysis and visualization
# Exp. 9-2 and 8-1 mRNA-seq analysis

# steps:
# 1. PCAs
  # 1b. all groups w/o neg.pulldowns
  # 1c. Aldh1l1.pos.pulldown (n=3/group) vs Aldh1l1.RiboTag.pos.pulldown (n=3/group)
# 2. Correlations
  # a. Aldh1l1.bulk.cortex vs Aldh1l1.RiboTag.bulk.cortex  (n=2/group)
  # b. Aldh1l1.RiboTag.bulk.cortex (n=2/group) vs Aldh1l1.RiboTag.pos.pulldown (n=3/group)
  # c. Aldh1l1.pos.pulldown (n=3/group) vs Aldh1l1.RiboTag.pos.pulldown (n=3/group)
# 3. Diffex volcanoes & GSEA
  # a. Aldh1l1.bulk.cortex vs Aldh1l1.RiboTag.bulk.cortex  (n=2/group)
  # b. Aldh1l1.RiboTag.bulk.cortex (n=2/group) vs Aldh1l1.RiboTag.pos.pulldown (n=3/group)
  # c. Aldh1l1.pos.pulldown (n=3/group) vs Aldh1l1.RiboTag.pos.pulldown (n=3/group)


################################################################################
## Step 1: Generate PCA plots
library(DESeq2)
library(ggplot2)

# 1a: all groups
# create DESeq data set matrix dds for all 20 samples
dds_all9_2_8_1 <- DESeqDataSetFromMatrix(countData = counts_clean_filtered_rows_9_2_8_1,
                                     colData = traits_9_2_8_1,
                                     design = ~Group)
dim(dds_all9_2_8_1)
# [1] 17054    31

# set the factor level
dds_all9_2_8_1$Group <- relevel(dds_all9_2_8_1$Group, ref = "Aldh1l1.bulk.cortex")

# run DESeq
dds_all9_2_8_1<- DESeq(dds_all9_2_8_1)

# create pdf output file
pdf(file = "1a-1. Exp. 9-2 and 8-1_RNA-seq_all groups_DESeq2_PCA_CR_06182024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

# perform DESeq2-basedPCA on rlog transformed data of all 20 samples
rldr <- rlog(dds_all9_2_8_1, blind = TRUE)
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

print(pca1)

dev.off()

################################################################################
# load library
library(limma)

# 1b: all groups w/o neg.pulldowns
pdf(file = "1b-1. Exp. 9-2 and 8-1_RNA-seq_all groups_nonegs_DESeq2_PCA_CR_06182024.pdf", width = 8, height = 8.5)
# perform DESeq2-basedPCA on rlog transformed data of 18 samples w/o neg pulldowns

rldr <- rlog(dds_all9_2_8_1_no_negs, blind = TRUE)
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

print(pca2)

dev.off()

###############################################################################
# CR didn't edit 06/20/2024
# perform limma:plotMDS-based PCA w/o neg pulldowns
library(limma)

pdf(file = "1b-2. Exp. 9-2 and 8-1_RNA-seq_all groups_nonegs_Limma_PCA_CR_06182024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- traits_9_2_8_1$Group[c(1:11, 14:24)]

pch_values <- ifelse(Grouping_vector == "control.bulk.cortex", 16,
                     ifelse(Grouping_vector == "control.LPS.bulk.cortex", 16,
                            ifelse(Grouping_vector == "Aldh1l1.bulk.cortex", 16,
                                   ifelse(Grouping_vector == "Aldh1l1.RiboTag.bulk.cortex", 16, 
                                                 ifelse(Grouping_vector == "Aldh1l1.pos.pulldown", 16, 
                                                        ifelse(Grouping_vector == "Aldh1l1.LPS.pos.pulldown", 16, 
                                                               ifelse(Grouping_vector == "RiboTag.bulk.cortex", 16, NA)))))))


pt_colors <- ifelse(Grouping_vector == "control.bulk.cortex", "magenta",
                    ifelse(Grouping_vector == "control.LPS.bulk.cortex", "pink4",
                           ifelse(Grouping_vector == "Aldh1l1.bulk.cortex", "coral",
                                  ifelse(Grouping_vector == "Aldh1l1.RiboTag.bulk.cortex", "green3", 
                                                ifelse(Grouping_vector == "Aldh1l1.pos.pulldown", "purple2", 
                                                       ifelse(Grouping_vector == "Aldh1l1.LPS.pos.pulldown", "darkturquoise",
                                                              ifelse(Grouping_vector == "RiboTag.bulk.cortex", "dodgerblue", NA)))))))

legend_groups <- c("control.bulk.cortex", "control.LPS.bulk.cortex", "Aldh1l1.bulk.cortex", "Aldh1l1.RiboTag.bulk.cortex", "Aldh1l1.pos.pulldown","Aldh1l1.LPS.pos.pulldown", "RiboTag.bulk.cortex")

plotMDS_1_4_allgroups_nonegs <- plotMDS(log2(normalized_counts_2), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for log2(mean mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("magenta", "pink4", "coral", "green3", "purple", "darkturquoise", "dodgerblue"), 
       pch = 16, title = "Transcriptome groups",cex=1.4)

dev.off()

###############################################################################
# DESeq2-based PCA for Aldh1l1.bulk.cortex, Aldh1l1.RiboTag.bulk.cortex, Aldh1l1.pos.pulldown and Aldh1l1.RiboTag.pos.pulldown

dim(counts_clean_filtered_rows_9_2_8_1)
# [1] 17054    31

dim(traits_9_2_8_1)
# [1] 31  2

# updated clean_counts
TurboID_Vs_RiboTag.bulk.cortexs_pulls_clean_counts <- counts_clean_filtered_rows_9_2_8_1[,c(4:6, 14:16, 24:26, 29:31)]

# updated traits
TurboID_Vs_RiboTag.bulk.cortexs_pulls_traits <- traits_9_2_8_1[c(4:6, 14:16, 24:26, 29:31),]

#                   Sample.ID..                        Group
# 21047FL-134-01-31   4-RNA 30L                Aldh1l1.bulk.cortex
# 21047FL-134-01-32   5-RNA 31B                Aldh1l1.bulk.cortex
# 21047FL-134-01-33   6-RNA 24L                Aldh1l1.bulk.cortex
# 21047FL-134-01-41  15-RNA 30L         Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B         Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L         Aldh1l1.pos.pulldown
# 23238R-04-03              CR3        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-04              CR4        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-05              CR5        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-08              CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09              CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10             CR10 Aldh1l1.RiboTag.pos.pulldown

# create DESeq data set matrix dds for all 20 samples
dds_TurboID_Vs_RiboTag.bulk.cortexs_pulls <- DESeqDataSetFromMatrix(countData = TurboID_Vs_RiboTag.bulk.cortexs_pulls_clean_counts,
                                         colData = TurboID_Vs_RiboTag.bulk.cortexs_pulls_traits,
                                         design = ~Group)
dim(dds_TurboID_Vs_RiboTag.bulk.cortexs_pulls)
# [1] 17054    12

# set the factor level
dds_TurboID_Vs_RiboTag.bulk.cortexs_pulls$Group <- relevel(dds_TurboID_Vs_RiboTag.bulk.cortexs_pulls$Group, ref = "Aldh1l1.bulk.cortex")

# run DESeq
dds_TurboID_Vs_RiboTag.bulk.cortexs_pulls <- DESeq(dds_TurboID_Vs_RiboTag.bulk.cortexs_pulls)

# create pdf output file
pdf(file = "1c-1. Exp. 9-2 and 8-1_RNA-seq_TurboID_Vs_RiboTag.bulk.cortexs_pulls_DESeq2_PCA_CR_06182024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

# perform DESeq2-basedPCA on rlog transformed data of all 20 samples
rldr <- rlog(dds_TurboID_Vs_RiboTag.bulk.cortexs_pulls, blind = TRUE)
plotPCA(rldr, intgroup = "Group") # plots based on normalized counts

pcaData <- plotPCA(rldr, intgroup = "Group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca3 <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
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

print(pca3)

dev.off()

################################################################################
# 1d Aldh1l1.pos.pulldown vs Aldh1l1.Ribotag.pos.pulldown
# subset DESeq2 norm counts to only include Aldh1l1 samples

#Updated clean_counts
Aldh1l1_clean_counts <- counts_clean_filtered_rows_9_2_8_1[,c(14:16, 29:31)]

#Updated traits
Aldh1l1_traits <- traits_9_2_8_1[c(14:16, 29:31),]

rownames(Aldh1l1_traits)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43" "23238R-04-08"      "23238R-04-09"     
# [6] "23238R-04-10" 

colnames(Aldh1l1_clean_counts)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43" "23238R-04-08"      "23238R-04-09"     
# [6] "23238R-04-10" 


#create DESeq data set (dds) matrix for the 3 Aldh1l1 and 3 Camk2a pos.pulldowns
Aldh1l1_dds <- DESeqDataSetFromMatrix(countData = Aldh1l1_clean_counts, #filtered, unnorm counts is the input here
                                              colData = Aldh1l1_traits,
                                              design = ~Group)
dim(Aldh1l1_dds)
# [1] 17054     6

# set the factor level
Aldh1l1_dds$Group <- relevel(Aldh1l1_dds$Group, ref = "Aldh1l1.RiboTag.pos.pulldown")

#Run DESeq
Aldh1l1_dds<- DESeq(Aldh1l1_dds)


# Create PCA plots
pdf(file = "1d-1. Exp. 9-2 and 8-1_RNA-seq_Aldh1l1 pulldown groups_DESeq2_PCA_CR_06182024.pdf", width = 8, height = 8.5)
# perform DESeq2-basedPCA on rlog transformed data of 16 samples w/o neg pulldowns

rldr <- rlog(Aldh1l1_dds, blind = TRUE)
plotPCA(rldr, intgroup = "Group") # plots based on normalized counts

pcaData <- plotPCA(rldr, intgroup = "Group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca4 <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
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

print(pca4)

dev.off()

################################################################################
# perform limma:plotMDS-based PCA of 18 samples w/o neg pulldowns
pdf(file = "1c-2. Exp. 9-2 and 8-1_RNA-seq_Aldh1l1 groups_Limma_PCA_CR_06182024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector_Aldh1l1 <- Aldh1l1_traits$Group

Aldh1l1_norm_counts <- normalized_counts_2[,c(4:11, 14:21)]

pch_values <- ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.bulk.cortex", 16,
                     ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.RiboTag.bulk.cortex", 16,
                            ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.pos.pulldown", 16,
                                   ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.LPS.pos.pulldown", 16, NA))))


pt_colors <- ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.bulk.cortex", "coral",
                    ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.RiboTag.bulk.cortex", "green3", 
                           ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.pos.pulldown", "purple", 
                                  ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.LPS.pos.pulldown",  "darkturquoise", NA))))

legend_groups <- c("Aldh1l1.bulk.cortex", "Aldh1l1.RiboTag.bulk.cortex", "Aldh1l1.pos.pulldown","Aldh1l1.LPS.pos.pulldown")

plotMDS_Aldh1l1<- plotMDS(log2(Aldh1l1_norm_counts), #top = 500, 
                          labels = NULL, pch = pch_values, col = pt_colors, 
                          cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                          xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)
mtext(side=3, text="MDS Plot for log2(mean mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("coral", "green3", "purple", "darkturquoise"), 
       pch = 16, title = "Transcriptome groups", cex = 1.4)

dev.off()

###############################################################################
# 1d. Aldh1l1.pos.pulldowns vs RiboTag.bulk.cortexs
# subset DESeq2 norm counts to only include Aldh1l1 and RiboTag.bulk.cortex samples

#Updated clean_counts
Aldh1l1_Vs_RiboTag.bulk.cortex_clean_counts <- counts_clean_filtered_rows_9_2_8_1[,c(14:16, 22:23)]

#Updated traits
Aldh1l1_Vs_RiboTag.bulk.cortex_traits <- traits_9_2_8_1[c(14:16, 22:23),]
#                    Sample.ID..                Group
# 21047FL-134-01-41  15-RNA 30L Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L Aldh1l1.pos.pulldown
# 23238R-04-01              CR1        RiboTag.bulk.cortex
# 23238R-04-02              CR2        RiboTag.bulk.cortex

colnames(Aldh1l1_Vs_RiboTag.bulk.cortex_clean_counts)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43" "23238R-04-01"      "23238R-04-02" 

#create DESeq data set (dds) matrix for the 3 Aldh1l1 and 3 Camk2a pos.pulldowns
Aldh1l1_Vs_RiboTag.bulk.cortex_dds <- DESeqDataSetFromMatrix(countData = Aldh1l1_Vs_RiboTag.bulk.cortex_clean_counts, #filtered, unnorm counts is the input here
                                              colData = Aldh1l1_Vs_RiboTag.bulk.cortex_traits,
                                              design = ~Group)
dim(Aldh1l1_Vs_RiboTag.bulk.cortex_dds)
# [1] 16825     6

# set the factor level
Aldh1l1_Vs_RiboTag.bulk.cortex_dds$Group <- relevel(Aldh1l1_Vs_RiboTag.bulk.cortex_dds$Group, ref = "RiboTag.bulk.cortex")

#Run DESeq
Aldh1l1_Vs_RiboTag.bulk.cortex_dds<- DESeq(Aldh1l1_Vs_RiboTag.bulk.cortex_dds)


# Create PCA plots
pdf(file = "1d-1. Exp. 9-2 and 8-1_RNA-seq_Aldh1l1 vs RiboTag.bulk.cortexs_DESeq2_PCA_CR_06182024.pdf", width = 8, height = 8.5)
# perform DESeq2-basedPCA on rlog transformed data of 16 samples w/o neg pulldowns

rldr <- rlog(Aldh1l1_Vs_RiboTag.bulk.cortex_dds, blind = TRUE)
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
pdf(file = "1d-2. Exp. 9-2 and 8-1_RNA-seq_Aldh1l1 vs Camk2a_Limma_PCA_CR_06182024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector_Aldh1l1VsCamk2a <- Aldh1l1VsCamk2a_traits$Group

Aldh1l1VsCamk2a_norm_counts <- normalized_counts_2[,c(12:14, 20:22)]

colnames(Aldh1l1VsCamk2a_norm_counts)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43" "22082R-12-01"      "22082R-12-02"     
# [6] "22082R-12-03"  

pch_values <- ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Aldh1l1.pos.pulldown", 16,
                     ifelse(Grouping_vector_Aldh1l1VsCamk2a== "RiboTag.bulk.cortex", 16, NA))


pt_colors <- ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Aldh1l1.pos.pulldown", "purple", 
                    ifelse(Grouping_vector_Aldh1l1VsCamk2a== "RiboTag.bulk.cortex",  "dodgerblue", NA))

legend_groups <- c("Aldh1l1.pos.pulldown", "RiboTag.bulk.cortex")

plotMDS_Aldh1l1<- plotMDS(log2(Aldh1l1VsCamk2a_norm_counts), #top = 500, 
                          labels = NULL, pch = pch_values, col = pt_colors, 
                          cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                          xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)
mtext(side=3, text="MDS Plot for log2(mean mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("purple", "dodgerblue"), 
       pch = 16, title = "Transcriptome groups", cex = 1.4)

dev.off()


################################################################################
################################################################################
# Step 2. Correlations of mean mRNA counts
  # a. Aldh1l1.bulk.cortex (n=3/group) vs Aldh1l1.RiboTag.bulk.cortex  (n=2/group)
  # b. Aldh1l1.RiboTag.bulk.cortex (n=2/group) vs Aldh1l1.RiboTag.pos.pulldown (n=3/group)
  # c. Aldh1l1.pos.pulldown (n=3/group) vs Aldh1l1.RiboTag.pos.pulldown (n=3/group)

# check for finite values
any(is.finite(normalized_counts_2))
# [1] TRUE zeroes are considered finite

dim(normalized_counts_2)
# [1] 17054    27 (4 neg pulldowns removed)

# CR: Because there are zeroes in the DESeq2 normalized counts data matrix, 
# the log2 transformation changes these values to infinite values (-ifn).
# A correlation coefficient cannot be calculated with infinite values, so
# I am going to add a "pseudocount" of + 1 to all the counts values to
# avoid taking a log of 0.

norm_counts_pseudo <- normalized_counts_2 + 1
any(is.finite(norm_counts_pseudo))

norm_counts_log <- log2(norm_counts_pseudo)
any(is.finite(norm_counts_log))

dim(norm_counts_log)
# [1] 17054    27 (4 neg pulldowns removed)

# load library to generate plot
library(ggplot2)

################################################################################
# 2a. Aldh1l1.bulk.cortex (n=3/group) vs Aldh1l1.RiboTag.bulk.cortex  (n=2/group)

pdf(file = "2a-c. Exp. 9-2 and 8-1_RNA-seq_counts correlations_CR_06182024.pdf", width = 6, height = 6, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# subset into groups of interest and take the row mean
Aldh1l1.bulk.cortex <- rowMeans(norm_counts_log[,c(4:6)])
Aldh1l1.RiboTag.bulk.cortex <- rowMeans(norm_counts_log[,c(22:24)])

# generate a linear correlation
correlation_2a <- cor(Aldh1l1.bulk.cortex, Aldh1l1.RiboTag.bulk.cortex)
correlation_2a
# [1] 0.9788304

# convert data to a data frame
df_2a <- data.frame(Aldh1l1.bulk.cortex = Aldh1l1.bulk.cortex, Aldh1l1.RiboTag.bulk.cortex = Aldh1l1.RiboTag.bulk.cortex)

# create correlation plot
correlation_plot_2a <- ggplot(df_2a, aes(x = Aldh1l1.bulk.cortex, y = Aldh1l1.RiboTag.bulk.cortex)) +
  geom_point(color = "black", fill = "grey", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Aldh1l1.bulk.cortex log2(mean mRNA counts)",
       y = "Aldh1l1.RiboTag.bulk.cortex log2(mean mRNA counts)") +
  ggtitle("2a. Aldh1l1.bulk.cortex vs Aldh1l1.RiboTag.bulk.cortex") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(Aldh1l1.bulk.cortex), y = min(Aldh1l1.RiboTag.bulk.cortex), 
           label = paste("Correlation coefficient:", round(correlation_2a, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2a$Aldh1l1.bulk.cortex), y = min(df_2a$Aldh1l1.RiboTag.bulk.cortex) - 2, 
           label = paste("\nNumber of genes:", nrow(df_2a)), 
           hjust = 0, vjust = 1, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2a)

################################################################################
# 2b.  Aldh1l1.RiboTag.bulk.cortex (n=2/group) vs Aldh1l1.RiboTag.pos.pulldown (n=3/group)
# subset into groups of interest and take the row mean
Aldh1l1.RiboTag.bulk.cortex <- rowMeans(norm_counts_log[,c(22:24)])
Aldh1l1.RiboTag.pos.pulldown <- rowMeans(norm_counts_log[,c(25:27)])

# convert data to a data frame
df_2b <- data.frame(Aldh1l1.RiboTag.bulk.cortex = Aldh1l1.RiboTag.bulk.cortex, Aldh1l1.RiboTag.pos.pulldown  = Aldh1l1.RiboTag.pos.pulldown)

# calculate correlation coefficient
correlation_2b <- cor(Aldh1l1.RiboTag.bulk.cortex, Aldh1l1.RiboTag.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others
correlation_2b
# [1] 0.9169184

correlation_plot_2b <- ggplot(df_2b, aes(x = Aldh1l1.RiboTag.bulk.cortex, y = Aldh1l1.RiboTag.pos.pulldown)) +
  geom_point(color = "black", fill = "red4", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Aldh1l1.RiboTag.bulk.cortex log2(mean mRNA counts)",
       y = "Aldh1l1.RiboTag.pos.pulldown log2(mean mRNA counts)") +
  ggtitle("2b.  Aldh1l1.Ribotag bulk.cortex vs pulldown") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(Aldh1l1.RiboTag.bulk.cortex), y = min(Aldh1l1.RiboTag.pos.pulldown), 
           label = paste("Correlation coefficient:", round(correlation_2b, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2b$Aldh1l1.RiboTag.bulk.cortex), y = min(df_2b$Aldh1l1.RiboTag.pos.pulldown) - 2, 
           label = paste("Number of genes:", nrow(df_2a)), # # of genes are the same for each df
           hjust = 0, vjust = 1, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2b)

################################################################################
# 2c.Aldh1l1.pos.pulldown (n=3/group) vs Aldh1l1.RiboTag.pos.pulldown (n=3/group)
# subset into groups of interest and take the row mean
Aldh1l1.pos.pulldown <- rowMeans(norm_counts_log[,c(12:14)])
Aldh1l1.RiboTag.pos.pulldown <- rowMeans(norm_counts_log[,c(25:27)])

# convert data to a data frame
df_2c <- data.frame(Aldh1l1.pos.pulldown = Aldh1l1.pos.pulldown, Aldh1l1.RiboTag.pos.pulldown = Aldh1l1.RiboTag.pos.pulldown)

# calculate correlation coefficient
correlation_2c <- cor(Aldh1l1.pos.pulldown, Aldh1l1.RiboTag.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others
correlation_2c
# [1] 0.9326028

correlation_plot_2c <- ggplot(df_2c, aes(x = Aldh1l1.pos.pulldown, y = Aldh1l1.RiboTag.pos.pulldown)) +
  geom_point(color = "black", fill = "red2", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Aldh1l1.pos.pulldown log2(mean mRNA counts)",
       y = "Aldh1l1.RiboTag log2(mean mRNA counts)") +
  ggtitle("2c. Aldh1l1 TurboID vs RiboTag pulldown RNA") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(Aldh1l1.pos.pulldown), y = min(Aldh1l1.RiboTag.pos.pulldown), 
           label = paste("Correlation coefficient:", round(correlation_2c, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2c$Aldh1l1.pos.pulldown), y = min(df_2c$Aldh1l1.RiboTag.pos.pulldown) - 2, 
           label = paste("Number of genes:", nrow(df_2a)), # # of genes are the same for each df
           hjust = 1, vjust = 0, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2c)

dev.off()


################################################################################
################################################################################
# Step 3 Part I: Diffex volcanoes

# C2: Astrocyte-RiboTag pulldown vs bulk.cortex
# C4: Astrocyte-RiboTag vs astrocyte-TurboID pulldown
################################################################################
#load Libraries
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)
library(ggplot2)
library(ggrepel)

################################################################################
# C2: Astrocyte-RiboTag pulldown vs bulk.cortex
pdf(file = "C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_06182024.pdf", height = 8, width = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

final_clean_counts <- counts_clean_filtered_rows_9_2_8_1
dim(final_clean_counts)
# [1] 17054    31

#Updated clean_counts for comparison 1
C2_clean_counts <- final_clean_counts[,c(24:26, 29:31)]

#Updated traits_9_2_8_1 for comparison 1
C2_traits_9_2_8_1 <- traits_9_2_8_1[c(24:26, 29:31),]
C2_traits_9_2_8_1
#               Sample.ID..                        Group
# 23238R-04-03         CR3        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-04         CR4        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-05         CR5        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-08         CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09         CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10        CR10 Aldh1l1.RiboTag.pos.pulldown

#Deseq2 data matrix is needed for downstream analysis
dds2 <- DESeqDataSetFromMatrix(countData = C2_clean_counts,
                               colData = C2_traits_9_2_8_1,
                               design = ~Group)   
dim(dds2)
# [1] 17054     6

dds2

# set the factor level
dds2$Group <- relevel(dds2$Group, ref = "Aldh1l1.RiboTag.bulk.cortex")

#Run DESeq 
dds2 <- DESeq(dds2)

#diff p adj value 
res05 <- results(dds2, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 17051 with nonzero bulk.cortex read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5895, 35%
# LFC < 0 (down)     : 5347, 31%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)

dim(res05)
# [1] 17054     6

#contrasts
resultsNames(dds2)
# [1] "Intercept"                                                  
# [2] "Group_Aldh1l1.RiboTag.pos.pulldown_vs_Aldh1l1.RiboTag.bulk.cortex"

#plot MA
#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res05, main = "Comparison 2: Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds2, blind=TRUE)
de2 <- as.data.frame(res05)
de2$padj[is.na(de2$padj)] <- 1

write.csv(de2, file = "C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_06182024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054") 

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_Up_06182024.csv")
write.csv(down, file = "C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_Down_06182024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_Up_06182024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_Down_06182024.csv", header =  TRUE, row.names = 1)

#bulk.cortex df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 5512    6

#volcano plot of diffex with padj df
de2 <- read.csv("C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_06182024.csv")

de2$diffexpressed <- "NS"
de2$diffexpressed[de2$log2FoldChange > 1.0 & de2$padj < 0.05] <- "UP"
de2$diffexpressed[de2$log2FoldChange < -1.0 & de2$padj < 0.05] <- "DOWN"

de2$delabel <- NA
de2$delabel[de2$diffexpressed != "NS"] <- de2$X[de2$diffexpressed != "NS"]

ggplot(data=de2, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("black", "grey", "red4")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(-5, 300)) +
  scale_x_continuous(name = "log2(difference) 
Aldh1l1.RiboTag.pulldown - Aldh1l1.RiboTag bulk.cortex RNA",
limits = c(-10, 10))

#Perform PCA on rlog transformed data
rldr <- rlog(dds2, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

################################################################################
################################################################################
# C4: Astrocyte-RiboTag pulldown vs Astrocyte-TurboID pulldown
#use base pdf function to extract plots into a single pdf
pdf(file = "C4_Aldh1l1.TurboID vs Aldh1l1.RiboTag pulldowns_rowmeans_17054_Diffex_06182024.pdf", height = 8, width = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

final_clean_counts <- counts_clean_filtered_rows_9_2_8_1
dim(final_clean_counts)
# [1] 17054    31

#Updated clean_counts for comparison 1
C4_clean_counts <- final_clean_counts[,c(14:16, 29:31)]

#Updated traits_9_2_8_1 for comparison 1
C4_traits_9_2_8_1 <- traits_9_2_8_1[c(14:16, 29:31),]
C4_traits_9_2_8_1
#                  Sample.ID..                        Group
# 21047FL-134-01-41  15-RNA 30L         Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B         Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L         Aldh1l1.pos.pulldown
# 23238R-04-08              CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09              CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10             CR10 Aldh1l1.RiboTag.pos.pulldown

#Deseq2 data matrix is needed for downstream analysis
dds4 <- DESeqDataSetFromMatrix(countData = C4_clean_counts,
                               colData = C4_traits_9_2_8_1,
                               design = ~Group)   
dim(dds4)
# [1] 17054     6

dds4

# set the factor level
dds4$Group <- relevel(dds4$Group, ref = "Aldh1l1.RiboTag.pos.pulldown")

#Run DESeq 
dds4 <- DESeq(dds4)

#diff p adj value 
res05 <- results(dds4, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 17054 with nonzero bulk.cortex read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 4519, 26%
# LFC < 0 (down)     : 4555, 27%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

dim(res05)
# [1] 17054     6

#contrasts
resultsNames(dds4)
# [1] "Intercept"                                                 
# [2] "Group_Aldh1l1.pos.pulldown_vs_Aldh1l1.RiboTag.pos.pulldown"

#plot MA
#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res05, main = "Comparison 4: Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds4, blind=TRUE)
de4 <- as.data.frame(res05)
de4$padj[is.na(de4$padj)] <- 1

write.csv(de4, file = "C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054_Diffex_06182024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054") 

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054_Diffex_Up_06182024.csv")
write.csv(down, file = "C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054_Diffex_Down_06182024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054_Diffex_Up_06182024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054_Diffex_Down_06182024.csv", header =  TRUE, row.names = 1)

#bulk.cortex df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 4362    6

#volcano plot of diffex with padj df
de4 <- read.csv("C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054_Diffex_06182024.csv")

de4$diffexpressed <- "NS"
de4$diffexpressed[de4$log2FoldChange > 1.0 & de4$padj < 0.05] <- "UP"
de4$diffexpressed[de4$log2FoldChange < -1.0 & de4$padj < 0.05] <- "DOWN"

de4$delabel <- NA
de4$delabel[de4$diffexpressed != "NS"] <- de4$X[de4$diffexpressed != "NS"]

ggplot(data=de4, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("red4", "grey", "purple4")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(-5, 200)) +
  scale_x_continuous(name = "log2(difference) 
Aldh1l1.TurboID - Aldh1l1.RiboTag pulldown RNA",
limits = c(-30, 30))

#Perform PCA on rlog transformed data
rldr <- rlog(dds1, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

################################################################################
################################################################################
# Step 3 Part II: Diffex GSEA
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
# Relevant diffex comparisons:
  # C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_06182024.csv
  # C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054_Diffex_06182024.csv

###############################################################################
# Comparsion 2: Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex
# organize DESeq2 diffex output file to mirror LFQ-MS ANOVAout file
C2_ANOVA<-read.csv("C2_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_06182024.csv", header=TRUE, row.names=1)
C2_ANOVAout<-C2_ANOVA[,c(4,1,5,2)]
colnames(C2_ANOVAout)[c(3,4)]<-c("Aldh1l1.RiboTag.pull-Aldh1l1.RiboTag.bulk.cortex","diff Aldh1l1.RiboTag.pull.Aldh1l1.RiboTag.bulk.cortex")
write.csv(C2_ANOVAout, file = "C2_RNA_ANOVAout_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_rowmeans_17054_Diffex_06182024.csv")

# GOpar
flip=c(0)
outFilename <- "C2_RNA_Aldh1l1.RiboTag.pulldown vs Aldh1l1.RiboTag.bulk.cortex_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C2_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

###############################################################################
# Comparison 4: Astrocyte-TurboID pulldown vs Astrocyte-RiboTag pulldown

# organize DESeq2 diffex output file to mirror LFQ-MS ANOVAout file
C4_ANOVA<-read.csv("C4_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054_Diffex_06182024.csv", header=TRUE, row.names=1)
C4_ANOVAout<-C4_ANOVA[,c(4,1,5,2)]
colnames(C4_ANOVAout)[c(3,4)]<-c("Aldh1l1.TurboID-Aldh1l1.RiboTag","diff Aldh1l1.TurboID.Aldh1l1.Ribotag")
write.csv(C4_ANOVAout, file = "C4_RNA_ANOVAout_Aldh1l1.pos.pulldown vs Aldh1l1.RiboTag.pos.pulldown_rowmeans_17054_Diffex_06182024.csv")

# GOpar
flip=c(0)
outFilename <- "C4_RNA_Aldh1l1.TurboID vs Aldh1l1.RiboTag pulldowns_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C4_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.


####################################################
# heatmap of bulk vs pulldown RNA
#Load necessary libraries
library(pheatmap)
library(dplyr)
library(grDevices)  # for colorRampPalette

# Load top 50 cell-type gene list from Barres lab
top50_cell_type_genes_Sloan <- read.csv("mouse_cell_markers_Sloan_03062023.csv", header = TRUE)

# Create a mapping of genes to cell types
gene_to_cell_type <- stack(top50_cell_type_genes_Sloan)
colnames(gene_to_cell_type) <- c("Gene", "CellType")

# Load expression data
final_clean_counts <- read.csv("1b.Exp. 9-2 and 8-1_RNA-seq_clean_counts_diffex filter-17054x31_CR_06182024.csv", row.names = 1, check.names=FALSE)

dim(final_clean_counts)
# [1] 17054    31

# Load traits data
traits_9_2_8_1 <- read.csv("Exp. 9-2 and 8-1 RNA traits.csv", row.names = 1)
dim(traits_9_2_8_1)
# [1] 31  2

################################################################################
# Comparison 2: Aldh1l1.RiboTag.pos.pulldown vs Aldh1l1.RiboTag.bulk.cortex
# subset counts data

# Subset counts data
Aldh1l1.RiboTag_pullVsBulk_clean_counts <- final_clean_counts[,c(24:26, 29:31)]

# Updated traits_9_2_8_1
Aldh1l1.RiboTag_pullVsBulk_traits <- traits_9_2_8_1[c(24:26, 29:31),]
Aldh1l1.RiboTag_pullVsBulk_traits
#              Sample.ID..                        Group
# 23238R-04-03         CR3        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-04         CR4        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-05         CR5        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-08         CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09         CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10        CR10 Aldh1l1.RiboTag.pos.pulldown

# Check if sample names match
if(!all(row.names(Aldh1l1.RiboTag_pullVsBulk_traits) %in% colnames(Aldh1l1.RiboTag_pullVsBulk_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1.RiboTag_pullVsBulk_traits$Group
names(group_info) <- row.names(Aldh1l1.RiboTag_pullVsBulk_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes <- gene_to_cell_type %>% filter(CellType == "Astrocyte") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1.RiboTag_pullVsBulk_filtered_data <- Aldh1l1.RiboTag_pullVsBulk_clean_counts[rownames(Aldh1l1.RiboTag_pullVsBulk_clean_counts) %in% astrocyte_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.RiboTag_pullVsBulk_filtered_data, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.RiboTag_pullVsBulk_filtered_data)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C2_Exp. 8-1_RNA-seq_Aldh1l1.RiboTag_pullVsBulk_Zhang_heatmap_07022024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
RiboTag_astrocyte_heatmap <- pheatmap(row_z_scores,
                                      cluster_rows = FALSE,
                                      cluster_cols = FALSE,
                                      annotation_col = annotation_col,
                                      show_rownames = TRUE,
                                      show_colnames = TRUE,
                                      color = custom_colors)
###############################################################################
# Filter for neuronal genes
neuronal_genes <- gene_to_cell_type %>% filter(CellType == "Neuron") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1.RiboTag_pullVsBulk_filtered_data <- Aldh1l1.RiboTag_pullVsBulk_clean_counts[rownames(Aldh1l1.RiboTag_pullVsBulk_clean_counts) %in% neuronal_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.RiboTag_pullVsBulk_filtered_data, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.RiboTag_pullVsBulk_filtered_data)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
RiboTag_neuron_heatmap <- pheatmap(row_z_scores,
                                   cluster_rows = FALSE,
                                   cluster_cols = FALSE,
                                   annotation_col = annotation_col,
                                   show_rownames = TRUE,
                                   show_colnames = TRUE,
                                   color = custom_colors)
dev.off()
################################################################################
################################################################################
# Comparison 4: Astrocyte-TurboID pulldown vs Astrocyte-RiboTag pulldown

# Subset counts data
Aldh1l1_TurboVsRibo_clean_counts <- final_clean_counts[,c(4:6, 29:31)]

# Updated traits_9_2_8_1
Aldh1l1_TurboVsRibo_traits <- traits_9_2_8_1[c(14:16, 29:31),]
#                   Sample.ID..                        Group
# 21047FL-134-01-41  15-RNA 30L         Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B         Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L         Aldh1l1.pos.pulldown
# 23238R-04-08              CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09              CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10             CR10 Aldh1l1.RiboTag.pos.pulldown

# Check if sample names match
if(!all(row.names(Aldh1l1_TurboVsRibo_traits) %in% colnames(Aldh1l1_TurboVsRibo_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1_TurboVsRibo_traits$Group
names(group_info) <- row.names(Aldh1l1_TurboVsRibo_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes <- gene_to_cell_type %>% filter(CellType == "Astrocyte") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1_TurboVsRibo_filtered_data <- Aldh1l1_TurboVsRibo_clean_counts[rownames(Aldh1l1_TurboVsRibo_clean_counts) %in% astrocyte_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_TurboVsRibo_filtered_data, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_TurboVsRibo_filtered_data)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C4_Exp. 8-1_RNA-seq_Aldh1l1_TurboVsRibo_pulls_Zhang_heatmap_07022024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_pulls_astrocyte_heatmap <- pheatmap(row_z_scores,
                                            cluster_rows = FALSE,
                                            cluster_cols = FALSE,
                                            annotation_col = annotation_col,
                                            show_rownames = TRUE,
                                            show_colnames = TRUE,
                                            color = custom_colors)
###############################################################################
# Filter for neuronal genes
neuronal_genes <- gene_to_cell_type %>% filter(CellType == "Neuron") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1_TurboVsRibo_filtered_data <- Aldh1l1_TurboVsRibo_clean_counts[rownames(Aldh1l1_TurboVsRibo_clean_counts) %in% neuronal_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_TurboVsRibo_filtered_data, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_TurboVsRibo_filtered_data)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_pulls_neuron_heatmap <- pheatmap(row_z_scores,
                                         cluster_rows = FALSE,
                                         cluster_cols = FALSE,
                                         annotation_col = annotation_col,
                                         show_rownames = TRUE,
                                         show_colnames = TRUE,
                                         color = custom_colors)
dev.off()

################################################################################
##################################################################################
# Filtering for cell-type markers from Zhang et al., 2014 and Sharma et al., 2015

# load Zhang and Sharma cell-type marker union list from Johnson et al. 2022 Nat Neurosci
Zhang_Sharma_markers <- read.csv("Johnson et al.. 2022 Nat Neurosci Sharma and Zhang marker union list_mouse_from Eric.csv", header = TRUE)

# Create a mapping of genes to cell types
gene_to_cell_type_2 <- stack(Zhang_Sharma_markers)
colnames(gene_to_cell_type_2) <- c("Gene", "CellType")

################################################################################
# Comparison 2: Aldh1l1.RiboTag.pos.pulldown vs Aldh1l1.RiboTag.bulk.cortex
# subset counts data

# Subset counts data
Aldh1l1.RiboTag_pullVsBulk_clean_counts <- final_clean_counts[,c(24:26, 29:31)]

# Updated traits_9_2_8_1
Aldh1l1.RiboTag_pullVsBulk_traits <- traits_9_2_8_1[c(24:26, 29:31),]
Aldh1l1.RiboTag_pullVsBulk_traits
#              Sample.ID..                        Group
# 23238R-04-03         CR3        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-04         CR4        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-05         CR5        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-08         CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09         CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10        CR10 Aldh1l1.RiboTag.pos.pulldown

# Check if sample names match
if(!all(row.names(Aldh1l1.RiboTag_pullVsBulk_traits) %in% colnames(Aldh1l1.RiboTag_pullVsBulk_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1.RiboTag_pullVsBulk_traits$Group
names(group_info) <- row.names(Aldh1l1.RiboTag_pullVsBulk_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes_2 <- gene_to_cell_type_2 %>% filter(CellType == "Astrocytes") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1.RiboTag_pullVsBulk_filtered_data_2 <- Aldh1l1.RiboTag_pullVsBulk_clean_counts[rownames(Aldh1l1.RiboTag_pullVsBulk_clean_counts) %in% astrocyte_genes_2, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.RiboTag_pullVsBulk_filtered_data_2, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.RiboTag_pullVsBulk_filtered_data_2)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

sum(is.na(row_z_scores))
# [1] 12
sum(is.nan(row_z_scores))
# [1] 12
sum(is.infinite(row_z_scores))
# [1] 0

# Identify rows with NA, NaN values
rows_with_na <- apply(row_z_scores, 1, function(x) any(is.na(x)))
rows_with_nan <- apply(row_z_scores, 1, function(x) any(is.nan(x)))

# Print the indices of rows with NA values
which(rows_with_na)
# Eif2s3y   Kdm5d 
#   320     525 
which(rows_with_nan)
# Eif2s3y   Kdm5d 
# 320     525

#rows_to_remove <- c(320, 525)
#row_z_scores_clean <- row_z_scores[-rows_to_remove, ]
row_z_scores_clean <- na.omit(row_z_scores)
rows_with_na <- apply(row_z_scores_clean, 1, function(x) any(is.na(x)))
which(rows_with_na)
# named integer(0)

pdf(file = "C2_Exp. 8-1_RNA-seq_Aldh1l1.RiboTag_pullVsBulk_Zhang and Sharma_heatmaps_07112024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
RiboTag_astrocyte_heatmap_2 <- pheatmap(row_z_scores_clean,
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE,
                                         annotation_col = annotation_col,
                                         show_rownames = FALSE,
                                         show_colnames = TRUE,
                                         color = custom_colors,
                                         main = "Top astrocytic markers from Zhang and Sharma")
###############################################################################
# Filter for neuronal genes
neuronal_genes_2 <- gene_to_cell_type_2 %>% filter(CellType == "Neurons") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1.RiboTag_pullVsBulk_filtered_data_2 <- Aldh1l1.RiboTag_pullVsBulk_clean_counts[rownames(Aldh1l1.RiboTag_pullVsBulk_clean_counts) %in% neuronal_genes_2, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.RiboTag_pullVsBulk_filtered_data_2, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.RiboTag_pullVsBulk_filtered_data_2)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
RiboTag_neuron_heatmap_2 <- pheatmap(row_z_scores,
                                     cluster_rows = TRUE,
                                     cluster_cols = TRUE,
                                     annotation_col = annotation_col,
                                     show_rownames = FALSE,
                                     show_colnames = TRUE,
                                     color = custom_colors,
                                     main = "Top neuronal markers from Zhang and Sharma")
dev.off()

################################################################################
###############################################################################
# Filter for neuronal genes
neuronal_genes <- gene_to_cell_type %>% filter(CellType == "Neuron") %>% pull(Gene)

# Filter counts for astrocyte genes
TurboVsRibo_bulk.cortexs_filtered_data <- TurboVsRibo_bulk.cortexs_clean_counts[rownames(TurboVsRibo_bulk.cortexs_clean_counts) %in% neuronal_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(TurboVsRibo_bulk.cortexs_filtered_data, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(TurboVsRibo_bulk.cortexs_filtered_data)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_bulk.cortexs_neuron_heatmap <- pheatmap(row_z_scores,
                                          cluster_rows = FALSE,
                                          cluster_cols = FALSE,
                                          annotation_col = annotation_col,
                                          show_rownames = TRUE,
                                          show_colnames = TRUE,
                                          color = custom_colors)
dev.off()

# Subset counts data
Aldh1l1_TurboVsRibo_clean_counts <- final_clean_counts[,c(14:16, 29:31)]

# Updated traits_9_2_8_1
Aldh1l1_TurboVsRibo_traits <- traits_9_2_8_1[c(14:16, 29:31),]
Aldh1l1_TurboVsRibo_traits
#                   Sample.ID..                        Group
# 21047FL-134-01-41  15-RNA 30L         Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B         Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L         Aldh1l1.pos.pulldown
# 23238R-04-08              CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09              CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10             CR10 Aldh1l1.RiboTag.pos.pulldown

# Check if sample names match
if(!all(row.names(Aldh1l1_TurboVsRibo_traits) %in% colnames(Aldh1l1_TurboVsRibo_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1_TurboVsRibo_traits$Group
names(group_info) <- row.names(Aldh1l1_TurboVsRibo_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes_2 <- gene_to_cell_type_2 %>% filter(CellType == "Astrocytes") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1_TurboVsRibo_filtered_data_2 <- Aldh1l1_TurboVsRibo_clean_counts[rownames(Aldh1l1_TurboVsRibo_clean_counts) %in% astrocyte_genes_2, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_TurboVsRibo_filtered_data_2, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_TurboVsRibo_filtered_data_2)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C4_Exp. 8-1_RNA-seq_Aldh1l1_TurboVsRibo_pulls_Zhang and Sharma_heatmaps_07112024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_pulls_astrocyte_heatmap_2 <- pheatmap(row_z_scores,
                                            cluster_rows = TRUE,
                                            cluster_cols = TRUE,
                                            annotation_col = annotation_col,
                                            show_rownames = FALSE,
                                            show_colnames = TRUE,
                                            color = custom_colors,
                                            main = "Top astrocytic markers from Zhang and Sharma")
###############################################################################
# Filter for neuronal genes
neuronal_genes_2 <- gene_to_cell_type_2 %>% filter(CellType == "Neurons") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1_TurboVsRibo_filtered_data_2 <- Aldh1l1_TurboVsRibo_clean_counts[rownames(Aldh1l1_TurboVsRibo_clean_counts) %in% neuronal_genes_2, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_TurboVsRibo_filtered_data_2, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_TurboVsRibo_filtered_data_2)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_pulls_neuron_heatmap_2 <- pheatmap(row_z_scores,
                                           cluster_rows = TRUE,
                                           cluster_cols = TRUE,
                                           annotation_col = annotation_col,
                                           show_rownames = FALSE,
                                           show_colnames = TRUE,
                                           color = custom_colors,
                                           main = "Top neuronal markers from Zhang and Sharma")
dev.off()

################################################################################
#################################################################################
# load Seyfried et al., 2017 Zhang only list
Zhang_markers <- read.csv("Seyfried et al., 2017 Cell Syst S4 Table (Zhang list).csv", header = TRUE)
gene_to_cell_type_3 <- stack(Zhang_markers)
colnames(gene_to_cell_type_3) <- c("Gene", "CellType")

# Load expression data
final_clean_counts <- read.csv("1b.Exp. 9-2 and 8-1_RNA-seq_clean_counts_diffex filter-17054x31_CR_06182024.csv", row.names = 1, check.names=FALSE)

dim(final_clean_counts)
# [1] 17054    31

# Load traits data
traits_9_2_8_1 <- read.csv("Exp. 9-2 and 8-1 RNA traits.csv", row.names = 1)
dim(traits_9_2_8_1)
# [1] 31  2

##################################################################################
##################################################################################
# Comparison 2: Aldh1l1.RiboTag.pos.pulldown vs Aldh1l1.RiboTag.bulk.cortex

# Subset counts data
Aldh1l1.RiboTag_pullVsBulk_clean_counts <- final_clean_counts[,c(24:26, 29:31)]

# Updated traits_9_2_8_1
Aldh1l1.RiboTag_pullVsBulk_traits <- traits_9_2_8_1[c(24:26, 29:31),]
Aldh1l1.RiboTag_pullVsBulk_traits
#              Sample.ID..                        Group
# 23238R-04-03         CR3        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-04         CR4        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-05         CR5        Aldh1l1.RiboTag.bulk.cortex
# 23238R-04-08         CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09         CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10        CR10 Aldh1l1.RiboTag.pos.pulldown

# Check if sample names match
if(!all(row.names(Aldh1l1.RiboTag_pullVsBulk_traits) %in% colnames(Aldh1l1.RiboTag_pullVsBulk_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1.RiboTag_pullVsBulk_traits$Group
names(group_info) <- row.names(Aldh1l1.RiboTag_pullVsBulk_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for Astrocytes genes
Astrocyte_genes_3 <- gene_to_cell_type_3 %>% filter(CellType == "Astrocytes") %>% pull(Gene)

# Filter counts for Astrocytes genes
Aldh1l1.RiboTag_pullVsBulk_filtered_data_3 <- Aldh1l1.RiboTag_pullVsBulk_clean_counts[rownames(Aldh1l1.RiboTag_pullVsBulk_clean_counts) %in% Astrocyte_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.RiboTag_pullVsBulk_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.RiboTag_pullVsBulk_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocytes", sep = " | ")
rownames(row_z_scores) <- new_row_names

sum(is.na(row_z_scores))
# [1] 12
sum(is.nan(row_z_scores))
# [1] 12
sum(is.infinite(row_z_scores))
# [1] 0

# Identify rows with NA, NaN values
rows_with_na <- apply(row_z_scores, 1, function(x) any(is.na(x)))
rows_with_nan <- apply(row_z_scores, 1, function(x) any(is.nan(x)))

# Print the indices of rows with NA values
which(rows_with_na)
# Eif2s3y   Kdm5d 
#   320     525 
which(rows_with_nan)
# Eif2s3y   Kdm5d 
# 320     525

#rows_to_remove <- c(320, 525)
#row_z_scores_clean <- row_z_scores[-rows_to_remove, ]
row_z_scores_clean <- na.omit(row_z_scores)
rows_with_na <- apply(row_z_scores_clean, 1, function(x) any(is.na(x)))
which(rows_with_na)
# named integer(0)

pdf(file = "C2_Exp. 8-1_RNA-seq_Aldh1l1.RiboTag_pullVsBulk_Zhang all_heatmap_07122024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
RiboTag_Astrocytes_heatmap_3 <- pheatmap(row_z_scores_clean,
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE,
                                         annotation_col = annotation_col,
                                         show_rownames = FALSE,
                                         show_colnames = TRUE,
                                         color = custom_colors,
                                         main = "Top astrocytic markers from Zhang")
###############################################################################
# Filter for neuronal genes
neuronal_genes_3 <- gene_to_cell_type_3 %>% filter(CellType == "Neurons") %>% pull(Gene)

# Filter counts for Astrocytes genes
Aldh1l1.RiboTag_pullVsBulk_filtered_data_3 <- Aldh1l1.RiboTag_pullVsBulk_clean_counts[rownames(Aldh1l1.RiboTag_pullVsBulk_clean_counts) %in% neuronal_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.RiboTag_pullVsBulk_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.RiboTag_pullVsBulk_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
RiboTag_neuron_heatmap_3 <- pheatmap(row_z_scores,
                                     cluster_rows = TRUE,
                                     cluster_cols = TRUE,
                                     annotation_col = annotation_col,
                                     show_rownames = FALSE,
                                     show_colnames = TRUE,
                                     color = custom_colors,
                                     main = "Top neuronal markers from Zhang")
dev.off()

################################################################################
# Comparison 4: Astrocyte-TurboID pulldown vs Astrocyte-RiboTag pulldown

# Subset counts data
Aldh1l1_TurboVsRibo_clean_counts <- final_clean_counts[,c(14:16, 29:31)]

# Updated traits_9_2_8_1
Aldh1l1_TurboVsRibo_traits <- traits_9_2_8_1[c(14:16, 29:31),]
Aldh1l1_TurboVsRibo_traits
#                   Sample.ID..                        Group
# 21047FL-134-01-41  15-RNA 30L         Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B         Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L         Aldh1l1.pos.pulldown
# 23238R-04-08              CR8 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-09              CR9 Aldh1l1.RiboTag.pos.pulldown
# 23238R-04-10             CR10 Aldh1l1.RiboTag.pos.pulldown

# Check if sample names match
if(!all(row.names(Aldh1l1_TurboVsRibo_traits) %in% colnames(Aldh1l1_TurboVsRibo_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1_TurboVsRibo_traits$Group
names(group_info) <- row.names(Aldh1l1_TurboVsRibo_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for Astrocytes genes
Astrocyte_genes_3 <- gene_to_cell_type_3 %>% filter(CellType == "Astrocytes") %>% pull(Gene)

# Filter counts for Astrocytes genes
Aldh1l1_TurboVsRibo_filtered_data_3 <- Aldh1l1_TurboVsRibo_clean_counts[rownames(Aldh1l1_TurboVsRibo_clean_counts) %in% Astrocyte_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_TurboVsRibo_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_TurboVsRibo_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocytes", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C4_Exp. 8-1_RNA-seq_Aldh1l1_TurboVsRibo_pulls_Zhang all_heatmaps_07122024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_pulls_Astrocytes_heatmap_3 <- pheatmap(row_z_scores,
                                               cluster_rows = TRUE,
                                               cluster_cols = TRUE,
                                               annotation_col = annotation_col,
                                               show_rownames = FALSE,
                                               show_colnames = TRUE,
                                               color = custom_colors,
                                               main = "Top astrocytic markers from Zhang")
###############################################################################
# Filter for neuronal genes
neuronal_genes_3 <- gene_to_cell_type_3 %>% filter(CellType == "Neurons") %>% pull(Gene)

# Filter counts for Astrocytes genes
Aldh1l1_TurboVsRibo_filtered_data_3 <- Aldh1l1_TurboVsRibo_clean_counts[rownames(Aldh1l1_TurboVsRibo_clean_counts) %in% neuronal_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_TurboVsRibo_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_TurboVsRibo_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_pulls_neuron_heatmap_3 <- pheatmap(row_z_scores,
                                           cluster_rows = TRUE,
                                           cluster_cols = TRUE,
                                           annotation_col = annotation_col,
                                           show_rownames = FALSE,
                                           show_colnames = TRUE,
                                           color = custom_colors,
                                           main = "Top neuronal markers from Zhang")
dev.off()

################################################################################
# end of data analysis and visualization