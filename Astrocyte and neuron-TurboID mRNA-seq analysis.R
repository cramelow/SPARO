########################################################################################################
# Code for performing Differential Gene Expression Analysis on Exp. 9-2 (Aldh1l1) and 10-1 (Camk2a)
# row means >=10 filter before dds and splicing
# Christina Ramelow, MS
# Date: 02/27/2024
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

rootdir <- setwd("/Users/christina/Desktop/TurboID dual-omics/Exp. 9-2 and 10-1/3. RNA-seq")

# load feature_counts.txt file
raw_counts <- read.table("Exp. 9-2 and 10-1, 9-2 and 10-1_feature_counts.txt", sep = "\t", header = TRUE, row.names = 1)
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
traits <- read.csv("Exp. 9-2 and 10-1, 9-2 and & 10-1 combined traits_RNA.csv",  header = TRUE, row.names = 1)
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

# splice the raw_counts_clean df to separate Exp. 9-2 and 10-1, 9-2 and 10-1
# Exp. 9-2 and 10-1 spliced
# raw_counts_clean_9_2_10_1 <- raw_counts_clean[,c(1:20)]
# dim(raw_counts_clean_9_2_10_1)
# [1] 35976    20

# colnames(raw_counts_clean_9_2_10_1)
# should be "21047FL-134-01-01" - "21047FL-134-01-21" CR: checked

# write .csv file for raw_counts_clean for Exp. 9-2 and 10-1
# write.csv(raw_counts_clean_9_2_10_1, file=paste0("0.Exp. 9-2 and 10-1_RNA-seq_raw_counts_clean_unfiltered-",dim(raw_counts_clean_9_2_10_1)[1],"x",dim(raw_counts_clean_9_2_10_1)[2],"_CR_02272024.csv"))

# Exp. 9-2 and 10-1 splced
raw_counts_clean_9_2_10_1 <- raw_counts_clean[,c(21:44)]
dim(raw_counts_clean_9_2_10_1)
# [1] 35976    24

colnames(raw_counts_clean_9_2_10_1)
# should be "21047FL-134-01-28" - "22082R-12-03"  CR: checked

# write .csv file for raw_counts_clean for Exp. 9-2 and 10-1
write.csv(raw_counts_clean_9_2_10_1, file=paste0("0.Exp. 9-2 and 10-1_RNA-seq_raw_counts_clean_unfiltered-",dim(raw_counts_clean_9_2_10_1)[1],"x",dim(raw_counts_clean_9_2_10_1)[2],"_CR_02272024.csv"))

################################################################################
################################################################################
## Step 1a: Filter out low counts using a row mean >/= 10 by group for Venn diagrams
# Exp. 9-2 and 10-1
# 1. subset the traits and counts dfs by sample group to create list of genes for Venn diagrams
# 2. calculate the row mean
# 3. filter for genes with a row mean >= 10
# 4. write a .csv file for each sample group
################################################################################
# subset the traits file to only include Exp. 9-2 and 10-1
traits_9_2_10_1 <- traits[c(21:44),]
dim(traits_9_2_10_1)
# [1] 24  2
################################################################################
# row mean >= 10 Aldh1l1.bulk
counts_Aldh1l1.bulk <- raw_counts_clean_9_2_10_1[,c(4:6)]
dim(counts_Aldh1l1.bulk)
# [1] 35976     3

traits_Aldh1l1.bulk  <- traits_9_2_10_1[c(4:6),]
traits_Aldh1l1.bulk
#                 Sample.ID..         Group
# 21047FL-134-01-31   4-RNA 30L Aldh1l1.bulk
# 21047FL-134-01-32   5-RNA 31B Aldh1l1.bulk
# 21047FL-134-01-33   6-RNA 24L Aldh1l1.bulk

rownames(traits_Aldh1l1.bulk)
# [1] "21047FL-134-01-31" "21047FL-134-01-32" "21047FL-134-01-33"

colnames(counts_Aldh1l1.bulk)
# [1] "21047FL-134-01-31" "21047FL-134-01-32" "21047FL-134-01-33"

# calculate row mean
Aldh1l1.bulk_row_mean <- rowMeans(counts_Aldh1l1.bulk)

# subset rows where row means >= 10
Aldh1l1.bulk_filtered_rows <- counts_Aldh1l1.bulk[Aldh1l1.bulk_row_mean >= 10, ]

dim(Aldh1l1.bulk_filtered_rows)
# [1] 16737     3

# write .csv file 
write.csv(Aldh1l1.bulk_filtered_rows, file=paste0("1a.Exp. 9-2 and 10-1_RNA-seq_Aldh1l1.bulk_Venn filter-",dim(Aldh1l1.bulk_filtered_rows)[1],"x",dim(Aldh1l1.bulk_filtered_rows)[2],"_CR_02272024.csv"))

################################################################################
# row mean >= 10 Aldh1l1.LPS.bulk
counts_Aldh1l1.LPS.bulk <- raw_counts_clean_9_2_10_1[,c(7:11)]
dim(counts_Aldh1l1.LPS.bulk)
# [1] 35976     5
traits_Aldh1l1.LPS.bulk <- traits_9_2_10_1[c(7:11),]
traits_Aldh1l1.LPS.bulk
#                    Sample.ID..             Group
# 21047FL-134-01-34    7-RNA 15L Aldh1l1.LPS.bulk
# 21047FL-134-01-35    8-RNA 16B Aldh1l1.LPS.bulk
# 21047FL-134-01-36    9-RNA 18R Aldh1l1.LPS.bulk
# 21047FL-134-01-37 10-RNA 21 2R Aldh1l1.LPS.bulk
# 21047FL-134-01-38   11-RNA 25B Aldh1l1.LPS.bulk

rownames(traits_Aldh1l1.LPS.bulk)
# [1] "21047FL-134-01-34" "21047FL-134-01-35" "21047FL-134-01-36" "21047FL-134-01-37" "21047FL-134-01-38"

colnames(counts_Aldh1l1.LPS.bulk)
# [1] "21047FL-134-01-34" "21047FL-134-01-35" "21047FL-134-01-36" "21047FL-134-01-37" "21047FL-134-01-38"

# calculate row mean
Aldh1l1.LPS.bulk_row_mean <- rowMeans(counts_Aldh1l1.LPS.bulk)

# subset rows where row means >= 10
Aldh1l1.LPS.bulk_filtered_rows <- counts_Aldh1l1.LPS.bulk[Aldh1l1.LPS.bulk_row_mean >= 10, ]

dim(Aldh1l1.LPS.bulk_filtered_rows)
# [1] 16475     5

# write .csv file 
write.csv(Aldh1l1.bulk_filtered_rows, file=paste0("1a.Exp. 9-2 and 10-1_RNA-seq_Aldh1l1.LPS.bulk_Venn filter-",dim(Aldh1l1.LPS.bulk_filtered_rows)[1],"x",dim(Aldh1l1.LPS.bulk_filtered_rows)[2],"_CR_02272024.csv"))

################################################################################
# row mean >= 10 Aldh1l1.pos.pulldown
counts_Aldh1l1.pos.pulldown <- raw_counts_clean_9_2_10_1[,c(14:16)]
dim(counts_Aldh1l1.pos.pulldown)
# [1] 35976     3
traits_Aldh1l1.pos.pulldown <- traits_9_2_10_1[c(14:16),]
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
write.csv(Aldh1l1.pos.pulldown_filtered_rows, file=paste0("1a.Exp. 9-2 and 10-1_RNA-seq_Aldh1l1.pos.pulldown_Venn filter-",dim(Aldh1l1.pos.pulldown_filtered_rows)[1],"x",dim(Aldh1l1.pos.pulldown_filtered_rows)[2],"_CR_02272024.csv"))
################################################################################
# row mean >= 10 Aldh1l1.LPS.pos.pulldown
counts_Aldh1l1.LPS.pos.pulldown <- raw_counts_clean_9_2_10_1[,c(17:21)]
dim(counts_Aldh1l1.LPS.pos.pulldown)
# [1] 35976     3

traits_Aldh1l1.LPS.pos.pulldown<- traits_9_2_10_1[c(18:20),]
traits_Aldh1l1.LPS.pos.pulldown
#                    Sample.ID..                    Group
# 21047FL-134-01-45   19-RNA 16B Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-46   20-RNA 18R Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-47 21-RNA 21 2R Aldh1l1.LPS.pos.pulldown

rownames(traits_Aldh1l1.LPS.pos.pulldown)
# [1] "21047FL-134-01-45" "21047FL-134-01-46" "21047FL-134-01-47"

colnames(counts_Aldh1l1.LPS.pos.pulldown)
# [1] "21047FL-134-01-45" "21047FL-134-01-46" "21047FL-134-01-47"

# calculate row mean
Aldh1l1.LPS.pos.pulldown_row_mean <- rowMeans(counts_Aldh1l1.LPS.pos.pulldown)

# subset rows where row means >= 10
Aldh1l1.LPS.pos.pulldown_filtered_rows <- counts_Aldh1l1.LPS.pos.pulldown[Aldh1l1.LPS.pos.pulldown_row_mean >= 10, ]

dim(Aldh1l1.LPS.pos.pulldown_filtered_rows)
# [1] 16473     5

# write .csv file 
write.csv(Aldh1l1.LPS.pos.pulldown_filtered_rows, file=paste0("1a.Exp. 9-2 and 10-1_RNA-seq_Aldh1l1.LPS.pos.pulldown_Venn filter-",dim(Aldh1l1.LPS.pos.pulldown_filtered_rows)[1],"x",dim(Aldh1l1.LPS.pos.pulldown_filtered_rows)[2],"_CR_02272024.csv"))

################################################################################
# row mean >= 10 Camk2a.pos.pulldown
counts_Camk2a.pos.pulldown <- raw_counts_clean_9_2_10_1[,c(22:24)]
dim(counts_Camk2a.pos.pulldown)
# [1] 35976     3

traits_Camk2a.pos.pulldown<- traits_9_2_10_1[c(22:24),]
traits_Camk2a.pos.pulldown
#              Sample.ID..               Group
# 22082R-12-01          12 Camk2a.pos.pulldown
# 22082R-12-02          13 Camk2a.pos.pulldown
# 22082R-12-03          14 Camk2a.pos.pulldown

rownames(traits_Camk2a.pos.pulldown)
# [1] "22082R-12-01" "22082R-12-02" "22082R-12-03"

colnames(counts_Camk2a.pos.pulldown)
# [1] "22082R-12-01" "22082R-12-02" "22082R-12-03"

# calculate row mean
Camk2a.pos.pulldown_row_mean <- rowMeans(counts_Camk2a.pos.pulldown)

# subset rows where row means >= 10
Camk2a.pos.pulldown_filtered_rows <- counts_Camk2a.pos.pulldown[Camk2a.pos.pulldown_row_mean >= 10, ]

dim(Camk2a.pos.pulldown_filtered_rows)
# [1] 17288     3

# write .csv file 
write.csv(Camk2a.pos.pulldown_filtered_rows, file=paste0("1a.Exp. 9-2 and 10-1_RNA-seq_Camk2a.pos.pulldown_Venn filter-",dim(Camk2a.pos.pulldown_filtered_rows)[1],"x",dim(Camk2a.pos.pulldown_filtered_rows)[2],"_CR_02272024.csv"))


################################################################################
################################################################################
# clean up filtered row dfs to only include gene lists

# Aldh1l1.bulk_filtered_rows
Aldh1l1.bulk_genes <- as.matrix(rownames(Aldh1l1.bulk_filtered_rows))
dim(Aldh1l1.bulk_genes)
# [1] 16737     1

# Aldh1l1.LPS.bulk_filtered_rows
Aldh1l1.LPS.bulk_genes <- as.matrix(rownames(Aldh1l1.LPS.bulk_filtered_rows))
dim(Aldh1l1.LPS.bulk_genes)
# [1] 16475     1

# Aldh1l1.pos.pulldown_filtered_rows
Aldh1l1.pos.pulldown_genes <- as.matrix(rownames(Aldh1l1.pos.pulldown_filtered_rows))
dim(Aldh1l1.pos.pulldown_genes)
# [1] 16632     1

# Aldh1l1.LPS.pos.pulldown_filtered_rows
Aldh1l1.LPS.pos.pulldown_genes <- as.matrix(rownames(Aldh1l1.LPS.pos.pulldown_filtered_rows))
dim(Aldh1l1.LPS.pos.pulldown_genes)
# [1] 16473     1

# Camk2a.pos.pulldown_filtered_rows
Camk2a.pos.pulldown_genes <- as.matrix(rownames(Camk2a.pos.pulldown_filtered_rows))
dim(Camk2a.pos.pulldown_genes)
# [1] 17288     1

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
pdf(file = "1a. Exp. 9-2 and 10-1__RNA-seq_Venns_CR_02272024.pdf", height = 11, width = 8.5, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# Aldh1l1.bulk vs Aldh1l1.LPS.bulk
biovenn1 <- draw.venn(Aldh1l1.bulk_genes, Aldh1l1.LPS.bulk_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Aldh1l1.bulk vs Aldh1l1.LPS.bulk cortex RNA", 
                      xtitle = "bulk", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "LPS bulk",
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
                      y_c = "turquoise4",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x bulk: 16737"
# [1] "y bulk: 16475"
# [1] "z bulk: 0"
# [1] "x only: 546"
# [1] "y only: 284"
# [1] "z only: 0"
# [1] "x-y bulk overlap: 16191"

# Aldh1l1.bulk vs Aldh1l1.pos.pulldown
biovenn2 <- draw.venn(Aldh1l1.bulk_genes, Aldh1l1.pos.pulldown_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Aldh1l1.bulk vs Aldh1l1.pos.pulldown RNA", 
                      xtitle = "bulk", 
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
                      y_c = "purple4",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)
# [1] "x bulk: 16737"
# [1] "y bulk: 16632"
# [1] "z bulk: 0"
# [1] "x only: 573"
# [1] "y only: 468"
# [1] "z only: 0"
# [1] "x-y bulk overlap: 16164"


# Aldh1l1.pos.pulldown vs Aldh1l1.LPS.pos.pulldown
biovenn3 <- draw.venn(Aldh1l1.pos.pulldown_genes, Aldh1l1.LPS.pos.pulldown_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Aldh1l1 vs Aldh1l1.LPS pulldown RNA", 
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
                      y_c = "turquoise2",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x bulk: 16632"
# [1] "y bulk: 16473"
# [1] "z bulk: 0"
# [1] "x only: 458"
# [1] "y only: 299"
# [1] "z only: 0"
# [1] "x-y bulk overlap: 16174"


# Aldh1l1.LPS.bulk vs Aldh1l1.LPS.pos.pulldown
biovenn4 <- draw.venn(Aldh1l1.LPS.bulk_genes, Aldh1l1.LPS.pos.pulldown_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Aldh1l1.LPS.bulk vs Aldh1l1.LPS.pulldown RNA", 
                      xtitle = "LPS bulk", 
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
                      x_c = "darkgrey",
                      y_c = "turquoise2",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x bulk: 16475"
# [1] "y bulk: 16473"
# [1] "z bulk: 0"
# [1] "x only: 469"
# [1] "y only: 467"
# [1] "z only: 0"
# [1] "x-y bulk overlap: 16006"

# Aldh1l1.pos.pulldown vs Camk2a.pos.pulldown
biovenn5 <- draw.venn(Aldh1l1.pos.pulldown_genes, Camk2a.pos.pulldown_genes, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Aldh1l1 vs Camk2a pulldown RNA", 
                      xtitle = "Aldh1l1", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "Camk2a",
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
                      y_c = "dodgerblue",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x bulk: 16632"
# [1] "y bulk: 17288"
# [1] "z bulk: 0"
# [1] "x only: 485"
# [1] "y only: 1141"
# [1] "z only: 0"
# [1] "x-y bulk overlap: 16147"


dev.off()

################################################################################
################################################################################
## Step 1b: Create row mean filter >=10 across all samples for downstream diffex analysis

# calculate row mean
counts_clean_row_mean <- rowMeans(raw_counts_clean_9_2_10_1)

# subset rows where row means >= 10
counts_clean_filtered_rows_9_2_10_1 <- raw_counts_clean_9_2_10_1[counts_clean_row_mean  >= 10, ]

dim(counts_clean_filtered_rows_9_2_10_1)
# [1] 16825    24

# write .csv file 
write.csv(counts_clean_filtered_rows_9_2_10_1, file=paste0("1b.Exp. 9-2 and 10-1_RNA-seq_clean_counts_diffex filter-",dim(counts_clean_filtered_rows_9_2_10_1)[1],"x",dim(counts_clean_filtered_rows_9_2_10_1)[2],"_CR_02272024.csv"))

################################################################################
################################################################################
## Step 2: DESeq2 normalization of counts based on median of ratios

# load library
library(DESeq2)

# create DESeq data set matrix (dds) for all 20 samples
dds_all9_2_10_1 <- DESeqDataSetFromMatrix(countData = counts_clean_filtered_rows_9_2_10_1,
                                     colData = traits_9_2_10_1,
                                     design = ~Group)
dim(dds_all9_2_10_1)
# [1] 16825    24

# set the factor level
dds_all9_2_10_1$Group <- relevel(dds_all9_2_10_1$Group, ref = "Aldh1l1.bulk")

# run DESeq
dds_all9_2_10_1<- DESeq(dds_all9_2_10_1)

View(counts(dds_all9_2_10_1))

dds <- dds_all9_2_10_1

dds <- estimateSizeFactors(dds)

sizeFactors(dds)
# 21047FL-134-01-28 21047FL-134-01-29 21047FL-134-01-30 21047FL-134-01-31 21047FL-134-01-32 21047FL-134-01-33 
# 0.9359811         0.8443147         0.8192482         0.8977261         1.0837245         1.2079097 
# 21047FL-134-01-34 21047FL-134-01-35 21047FL-134-01-36 21047FL-134-01-37 21047FL-134-01-38 21047FL-134-01-39 
# 1.0570472         0.6690118         0.8282205         1.0290639         0.9633570         1.2909067 
# 21047FL-134-01-40 21047FL-134-01-41 21047FL-134-01-42 21047FL-134-01-43 21047FL-134-01-44 21047FL-134-01-45 
# 1.6434555         1.0804643         1.1571509         0.9388293         1.1121398         0.7473962 
# 21047FL-134-01-46 21047FL-134-01-47 21047FL-134-01-48      22082R-12-01      22082R-12-02      22082R-12-03 
# 0.8843925         0.8813361         0.8654553         1.4092297         1.4089887         0.9357295 

normalized_counts_1 <- counts(dds, normalized=TRUE)

write.csv(normalized_counts_1, file=paste0("2a.Exp. 9-2 and 10-1_RNA-seq_clean_filtered_DESeq2 norm counts_diffex filter-",dim(normalized_counts_1)[1],"x",dim(normalized_counts_1)[2],"_CR_02272024.csv"))

################################################################################
# updated clean_counts with neg.pulldowns removed
all9_2_10_1_no_negs_clean_counts <- counts_clean_filtered_rows_9_2_10_1[,c(21:31, 34:44)]

# udated traits with neg.pulldowns removed
all9_2_10_1_no_negs_traits <- traits[c(21:31, 34:44),]
all9_2_10_1_no_negs_traits


# create DESeq data set matrix dds for all 22 samples
dds_all9_2_10_1_no_negs <- DESeqDataSetFromMatrix(countData = all9_2_10_1_no_negs_clean_counts,
                                             colData = all9_2_10_1_no_negs_traits,
                                             design = ~Group)
dim(dds_all9_2_10_1_no_negs)
# [1] 16825    22

# set the factor level
dds_all9_2_10_1_no_negs$Group <- relevel(dds_all9_2_10_1_no_negs$Group, ref = "Aldh1l1.bulk")

# run DESeq
dds_all9_2_10_1_no_negs<- DESeq(dds_all9_2_10_1_no_negs)

View(counts(dds_all9_2_10_1_no_negs))

dds <- dds_all9_2_10_1_no_negs

dds <- estimateSizeFactors(dds)

sizeFactors(dds)
# 21047FL-134-01-28 21047FL-134-01-29 21047FL-134-01-30 21047FL-134-01-31 21047FL-134-01-32 21047FL-134-01-33 
# 0.9616844         0.8733660         0.8445953         0.9296244         1.1213567         1.2470030 
# 21047FL-134-01-34 21047FL-134-01-35 21047FL-134-01-36 21047FL-134-01-37 21047FL-134-01-38 21047FL-134-01-41 
# 1.0909104         0.6886373         0.8556927         1.0648926         0.9960627         1.1176854 
# 21047FL-134-01-42 21047FL-134-01-43 21047FL-134-01-44 21047FL-134-01-45 21047FL-134-01-46 21047FL-134-01-47 
# 1.1920365         0.9718914         1.1436330         0.7732927         0.9139097         0.9126985 
# 21047FL-134-01-48      22082R-12-01      22082R-12-02      22082R-12-03 
# 0.8953743         1.4549212         1.4494015         0.9647122 

normalized_counts_2 <- counts(dds, normalized=TRUE)

write.csv(normalized_counts_2, file=paste0("2a.Exp. 9-2 and 10-1_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-",dim(normalized_counts_2)[1],"x",dim(normalized_counts_2)[2],"_CR_02272024.csv"))

## end of data processing and clean up ##

################################################################################
################################################################################
## Part II: Data analysis and visualization
# Exp. 9-2 and 10-1 mRNA-seq analysis

# steps:
# 1. PCAs
  # 1a. all groups
  # 1b. all groups w/o neg.pulldowns
  # 1c. Aldh1l1.bulk +/- LPS (n=3-5/group) vs Aldh1l1.pos.pulldown (n=3-5/group)
# 2. Correlations
  # a. Aldh1l1.LPS.bulk vs Aldh1l1.bulk
  # b. Aldh1l1.bulk vs Aldh1l1.pos.pulldown
  # c. Aldh1l1.LPS.bulk vs Aldh1l1.LPS.pos.pulldown
  # d. Aldh1l1.pos.pulldown vs Aldh1l1.LPS.pos.pulldown
# 3. Diffex volcanoes & GSEA
  # a. Aldh1l1.LPS.bulk vs Aldh1l1.bulk
  # b. Aldh1l1.pos.pulldown vs Aldh1l1.bulk
  # c. Aldh1l1.LPS.bulk vs Aldh1l1.LPS.pos.pulldown
  # d. Aldh1l1.pos.pulldown vs Aldh1l1.LPS.pos.pulldown
  # e. Aldh1l1.pos.pulldown vs Camk2a.pos.pulldown

################################################################################
## Step 1: Generate PCA plots
library(DESeq2)
library(ggplot2)

# 1a: all groups
# create DESeq data set matrix dds for all 20 samples
dds_all9_2_10_1 <- DESeqDataSetFromMatrix(countData = counts_clean_filtered_rows_9_2_10_1,
                                     colData = traits_9_2_10_1,
                                     design = ~Group)
dim(dds_all9_2_10_1)
# [1] 16825    24

# set the factor level
dds_all9_2_10_1$Group <- relevel(dds_all9_2_10_1$Group, ref = "Aldh1l1.bulk")

# run DESeq
dds_all9_2_10_1<- DESeq(dds_all9_2_10_1)

# create pdf output file
pdf(file = "1a-1. Exp. 9-2 and 10-1_RNA-seq_all groups_DESeq2_PCA_CR_02272024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

# perform DESeq2-basedPCA on rlog transformed data of all 20 samples
rldr <- rlog(dds_all9_2_10_1, blind = TRUE)
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

pdf(file = "1a-2. Exp. 9-2 and 10-1_RNA-seq_all groups_Limma_PCA_CR_02272024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- traits_9_2_10_1$Group

pch_values <- ifelse(Grouping_vector == "control.bulk", 16,
                     ifelse(Grouping_vector == "control.LPS.bulk", 16,
                            ifelse(Grouping_vector == "Aldh1l1.bulk", 16,
                                   ifelse(Grouping_vector == "Aldh1l1.LPS.bulk", 16, 
                                          ifelse(Grouping_vector == "Aldh1l1.neg.pulldown", 16, 
                                                 ifelse(Grouping_vector == "Aldh1l1.pos.pulldown", 16, 
                                                        ifelse(Grouping_vector == "Aldh1l1.LPS.pos.pulldown", 16, 
                                                               ifelse(Grouping_vector == "Camk2a.pos.pulldown", 16, NA))))))))


pt_colors <- ifelse(Grouping_vector == "control.bulk", "magenta",
                    ifelse(Grouping_vector == "control.LPS.bulk", "pink4",
                           ifelse(Grouping_vector == "Aldh1l1.bulk", "coral",
                                  ifelse(Grouping_vector == "Aldh1l1.LPS.bulk", "green3", 
                                       ifelse(Grouping_vector == "Aldh1l1.neg.pulldown", "grey", 
                                         ifelse(Grouping_vector == "Aldh1l1.pos.pulldown", "purple2", 
                                                ifelse(Grouping_vector == "Aldh1l1.LPS.pos.pulldown", "darkturquoise",
                                                       ifelse(Grouping_vector == "Camk2a.pos.pulldown", "dodgerblue", NA))))))))
  
legend_groups <- c("control.bulk", "control.LPS.bulk", "Aldh1l1.bulk", "Aldh1l1.LPS.bulk", "Aldh1l1.neg.pulldown", "Aldh1l1.pos.pulldown","Aldh1l1.LPS.pos.pulldown", "Camk2a.pos.pulldown")

plotMDS_1_4_allgroups_nonegs <- plotMDS(log2(normalized_counts_1), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for log2(mean mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("magenta", "pink4", "coral", "green3", "grey", "purple", "darkturquoise", "dodgerblue"), 
       pch = 16, title = "Transcriptome groups",cex=1.4)

dev.off()


################################################################################
# 1b: all groups w/o neg.pulldowns
pdf(file = "1b-1. Exp. 9-2 and 10-1_RNA-seq_all groups_nonegs_DESeq2_PCA_CR_02272024.pdf", width = 8, height = 8.5)
# perform DESeq2-basedPCA on rlog transformed data of 18 samples w/o neg pulldowns

rldr <- rlog(dds_all9_2_10_1_no_negs, blind = TRUE)
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

pdf(file = "1b-2. Exp. 9-2 and 10-1_RNA-seq_all groups_nonegs_Limma_PCA_CR_02272024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- traits_9_2_10_1$Group[c(1:11, 14:24)]

pch_values <- ifelse(Grouping_vector == "control.bulk", 16,
                     ifelse(Grouping_vector == "control.LPS.bulk", 16,
                            ifelse(Grouping_vector == "Aldh1l1.bulk", 16,
                                   ifelse(Grouping_vector == "Aldh1l1.LPS.bulk", 16, 
                                                 ifelse(Grouping_vector == "Aldh1l1.pos.pulldown", 16, 
                                                        ifelse(Grouping_vector == "Aldh1l1.LPS.pos.pulldown", 16, 
                                                               ifelse(Grouping_vector == "Camk2a.pos.pulldown", 16, NA)))))))


pt_colors <- ifelse(Grouping_vector == "control.bulk", "magenta",
                    ifelse(Grouping_vector == "control.LPS.bulk", "pink4",
                           ifelse(Grouping_vector == "Aldh1l1.bulk", "coral",
                                  ifelse(Grouping_vector == "Aldh1l1.LPS.bulk", "green3", 
                                                ifelse(Grouping_vector == "Aldh1l1.pos.pulldown", "purple2", 
                                                       ifelse(Grouping_vector == "Aldh1l1.LPS.pos.pulldown", "darkturquoise",
                                                              ifelse(Grouping_vector == "Camk2a.pos.pulldown", "dodgerblue", NA)))))))

legend_groups <- c("control.bulk", "control.LPS.bulk", "Aldh1l1.bulk", "Aldh1l1.LPS.bulk", "Aldh1l1.pos.pulldown","Aldh1l1.LPS.pos.pulldown", "Camk2a.pos.pulldown")

plotMDS_1_4_allgroups_nonegs <- plotMDS(log2(normalized_counts_2), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for log2(mean mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("magenta", "pink4", "coral", "green3", "purple", "darkturquoise", "dodgerblue"), 
       pch = 16, title = "Transcriptome groups",cex=1.4)

dev.off()
################################################################################
# 1c Aldh1l1.bulk +/- LPS vs Aldh1l1.pos.pulldown +/- LPS
# subset DESeq2 norm counts to only include Aldh1l1 samples

#Updated clean_counts
Aldh1l1_clean_counts <- counts_clean_filtered_rows_9_2_10_1[,c(4:11, 14:21)]

#Updated traits
Aldh1l1_traits <- traits_9_2_10_1[c(4:11, 14:21),]

colnames(Aldh1l1_clean_counts)
# [1] "21047FL-134-01-31" "21047FL-134-01-32" "21047FL-134-01-33" "21047FL-134-01-34" "21047FL-134-01-35"
# [6] "21047FL-134-01-36" "21047FL-134-01-37" "21047FL-134-01-38" "21047FL-134-01-41" "21047FL-134-01-42"
# [11] "21047FL-134-01-43" "21047FL-134-01-44" "21047FL-134-01-45" "21047FL-134-01-46" "21047FL-134-01-47"
# [16] "21047FL-134-01-48"

#create DESeq data set (dds) matrix for the 3 Aldh1l1 and 3 Camk2a pos.pulldowns
Aldh1l1_dds <- DESeqDataSetFromMatrix(countData = Aldh1l1_clean_counts, #filtered, unnorm counts is the input here
                                              colData = Aldh1l1_traits,
                                              design = ~Group)
dim(Aldh1l1_dds)
# [1] 16825     16

# set the factor level
Aldh1l1_dds$Group <- relevel(Aldh1l1_dds$Group, ref = "Aldh1l1.bulk")

#Run DESeq
Aldh1l1_dds<- DESeq(Aldh1l1_dds)


# Create PCA plots
pdf(file = "1c-1. Exp. 9-2 and 10-1_RNA-seq_Aldh1l1 groups_DESeq2_PCA_CR_02272024.pdf", width = 8, height = 8.5)
# perform DESeq2-basedPCA on rlog transformed data of 16 samples w/o neg pulldowns

rldr <- rlog(Aldh1l1_dds, blind = TRUE)
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
pdf(file = "1c-2. Exp. 9-2 and 10-1_RNA-seq_Aldh1l1 groups_Limma_PCA_CR_02272024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector_Aldh1l1 <- Aldh1l1_traits$Group

Aldh1l1_norm_counts <- normalized_counts_2[,c(4:11, 14:21)]

pch_values <- ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.bulk", 16,
                     ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.LPS.bulk", 16,
                            ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.pos.pulldown", 16,
                                   ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.LPS.pos.pulldown", 16, NA))))


pt_colors <- ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.bulk", "coral",
                    ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.LPS.bulk", "green3", 
                           ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.pos.pulldown", "purple", 
                                  ifelse(Grouping_vector_Aldh1l1== "Aldh1l1.LPS.pos.pulldown",  "darkturquoise", NA))))

legend_groups <- c("Aldh1l1.bulk", "Aldh1l1.LPS.bulk", "Aldh1l1.pos.pulldown","Aldh1l1.LPS.pos.pulldown")

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
# 1d. Aldh1l1.pos.pulldowns vs Camk2a.pos.pulldowns
# subset DESeq2 norm counts to only include Aldh1l1 and Camk2a samples

#Updated clean_counts
Aldh1l1VsCamk2a_clean_counts <- counts_clean_filtered_rows_9_2_10_1[,c(14:16, 22:24)]

#Updated traits
Aldh1l1VsCamk2a_traits <- traits_9_2_10_1[c(14:16, 22:24),]
#                   Sample.ID..                Group
# 21047FL-134-01-41  15-RNA 30L Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L Aldh1l1.pos.pulldown
# 22082R-12-01               12  Camk2a.pos.pulldown
# 22082R-12-02               13  Camk2a.pos.pulldown
# 22082R-12-03               14  Camk2a.pos.pulldown

colnames(Aldh1l1VsCamk2a_clean_counts)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43" "22082R-12-01"      "22082R-12-02"     
# [6] "22082R-12-03"  

#create DESeq data set (dds) matrix for the 3 Aldh1l1 and 3 Camk2a pos.pulldowns
Aldh1l1VsCamk2a_dds <- DESeqDataSetFromMatrix(countData = Aldh1l1VsCamk2a_clean_counts, #filtered, unnorm counts is the input here
                                              colData = Aldh1l1VsCamk2a_traits,
                                              design = ~Group)
dim(Aldh1l1VsCamk2a_dds)
# [1] 16825     6

# set the factor level
Aldh1l1VsCamk2a_dds$Group <- relevel(Aldh1l1VsCamk2a_dds$Group, ref = "Camk2a.pos.pulldown")

#Run DESeq
Aldh1l1VsCamk2a_dds<- DESeq(Aldh1l1VsCamk2a_dds)


# Create PCA plots
pdf(file = "1d-1. Exp. 9-2 and 10-1_RNA-seq_Aldh1l1 vs Camk2a_DESeq2_PCA_CR_02272024.pdf", width = 8, height = 8.5)
# perform DESeq2-basedPCA on rlog transformed data of 16 samples w/o neg pulldowns

rldr <- rlog(Aldh1l1VsCamk2a_dds, blind = TRUE)
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
################################################################################
# CR edited June 2024 to fix plot dimension within the pdf function
pdf(file = "1d-2. Exp. 9-2 and 10-1_RNA-seq_Aldh1l1 vs Camk2a_Limma_PCA_CR_02272024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

Grouping_vector_Aldh1l1VsCamk2a <- Aldh1l1VsCamk2a_traits$Group

Aldh1l1VsCamk2a_norm_counts <- normalized_counts_2[,c(12:14, 20:22)]

colnames(Aldh1l1VsCamk2a_norm_counts)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43" "22082R-12-01"      "22082R-12-02"     
# [6] "22082R-12-03"  

pch_values <- ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Aldh1l1.pos.pulldown", 16,
                     ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Camk2a.pos.pulldown", 16, NA))


pt_colors <- ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Aldh1l1.pos.pulldown", "purple", 
                    ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Camk2a.pos.pulldown",  "dodgerblue", NA))

legend_groups <- c("Aldh1l1.pos.pulldown", "Camk2a.pos.pulldown")

plotMDS_Aldh1l1<- plotMDS(log2(Aldh1l1VsCamk2a_norm_counts), #top = 500, 
                          labels = NULL, pch = pch_values, col = pt_colors, 
                          cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                          xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)
mtext(side=3, text="MDS Plot for log2(mean mRNA counts)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("purple", "dodgerblue"), 
       pch = 16, title = "Transcriptome groups", cex = 1.4)

dev.off()
###############################################################################
# perform limma:plotMDS-based PCA of 18 samples w/o neg pulldowns
pdf(file = "1d-2. Exp. 9-2 and 10-1_RNA-seq_Aldh1l1 vs Camk2a_Limma_PCA_CR_02272024.pdf", width = 8, height = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector_Aldh1l1VsCamk2a <- Aldh1l1VsCamk2a_traits$Group

Aldh1l1VsCamk2a_norm_counts <- normalized_counts_2[,c(12:14, 20:22)]

colnames(Aldh1l1VsCamk2a_norm_counts)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43" "22082R-12-01"      "22082R-12-02"     
# [6] "22082R-12-03"  

pch_values <- ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Aldh1l1.pos.pulldown", 16,
                     ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Camk2a.pos.pulldown", 16, NA))


pt_colors <- ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Aldh1l1.pos.pulldown", "purple", 
                    ifelse(Grouping_vector_Aldh1l1VsCamk2a== "Camk2a.pos.pulldown",  "dodgerblue", NA))

legend_groups <- c("Aldh1l1.pos.pulldown", "Camk2a.pos.pulldown")

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
# a. Aldh1l1.LPS.bulk vs Aldh1l1.bulk
# b. Aldh1l1.bulk vs Aldh1l1.pos.pulldown
# c. Aldh1l1.LPS.bulk vs Aldh1l1.LPS.pos.pulldown
# d. Aldh1l1.LPS.pos.pulldown vs Aldh1l1.pos.pulldown

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
# 2a. Aldh1l1.bulk vs Aldh1l1.LPS.bulk

pdf(file = "2a-d. Exp. 9-2 and 10-1_RNA-seq_counts correlations_CR_02272024.pdf", width = 6, height = 6, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# subset into groups of interest and take the row mean
Aldh1l1.bulk <- rowMeans(norm_counts_log[,c(4:6)])
Aldh1l1.LPS.bulk <- rowMeans(norm_counts_log[,c(7:11)])

# generate a linear correlation
correlation_2a <- cor(Aldh1l1.bulk, Aldh1l1.LPS.bulk)
correlation_2a
# [1] 0.9910909

# convert data to a data frame
df_2a <- data.frame(Aldh1l1.bulk = Aldh1l1.bulk, Aldh1l1.LPS.bulk = Aldh1l1.LPS.bulk)

# create correlation plot
correlation_plot_2a <- ggplot(df_2a, aes(x = Aldh1l1.bulk, y = Aldh1l1.LPS.bulk)) +
  geom_point(color = "black", fill = "turquoise4", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Aldh1l1.bulk log2(mean mRNA counts)",
       y = "Aldh1l1.LPS.bulk log2(mean mRNA counts)") +
  ggtitle("2a. Aldh1l1.bulk vs Aldh1l1.LPS.bulk") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(Aldh1l1.bulk), y = min(Aldh1l1.LPS.bulk), 
           label = paste("Correlation coefficient:", round(correlation_2a, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2a$Aldh1l1.pos.pulldown), y = min(df_2a$Aldh1l1.LPS.pos.pulldown) - 2, 
           label = paste("\nNumber of genes:", nrow(df_2a)), 
           hjust = 0, vjust = 1, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2a)

################################################################################
# 2b.  Aldh1l1.bulk vs Aldh1l1.pos.pulldown
# subset into groups of interest and take the row mean
Aldh1l1.bulk <- rowMeans(norm_counts_log[,c(4:6)])
Aldh1l1.pos.pulldown <- rowMeans(norm_counts_log[,c(12:14)])

# convert data to a data frame
df_2b <- data.frame(Aldh1l1.bulk = Aldh1l1.bulk, Aldh1l1.pos.pulldown = Aldh1l1.pos.pulldown)

# calculate correlation coefficient
correlation_2b <- cor(Aldh1l1.bulk, Aldh1l1.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others
correlation_2b
# [1] 0.9717215

correlation_plot_2b <- ggplot(df_2b, aes(x = Aldh1l1.bulk, y = Aldh1l1.pos.pulldown)) +
  geom_point(color = "black", fill = "turquoise4", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Aldh1l1.bulk log2(mean mRNA counts)",
       y = "Aldh1l1.pos.pulldown log2(mean mRNA counts)") +
  ggtitle("2b.  Aldh1l1.bulk vs Aldh1l1.pos.pulldown") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(Aldh1l1.bulk), y = min(Aldh1l1.pos.pulldown), 
           label = paste("Correlation coefficient:", round(correlation_2b, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2b$Aldh1l1.pos.pulldown), y = min(df_2b$Aldh1l1.LPS.pos.pulldown) - 2, 
           label = paste("Number of genes:", nrow(df_2a)), # # of genes are the same for each df
           hjust = 0, vjust = 1, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2b)

################################################################################
# 2c. Aldh1l1.LPS.bulk vs Aldh1l1.LPS.pos.pulldown
# subset into groups of interest and take the row mean
Aldh1l1.LPS.bulk <- rowMeans(norm_counts_log[,c(7:11)])
Aldh1l1.LPS.pos.pulldown<- rowMeans(norm_counts_log[,c(15:19)])

# convert data to a data frame
df_2c <- data.frame(Aldh1l1.LPS.bulk = Aldh1l1.LPS.bulk, Aldh1l1.LPS.pos.pulldown= Aldh1l1.LPS.pos.pulldown)

# calculate correlation coefficient
correlation_2c <- cor(Aldh1l1.LPS.bulk, Aldh1l1.LPS.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others
correlation_2c
# [1] 0.9701249

correlation_plot_2c <- ggplot(df_2c, aes(x = Aldh1l1.LPS.bulk, y = Aldh1l1.LPS.pos.pulldown)) +
  geom_point(color = "black", fill = "turquoise4", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Aldh1l1.LPS.bulk log2(mean mRNA counts)",
       y = "Aldh1l1LPS pulldown log2(mean mRNA counts)") +
  ggtitle("2c. Aldh1l1.LPS.bulk vs Aldh1l1.LPS.pos.pulldown") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(Aldh1l1.LPS.bulk), y = min(Aldh1l1.LPS.pos.pulldown), 
           label = paste("Correlation coefficient:", round(correlation_2c, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2c$Aldh1l1.pos.pulldown), y = min(df_2c$Aldh1l1.LPS.pos.pulldown) - 2, 
           label = paste("Number of genes:", nrow(df_2a)), # # of genes are the same for each df
           hjust = 1, vjust = 0, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2c)

################################################################################
# 2d. Aldh1l1.LPS.pos.pulldownvs Aldh1l1.pos.pulldown
# subset into groups of interest and take the row mean
Aldh1l1.pos.pulldown <- rowMeans(norm_counts_log[,c(12:14)])
Aldh1l1.LPS.pos.pulldown<- rowMeans(norm_counts_log[,c(15:19)])

# convert data to a data frame
df_2d <- data.frame(Aldh1l1.pos.pulldown = Aldh1l1.pos.pulldown, Aldh1l1.LPS.pos.pulldown = Aldh1l1.LPS.pos.pulldown)

# calculate correlation coefficient
correlation_2d <- cor(Aldh1l1.pos.pulldown, Aldh1l1.LPS.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others
# [1] 0.9906711

correlation_plot_2d <- ggplot(df_2d, aes(x = Aldh1l1.pos.pulldown, y = Aldh1l1.LPS.pos.pulldown)) +
  geom_point(color = "black", fill = "turquoise4", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Aldh1l1.pos.pulldown log2(mean mRNA counts)",
       y = "Aldh1l1LPS pulldown log2(mean mRNA counts)") +
  ggtitle("2d. Aldh1l1.LPS.pos.pulldown vs Aldh1l1.pos.pulldown") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add black x and y axes
    axis.ticks = element_line(color = "black"),  # Add ticks to the axes
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  annotate("text", x = max(Aldh1l1.pos.pulldown), y = min(Aldh1l1.LPS.pos.pulldown), 
           label = paste("Correlation coefficient:", round(correlation_2d, 2)), 
           hjust = 1, vjust = 0, color = "black") +
  annotate("text", x = min(df_2d$Aldh1l1.pos.pulldown), y = min(df_2d$Aldh1l1.LPS.pos.pulldown) - 2, 
           label = paste("Number of genes:", nrow(df_2a)), # # of genes are the same for each df
           hjust = 1, vjust = 0, color = "black") +  # Adjusted y coordinate
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20))

print(correlation_plot_2d)

dev.off()


################################################################################
################################################################################
# Step 3 Part I: Diffex volcanoes

# C1: Cre only bulk cortex RNA vs LPS bulk cortex RNA
# C2: Cre only bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA
# C3: Control LPS bulk cortex RNA vs Astrocyte-TurboID LPS bulk cortex RNA
# C4: Astrocyte-TurboID bulk cortex RNA vs Astrocyte-TurboID LPS bulk cortex RNA
# C5: Astrocyte-TurboID bulk cortex RNA vs Astrocyte-TurboID pulldown RNA
# C6: Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID LPS pulldown RNA
# C8: Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID LPS pulldown RNA
# C11: Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA

#load Libraries
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)
library(ggplot2)
library(ggrepel)

#use base pdf function to extract plots into a single pdf
pdf(file = "C1_Cre only bulk cortexvs Control.LPS.bulk cortex RNA_rowmeans_16825_Diffex_02272024.pdf", height = 8, width = 8.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

final_clean_counts <- counts_clean_filtered_rows_9_2_10_1

#Updated clean_counts for comparison 1
C1_clean_counts <- final_clean_counts[,c(1:3)]

#Updated traits_9_2_10_1 for comparison 1
C1_traits_9_2_10_1 <- traits_9_2_10_1[c(1:3),]
C1_traits_9_2_10_1
#                   Sample.ID..             Group
# 21047FL-134-01-28 1-RNA 26 2R     control.bulk
# 21047FL-134-01-29   2-RNA 19L control.LPS.bulk
# 21047FL-134-01-30   3-RNA 14R control.LPS.bulk

#Deseq2 data matrix is needed for downstream analysis
dds1 <- DESeqDataSetFromMatrix(countData = C1_clean_counts,
                               colData = C1_traits_9_2_10_1,
                               design = ~Group)   
dim(dds1)
#[1] 16825     3

dds1

# set the factor level
dds1$Group <- relevel(dds1$Group, ref = "control.bulk")

#Run DESeq 
dds1 <- DESeq(dds1)

#diff p adj value 
res05 <- results(dds1, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 16822 with nonzero bulk read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 67, 0.4%
# LFC < 0 (down)     : 28, 0.17%
# outliers [1]       : 0, 0%
# low counts [2]     : 1957, 12%
# (mean count < 18)

dim(res05)
#[1] 16825     6

#contrasts
resultsNames(dds1)
#[1] "Intercept"                                "Group_control.LPS.bulk.RNA_vs_control.bulk.RNA"

#plot MA
#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res05, main = "Comparison 1: Control LPS bulk cortex RNA vs Control bulk")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds1, blind=TRUE)
de1 <- as.data.frame(res05)
de1$padj[is.na(de1$padj)] <- 1

write.csv(de1, file = "C1_Control LPS bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_02272024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C1_Control LPS bulk cortex RNA vs control.bulk_rowmeans_16825") 

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C1_Control LPS bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_Up_02272024.csv")
write.csv(down, file = "C1_Control LPS bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_Down_02272024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C1_Control LPS bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_Up_02272024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C1_Control LPS bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_Down_02272024.csv", header =  TRUE, row.names = 1)

#bulk df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
#[1] 72  6

#volcano plot of diffex with padj df
de1 <- read.csv("C1_Control LPS bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_02272024.csv")

de1$diffexpressed <- "NS"
de1$diffexpressed[de1$log2FoldChange > 1.0 & de1$padj < 0.05] <- "UP"
de1$diffexpressed[de1$log2FoldChange < -1.0 & de1$padj < 0.05] <- "DOWN"

de1$delabel <- NA
de1$delabel[de1$diffexpressed != "NS"] <- de1$X[de1$diffexpressed != "NS"]

ggplot(data=de1, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C1_Control LPS bulk cortex RNA vs control.bulk_rowmeans_16825") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(-5, 25)) +
  scale_x_continuous(name = "log2(difference) 
Control LPS bulk cortex RNA - control bulk",
limits = c(-10, 10))

#Perform PCA on rlog transformed data
rldr <- rlog(dds1, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

######################################################################################################################################
######################################################################################################################################
#Comparison 2 - C2 signifies the below mentioned comparison!
#2. Astrocyte-TurboID bulk cortex RNA vs control.bulk
#Hypothesis: No DEGs between groups.

#use base pdf function to extract plots into a single pdf
pdf(file = "C2_Control bulk vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.pdf", height = 8, width = 8.5)

#Updated clean_counts for comparison 2
C2_clean_counts <- final_clean_counts[,c(1, 4:6)]

#Updated traits_9_2_10_1 for comparison 2
C2_traits_9_2_10_1 <- traits_9_2_10_1[c(1, 4:6),]
C2_traits_9_2_10_1
#                   Sample.ID..         Group
# 21047FL-134-01-28 1-RNA 26 2R control.bulk
# 21047FL-134-01-31   4-RNA 30L Aldh1l1.bulk
# 21047FL-134-01-32   5-RNA 31B Aldh1l1.bulk
# 21047FL-134-01-33   6-RNA 24L Aldh1l1.bulk

#Deseq2 data matrix
dds2 <- DESeqDataSetFromMatrix(countData = C2_clean_counts,
                               colData = C2_traits_9_2_10_1,
                               design = ~Group)   

dim(dds2)
#[1] 16825     4

dds2

# set the factor level
dds2$Group <- relevel(dds2$Group, ref = "control.bulk")

#Run DESeq
dds2 <- DESeq(dds2)

#diff p adj value 
res05 <- results(dds2, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 16825 with nonzero bulk read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 119, 0.71%
# LFC < 0 (down)     : 154, 0.92%
# outliers [1]       : 0, 0%
# low counts [2]     : 3588, 21%
# (mean count < 45)

dim(res05)
#[1] 16825     6

#contrasts
resultsNames(dds2)
#[1] "Intercept"                            "Group_Aldh1l1.bulk_vs_control.bulk"

#plot MA
plotMA(res05, main = "Comparison 2_Astrocyte-TurboID bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds2, blind=TRUE)
de2 <- as.data.frame(res05)
de2$padj[is.na(de2$padj)] <- 1

write.csv(de2, file = "C2_Astrocyte-TurboID bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_02272024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C2_Astrocyte-TurboID bulk cortex RNA vs control.bulk_rowmeans_16825") 

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C2_Astrocyte-TurboID bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_Up_02272024.csv")
write.csv(down, file = "C2_Astrocyte-TurboID bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_Down_02272024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C2_Astrocyte-TurboID bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_Up_02272024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C2_Astrocyte-TurboID bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_Down_02272024.csv", header =  TRUE, row.names = 1)

#bulk df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
#[1] 88   6

#volcano plot of diffex with padj df
de2 <- read.csv("C2_Astrocyte-TurboID bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex_02272024.csv")

de2$diffexpressed <- "NS"
de2$diffexpressed[de2$log2FoldChange > 1.0 & de2$padj < 0.05] <- "UP"
de2$diffexpressed[de2$log2FoldChange < -1.0 & de2$padj < 0.05] <- "DOWN"

de2$delabel <- NA
de2$delabel[de2$diffexpressed != "NS"] <- de2$X[de2$diffexpressed != "NS"]

ggplot(data=de2, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C2_Astrocyte-TurboID bulk cortex RNA vs control.bulk_rowmeans_16825_Diffex") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(-5, 25)) +
  scale_x_continuous(name = "log2(difference) 
Astrocyte-TurboID bulk cortex RNA - control bulk",
limits = c(-5, 5))

#Perform PCA on rlog transformed data 
rldr <- rlog(dds2, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

######################################################################################################################################
######################################################################################################################################
#Comparison 3 - C3 signifies the below mentioned comparison!
#3. Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA
#Hypothesis: No DEGs between groups.

#use base pdf function to extract plots into a single pdf
pdf(file = "C3_Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825_Diffex_02272024.pdf", height = 8, width = 8.5)

#Updated clean_counts for comparison 3
C3_clean_counts <- final_clean_counts[,c(2:3, 7:11)]

#Updated traits_9_2_10_1 for comparison 3
C3_traits_9_2_10_1 <- traits_9_2_10_1[c(2:3, 7:11),]
C3_traits_9_2_10_1
#                    Sample.ID..             Group
# 21047FL-134-01-29    2-RNA 19L control.LPS.bulk
# 21047FL-134-01-30    3-RNA 14R control.LPS.bulk
# 21047FL-134-01-34    7-RNA 15L Aldh1l1.LPS.bulk
# 21047FL-134-01-35    8-RNA 16B Aldh1l1.LPS.bulk
# 21047FL-134-01-36    9-RNA 18R Aldh1l1.LPS.bulk
# 21047FL-134-01-37 10-RNA 21 2R Aldh1l1.LPS.bulk
# 21047FL-134-01-38   11-RNA 25B Aldh1l1.LPS.bulk


#Deseq2 data matrix
dds3 <- DESeqDataSetFromMatrix(countData = C3_clean_counts,
                               colData = C3_traits_9_2_10_1,
                               design = ~Group)   

dim(dds3)
#[1] 16825     7

dds3

# set the factor level
dds3$Group <- relevel(dds3$Group, ref = "control.LPS.bulk")

#Run DESeq
dds3 <- DESeq(dds3)

#diff p adj value 
res05 <- results(dds3, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 16825 with nonzero bulk read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 8, 0.048%
# LFC < 0 (down)     : 21, 0.12%
# outliers [1]       : 4, 0.024%
# low counts [2]     : 1958, 12%
# (mean count < 19)

dim(res05)
#[1] 16825     6

#contrasts
resultsNames(dds3)
#[1] "Intercept"                                    "Group_astro.CIBOP.LPS.bulk.RNA_vs_control.LPS.bulk.RNA"

#plot MA
plotMA(res05, main = "Comparison 3: Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825_Diffex")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds3, blind=TRUE)
de3 <- as.data.frame(res05)
de3$padj[is.na(de3$padj)] <- 1

write.csv(de3, file = "C3_Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=FALSE, main = "Top 50 padj C3_Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C3_Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825_Diffex_Up_02272024.csv")
write.csv(down, file = "C3_Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825_Diffex_Down_02272024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C3_Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825_Diffex_Up_02272024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C3_Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825_Diffex_Down_02272024.csv", header =  TRUE, row.names = 1)

#bulk df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 7 6

#volcano plot of diffex with padj df
de3 <- read.csv("C3_Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv")

de3$diffexpressed <- "NS"
de3$diffexpressed[de3$log2FoldChange > 1.0 & de3$padj < 0.05] <- "UP"
de3$diffexpressed[de3$log2FoldChange < -1.0 & de3$padj < 0.05] <- "DOWN"

de3$delabel <- NA
de3$delabel[de3$diffexpressed != "NS"] <- de3$X[de3$diffexpressed != "NS"]

ggplot(data=de3, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C3_Astrocyte-TurboID LPS bulk cortex RNA vs Control LPS bulk cortex RNA_rowmeans_16825_Diffex") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(0, 2.5)) +
  scale_x_continuous(name = "log2(difference) 
Astrocyte-TurboID LPS bulk cortex RNA - Control LPS bulk cortex RNA",
limits = c(-5, 5))

#Perform PCA on rlog transformed data 
rldr <- rlog(dds3, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

######################################################################################################################################
######################################################################################################################################
#Comparison 4 - C4 signifies the below mentioned comparison!
#4. Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA
#Hypothesis: DEGs in Astrocyte-TurboID LPS bulk cortex RNA not seen in Astrocyte-TurboID bulk cortex RNA.

#use base pdf function to extract plots into a single pdf
pdf(file = "C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.pdf", height = 8, width = 8.5)

#Updated clean_counts for comparison 4
C4_clean_counts <- final_clean_counts[,c(4:6, 7:11)]

#Updated traits_9_2_10_1 for comparison 4
C4_traits_9_2_10_1 <- traits_9_2_10_1[c(4:6, 7:11),]
C4_traits_9_2_10_1
#                    Sample.ID..             Group
# 21047FL-134-01-31    4-RNA 30L     Aldh1l1.bulk
# 21047FL-134-01-32    5-RNA 31B     Aldh1l1.bulk
# 21047FL-134-01-33    6-RNA 24L     Aldh1l1.bulk
# 21047FL-134-01-34    7-RNA 15L Aldh1l1.LPS.bulk
# 21047FL-134-01-35    8-RNA 16B Aldh1l1.LPS.bulk
# 21047FL-134-01-36    9-RNA 18R Aldh1l1.LPS.bulk
# 21047FL-134-01-37 10-RNA 21 2R Aldh1l1.LPS.bulk
# 21047FL-134-01-38   11-RNA 25B Aldh1l1.LPS.bulk

#Deseq2 data matrix
dds4 <- DESeqDataSetFromMatrix(countData = C4_clean_counts,
                               colData = C4_traits_9_2_10_1,
                               design = ~Group)   

dim(dds4)
#[1] 16825     8

dds4

# set the factor level
dds4$Group <- relevel(dds4$Group, ref = "Aldh1l1.bulk")

#Run DESeq
dds4 <- DESeq(dds4)

#diff p adj value 
res05 <- results(dds4, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 16825 with nonzero bulk read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 51, 0.3%
# LFC < 0 (down)     : 38, 0.23%
# outliers [1]       : 12, 0.071%
# low counts [2]     : 978, 5.8%
# (mean count < 13)

dim(res05)
#[1] 16825     6

#contrasts
resultsNames(dds4)
#[1] "Intercept"                                "Group_Aldh1l1.LPS.bulk_vs_Aldh1l1.bulk"

#plot MA
plotMA(res05, main = "Comparison 4: Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds4, blind=TRUE)
de4 <- as.data.frame(res05)
de4$padj[is.na(de4$padj)] <- 1

write.csv(de4, file = "C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=FALSE, main = "Top 50 padj C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_Up_02272024.csv")
write.csv(down, file = "C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_Down_02272024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_Up_02272024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_Down_02272024.csv", header =  TRUE, row.names = 1)

#bulk df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
#[1] 17  6

#volcano plot of diffex with padj df
de4 <- read.csv("C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv")

de4$diffexpressed <- "NS"
de4$diffexpressed[de4$log2FoldChange > 1.0 & de4$padj < 0.05] <- "UP"
de4$diffexpressed[de4$log2FoldChange < -1.0 & de4$padj < 0.05] <- "DOWN"

de4$delabel <- NA
de4$delabel[de4$diffexpressed != "NS"] <- de4$X[de4$diffexpressed != "NS"]

ggplot(data=de4, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(0, 5)) +
  scale_x_continuous(name = "log2(difference) 
Astrocyte-TurboID LPS bulk cortex RNA - Astrocyte-TurboID bulk cortex RNA",
limits = c(-5, 5))

#Perform PCA on rlog transformed data 
rldr <- rlog(dds4, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

#################################################################################
#################################################################################
#Comparison 5 - C5 signifies the below mentioned comparison!
#5. Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA
#Hypothesis: Astrocyte genes up in the pulldown group.

#use base pdf function to extract plots into a single pdf
pdf(file = "C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.pdf", height = 8, width = 8.5)

#Updated clean_counts for comparison 5
C5_clean_counts <- final_clean_counts[,c(4:6, 14:16)]

#Updated traits_9_2_10_1 for comparison 5
C5_traits_9_2_10_1 <- traits_9_2_10_1[c(4:6, 14:16),]
C5_traits_9_2_10_1
#                   Sample.ID..                Group
# 21047FL-134-01-31   4-RNA 30L        Aldh1l1.bulk
# 21047FL-134-01-32   5-RNA 31B        Aldh1l1.bulk
# 21047FL-134-01-33   6-RNA 24L        Aldh1l1.bulk
# 21047FL-134-01-41  15-RNA 30L Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L Aldh1l1.pos.pulldown

#Deseq2 data matrix
dds5 <- DESeqDataSetFromMatrix(countData = C5_clean_counts,
                               colData = C5_traits_9_2_10_1,
                               design = ~Group)   
dim(dds5)
#[1] 16825     6

dds5

# set the factor level
dds5$Group <- relevel(dds5$Group, ref = "Aldh1l1.bulk")

#Run DESeq
dds5 <- DESeq(dds5)

#diff p adj value 
res05 <- results(dds5, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 16825 with nonzero bulk read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1904, 11%
# LFC < 0 (down)     : 1820, 11%
# outliers [1]       : 0, 0%
# low counts [2]     : 979, 5.8%
# (mean count < 14)

dim(res05)
#[1] 16825     6

#contrasts
resultsNames(dds5)
# [1] "Intercept"                                   "Group_Aldh1l1.pos.pulldown_vs_Aldh1l1.bulk"

#plot MA
plotMA(res05, main = "C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds5, blind=TRUE)
de5<- as.data.frame(res05)
de5$padj[is.na(de5$padj)] <- 1

write.csv(de5, file = "C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_Up_02272024.csv")
write.csv(down, file = "C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_Down_02272024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_Up_02272024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_Down_02272024.csv", header =  TRUE, row.names = 1)

#bulk df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 961   6

#volcano plot of diffex with padj df
de5 <- read.csv("C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv")

de5$diffexpressed <- "NS"
de5$diffexpressed[de5$log2FoldChange > 1.0 & de5$padj < 0.05] <- "UP"
de5$diffexpressed[de5$log2FoldChange < -1.0 & de5$padj < 0.05] <- "DOWN"

de5$delabel <- NA
de5$delabel[de5$diffexpressed != "NS"] <- de5$X[de5$diffexpressed != "NS"]

ggplot(data=de5, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(-5, 50)) +
  scale_x_continuous(name = "log2(difference) 
Astrocyte-TurboID pulldown RNA - Astrocyte-TurboID bulk cortex RNA",
limits = c(-10, 10))

#Perform PCA on rlog transformed data 
rldr <- rlog(dds5, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

#################################################################################
#################################################################################
#Comparison 6 - C6 signifies the below mentioned comparison!
#6. Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA
#Hypothesis: DEGs in Astrocyte-TurboID LPS pulldown not seen in the Astrocyte-TurboID pulldown RNA.

#use base pdf function to extract plots into a single pdf
pdf(file = "C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024.pdf", height = 8, width = 8.5)

#Updated clean_counts for comparison 6
C6_clean_counts <- final_clean_counts[,c(14:21)]

#Updated traits_9_2_10_1 for comparison 6
C6_traits_9_2_10_1 <- traits_9_2_10_1[c(14:21),]
C6_traits_9_2_10_1
#                    Sample.ID..                    Group
# 21047FL-134-01-41   15-RNA 30L     Aldh1l1.pos.pulldown
# 21047FL-134-01-42   16-RNA 31B     Aldh1l1.pos.pulldown
# 21047FL-134-01-43   17-RNA 24L     Aldh1l1.pos.pulldown
# 21047FL-134-01-44   18-RNA 15L Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-45   19-RNA 16B Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-46   20-RNA 18R Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-47 21-RNA 21 2R Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-48   22-RNA 25B Aldh1l1.LPS.pos.pulldown

#Deseq2 data matrix
dds6 <- DESeqDataSetFromMatrix(countData = C6_clean_counts,
                               colData = C6_traits_9_2_10_1,
                               design = ~Group)   
dim(dds6)
# [1] 16825     8

dds6

# set the factor level
dds6$Group <- relevel(dds6$Group, ref = "Aldh1l1.pos.pulldown")

#Run DESeq
dds6 <- DESeq(dds6)

#diff p adj value 
res05 <- results(dds6, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 16825 with nonzero bulk read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 11, 0.065%
# LFC < 0 (down)     : 7, 0.042%
# outliers [1]       : 9, 0.053%
# low counts [2]     : 2935, 17%
# (mean count < 31)

dim(res05)
# [1] 16825     6

#contrasts
resultsNames(dds6)
# [1] "Intercept"                                              "Group_Aldh1l1.LPS.pos.pulldown_vs_Aldh1l1.pos.pulldown"

#plot MA
plotMA(res05, main = "C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds6, blind=TRUE)
de6 <- as.data.frame(res05)
de6$padj[is.na(de6$padj)] <- 1

write.csv(de6, file = "C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_02272024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_022720248")

# heatmap with all genes, row names excluded
var <- order(de6$padj, decreasing=F)
matrix_all <- assay(rld)
matrix_all <- matrix_all - rowMeans(matrix_all)

matrixx <- matrix_all
#thresh_val <- min(abs(max(matrix)),abs(min(matrix)))
thresh_val <- 0.5
matrixx[matrixx < -thresh_val] <- -thresh_val
matrixx[matrixx > thresh_val] <- thresh_val

annotation_data <- as.data.frame(colData(rld)[c("Group")])
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = F,cluster_cols=TRUE, main = "All genes_C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_02272024")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_Up_02272024.csv")
write.csv(down, file = "C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_Down_02272024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_Up_02272024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_Down_02272024.csv", header =  TRUE, row.names = 1)

#bulk df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 7 6

#volcano plot of diffex with padj df
de6 <- read.csv("C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_02272024.csv")

de6$diffexpressed <- "NS"
de6$diffexpressed[de6$log2FoldChange > 1.0 & de6$padj < 0.05] <- "UP"
de6$diffexpressed[de6$log2FoldChange < -1.0 & de6$padj < 0.05] <- "DOWN"

de6$delabel <- NA
de6$delabel[de6$diffexpressed != "NS"] <- de6$X[de6$diffexpressed != "NS"]

ggplot(data=de6, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_02272024") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(0, 10)) +
  scale_x_continuous(name = "log2(difference) 
Astrocyte-TurboID LPS pulldown RNA - Astrocyte-TurboID pulldown RNA",
limits = c(-10, 10))

#Perform PCA on rlog transformed data 
rldr <- rlog(dds6, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

#################################################################################
#################################################################################
#Comparison 8 - C8 signifies the below mentioned comparison!
#8. Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA
#Hypothesis: No DEGs between groups.

#use base pdf function to extract plots into a single pdf
pdf(file = "C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_02272024.pdf", height = 8, width = 8.5)

#Updated clean_counts for comparison 8
C8_clean_counts <- final_clean_counts[,c(7:11, 17:21)]

#Updated traits_9_2_10_1 for comparison 8
C8_traits_9_2_10_1 <- traits_9_2_10_1[c(7:11, 17:21),]
C8_traits_9_2_10_1
#                    Sample.ID..                    Group
# 21047FL-134-01-34    7-RNA 15L        Aldh1l1.LPS.bulk
# 21047FL-134-01-35    8-RNA 16B        Aldh1l1.LPS.bulk
# 21047FL-134-01-36    9-RNA 18R        Aldh1l1.LPS.bulk
# 21047FL-134-01-37 10-RNA 21 2R        Aldh1l1.LPS.bulk
# 21047FL-134-01-38   11-RNA 25B        Aldh1l1.LPS.bulk
# 21047FL-134-01-44   18-RNA 15L Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-45   19-RNA 16B Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-46   20-RNA 18R Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-47 21-RNA 21 2R Aldh1l1.LPS.pos.pulldown
# 21047FL-134-01-48   22-RNA 25B Aldh1l1.LPS.pos.pulldown

#Deseq2 data matrix
dds8 <- DESeqDataSetFromMatrix(countData = C8_clean_counts,
                               colData = C8_traits_9_2_10_1,
                               design = ~Group)   
dim(dds8)
# [1] 16825    10

dds8

# set the factor level
dds8$Group <- relevel(dds8$Group, ref = "Aldh1l1.LPS.bulk")

#Run DESeq
dds8 <- DESeq(dds8)

#diff p adj value 
res05 <- results(dds8, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 16825 with nonzero bulk read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2622, 16%
# LFC < 0 (down)     : 2302, 14%
# outliers [1]       : 16, 0.095%
# low counts [2]     : 0, 0%
# (mean count < 3)

dim(res05)
# [1] 16825     6

#contrasts
resultsNames(dds8)
#[1] "Intercept"                                      "Group_Aldh1l1.LPS.pos.pulldown_vs_Aldh1l1.LPS.bulk"

#plot MA
plotMA(res05, main = "C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds8, blind=TRUE)
de8<- as.data.frame(res05)
de8$padj[is.na(de8$padj)] <- 1

write.csv(de8, file = "C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv")

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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_Up_02272024.csv")
write.csv(down, file = "C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_Down_02272024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_Up_02272024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_Down_02272024.csv", header =  TRUE, row.names = 1)

#bulk df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 1227    6

#volcano plot of diffex with padj df
de8 <- read.csv("C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv")

de8$diffexpressed <- "NS"
de8$diffexpressed[de8$log2FoldChange > 1.0 & de8$padj < 0.05] <- "UP"
de8$diffexpressed[de8$log2FoldChange < -1.0 & de8$padj < 0.05] <- "DOWN"

de8$delabel <- NA
de8$delabel[de8$diffexpressed != "NS"] <- de8$X[de8$diffexpressed != "NS"]

ggplot(data=de8, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(0, 50)) +
  scale_x_continuous(name = "log2(difference) 
Astrocyte-TurboID LPS bulk cortex RNA - Astrocyte-TurboID LPS pulldown RNA",
limits = c(-5, 5))

#Perform PCA on rlog transformed data 
rldr <- rlog(dds8, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()
###############################################################################
################################################################################
# 11. Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA
#Astrocyte enriched genes in the Astrocyte-TurboID pulldown and excitatory neuron-enriched genes in the 
#Neuron-TurboID pulldown RNA.

#use base pdf function to extract plots into a single pdf
pdf(file = "C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024.pdf", height = 8, width = 8.5)

#Updated clean_counts for comparison 8
C11_clean_counts <- final_clean_counts[,c(14:16, 22:24)]

#Updated traits_9_2_10_1 for comparison 8
C11_traits_9_2_10_1 <- traits_9_2_10_1[c(14:16, 22:24),]
C11_traits_9_2_10_1
#                   Sample.ID..                Group
# 21047FL-134-01-41  15-RNA 30L Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L Aldh1l1.pos.pulldown
# 22082R-12-01               12  Camk2a.pos.pulldown
# 22082R-12-02               13  Camk2a.pos.pulldown
# 22082R-12-03               14  Camk2a.pos.pulldown

#Deseq2 data matrix
dds11 <- DESeqDataSetFromMatrix(countData = C11_clean_counts,
                               colData = C11_traits_9_2_10_1,
                               design = ~Group)   
dim(dds11)
# [1] 16825     6

dds11

# set the factor level
dds11$Group <- relevel(dds11$Group, ref = "Camk2a.pos.pulldown")

#Run DESeq
dds11 <- DESeq(dds11)

#diff p adj value 
res05 <- results(dds11, alpha = 0.05)
res05 #table w statistics

#diff ex analysis
summary(res05)
# out of 16825 with nonzero bulk read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3075, 18%
# LFC < 0 (down)     : 2253, 13%
# outliers [1]       : 26, 0.15%
# low counts [2]     : 327, 1.9%
# (mean count < 11)

dim(res05)
# [1] 16825     6

#contrasts
resultsNames(dds11)
# [1] "Intercept"                                         "Group_Aldh1l1.pos.pulldown_vs_Camk2a.pos.pulldown"

#plot MA
plotMA(res05, main = "C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024")

#use the log transform on the data set
#The above will plot the most varied genes across all conditions
rld <- rlog(dds11, blind=TRUE)
de11<- as.data.frame(res05)
de11$padj[is.na(de11$padj)] <- 1

write.csv(de11, file = "C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024.csv")

var <- order(de11$padj, decreasing=F)
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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,show_colnames = T, cluster_cols=TRUE, main = "Top 50 padj C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024")

# heatmap with all genes, row names excluded
var <- order(de11$padj, decreasing=F)
matrix_all <- assay(rld)
matrix_all <- matrix_all - rowMeans(matrix_all)

matrixx <- matrix_all
#thresh_val <- min(abs(max(matrix)),abs(min(matrix)))
thresh_val <- 0.5
matrixx[matrixx < -thresh_val] <- -thresh_val
matrixx[matrixx > thresh_val] <- thresh_val

annotation_data <- as.data.frame(colData(rld)[c("Group")])
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = F,cluster_cols=TRUE, main = "All genes_C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024")

#up and down list created to pull out diff ex genes
up <- res05[which(res05$log2FoldChange > 1 & res05$padj < .05),]
down <- res05[which(res05$log2FoldChange < -1 & res05$padj < .05),]

write.csv(up, file = "C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_Up_02272024.csv")
write.csv(down, file = "C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_Down_02272024.csv")

#The counts here need to match with the one counter above by the model!
up_genes <- read.csv("C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_Up_02272024.csv", header = TRUE, row.names = 1)
down_genes <- read.csv("C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_Down_02272024.csv", header =  TRUE, row.names = 1)

#bulk df of up and down genes based on a criteria that were extracted! 
up_down_genes <- rbind(up_genes,down_genes)

dim(up_down_genes)
# [1] 2281    6

#volcano plot of diffex with padj df
de11 <- read.csv("C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024.csv")

de11$diffexpressed <- "NS"
de11$diffexpressed[de11$log2FoldChange > 1.0 & de11$padj < 0.05] <- "UP"
de11$diffexpressed[de11$log2FoldChange < -1.0 & de11$padj < 0.05] <- "DOWN"

de11$delabel <- NA
de11$delabel[de11$diffexpressed != "NS"] <- de11$X[de11$diffexpressed != "NS"]

ggplot(data=de11, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  ggtitle(label = "C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024") +
  theme_minimal() +
  geom_text_repel(segment.colour = NA, max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), linetype = "longdash", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "longdash", col="black") +
  theme(legend.position="none") +
  scale_y_continuous(name = "-log10(padj)",
                     limits = c(0, 50)) +
  scale_x_continuous(name = "log2(difference) 
Astrocyte-TurboID pulldown RNA - Neuron-TurboID pulldown RNA",
# limits = c(-5, 5)
  )

#Perform PCA on rlog transformed data 
rldr <- rlog(dds11, blind = TRUE)
plotPCA(rldr, intgroup = "Group")

dev.off()

################################################################################
################################################################################
# Step 3 Part II: Diffex GSEA

# C4. Aldh1l1.LPS.bulk vs Aldh1l1.bulk
# C5. Aldh1l1.pos.pulldown vs Aldh1l1.bulk
# C6. Aldh1l1.LPS.pos.pulldownvs Aldh1l1.pos.pulldown
# C8. Aldh1l1.LPS.pos.pulldownvs Aldh1l1.LPS.bulk

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
  # C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv
  # C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv
  # C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_02272024.csv
  # C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv
  # C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024.csv

###############################################################################
# Comparison 4: Aldh1l1.LPS.bulk vs Aldh1l1.bulk

# organize DESeq2 diffex output file to mirror LFQ-MS ANOVAout file
C4_ANOVA<-read.csv("C4_Astrocyte-TurboID LPS bulk cortex RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv", header=TRUE, row.names=1)
C4_ANOVAout<-C4_ANOVA[,c(4,1,5,2)]
colnames(C4_ANOVAout)[c(3,4)]<-c("Aldh1l1.LPS.bulk-Aldh1l1.bulk","diff Aldh1l1.LPS.bulk.Aldh1l1.bulk")
write.csv(C4_ANOVAout, file = "C4_RNA_ANOVAout_Aldh1l1.LPS.bulk vs Aldh1l1.bulk_rowmeans_16825_Diffex_02272024.csv")

# GOpar
flip=c(0)
outFilename <- "C4_RNA_Aldh1l1.LPS.bulk_vs_Aldh1l1.bulk_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C4_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.


###############################################################################
# Comparison 5: Aldh1l1.pos.pulldown vs Aldh1l1.bulk

C5_ANOVA<-read.csv("C5_Astrocyte-TurboID pulldown RNA vs Astrocyte-TurboID bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv", header=TRUE, row.names=1)
C5_ANOVAout<-C5_ANOVA[,c(4,1,5,2)]
colnames(C5_ANOVAout)[c(3,4)]<-c("Aldh1l1.pulldown-Aldh1l1.bulk","diff Aldh1l1.pulldown.Aldh1l1.bulk")
write.csv(C5_ANOVAout, file = "C5_RNA_ANOVAout_Aldh1l1.pos.pulldown vs Aldh1l1.bulk_rowmeans_16825_Diffex_02272024.csv" )

# GOpar
flip=c(0)
outFilename <- "C5_RNA_Aldh1l1.pulldown_vs_Aldh1l1.bulk_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C5_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

###############################################################################
# Comparison 6: Aldh1l1LPS pulldown bs Aldh1l1.pos.pulldown

C6_ANOVA<-read.csv("C6_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID pulldown RNA_16825_Diffex_02272024.csv", header=TRUE, row.names=1)
C6_ANOVAout<-C6_ANOVA[,c(4,1,5,2)]
colnames(C6_ANOVAout)[c(3,4)]<-c("Aldh1l1.LPS.pulldown-Aldh1l1.pulldown","diff Aldh1l1.LPS.pulldown.Aldh1l1.pulldown")
write.csv(C6_ANOVAout, file = "C6_RNA_ANOVAout_Aldh1l1.LPS.pos.pulldown vs Aldh1l1.pos.pulldown_rowmeans_16825_Diffex_02272024.csv")

# GOpar
flip=c(0)
outFilename <- "C6_RNA_Aldh1l1.LPS.pulldown_vs_Aldh1l1.pulldown_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C6_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

###############################################################################
# Comparison 8: Aldh1l1.LPS.pos.pulldownvs Aldh1l1.LPS.bulk

C8_ANOVA<-read.csv("C8_Astrocyte-TurboID LPS pulldown RNA vs Astrocyte-TurboID LPS bulk cortex RNA_rowmeans_16825_Diffex_02272024.csv", header=TRUE, row.names=1)
C8_ANOVAout<-C8_ANOVA[,c(4,1,5,2)]
colnames(C8_ANOVAout)[c(3,4)]<-c("Aldh1l1.LPS.pulldown-Aldh1l1.LPS.bulk","diff Aldh1l1.LPS.pulldown.Aldh1l1.LPS.bulk")
write.csv(C8_ANOVAout, file = "C8_RNA_ANOVAout_Aldh1l1LPS pulldown vs Aldh1l1.LPS.bulk_rowmeans_16825_Diffex_02272024.csv" )

# GOpar
flip=c(0)
outFilename <- "C8_RNA_Aldh1l1.LPS.pos.pulldown_vs_Aldh1l1.LPS.bulk_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C8_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

###############################################################################
# Comparison 11: Aldh1l1.pos.pulldown vs Camk2a.pos.pulldown

#C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024.csv

C11_ANOVA<-read.csv("C11_Astrocyte-TurboID pulldown RNA vs Neuron-TurboID pulldown RNA_rowmeans_16825_Diffex_02272024.csv", header=TRUE, row.names=1)
C11_ANOVAout<-C11_ANOVA[,c(4,1,5,2)]
colnames(C11_ANOVAout)[c(3,4)]<-c("Aldh1l1.pulldown-Camk2a.pulldown","diff Aldh1l1.pulldown.Camk2a.pulldown")
write.csv(C11_ANOVAout, file = "C11_RNA_ANOVAout_Aldh1l1.pos.pulldown vs Camk2a.pos.pulldown_rowmeans_16825_Diffex_02272024.csv" )

# GOpar
flip=c(0)
outFilename <- "C11_RNA_Aldh1l1.pos.pulldown_vs_Camk2a.pos.pulldown_pval_DEGs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(C11_ANOVAout)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.

###############################################################################
# Step 4: Heatmaps of neuronal- and astrocytic filtered genes from  Zhang et al., 2014 (Barres lab)

# Comparisons
# C5. Aldh1l1.pulldown (astrocyte-TurboID) vs bulk cortex
# C11. Aldh1l1.pulldown (astrocyte-TurboID) vs Camk2a.pulldown (neuron-TurboID
# C12. Camk2a.pulldown (neuron-TurboID) vs bulk cortex

#Load necessary libraries
library(pheatmap)
library(dplyr)
library(grDevices)  # for colorRampPalette

# Load top 50 cell-type gene list from Barres lab
top50_cell_type_genes_Sloan <- read.csv("mouse_cell_markers_Sloan_03062023.csv", header = TRUE)

# Create a mapping of genes to cell types
gene_to_cell_type <- stack(top50_cell_type_genes_Sloan)
colnames(gene_to_cell_type) <- c("Gene", "CellType")

# Filtering for cell-type markers from Zhang et al., 2014 and Sharma et al., 2015 combined lists
# load Zhang and Sharma cell-type marker union list from Johnson et al. 2022 Nat Neurosci
Zhang_Sharma_markers <- read.csv("Johnson et al.. 2022 Nat Neurosci Sharma and Zhang marker union list_mouse_from Eric.csv", header = TRUE)

# Create a mapping of genes to cell types
gene_to_cell_type_2 <- stack(Zhang_Sharma_markers)
colnames(gene_to_cell_type_2) <- c("Gene", "CellType")

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

###############################################################################
# Comparison 5: Aldh1l1.pulldown vs bulk cortex
# Subset counts data
Aldh1l1.TurboID_pullVsBulk_clean_counts <- final_clean_counts[,c(4:6,14:16)]

# Updated traits_9_2_8_1
Aldh1l1.TurboID_pullVsBulk_traits <- traits_9_2_8_1[c(4:6,14:16),]

# Check if sample names match
if(!all(row.names(Aldh1l1.TurboID_pullVsBulk_traits) %in% colnames(Aldh1l1.TurboID_pullVsBulk_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1.TurboID_pullVsBulk_traits$Group
names(group_info) <- row.names(Aldh1l1.TurboID_pullVsBulk_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes <- gene_to_cell_type %>% filter(CellType == "Astrocyte") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1.TurboID_pullVsBulk_filtered_data <- Aldh1l1.TurboID_pullVsBulk_clean_counts[rownames(Aldh1l1.TurboID_pullVsBulk_clean_counts) %in% astrocyte_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.TurboID_pullVsBulk_filtered_data, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.TurboID_pullVsBulk_filtered_data)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C5_Exp. 9-2_RNA-seq_Aldh1l1 vs bulk_Zhang_heatmap_07022024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
astrocyte_heatmap <- pheatmap(row_z_scores,
                              cluster_rows = FALSE,
                              cluster_cols = FALSE,
                              annotation_col = annotation_col,
                              show_rownames = TRUE,
                              show_colnames = TRUE,
                              color = custom_colors)
###############################################################################
# Filter for neuronal genes
neuronal_genes <- gene_to_cell_type %>% filter(CellType == "Neuron") %>% pull(Gene)

# Filter counts for neuronal genes
Aldh1l1.TurboID_pullVsBulk_filtered_data <- Aldh1l1.TurboID_pullVsBulk_clean_counts[rownames(Aldh1l1.TurboID_pullVsBulk_clean_counts) %in% neuronal_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.TurboID_pullVsBulk_filtered_data, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.TurboID_pullVsBulk_filtered_data)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
neuron_heatmap <- pheatmap(row_z_scores,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           annotation_col = annotation_col,
                           show_rownames = TRUE,
                           show_colnames = TRUE,
                           color = custom_colors)

dev.off()

###############################################################################
# Comparison 11: Aldh1l1.pulldown vs Camk2a.pulldown
# Updated traits_9_2_10_1
traits_9_2_10_1 <- traits[c(21:31, 34:44),]
dim(traits_9_2_10_1)
# [1] 22  2

traits_9_2_10_1

# load Exp. 9-2 and 10-1 norm counts matrix
norm_counts_9_2_10_1 <- read.csv("2a.Exp. 9-2 and 10-1_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-16825x22_CR_02272024.csv", check.names = FALSE, header = TRUE, row.names = 1)

dim(norm_counts_9_2_10_1)
# [1] 16825    22

final_clean_counts <- norm_counts_9_2_10_1
dim(final_clean_counts)
# [1] 16825    22

# subset counts data
Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts <- final_clean_counts[,c(12:14, 20:22)]
dim(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts)
# [1] 16825     6
colnames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43"
# [4] "22082R-12-01"      "22082R-12-02"      "22082R-12-03" 

# subset traits 
Aldh1l1_Vs_Camk2a_traits <- traits_9_2_10_1[c(12:14, 20:22),]
Aldh1l1_Vs_Camk2a_traits
#                   Sample.ID..                Group
# 21047FL-134-01-41  15-RNA 30L Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L Aldh1l1.pos.pulldown
# 22082R-12-01               12  Camk2a.pos.pulldown
# 22082R-12-02               13  Camk2a.pos.pulldown
# 22082R-12-03               14  Camk2a.pos.pulldown

# Check if sample names match
if(!all(row.names(Aldh1l1_Vs_Camk2a_traits) %in% colnames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1_Vs_Camk2a_traits$Group
names(group_info) <- row.names(Aldh1l1_Vs_Camk2a_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes <- gene_to_cell_type %>% filter(CellType == "Astrocyte") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1_Vs_Camk2a_filtered_data <- Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts[rownames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts) %in% astrocyte_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_Vs_Camk2a_filtered_data, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_Vs_Camk2a_filtered_data)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C11_Exp. 9-2_RNA-seq_Aldh1l1 vs Camk2a_Zhang_heatmap_07022024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_Vs_Camk2a_pulls_astrocyte_heatmap <- pheatmap(row_z_scores,
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
Aldh1l1_Vs_Camk2a_filtered_data <- Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts[rownames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts) %in% neuronal_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_Vs_Camk2a_filtered_data, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_Vs_Camk2a_filtered_data)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_Vs_Camk2a_neuron_heatmap <- pheatmap(row_z_scores,
                                             cluster_rows = FALSE,
                                             cluster_cols = FALSE,
                                             annotation_col = annotation_col,
                                             show_rownames = TRUE,
                                             show_colnames = TRUE,
                                             color = custom_colors)
dev.off()
################################################################################
# Comparison 11: Aldh1l1.pulldown vs Camk2a.pulldown
# Updated traits_9_2_10_1
traits_9_2_10_1 <- traits[c(21:31, 34:44),]
dim(traits_9_2_10_1)
# [1] 22  2

traits_9_2_10_1

# load Exp. 9-2 and 10-1 norm counts matrix
norm_counts_9_2_10_1 <- read.csv("2a.Exp. 9-2 and 10-1_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-16825x22_CR_02272024.csv", check.names = FALSE, header = TRUE, row.names = 1)

dim(norm_counts_9_2_10_1)
# [1] 16825    22

final_clean_counts <- norm_counts_9_2_10_1
dim(final_clean_counts)
# [1] 16825    22

# subset counts data
Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts <- final_clean_counts[,c(12:14, 20:22)]
dim(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts)
# [1] 16825     6
colnames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43"
# [4] "22082R-12-01"      "22082R-12-02"      "22082R-12-03" 

# subset traits 
Aldh1l1_Vs_Camk2a_traits <- traits_9_2_10_1[c(12:14, 20:22),]
Aldh1l1_Vs_Camk2a_traits
#                   Sample.ID..                Group
# 21047FL-134-01-41  15-RNA 30L Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L Aldh1l1.pos.pulldown
# 22082R-12-01               12  Camk2a.pos.pulldown
# 22082R-12-02               13  Camk2a.pos.pulldown
# 22082R-12-03               14  Camk2a.pos.pulldown

# Check if sample names match
if(!all(row.names(Aldh1l1_Vs_Camk2a_traits) %in% colnames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1_Vs_Camk2a_traits$Group
names(group_info) <- row.names(Aldh1l1_Vs_Camk2a_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes_2 <- gene_to_cell_type_2 %>% filter(CellType == "Astrocytes") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1_Vs_Camk2a_filtered_data_2 <- Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts[rownames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts) %in% astrocyte_genes_2, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_Vs_Camk2a_filtered_data_2, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_Vs_Camk2a_filtered_data_2)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C11_Exp. 9-2_RNA-seq_Aldh1l1 vs Camk2a_Zhang and Sharma_heatmaps_07112024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_Vs_Camk2a_pulls_astrocyte_heatmap_2 <- pheatmap(row_z_scores,
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
Aldh1l1_Vs_Camk2a_filtered_data_2 <- Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts[rownames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts) %in% neuronal_genes_2, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_Vs_Camk2a_filtered_data_2, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_Vs_Camk2a_filtered_data_2)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_Vs_Camk2a_neuron_heatmap_2 <- pheatmap(row_z_scores,
                                               cluster_rows = TRUE,
                                               cluster_cols = TRUE,
                                               annotation_col = annotation_col,
                                               show_rownames = FALSE,
                                               show_colnames = TRUE,
                                               color = custom_colors,
                                               main = "Top neuronal markers from Zhang and Sharma")
dev.off()
###############################################################################
# Comparison 12: Camk2a.pulldown vs bulk cortex
traits_9_2_10_1 <- traits[c(21:31, 34:44),]
dim(traits_9_2_10_1)
# [1] 22  2

traits_9_2_10_1

# load Exp. 9-2 and 10-1 norm counts matrix
norm_counts_9_2_10_1 <- read.csv("2a.Exp. 9-2 and 10-1_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-16825x22_CR_02272024.csv", check.names = FALSE, header = TRUE, row.names = 1)

dim(norm_counts_9_2_10_1)
# [1] 16825    22

final_clean_counts <- norm_counts_9_2_10_1
dim(final_clean_counts)
# [1] 16825    22

# subset counts data
Camk2a_Vs_Bulk_clean_counts <- final_clean_counts[,c(4:6, 20:22)]
dim(Camk2a_Vs_Bulk_clean_counts)
# [1] 16825     6
colnames(Camk2a_Vs_Bulk_clean_counts)
# [1] "22082R-12-01"      "22082R-12-02"      "22082R-12-03"      "21047FL-134-01-31" "21047FL-134-01-32"
# [6] "21047FL-134-01-33"

# subset traits 
Camk2a_Vs_Bulk_traits <- traits_9_2_10_1[c(4:6, 20:22),]
Camk2a_Vs_Bulk_traits
#               Sample.ID..               Group

# 21047FL-134-01-31   4-RNA 30L       Aldh1l1.bulk
# 21047FL-134-01-32   5-RNA 31B       Aldh1l1.bulk
# 21047FL-134-01-33   6-RNA 24L       Aldh1l1.bulk
# 22082R-12-01               12 Camk2a.pos.pulldown
# 22082R-12-02               13 Camk2a.pos.pulldown
# 22082R-12-03               14 Camk2a.pos.pulldown

# Check if sample names match
if(!all(row.names(Camk2a_Vs_Bulk_traits) %in% colnames(Camk2a_Vs_Bulk_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Camk2a_Vs_Bulk_traits$Group
names(group_info) <- row.names(Camk2a_Vs_Bulk_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes <- gene_to_cell_type %>% filter(CellType == "Astrocyte") %>% pull(Gene)

# Filter counts for astrocyte genes
Camk2a_Vs_Bulk_filtered_data <- Camk2a_Vs_Bulk_clean_counts[rownames(Camk2a_Vs_Bulk_clean_counts) %in% astrocyte_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Camk2a_Vs_Bulk_filtered_data , 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Camk2a_Vs_Bulk_filtered_data )

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C12_Exp. 9-2_RNA-seq_Camk2a vs bulk_Zhang_heatmap_07022024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_Vs_Camk2a_pulls_astrocyte_heatmap <- pheatmap(row_z_scores,
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
Camk2a_Vs_Bulk_filtered_data <- Camk2a_Vs_Bulk_clean_counts[rownames(Camk2a_Vs_Bulk_clean_counts) %in% neuronal_genes, ]

# Calculate row z-scores
row_z_scores <- t(apply(Camk2a_Vs_Bulk_filtered_data , 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Camk2a_Vs_Bulk_filtered_data )

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Camk2a_Vs_Bulk_neuron_heatmap <- pheatmap(row_z_scores,
                                          cluster_rows = FALSE,
                                          cluster_cols = FALSE,
                                          annotation_col = annotation_col,
                                          show_rownames = TRUE,
                                          show_colnames = TRUE,
                                          color = custom_colors)

dev.off()
################################################################################
# Comparison 12: Camk2a.pulldown vs bulk cortex
traits_9_2_10_1 <- traits[c(21:31, 34:44),]
dim(traits_9_2_10_1)
# [1] 22  2

traits_9_2_10_1

# load Exp. 9-2 and 10-1 norm counts matrix
norm_counts_9_2_10_1 <- read.csv("2a.Exp. 9-2 and 10-1_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-16825x22_CR_02272024.csv", check.names = FALSE, header = TRUE, row.names = 1)

dim(norm_counts_9_2_10_1)
# [1] 16825    22

final_clean_counts <- norm_counts_9_2_10_1
dim(final_clean_counts)
# [1] 16825    22

# subset counts data
Camk2a_Vs_Bulk_clean_counts <- final_clean_counts[,c(4:6, 20:22)]
dim(Camk2a_Vs_Bulk_clean_counts)
# [1] 16825     6
colnames(Camk2a_Vs_Bulk_clean_counts)
# [1] "22082R-12-01"      "22082R-12-02"      "22082R-12-03"      "21047FL-134-01-31" "21047FL-134-01-32"
# [6] "21047FL-134-01-33"

# subset traits 
Camk2a_Vs_Bulk_traits <- traits_9_2_10_1[c(4:6, 20:22),]
Camk2a_Vs_Bulk_traits
#               Sample.ID..               Group
# 21047FL-134-01-31   4-RNA 30L       Aldh1l1.bulk
# 21047FL-134-01-32   5-RNA 31B       Aldh1l1.bulk
# 21047FL-134-01-33   6-RNA 24L       Aldh1l1.bulk
# 22082R-12-01               12 Camk2a.pos.pulldown
# 22082R-12-02               13 Camk2a.pos.pulldown
# 22082R-12-03               14 Camk2a.pos.pulldown


# Check if sample names match
if(!all(row.names(Camk2a_Vs_Bulk_traits) %in% colnames(Camk2a_Vs_Bulk_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Camk2a_Vs_Bulk_traits$Group
names(group_info) <- row.names(Camk2a_Vs_Bulk_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes_2 <- gene_to_cell_type_2 %>% filter(CellType == "Astrocytes") %>% pull(Gene)

# Filter counts for astrocyte genes
Camk2a_Vs_Bulk_filtered_data_2 <- Camk2a_Vs_Bulk_clean_counts[rownames(Camk2a_Vs_Bulk_clean_counts) %in% astrocyte_genes_2, ]

# Calculate row z-scores
row_z_scores <- t(apply(Camk2a_Vs_Bulk_filtered_data_2, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Camk2a_Vs_Bulk_filtered_data_2)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C12_Exp. 9-2_RNA-seq_Camk2a_Vs_Bulk_Zhang and Sharma_heatmaps_07112024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Camk2a_Vs_Bulk_astrocyte_heatmap_2 <- pheatmap(row_z_scores,
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
Camk2a_Vs_Bulk_filtered_data_2 <- Camk2a_Vs_Bulk_clean_counts[rownames(Camk2a_Vs_Bulk_clean_counts) %in% neuronal_genes_2, ]

# Calculate row z-scores
row_z_scores <- t(apply(Camk2a_Vs_Bulk_filtered_data_2, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Camk2a_Vs_Bulk_filtered_data_2)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Camk2a_Vs_Bulk_neuron_heatmap_2 <- pheatmap(row_z_scores,
                                            cluster_rows = TRUE,
                                            cluster_cols = TRUE,
                                            annotation_col = annotation_col,
                                            show_rownames = FALSE,
                                            show_colnames = TRUE,
                                            color = custom_colors,
                                            main = "Top neuronal markers from Zhang and Sharma")
dev.off()
################################################################################
# Seyfried et al., 2017 Zhang only list filtering

# Comparison 5: Aldh1l1.pulldown vs bulk cortex
# Subset counts data
Aldh1l1.TurboID_pullVsBulk_clean_counts <- final_clean_counts[,c(4:6,14:16)]

# Updated traits_9_2_8_1
Aldh1l1.TurboID_pullVsBulk_traits <- traits_9_2_8_1[c(4:6,14:16),]

# Check if sample names match
if(!all(row.names(Aldh1l1.TurboID_pullVsBulk_traits) %in% colnames(Aldh1l1.TurboID_pullVsBulk_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1.TurboID_pullVsBulk_traits$Group
names(group_info) <- row.names(Aldh1l1.TurboID_pullVsBulk_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes_3 <- gene_to_cell_type_3 %>% filter(CellType == "Astrocytes") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1.TurboID_pullVsBulk_filtered_data_3 <- Aldh1l1.TurboID_pullVsBulk_clean_counts[rownames(Aldh1l1.TurboID_pullVsBulk_clean_counts) %in% astrocyte_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.TurboID_pullVsBulk_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.TurboID_pullVsBulk_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C5_Exp. 9-2_RNA-seq_Aldh1l1 vs bulk_Zhang_heatmap_07152024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
astrocyte_heatmap_3 <- pheatmap(row_z_scores,
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

# Filter counts for neuronal genes
Aldh1l1.TurboID_pullVsBulk_filtered_data_3 <- Aldh1l1.TurboID_pullVsBulk_clean_counts[rownames(Aldh1l1.TurboID_pullVsBulk_clean_counts) %in% neuronal_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1.TurboID_pullVsBulk_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1.TurboID_pullVsBulk_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
neuron_heatmap_3 <- pheatmap(row_z_scores,
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             annotation_col = annotation_col,
                             show_rownames = FALSE,
                             show_colnames = TRUE,
                             color = custom_colors,
                             main = "Top neuronal markers from Zhang")

dev.off()

###############################################################################
# Comparison 11: Aldh1l1.pulldown vs Camk2a.pulldown
# Updated traits_9_2_10_1
traits_9_2_10_1 <- traits[c(21:31, 34:44),]
dim(traits_9_2_10_1)
# [1] 22  2

traits_9_2_10_1

# load Exp. 9-2 and 10-1 norm counts matrix
norm_counts_9_2_10_1 <- read.csv("2a.Exp. 9-2 and 10-1_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-16825x22_CR_02272024.csv", check.names = FALSE, header = TRUE, row.names = 1)

dim(norm_counts_9_2_10_1)
# [1] 16825    22

final_clean_counts <- norm_counts_9_2_10_1
dim(final_clean_counts)
# [1] 16825    22

# subset counts data
Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts <- final_clean_counts[,c(12:14, 20:22)]
dim(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts)
# [1] 16825     6
colnames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts)
# [1] "21047FL-134-01-41" "21047FL-134-01-42" "21047FL-134-01-43"
# [4] "22082R-12-01"      "22082R-12-02"      "22082R-12-03" 

# subset traits 
Aldh1l1_Vs_Camk2a_traits <- traits_9_2_10_1[c(12:14, 20:22),]
Aldh1l1_Vs_Camk2a_traits
#                   Sample.ID..                Group
# 21047FL-134-01-41  15-RNA 30L Aldh1l1.pos.pulldown
# 21047FL-134-01-42  16-RNA 31B Aldh1l1.pos.pulldown
# 21047FL-134-01-43  17-RNA 24L Aldh1l1.pos.pulldown
# 22082R-12-01               12  Camk2a.pos.pulldown
# 22082R-12-02               13  Camk2a.pos.pulldown
# 22082R-12-03               14  Camk2a.pos.pulldown

# Check if sample names match
if(!all(row.names(Aldh1l1_Vs_Camk2a_traits) %in% colnames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Aldh1l1_Vs_Camk2a_traits$Group
names(group_info) <- row.names(Aldh1l1_Vs_Camk2a_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes_3 <- gene_to_cell_type_3 %>% filter(CellType == "Astrocytes") %>% pull(Gene)

# Filter counts for astrocyte genes
Aldh1l1_Vs_Camk2a_filtered_data_3 <- Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts[rownames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts) %in% astrocyte_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_Vs_Camk2a_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_Vs_Camk2a_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C11_Exp. 9-2_RNA-seq_Aldh1l1 vs Camk2a_Zhang_heatmap_07152024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_Vs_Camk2a_pulls_astrocyte_heatmap_3 <- pheatmap(row_z_scores,
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

# Filter counts for astrocyte genes
Aldh1l1_Vs_Camk2a_filtered_data_3 <- Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts[rownames(Aldh1l1_Vs_Camk2a.TurboIDs_clean_counts) %in% neuronal_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Aldh1l1_Vs_Camk2a_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Aldh1l1_Vs_Camk2a_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Aldh1l1_Vs_Camk2a_neuron_heatmap_3 <- pheatmap(row_z_scores,
                                               cluster_rows = TRUE,
                                               cluster_cols = TRUE,
                                               annotation_col = annotation_col,
                                               show_rownames = FALSE,
                                               show_colnames = TRUE,
                                               color = custom_colors,
                                               main = "Top neuronal markers from Zhang")
dev.off()
################################################################################
# Comparison 12: Camk2a.pulldown vs bulk cortex
traits_9_2_10_1 <- traits[c(21:31, 34:44),]
dim(traits_9_2_10_1)
# [1] 22  2

traits_9_2_10_1

# load Exp. 9-2 and 10-1 norm counts matrix
norm_counts_9_2_10_1 <- read.csv("2a.Exp. 9-2 and 10-1_RNA-seq_clean_filtered_DESeq2 norm counts_no negs_diffex filter-16825x22_CR_02272024.csv", check.names = FALSE, header = TRUE, row.names = 1)

dim(norm_counts_9_2_10_1)
# [1] 16825    22

final_clean_counts <- norm_counts_9_2_10_1
dim(final_clean_counts)
# [1] 16825    22

# subset counts data
Camk2a_Vs_Bulk_clean_counts <- final_clean_counts[,c(4:6, 20:22)]
dim(Camk2a_Vs_Bulk_clean_counts)
# [1] 16825     6
colnames(Camk2a_Vs_Bulk_clean_counts)
# [1] "22082R-12-01"      "22082R-12-02"      "22082R-12-03"      "21047FL-134-01-31" "21047FL-134-01-32"
# [6] "21047FL-134-01-33"

# subset traits 
Camk2a_Vs_Bulk_traits <- traits_9_2_10_1[c(4:6, 20:22),]
Camk2a_Vs_Bulk_traits
#               Sample.ID..               Group
# 21047FL-134-01-31   4-RNA 30L       Aldh1l1.bulk
# 21047FL-134-01-32   5-RNA 31B       Aldh1l1.bulk
# 21047FL-134-01-33   6-RNA 24L       Aldh1l1.bulk
# 22082R-12-01               12 Camk2a.pos.pulldown
# 22082R-12-02               13 Camk2a.pos.pulldown
# 22082R-12-03               14 Camk2a.pos.pulldown


# Check if sample names match
if(!all(row.names(Camk2a_Vs_Bulk_traits) %in% colnames(Camk2a_Vs_Bulk_clean_counts))) {
  stop("Sample names in traits file do not match column names in expression data.")
}

group_info <- Camk2a_Vs_Bulk_traits$Group
names(group_info) <- row.names(Camk2a_Vs_Bulk_traits)

annotation_col <- data.frame(Group = group_info)
rownames(annotation_col) <- names(group_info)

# Filter for astrocyte genes
astrocyte_genes_3 <- gene_to_cell_type_3 %>% filter(CellType == "Astrocytes") %>% pull(Gene)

# Filter counts for astrocyte genes
Camk2a_Vs_Bulk_filtered_data_3 <- Camk2a_Vs_Bulk_clean_counts[rownames(Camk2a_Vs_Bulk_clean_counts) %in% astrocyte_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Camk2a_Vs_Bulk_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Camk2a_Vs_Bulk_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Astrocyte", sep = " | ")
rownames(row_z_scores) <- new_row_names

pdf(file = "C12_Exp. 9-2_RNA-seq_Camk2a_Vs_Bulk_Zhang_heatmaps_07112024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Camk2a_Vs_Bulk_astrocyte_heatmap_3 <- pheatmap(row_z_scores,
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

# Filter counts for astrocyte genes
Camk2a_Vs_Bulk_filtered_data_3 <- Camk2a_Vs_Bulk_clean_counts[rownames(Camk2a_Vs_Bulk_clean_counts) %in% neuronal_genes_3, ]

# Calculate row z-scores
row_z_scores <- t(apply(Camk2a_Vs_Bulk_filtered_data_3, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
colnames(row_z_scores) <- colnames(Camk2a_Vs_Bulk_filtered_data_3)

# Adjust row names to include cell type information
new_row_names <- paste(rownames(row_z_scores), "Neuron", sep = " | ")
rownames(row_z_scores) <- new_row_names

# Define a custom color palette
custom_colors <- colorRampPalette(c("lightblue", "white", "darkblue"))(100)

# Create the heatmap with annotation and custom colors
Camk2a_Vs_Bulk_neuron_heatmap_3 <- pheatmap(row_z_scores,
                                            cluster_rows = TRUE,
                                            cluster_cols = TRUE,
                                            annotation_col = annotation_col,
                                            show_rownames = FALSE,
                                            show_colnames = TRUE,
                                            color = custom_colors,
                                            main = "Top neuronal markers from Zhang")
dev.off()

#################################################################################

# end of data analysis and visualization
