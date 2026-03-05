#####################
# Exp. 13-1 HEK393 SR64 vs SR71 LV infection (3 uL) DIS LFQ-MS
# Spectronaut DIA search with human UP000005640_9606 Uniprot ID & TurboID-V5
# Christina Ramelow, PhD, MS
# Wrapper script from Sarah Shapley, PhD
# August 30th, 2025
#####################

### STEP 1 SCRIPT ###
###################################################################################################################################################
# ## Goals ###
# 0. Set up your environment
# 1. Read in protein LFQ-MS Spectronaut search
#       Load in dataset, load in traits file, Log2 transform data
# 2. Perform >50% missingness per experimental group, outer_join by common protein IDs
# 3. Impute missing values via Perseus-style imputation
#       Impute missing values so are able to do statistical analyses
###################################################################################################################################################
#######################################################################
# 0. Set up your environment
#######################################################################
#Set up your environment here
options(stringsAsFactors=FALSE)

rootdir="/Users/christina/Desktop/Exp. 13-1 HEK SR64 and SR71 DIA-MS analysis_August 2025"
datadir=data

setwd(rootdir)
library(tidyverse)
########################################################################
# 1. Read in protein groups tsv output from Spectronout
# Based off Seyfried Lab Dammer "Loader-MaxQuant SummedIntensity (proteinGroups.txt).R" script
########################################################################
##########################################
# Load required package
library(data.table)

# Set file path
input_file.prot <- "20250830_092247_Exp. 13-1 HEK SR64 vs SR71 DIA-MS_08272025_Peptide Quant Report_CR.tsv"
# semi-tryptic, filtered: "20250628_091818_Bader CSF 210 singleShot 3Cohorts (semitryptic) DDA 48 frac library search_Report.tsv"
# full-tryptic, "20250624_144145_SNA-fSN1.1_CLib1.1 with deamidation_SNA-SER1.1 BGS default_all final raw files Zet Maj Dei_Report.tsv"  # Replace with your file path

# Read in the data efficiently
df.prot <- fread(input_file.prot, sep = "\t", header = TRUE, data.table = FALSE)
dim(df.prot)
# [1] 999971     18

colnames(df.prot)
#  [1] "R.Condition"                 "R.FileName"                  "R.Replicate"                 "PG.FastaHeaders"            
# [5] "PG.Genes"                    "PG.ProteinAccessions"        "PG.ProteinDescriptions"      "PG.ProteinGroups"           
# [9] "PG.ProteinNames"             "PG.UniProtIds"               "PG.Coverage"                 "PG.IsSingleHit"             
# [13] "PG.Qvalue"                   "PG.RunEvidenceCount"         "PG.Quantity"                 "EG.PrecursorId"             
# [17] "EG.Qvalue"                   "EG.TotalQuantity (Settings)"


library(dplyr)
library(tidyr)
library(stringr)

# Step 1: Subset the minimal data
df.prot.minimal <- df.prot[, c("R.FileName", "PG.Quantity", "PG.ProteinAccessions", "PG.Genes", "PG.UniProtIds")]

# Step 2: Create UniqueID: first gene symbol | first UniProt ID
first_gene <- str_split_fixed(df.prot.minimal$PG.Genes, ";", 2)[,1]
first_uniprot <- str_split_fixed(df.prot.minimal$PG.UniProtIds, ";", 2)[,1]
df.prot.minimal$UniqueID <- paste0(first_gene, "|", first_uniprot)

# Step 3: Remove duplicates per sample/protein combo
df.prot.minimal$keepUnique <- paste0(df.prot.minimal$R.FileName, "|", df.prot.minimal$PG.ProteinAccessions)
df.prot.minimal <- df.prot.minimal[!duplicated(df.prot.minimal$keepUnique), ]

# Step 4: Pivot to wide format using UniqueID as row identifier
df_wide.prot <- pivot_wider(
  data = df.prot.minimal,
  names_from = R.FileName,
  values_from = PG.Quantity,
  id_cols = UniqueID
)

# Step 5: Set rownames
df_wide.prot <- as.data.frame(df_wide.prot)
rownames(df_wide.prot) <- df_wide.prot$UniqueID
df_wide.prot$UniqueID <- NULL  # optionally remove column now that it's rownames


dim(df_wide.prot)
# [1] 8473   12

write.csv(df_wide.prot, file=paste0("Exp. 13-1 LFQ-MS-DIA_Spectonaut_raw intensities_uniqueIDs-",dim(df_wide.prot)[1],"x",dim(df_wide.prot)[2],".csv_08302025.csv"))

################################################################################
## Load traits (numericMeta)
numericMeta<-read.csv("Exp. 13-1 traits_protein.csv", header = T)
rownames(numericMeta) <- numericMeta[,1]

dat.sumInt<-df_wide.prot[,match(rownames(numericMeta),colnames(df_wide.prot))]

dat.sumInt.log<-log2(dat.sumInt)
#dat.sumInt.log[!is.finite(dat.sumInt.log)]<-NA

dim(dat.sumInt.log)
# [1] 8473   12

## Finalize Grouping of Samples for t.test
Grouping<-numericMeta$Group.simple
Grouping
#  [1] "Global"      "Global"      "Global"      "neg_pulldn"  "neg_pulldn"  "neg_pulldn"  "SR64_pulldn" "SR64_pulldn"
# [9] "SR64_pulldn" "SR71_pulldn" "SR71_pulldn" "SR71_pulldn"

################################################################################
## Part 1: Data processing and clean up
# steps:
# 2. Normalization
# a. pos.pulldown - neg.pulldown
# b. totals and (pos.pulldown - neg.pulldown) median zero-centered
################################################################################
# background subtraction

## From the above, select your group names for negative (.neg) and positive (.pos) in paired order to be used in the for loop within the apply below. Total samples can be left out.
myGroups.neg=c("neg_pulldn")
myGroups.pos=c("SR64_pulldn","SR71_pulldn")

dat.sumInt.bgSubtr <- t(apply(dat.sumInt,1,function(x) {
  for (i in 1:length(myGroups.pos)) {
    bkgr<-mean(x[which(Grouping==myGroups.neg[i])], na.rm=T)
    bkgr[!is.finite(bkgr)]<- 0  #** when no value in the background was present, we can't subtract anything (without first imputing) -- saves 31000 values from NaN status!
    #unfortunately, we cannot subtract background for any value in a row if one group is missed! -- keep NaN values in and fix by replacing the whole row back below.
    newRow=x
    newRow[which(Grouping==myGroups.pos[i])] <- x[which(Grouping==myGroups.pos[i])] - bkgr
  }
  newRow }))


## how many background subtracted values went <0 ? (should be a very small fraction of the large number)
length(which(dat.sumInt.bgSubtr<0))
# [1] 0

# we can set these to 0, for imputation later  (note: log2(0)=-Inf)
dat.sumInt.bgSubtr[dat.sumInt.bgSubtr<0]<- 0
#dat.sumInt.bgSubtr[!is.finite(dat.sumInt.bgSubtr)]<- NA

## find rows for which no background(s) were available -- they contain NaN -- and put back the unsubtracted values for the whole row in each case
length(which(is.nan(dat.sumInt.bgSubtr)))
# [1] 58

sum(unlist(apply(dat.sumInt.bgSubtr,1,function(x) if(length(which(is.nan(x)))>0) 1)))
# [1] 53

#write.csv(dat.sumInt.bgSubtr, file = "Exp. 13-1_dat.sumInt.bgSubtr_check.csv")

# perform cleanup, returning the unsubtracted matrix rows where there is at least 1 NaN
dat.sumInt.bgSubtr2<-t(sapply(rownames(dat.sumInt.bgSubtr),function(x) {
  if(length(which(is.nan(dat.sumInt.bgSubtr[x,])))>0) { 
    #cat(paste0(x,"\n"))
    dat.sumInt[x,]
  } else {
    dat.sumInt.bgSubtr[x,]
  }
}))

#write.csv(dat.sumInt.bgSubtr2, file = "Exp. 13-1_dat.sumInt.bgSubtr2_check.csv")
# there are NaN values here for some reason

length(which(is.nan(dat.sumInt.bgSubtr2))) # should always be 0
# Error in is.nan(dat.sumInt.bgSubtr2) : 
# default method not implemented for type 'list'

# CR troubleshooting for error on lone 165

dat.sumInt.bgSubtr2 <- t(sapply(rownames(dat.sumInt.bgSubtr), function(x) {
  # ensure we get a numeric vector, not a list
  row_vals <- as.numeric(dat.sumInt.bgSubtr[x, , drop = FALSE])
  
  # restore original row if there was NaN
  if (any(is.nan(row_vals))) {
    row_vals <- as.numeric(dat.sumInt[x, , drop = FALSE])
  }
  
  # replace any NaN with NA (preserve missing data)
  row_vals[is.nan(row_vals)] <- NA
  
  row_vals
}))

dat.sumInt.bgSubtr2 <- as.matrix(dat.sumInt.bgSubtr2)

# Verify
sum(is.nan(dat.sumInt.bgSubtr2))  # should now be 0
# [1] 0
sum(is.na(dat.sumInt.bgSubtr2))   # counts original NAs
# [1] 20970

row.names(dat.sumInt.bgSubtr2)
#    [1] "IGLV4-60|A0A075B6I1"    "C1orf232|A0A0U1RR37"    "TAF11L11|A0A1W2PQ09"    "MCTS2|A0A3B3IRV3"  

colnames(dat.sumInt.bgSubtr2) <- colnames(dat.sumInt.bgSubtr)
colnames(dat.sumInt.bgSubtr2)
# [1] "Input_1" "Input_2" "Input_3" "IP_1"    "IP_2"    "IP_3"    "IP_4"    "IP_5"    "IP_6"    "IP_7"    "IP_8"   
# [12] "IP_9" 

dim(dat.sumInt.bgSubtr2)
# [1] 8473   12

dat.sumInt.bgSubtr.log<-log2(dat.sumInt.bgSubtr2)
dat.sumInt.bgSubtr.log[!is.finite(dat.sumInt.bgSubtr.log)]<-NA  # handles log2(0)

# subset to remove neg.pulldowns
dat.sumInt.bgSubtr.nonegs.log <- dat.sumInt.bgSubtr.log[,which(!Grouping %in% myGroups.neg)]

dim(dat.sumInt.bgSubtr.nonegs.log)
# [1] 8473    9

colnames(dat.sumInt.bgSubtr.nonegs.log) <- c("Input_1", "Input_2", "Input_3",
                                            "IP_4", "IP_5", "IP_6", 
                                            "IP_7", "IP_8", "IP_9")

#################################################################################

# c. median samplewise zero centered normalization for pulldown vs global diffabund analysis
## Simple median samplewise subtraction will center the distributions all at 0.

class(dat.sumInt.bgSubtr.nonegs.log)
# [1] "matrix" "array" 

colnames(dat.sumInt.bgSubtr.nonegs.log)
# [1] "Input_1" "Input_2" "Input_3" "IP_4"    "IP_5"    "IP_6"    "IP_7"    "IP_8"    "IP_9"  

# Calculate median sample-wise
sampleMedianLFQ <- apply(dat.sumInt.bgSubtr.nonegs.log, 2, function(x) median(x, na.rm = TRUE))

# Calculate the mean of sample medians
meanMedian <- mean(sampleMedianLFQ)

# Apply median sample-wise subtraction
#dat.sumInt.log.meanMedian <-sapply(colnames(dat.sumInt.bgSubtr.nonegs.log),function(x) dat.sumInt.bgSubtr.nonegs.log[,x] -sampleMedianLFQ[x] )  # keeps values positive: -(sampleMedianLFQ[x]-meanMedian) )

# Subtract column medians from each column
dat.sumInt.log.meanMedian <- sweep(
  dat.sumInt.bgSubtr.nonegs.log, 2, sampleMedianLFQ, FUN = "-"
)

colnames(dat.sumInt.log.meanMedian)
# NULL

colnames(dat.sumInt.log.meanMedian) <- colnames(dat.sumInt.bgSubtr.nonegs.log)
rownames(dat.sumInt.log.meanMedian) <- rownames(dat.sumInt.bgSubtr.nonegs.log)

#################################################################################
# TurboID normalization

# 1. Extract the TurboID row from the log2 data
turboid_row <- dat.sumInt.bgSubtr.nonegs.log[turboid_row_index, , drop = FALSE]  # keeps as 1-row matrix
rownames(turboid_row) <- "|TurboID"

# 2. Identify positive pulldown columns
posPulldownCols <- which(Grouping_2 %in% myGroups.pos)
nonPosCols <- setdiff(seq_along(Grouping_2), posPulldownCols)

# 3. Extract TurboID values for positive pulldowns and ensure numeric
turboid_values_pos <- as.numeric(turboid_row[1, posPulldownCols])
turboid_values_pos
# [1] 16.66266 16.66829 16.85558 16.74465 16.77900 16.88203

# 4. Compute the mean TurboID signal across positive pulldowns
turboid_mean_pos <- mean(turboid_values_pos, na.rm = TRUE)
turboid_mean_pos
# [1] 16.76537

# 5. Compute per-column shift for positive pulldowns
turboid_shifts <- rep(0, ncol(dat.sumInt.bgSubtr.nonegs.log))  # initialize all zeros
turboid_shifts[posPulldownCols] <- turboid_mean_pos - turboid_values_pos

# 6. Apply shifts column-wise to the positive pulldown columns only
dat.sumInt.turboNorm.log2 <- dat.sumInt.bgSubtr.nonegs.log  # start with original
dat.sumInt.turboNorm.log2[, posPulldownCols] <- sweep(
  dat.sumInt.bgSubtr.nonegs.log[, posPulldownCols, drop = FALSE], 2,
  turboid_shifts[posPulldownCols], "+"
)

# 7. QC: Check that non-positive columns remain unchanged
column_shift_check <- colMeans(dat.sumInt.turboNorm.log2, na.rm = TRUE) -
  colMeans(dat.sumInt.bgSubtr.nonegs.log, na.rm = TRUE)

column_shift_check  # Positive pulldown columns should match turboid_shifts; others = 0
#     Input_1     Input_2     Input_3        IP_4        IP_5        IP_6        IP_7        IP_8        IP_9 
# 0.00000000  0.00000000  0.00000000  0.10271197  0.09708120 -0.09020991  0.02071996 -0.01363561 -0.11666762 

#################################################################################
## QC Plot - step 2 - Unnorm Sum Intensity vs. Enforced Equal Loading Assumption WITHIN TREATMENT GROUP (samplewise summed intensity normalization) -- not MaxLFQ norm
pdf(file="2. Exp. 13-1 QC-runorder-SummedIntensity.Vs.meanMedian zero-centered norm_08302025.pdf",width=21,height=14)
par(mar=c(12,5,3,3))
par(cex.axis=1.5)
par(cex.lab=1.5)
par(cex.main=3)

layout(matrix(c(1,2, 3,4, 5,6), nrow = 3, ncol = 2, byrow=TRUE),
       heights = c(0.95,0.95,0.95), # Heights of the rows
       widths = c(4.4,1)) # Widths of the columns

# colors for samples by Group.simple
colvec=WGCNA::labels2colors(as.numeric(factor(Grouping)))
colvec.legend=WGCNA::labels2colors(as.numeric(factor(levels(factor(Grouping)))))

# boxplot 1: a. no normalization
boxplot(dat.sumInt.log, ylab=bquote("log"[2]~"Spectronaut Summed Intensity"), main="a. Summed Intensity Protein Data (NO Normalization)", col=colvec, las=2)
par(mar=c(12,0,3,0))
plot.new()
legend("topleft", legend=levels(factor(Grouping)), fill=colvec.legend, cex=2.2)
par(mar=c(12,5,3,3))

# boxplot 2: b. neg.pulldown background subtraction
boxplot(dat.sumInt.bgSubtr.log[,which(Grouping %in% myGroups.pos)], ylab=bquote("log"[2]~"Spectronaut Summed Intensity"), main="b. Background Subtracted (Pulldowns Only)", col=colvec[which(Grouping %in% myGroups.pos)], las=2)
plot.new()

# adjust "Grouping" to remove neg.pulldowns
traits_nonegs <- numericMeta[which(!Grouping %in% myGroups.neg),]

Grouping_2<-traits_nonegs$Group.simple
Grouping_2
# [1] "Global"      "Global"      "Global"      "SR64_pulldn" "SR64_pulldn" "SR64_pulldn" "SR71_pulldn" "SR71_pulldn"
# [9] "SR71_pulldn"

# colors for samples by Group.simple
#colvec_2=WGCNA::labels2colors(as.numeric(factor(Grouping_2)))
#colvec.legend=WGCNA::labels2colors(as.numeric(factor(levels(factor(Grouping_2)))))
colvec_2=colvec[which(!Grouping %in% myGroups.neg)]

# boxplot 3. c. median samplewise zero centered
boxplot(dat.sumInt.log.meanMedian, ylab=bquote("log"[2]~"Spectronaut Summed Intensity"), main="c. Median Samplewise (Zero-)Centered Protein Data", col=colvec_2, las=2)
plot.new()

# boxplot 4. d. TurboID abundance normalization
boxplot(dat.sumInt.turboNorm.log2, ylab=bquote("log"[2]~"Spectronaut Summed Intensity"), main="d. TurboID normalization (pulldowns only) Protein Data", col=colvec_2, las=2)
plot.new()

dev.off()

################################################################################
# 3. Determine which rows to keep to control missing values
#    (Enforce missingness across samples with logical criteria)
###############################################################################
###############################################################################
###############################################################################
# missingness threshold for diffabund analysis: keep a protein if present in 2/3 samples/pos.pulldown group

dat.sumInt.log.norm <- dat.sumInt.log.meanMedian
dim(dat.sumInt.log.norm)
# [1] 8473    9

# the input df does not have duplicate gene IDs collapsed
# the input df is dat.sumInt.log.meanMedian
rownames.toKeep<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="SR64_pulldn")]))) >= 2 |
     length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="SR71_pulldn")]))) >= 2) x }))

length(rownames.toKeep)
# [1] 7333

#dat.sumInt.log.norm <- dat.sumInt.turboNorm.log2
#dim(dat.sumInt.log.norm)
# [1] 8473    9

# the input df does not have duplicate gene IDs collapsed
# the input df is dat.sumInt.turboNorm.log2
#rownames.toKeep<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
#  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="SR64_pulldn")]))) >= 2 |
#     length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="SR71_pulldn")]))) >= 2) x }))

#length(rownames.toKeep)
# [1] 7333

# CR: because the TurboID abundance is similar across all positive pulldown samples,
# I decided to moved forward with the meanMedian zero-centered norm approach
# since we will be comparing pulldown vs global samples per TurboID enzyme type (i.e.. SR64 vs SR71)

# missingness threshold for Venn diagrams: keep a protein if present in 2/3 samples/group
dat.sumInt.log.norm <- dat.sumInt.log.meanMedian
dim(dat.sumInt.log.norm)
# [1] 8473    9

# the input df does not have duplicate gene IDs collapsed
# Global filter
rownames.toKeep_Global <-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="Global")]))) >= 2 ) x }))

length(rownames.toKeep_Global)
# [1] 8071

cleanDat_Global_filtered<-dat.sumInt.log.norm[rownames.toKeep_Global,]

write.csv(cleanDat_Global_filtered, file=paste0("Exp. 13-1_DIA LFQ-MS_HEK.global_Venn filtered-",dim(cleanDat_Global_filtered)[1],"x",dim(cleanDat_Global_filtered)[2],"_08302025.csv"))


# SR64.pulldn filter
rownames.toKeep_SR64_pulldn<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="SR64_pulldn")]))) >= 2 ) x }))

length(rownames.toKeep_SR64_pulldn)
# [1] 6278

cleanDat_SR64_pulldn_filtered <- dat.sumInt.log.norm[rownames.toKeep_SR64_pulldn,]

write.csv(cleanDat_SR64_pulldn_filtered, file=paste0("Exp. 13-1_DIA LFQ-MS_SR64_pulldn_Venn filtered-",dim(cleanDat_SR64_pulldn_filtered)[1],"x",dim(cleanDat_SR64_pulldn_filtered)[2],"_08302025.csv"))

# SR71.pulldn filter
rownames.toKeep_SR71_pulldn<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="SR71_pulldn")]))) >= 2 ) x }))

length(rownames.toKeep_SR71_pulldn)
# [1] 7265

cleanDat_SR71_pulldn_filtered <- dat.sumInt.log.norm[rownames.toKeep_SR71_pulldn,]

write.csv(cleanDat_SR71_pulldn_filtered, file=paste0("Exp. 13-1_DIA LFQ-MS_SR71_pulldn_Venn filtered-",dim(cleanDat_SR71_pulldn_filtered)[1],"x",dim(cleanDat_SR71_pulldn_filtered)[2],"_08302025.csv"))

# file names:
# Exp. 13-1_DIA LFQ-MS_HEK.global_Venn filtered-8071x9_08302025.csv
# Exp. 13-1_DIA LFQ-MS_SR64_pulldn_Venn filtered-6278x9_08302025.csv
# Exp. 13-1_DIA LFQ-MS_SR71_pulldn_Venn filtered-7265x9_08302025.csv

#Load intensity files
global_protein <- read.csv("Exp. 13-1_DIA LFQ-MS_HEK.global_Venn filtered-8071x9_08302025.csv", header = TRUE)
SR64_pulldn_protein <- read.csv("Exp. 13-1_DIA LFQ-MS_SR64_pulldn_Venn filtered-6278x9_08302025.csv", header = TRUE)
SR71_pulldn_protein <- read.csv("Exp. 13-1_DIA LFQ-MS_SR71_pulldn_Venn filtered-7265x9_08302025.csv", header = TRUE)


## Subset the intensity matrices to extract only the protein ID columns to create Venn diagrams
global_protein_IDs  <- data.frame(global_protein[,1])
colnames(global_protein_IDs )
colnames(global_protein_IDs ) <- c("gene_protein_IDs")
library(tidyr)
global_protein_IDs_sep <- separate(global_protein_IDs , gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

SR64_pulldn_protein_IDs  <- data.frame(SR64_pulldn_protein[,1])
colnames(SR64_pulldn_protein_IDs )
colnames(SR64_pulldn_protein_IDs ) <- c("gene_protein_IDs")
library(tidyr)
SR64_pulldn_protein_IDs_sep <- separate(SR64_pulldn_protein_IDs , gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

SR71_pulldn_protein_IDs  <- data.frame(SR71_pulldn_protein[,1])
colnames(SR71_pulldn_protein_IDs )
colnames(SR71_pulldn_protein_IDs ) <- c("gene_protein_IDs")
library(tidyr)
SR71_pulldn_protein_IDs_sep <- separate(SR71_pulldn_protein_IDs , gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

global_protein_list <- global_protein_IDs_sep [,2]
SR64_pulldn_protein_list <- SR64_pulldn_protein_IDs_sep [,2]
SR71_pulldn_protein_list <-  SR71_pulldn_protein_IDs_sep [,2]
list_z <- NULL


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
pdf(file = "Exp. 13-1_DIA LFQ-MS_protein_Venns_CR_08302025.pdf", height = 11, width = 8.5, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# total (x) vs LPS.total (y)
biovenn1 <- draw.venn(global_protein_list, SR64_pulldn_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Global vs SR64_pulldn protein", 
                      xtitle = "Global", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "SR64_pulldn",
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
                      y_c = "darkorange3",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x total: 8071"
# [1] "y total: 6278"
# [1] "z total: 0"
# [1] "x only: 1985"
# [1] "y only: 192"
# [1] "z only: 0"
# [1] "x-y total overlap: 6086"


# total (x) vs BV2T pulldown (y)
biovenn2 <- draw.venn(global_protein_list, SR71_pulldn_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Global vs SR71_pulldn protein", 
                      xtitle = "Globall", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "SR71_pulldn",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "brown3",
                      y_c = "azure",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)
# [1] "x total: 8071"
# [1] "y total: 7265"
# [1] "z total: 0"
# [1] "x only: 1051"
# [1] "y only: 245"
# [1] "z only: 0"
# [1] "x-y total overlap: 7020"

# BV2T.pos.pulldown (x) vs SR64_pulldn (y)
biovenn3 <- draw.venn(SR64_pulldn_protein_list, SR71_pulldn_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="SR64 vs SR71 pulldown protein", 
                      xtitle = "SR64_pulldn", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "SR71_pulldn",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "azure3",
                      y_c = "brown3",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# 1] "x total: 6278"
# [1] "y total: 7265"
# [1] "z total: 0"
# [1] "x only: 68"
# [1] "y only: 1055"
# [1] "z only: 0"
# [1] "x-y total overlap: 6210"

dev.off()




##################################################################################
# 4. Impute missing values so we are able to do better-powered statistical analyses
###################################################################################
# Perseus Style imputation in R -- by Eric Dammer  ################################
#
# Assumes log2(LFQ abundance) normal distribution, and noise level -1.8SD from mean
# imputes according to a normal distribution +/-0.3SD from noise level
#
# -includes intelligent background IgG subtraction for IP Intensity data
###################################################################################
library(stats)
library("doParallel")
library("snow")
parallelThreads=7
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)

## Define imputation function
imputePerseusStyle <- function(cleanDat.unfiltered.withNA, ...) {
  for (i in 1:ncol(cleanDat.unfiltered.withNA)) { cleanDat.unfiltered.withNA[,i]<-as.numeric(cleanDat.unfiltered.withNA[,i]) }
  # set 0 values (missing data from macquant non-log-transformed data to 0)
  cleanDat.unfiltered.withNA[cleanDat.unfiltered.withNA==0]<-NA
  #cleanDat.unfiltered.withNA<-log2(cleanDat.unfiltered.withNA) #will also populate NAs from 0 values 
  
  #Create a matrix with NA positions from cleanDat.unfiltered.withNA (you can multiply this matrix by imputed exprMat to replace imputed values if desired later)
  NAmatrix<-cleanDat.unfiltered.withNA
  NAmatrix[!is.na(cleanDat.unfiltered.withNA)] <- 1
  
  #Calculate imputation parameters by column
  masterAvg=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  masterSD=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  downshift=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  noisePlusMinus=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  noiseAvg=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  for (n in 1:ncol(cleanDat.unfiltered.withNA)) {
    masterAvg[n] <- mean(cleanDat.unfiltered.withNA[!is.na(cleanDat.unfiltered.withNA[,n]),n])
    masterSD[n] <- sd(cleanDat.unfiltered.withNA[!is.na(cleanDat.unfiltered.withNA[,n]),n])
    downshift[n] <- 1.8*masterSD[n]
    noisePlusMinus[n] <- 0.3*masterSD[n]
    noiseAvg[n] <- masterAvg[n]-downshift[n]
  }
  
  #create a row of imputation parameters for each column (sample) we will be imputing NAs in cleanDat.unfiltered.withNA
  rbind(masterAvg,masterSD,downshift,noisePlusMinus,noiseAvg)
  
  #Impute Noise Level Signals, Perseus Style, Column-wise, each column considered Independently - all columns!
  exprMatTemp <- cleanDat.unfiltered.withNA
  exprMatImp <- foreach (n=1:ncol(cleanDat.unfiltered.withNA), .combine=cbind) %dopar% {
    avgVec <- cleanDat.unfiltered.withNA[,n] #rowMeans(cleanDat.unfiltered.withNA,na.rm=TRUE)
    seed=0
    set.seed(seed+3)
    NAvec<-which(is.na(exprMatTemp[,n]))
    noise<-noiseAvg[n]
    # random imputed values falling within the left-tail distribution
    randVec<-runif(sum(is.na(exprMatTemp[,n]))+1,noise-noisePlusMinus[n],noise+noisePlusMinus[n])
    # sample from this random distribution
    ImputeVec <- exprMatTemp[NAvec,n] <- sample(randVec,sum(is.na(cleanDat.unfiltered.withNA[,n])),replace=FALSE)
    exprMatTemp[,n]
  }
  colnames(exprMatImp)<-colnames(cleanDat.unfiltered.withNA)
  rownames(exprMatImp)<-rownames(cleanDat.unfiltered.withNA)
  return(exprMatImp)
}

# a. non-background and no normalization
cleanDat.unfiltered.imputed<-imputePerseusStyle(dat.sumInt.log) #full data with imputation Perseus Style Column-wise Independent - all columns
head(cleanDat.unfiltered.imputed)

# write output with imputed values (log2-transformed), cleanDat.unfiltered.imputed
write.csv(cleanDat.unfiltered.imputed, file=paste0("4a.cleanDat.nobackgsub.nonorm.unfiltered.imputed-",dim(cleanDat.unfiltered.imputed)[1],"x",dim(cleanDat.unfiltered.imputed)[2],".csv_Perseus-Imputation_083020254.csv"))

dim(cleanDat.unfiltered.imputed)
# [1] 8473   12

# now includes imputed data points -- and we have just applied the filtering to keep only rows meeting criteria from step 3
cleanDat<-cleanDat.filtered.imputed<-cleanDat.unfiltered.imputed[rownames.toKeep,]
write.csv(cleanDat.filtered.imputed, file=paste0("4a.cleanDat.nobackgsub.nonorm.filtered.imputed-",dim(cleanDat.filtered.imputed)[1],"x",dim(cleanDat.filtered.imputed)[2],".csv_Perseus-Imputation_083020254.csv"))

dim(cleanDat)
# [1] 7333   12

# b. negative pulldown background subtraction
cleanDat.bgSubtr.unfiltered.imputed<-imputePerseusStyle(dat.sumInt.bgSubtr.log) #full data with imputation Perseus Style Column-wise Independent - all columns
head(cleanDat.bgSubtr.unfiltered.imputed)

# write output with imputed values (log2-transformed), cleanDat.unfiltered.imputed
write.csv(cleanDat.bgSubtr.unfiltered.imputed, file=paste0("4b.cleanDat.bgSubtr.unfiltered.imputed-",dim(cleanDat.bgSubtr.unfiltered.imputed)[1],"x",dim(cleanDat.bgSubtr.unfiltered.imputed)[2],".csv_Perseus-Imputation_08302025.csv"))

dim(cleanDat.bgSubtr.unfiltered.imputed)
# [1] 8473   12

# now includes imputed data points -- and we have just applied the filtering to keep only rows meeting criteria from step 3
cleanDat<-cleanDat.bgSubtr.filtered.imputed<-cleanDat.bgSubtr.unfiltered.imputed[rownames.toKeep,]
write.csv(cleanDat.bgSubtr.filtered.imputed, file=paste0("4b.cleanDat.bgSubtr.filtered.imputed-",dim(cleanDat.bgSubtr.filtered.imputed)[1],"x",dim(cleanDat.bgSubtr.filtered.imputed)[2],".csv_Perseus-Imputation_08302025.csv"))

dim(cleanDat)
# [1] 7333   12

# c. meanMedian zero-centered norm
cleanDat.bgSubtr.0centered.unfiltered.imputed<-imputePerseusStyle(dat.sumInt.log.meanMedian) #full data with imputation Perseus Style Column-wise Independent - all columns
head(cleanDat.bgSubtr.0centered.unfiltered.imputed)

# write output with imputed values (log2-transformed), cleanDat.unfiltered.imputed
write.csv(cleanDat.bgSubtr.0centered.unfiltered.imputed, file=paste0("4c.cleanDat.bgSubtr.0centered.unfiltered.imputed-",dim(cleanDat.bgSubtr.0centered.unfiltered.imputed)[1],"x",dim(cleanDat.bgSubtr.0centered.unfiltered.imputed)[2],".csv_Perseus-Imputation_08302025.csv"))

dim(cleanDat.bgSubtr.0centered.unfiltered.imputed)
# [1] 8473   9

# now includes imputed data points -- and we have just applied the filtering to keep only rows meeting criteria from step 3
cleanDat<-cleanDat.bgSubtr.0centered.filtered.imputed<-cleanDat.bgSubtr.0centered.unfiltered.imputed[rownames.toKeep,]
write.csv(cleanDat.bgSubtr.0centered.filtered.imputed, file=paste0("4c.cleanDat.bgSubtr.0centered.filtered.imputed-",dim(cleanDat.bgSubtr.0centered.filtered.imputed)[1],"x",dim(cleanDat.bgSubtr.0centered.filtered.imputed)[2],".csv_Perseus-Imputation_08302025.csv"))

dim(cleanDat)
# [1] 7333   9

# d. TurboID norm
cleanDat.turboNorm.unfiltered.imputed<-imputePerseusStyle(dat.sumInt.turboNorm.log2) #full data with imputation Perseus Style Column-wise Independent - all columns
head(cleanDat.turboNorm.unfiltered.imputed)

# write output with imputed values (log2-transformed), cleanDat.unfiltered.imputed
write.csv(cleanDat.turboNorm.unfiltered.imputed, file=paste0("4d.cleanDat.TurboNorm.unfiltered.imputed-",dim(cleanDat.turboNorm.unfiltered.imputed)[1],"x",dim(cleanDat.turboNorm.unfiltered.imputed)[2],".csv_Perseus-Imputation_08302025.csv"))

dim(cleanDat.turboNorm.unfiltered.imputed)
# [1] 8473    9

# now includes imputed data points -- and we have just applied the filtering to keep only rows meeting criteria from step 3
cleanDat<-cleanDat.turboNorm.filtered.imputed<-cleanDat.turboNorm.unfiltered.imputed[rownames.toKeep,]
write.csv(cleanDat.turboNorm.filtered.imputed, file=paste0("4d.cleanDat.TurboNorm.filtered.imputed-",dim(cleanDat.turboNorm.filtered.imputed)[1],"x",dim(cleanDat.turboNorm.filtered.imputed)[2],".csv_Perseus-Imputation_08302025.csv"))

dim(cleanDat)
# [1] 7333    9

################################################################################
################################################################################
################################################################################
## Part 1: Data processing and clean up
# steps:
# 3. Collapse duplicate gene IDs 
# a. Keep a protein isoform with the most variance
###############################################################################
###############################################################################
# a. normalized, filtered and imputed intensity matrix with total and pos.pulldowns as input
# same as BV2 data e.g., 5a. cleanDat.bgSubtr.filtered.imputed.collapsed-1436x12.csv_Perseus-Imputation_02132024.csv
cleanDat <- cleanDat.bgSubtr.0centered.filtered.imputed
dim(cleanDat)
# [1] 7333    9

# check to see if any gene symbols are missing, should be 0
badRows<-which(grepl("^\\|",rownames(cleanDat)))
rownames(cleanDat)[badRows]
# [1] "|Q6ZSR9"     "|TurboID"    "|A0A0J9YY99" "|A0AAG2SWI5"

# add missing gene symbols
# "AAK1_predicted|Q6ZSR9"
# "TurboID|TurboID"
# "Ig-like domain-containing protein|A0A0J9YY99"
# "Ovostatin homolog 2-like|A0AAG2SWI5"

# Your lookup table of missing mappings
fix_map <- c(
  "|Q6ZSR9"     = "AAK1_predicted|Q6ZSR9",
  "|TurboID"    = "TurboID|TurboID",
  "|A0A0J9YY99" = "Ig-like domain-containing protein|A0A0J9YY99",
  "|A0AAG2SWI5" = "Ovostatin homolog 2-like|A0AAG2SWI5"
)

# Apply replacements
rownames(cleanDat) <- ifelse(rownames(cleanDat) %in% names(fix_map),
                             fix_map[rownames(cleanDat)],
                             rownames(cleanDat))

# Verify replacements
rownames(cleanDat)[badRows]
# [1] "AAK1_predicted|Q6ZSR9"                        "TurboID|TurboID"                             
# [3] "Ig-like domain-containing protein|A0A0J9YY99" "Ovostatin homolog 2-like|A0AAG2SWI5"  

# The norm intensity matrix row names are labeled gene symbol|protein ID, so we
# need to extract the first part of each row name to obtain gene symbols and put into a df
symbols <- as.data.frame(do.call(rbind, strsplit(rownames(cleanDat), "[|]")))[, 1]
symbols <- as.data.frame(do.call(rbind, strsplit(symbols, "[;]")))[, 1]

# Load library
library(WGCNA)

# Collapse rows based on variance
collapsed_data <- collapseRows(cleanDat, symbols, rownames(cleanDat), method = "MaxMean")

cleanDat.collapsed<-collapsed_data$datETcollapsed
dim(cleanDat.collapsed)
# [1] 7332    9

head(cleanDat.collapsed)
#          Input_1    Input_2     Input_3       IP_4       IP_5       IP_6       IP_7       IP_8       IP_9
# A2M    0.2241254 -0.2219995  0.67672712 -0.2398390 -0.1565929 -0.2104770  1.4668954  1.3415521  1.7029624
# A2ML1 -2.0515504 -1.6933404 -1.83892097  1.6800395  1.9649501  1.5768954  1.3397025  0.7151134  1.9992630


write.csv(cleanDat.collapsed, file=paste0("5a. cleanDat.bgSubtr.0center.filtered.imputed.collapsed-",dim(cleanDat.collapsed)[1],"x",dim(cleanDat.collapsed)[2],".csv_Perseus-Imputation_08302025.csv"))

# end of data processing and clean up pipeline

#################################################################################
# dat.sumInt.turboNorm.log2 as input for collapsed df for concord/discord analysis

dat.sumInt.log.norm <- dat.sumInt.turboNorm.log2
dim(dat.sumInt.log.norm)
# [1] 8473    9

# the input df does not have duplicate gene IDs collapsed
# the input df is dat.sumInt.turboNorm.log2
rownames.toKeep<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="SR64_pulldn")]))) >= 2 |
     length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="SR71_pulldn")]))) >= 2) x }))

length(rownames.toKeep)
# [1] 7333

##################################################################################
# 4. Impute missing values so we are able to do better-powered statistical analyses
###################################################################################
# Perseus Style imputation in R -- by Eric Dammer  ################################
#
# Assumes log2(LFQ abundance) normal distribution, and noise level -1.8SD from mean
# imputes according to a normal distribution +/-0.3SD from noise level
#
# -includes intelligent background IgG subtraction for IP Intensity data
###################################################################################
library(stats)
library("doParallel")
library("snow")
parallelThreads=7
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)

## Define imputation function
imputePerseusStyle <- function(cleanDat.unfiltered.withNA, ...) {
  for (i in 1:ncol(cleanDat.unfiltered.withNA)) { cleanDat.unfiltered.withNA[,i]<-as.numeric(cleanDat.unfiltered.withNA[,i]) }
  # set 0 values (missing data from macquant non-log-transformed data to 0)
  cleanDat.unfiltered.withNA[cleanDat.unfiltered.withNA==0]<-NA
  #cleanDat.unfiltered.withNA<-log2(cleanDat.unfiltered.withNA) #will also populate NAs from 0 values 
  
  #Create a matrix with NA positions from cleanDat.unfiltered.withNA (you can multiply this matrix by imputed exprMat to replace imputed values if desired later)
  NAmatrix<-cleanDat.unfiltered.withNA
  NAmatrix[!is.na(cleanDat.unfiltered.withNA)] <- 1
  
  #Calculate imputation parameters by column
  masterAvg=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  masterSD=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  downshift=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  noisePlusMinus=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  noiseAvg=as.vector(rep(0,ncol(cleanDat.unfiltered.withNA)))
  for (n in 1:ncol(cleanDat.unfiltered.withNA)) {
    masterAvg[n] <- mean(cleanDat.unfiltered.withNA[!is.na(cleanDat.unfiltered.withNA[,n]),n])
    masterSD[n] <- sd(cleanDat.unfiltered.withNA[!is.na(cleanDat.unfiltered.withNA[,n]),n])
    downshift[n] <- 1.8*masterSD[n]
    noisePlusMinus[n] <- 0.3*masterSD[n]
    noiseAvg[n] <- masterAvg[n]-downshift[n]
  }
  
  #create a row of imputation parameters for each column (sample) we will be imputing NAs in cleanDat.unfiltered.withNA
  rbind(masterAvg,masterSD,downshift,noisePlusMinus,noiseAvg)
  
  #Impute Noise Level Signals, Perseus Style, Column-wise, each column considered Independently - all columns!
  exprMatTemp <- cleanDat.unfiltered.withNA
  exprMatImp <- foreach (n=1:ncol(cleanDat.unfiltered.withNA), .combine=cbind) %dopar% {
    avgVec <- cleanDat.unfiltered.withNA[,n] #rowMeans(cleanDat.unfiltered.withNA,na.rm=TRUE)
    seed=0
    set.seed(seed+3)
    NAvec<-which(is.na(exprMatTemp[,n]))
    noise<-noiseAvg[n]
    # random imputed values falling within the left-tail distribution
    randVec<-runif(sum(is.na(exprMatTemp[,n]))+1,noise-noisePlusMinus[n],noise+noisePlusMinus[n])
    # sample from this random distribution
    ImputeVec <- exprMatTemp[NAvec,n] <- sample(randVec,sum(is.na(cleanDat.unfiltered.withNA[,n])),replace=FALSE)
    exprMatTemp[,n]
  }
  colnames(exprMatImp)<-colnames(cleanDat.unfiltered.withNA)
  rownames(exprMatImp)<-rownames(cleanDat.unfiltered.withNA)
  return(exprMatImp)
}

# d. TurboID norm
cleanDat.bgSubtr.turboNorm.unfiltered.imputed<-imputePerseusStyle(dat.sumInt.turboNorm.log2) #full data with imputation Perseus Style Column-wise Independent - all columns
head(cleanDat.bgSubtr.turboNorm.unfiltered.imputed)

# write output with imputed values (log2-transformed), cleanDat.unfiltered.imputed
write.csv(cleanDat.bgSubtr.turboNorm.unfiltered.imputed, file=paste0("4d.cleanDat.bgSubtr.turboNorm.unfiltered.imputed-",dim(cleanDat.bgSubtr.turboNorm.unfiltered.imputed)[1],"x",dim(cleanDat.bgSubtr.turboNorm.unfiltered.imputed)[2],"_Perseus-Imputation_08302025.csv"))

dim(cleanDat.bgSubtr.turboNorm.unfiltered.imputed)
# [1] 8473    9

# now includes imputed data points -- and we have just applied the filtering to keep only rows meeting criteria from step 3
cleanDat<-cleanDat.bgSubtr.turboNorm.filtered.imputed<-cleanDat.bgSubtr.turboNorm.unfiltered.imputed[rownames.toKeep,]
write.csv(cleanDat.bgSubtr.turboNorm.filtered.imputed, file=paste0("4d.cleanDat.bgSubtr.turboNorm.filtered.imputed-",dim(cleanDat.bgSubtr.turboNorm.filtered.imputed)[1],"x",dim(cleanDat.bgSubtr.turboNorm.filtered.imputed)[2],"_Perseus-Imputation_08302025.csv"))

dim(cleanDat)
# [1] 7333    9

################################################################################
################################################################################
################################################################################
## Part 1: Data processing and clean up
# steps:
# 3. Collapse duplicate gene IDs 
# a. Keep a protein isoform with the most variance
###############################################################################
###############################################################################
# a. normalized, filtered and imputed intensity matrix with total and pos.pulldowns as input
# same as BV2 data e.g., 5a. cleanDat.bgSubtr.filtered.imputed.collapsed-1436x12.csv_Perseus-Imputation_02132024.csv
#cleanDat <- cleanDat.bgSubtr.turboNorm.filtered.imputed
#dim(cleanDat)
# [1] 7333    9

# check to see if any gene symbols are missing, should be 0
badRows<-which(grepl("^\\|",rownames(cleanDat)))
rownames(cleanDat)[badRows]
# [1] "|Q6ZSR9"     "|TurboID"    "|A0A0J9YY99"

# add missing gene symbols
# "AAK1_predicted|Q6ZSR9"
# "TurboID|TurboID"
# "Ig-like domain-containing protein|A0A0J9YY99"
# "Ovostatin homolog 2-like|A0AAG2SWI5"

# Your lookup table of missing mappings
fix_map <- c(
  "|Q6ZSR9"     = "AAK1_predicted|Q6ZSR9",
  "|TurboID"    = "TurboID|TurboID",
  "|A0A0J9YY99" = "Ig-like domain-containing protein|A0A0J9YY99"
)

# Apply replacements
rownames(cleanDat) <- ifelse(rownames(cleanDat) %in% names(fix_map),
                             fix_map[rownames(cleanDat)],
                             rownames(cleanDat))

# Verify replacements
rownames(cleanDat)[badRows]
# [1] "AAK1_predicted|Q6ZSR9"                        "TurboID|TurboID"                             
# [3] "Ig-like domain-containing protein|A0A0J9YY99"

# The norm intensity matrix row names are labeled gene symbol|protein ID, so we
# need to extract the first part of each row name to obtain gene symbols and put into a df
symbols <- as.data.frame(do.call(rbind, strsplit(rownames(cleanDat), "[|]")))[, 1]
symbols <- as.data.frame(do.call(rbind, strsplit(symbols, "[;]")))[, 1]

# Load library
library(WGCNA)

# Collapse rows based on variance
collapsed_data <- collapseRows(cleanDat, symbols, rownames(cleanDat), method = "MaxMean")

cleanDat.collapsed<-collapsed_data$datETcollapsed
dim(cleanDat.collapsed)
# [1] 7332    9

head(cleanDat.collapsed)
#        Input_1  Input_2  Input_3     IP_4     IP_5     IP_6     IP_7     IP_8     IP_9
# A2M   15.89629 15.29132 15.95043 12.81232 13.12448 12.92067 15.19241 15.07876 15.40791
# A2ML1 13.62061 13.81998 13.43478 14.73220 15.24602 14.70804 15.06522 14.45232 15.70421

write.csv(cleanDat.collapsed, file=paste0("5d. cleanDat.bgSubtr.turboNorm.filtered.imputed.collapsed-",dim(cleanDat.collapsed)[1],"x",dim(cleanDat.collapsed)[2],"_Perseus-Imputation_09012025.csv"))

# end of data processing and clean up pipeline

#################################################################################
#################################################################################
#################################################################################

## Part II: Data analysis and visualization
# 1. PCAs
# 2. log2(intensity) correlations
# 3. Volcanoes
# 4. GOparallel GSEA on volcano-defined DEP lists 

################################################################################
################################################################################
################################################################################
## Part 2: Data analysis and visualization
# 1. PCAs of uncollapsed data
# change the cleanDat to be the normalized intensity matrix
cleanDat <- dat.sumInt.log.norm
dim(cleanDat)
# [1] 8473    9

head(cleanDat)
#                        Input_1    Input_2    Input_3       IP_4      IP_5       IP_6       IP_7       IP_8       IP_9
# IGLV4-60|A0A075B6I1 -1.7599127 -0.7577473  0.4910774         NA        NA         NA -2.7413206         NA         NA
# C1orf232|A0A0U1RR37 -3.1078571 -3.3529860 -3.4335576         NA        NA         NA         NA         NA         NA
# TAF11L11|A0A1W2PQ09 -2.2083268 -1.7193634 -2.3181059         NA        NA         NA -0.9435292 -1.7597154 -0.8521431

# load library
library(limma)

pdf(file = "1a. Exp. 13-1_DIA LFQ-MS_all groups_nonegs_uncollapsed_Limma_PCA_CR_08302025.pdf", width = 8, height = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- numericMeta_nonegs$Group.simple
# [1] "Global"      "Global"      "Global"      "SR64_pulldn" "SR64_pulldn" "SR64_pulldn" "SR71_pulldn" "SR71_pulldn"
# [9] "SR71_pulldn"

pch_values <- ifelse(Grouping_vector == "Global", 16,
                     ifelse(Grouping_vector == "SR64_pulldn", 16,
                            ifelse(Grouping_vector == "SR71_pulldn", 16, NA)))


pt_colors <- ifelse(Grouping_vector == "Global", "green3",
                    ifelse(Grouping_vector == "SR64_pulldn", "orange",
                           ifelse(Grouping_vector == "SR71_pulldn", "dodgerblue",NA)))

legend_groups <- c("Global", "SR64_pulldn", "SR71_pulldn")

plotMDS_allgroups_nonegs <- plotMDS((cleanDat), #top = 500, 
                                    labels = NULL, pch = pch_values, col = pt_colors, 
                                    cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                    xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for uncollapsed log2(intensity)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("green3", "orange", "dodgerblue"), 
       pch = 16, title = "Proteome groups",cex=1.4)

################################################################################
pdf(file = "1b. Exp. 13-1_DIA LFQ-MS_SR64 vs SR71 pulldn_uncollapsed_Limma_PCA_CR_08302025.pdf", width = 8, height = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

pch_values <- ifelse(Grouping_vector == "SR64_pulldn", 16,
                     ifelse(Grouping_vector == "SR71_pulldn", 16, NA))


pt_colors <- ifelse(Grouping_vector == "SR64_pulldn", "orange",
                    ifelse(Grouping_vector == "SR71_pulldn", "dodgerblue",NA))

legend_groups <- c("SR64_pulldn", "SR71_pulldn")

plotMDS_SR64VsSR71 <- plotMDS((cleanDat), #top = 500, 
                              labels = NULL, pch = pch_values, col = pt_colors, 
                              cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                              xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for uncollapsed log2(intensity)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("orange", "dodgerblue"), 
       pch = 16, title = "Proteome groups",cex=1.4)

dev.off()

################################################################################

################################################################################
## Part 2: Data analysis and visualization
# 2. log2(intensity) correlations of uncollapsed df

# Step 2. Correlations of protein abundances
# Comparisons
# a. HEK.SR64.pulldn vs HEK.global
# b. HEK.SR71.global vs HEK.global
# c. HEK.SR64.pulldn vs HEK.SR71.pulldn

# load library to generate plot
library(ggplot2)

# change the cleanDat to be the normalized, filtered and imputed intensity matrix
cleanDat <- read.csv("4c.cleanDat.bgSubtr.0centered.filtered.imputed-7333x9.csv_Perseus-Imputation_08302025.csv", header = TRUE, row.names = 1)
dim(cleanDat)
# [1] 7333    9

################################################################################
# a. HEK.SR64.pulldn vs HEK.global

pdf(file = "2a-c. Exp. 13-1_LFQ-MS_intensity correlations_uncollapsed_0centered_CR_08302025.pdf", width = 6, height = 6)
par(mfrow= c(2,1)) # centers plot on single page (I think)

# subset into groups of interest and take the row mean
global <- rowMeans(cleanDat[,c(1:3)])
length(global)
# [1] 7333

SR64_pulldn <- rowMeans(cleanDat[,c(4:6)])
length(SR64_pulldn)
# [1] 7333

# generate a linear correlation
correlation_2a <- cor(global, SR64_pulldn)
correlation_2a
# [1] 0.5808916

correlation_2a_test <- cor.test(global, SR64_pulldn, method = "pearson")
correlation_2a_test
# 	Pearson's product-moment correlation

# data:  global and SR64_pulldn
# t = 61.102, df = 7331, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.5655178 0.5958538
# sample estimates:
#   cor 
# 0.5808875 

# convert data to a data frame
df_2a <- data.frame(global = global, SR64_pulldn = SR64_pulldn)

# create the correlation plot
correlation_plot_2a <- ggplot(df_2a, aes(x = global, y = SR64_pulldn)) +
  geom_point(color = "black", fill = "darkorange", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "HEK.global log2(intensity)",
       y = "HEK.SR64_pulldn log2(intensity)",
       title = "Linear correlation for uncollapsed intensities_0centerNorm",
       caption = paste("Correlation coefficient:", round(correlation_2a, 2),
                       "\nNumber of proteins:", nrow(df_2a))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 1, vjust = 1))  # Align the caption to the bottom right
#) +
#scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
#scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

# Print the plot
print(correlation_plot_2a)

################################################################################
# b. HEK.SR71.global vs HEK.global

# subset into groups of interest and take the row mean
global <- rowMeans(cleanDat[,c(1:3)])
length(global)
# [1] 7333

SR71_pulldn <- rowMeans(cleanDat[,c(7:9)])
length(SR71_pulldn)
# [1] 7333

# generate a linear correlation
correlation_2b <- cor(global, SR71_pulldn)
correlation_2b
# [1] 0.6821919

correlation_2b_test <- cor.test(global, SR71_pulldn, method = "pearson")
correlation_2b_test
# 	Pearson's product-moment correlation

# data:  global and SR71_pulldn
# t = 79.795, df = 7331, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.6693357 0.6938406
# sample estimates:
#   cor 
# 0.6817794 

# convert data to a data frame
df_2b <- data.frame(global = global, SR71_pulldn = SR71_pulldn)

# create the correlation plot
correlation_plot_2b <- ggplot(df_2b, aes(x = global, y = SR71_pulldn)) +
  geom_point(color = "black", fill = "darkorange3", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "HEK.global log2(intensity)",
       y = "HEK.SR71_pulldn log2(intensity)",
       title = "Linear correlation for uncollapsed intensities_0centerNorm",
       caption = paste("Correlation coefficient:", round(correlation_2b, 2),
                       "\nNumber of proteins:", nrow(df_2b))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 1, vjust = 1))  # Align the caption to the bottom right
#) +
#scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
#scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

# Print the plot
print(correlation_plot_2b)

################################################################################
# c. HEK.SR71.pulldn vs HEK.SR71.pulldn

# subset into groups of interest and take the row mean
SR64_pulldn <- rowMeans(cleanDat[,c(4:6)])
length(SR64_pulldn)
# [1] 7333

SR71_pulldn <- rowMeans(cleanDat[,c(7:9)])
length(SR71_pulldn)
# [1] 7333

# generate a linear correlation
correlation_2c <- cor(SR64_pulldn, SR71_pulldn)
correlation_2c
# [1] 0.9264055

correlation_2c_test <- cor.test(SR64_pulldn, SR71_pulldn, method = "pearson")
correlation_2c_test
# 	Pearson's product-moment correlation

# data:  SR64_pulldn and SR71_pulldn
# t = 210.76, df = 7331, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9231538 0.9296415
# sample estimates:
#  cor 
# 0.9264664 

# convert data to a data frame
df_2c <- data.frame(SR64_pulldn = SR64_pulldn, SR71_pulldn = SR71_pulldn)

# create the correlation plot
correlation_plot_2c <- ggplot(df_2c, aes(x = SR64_pulldn, y = SR71_pulldn)) +
  geom_point(color = "black", fill = "darkorange3", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "HEK.SR64_pulldn log2(intensity)",
       y = "HEK.SR71_pulldn log2(intensity)",
       title = "Linear correlation for uncollapsed intensities_0centerNorm",
       caption = paste("Correlation coefficient:", round(correlation_2c, 2),
                       "\nNumber of proteins:", nrow(df_2c))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 1, vjust = 1))  # Align the caption to the bottom right
#) +
#scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
#scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

# Print the plot
print(correlation_plot_2c)

dev.off()
################################################################################
################################################################################
# corr plots with TurboID norm df as input
# change the cleanDat to be the normalized, filtered and imputed intensity matrix
cleanDat <- read.csv("4d.cleanDat.TurboNorm.filtered.imputed-7333x9.csv_Perseus-Imputation_08302025.csv", header = TRUE, row.names = 1)
dim(cleanDat)
# [1] 7333    9

################################################################################
# a. HEK.SR64.pulldn vs HEK.global

pdf(file = "2a-c. Exp. 13-1_LFQ-MS_intensity correlations_uncollapsed_TurboNorm_CR_08302025.pdf", width = 6, height = 6)
par(mfrow= c(2,1)) # centers plot on single page (I think)

# subset into groups of interest and take the row mean
global <- rowMeans(cleanDat[,c(1:3)])
length(global)
# [1] 7333

SR64_pulldn <- rowMeans(cleanDat[,c(4:6)])
length(SR64_pulldn)
# [1] 7333

# generate a linear correlation
correlation_2a_Turbo <- cor(global, SR64_pulldn)
# [1] 0.5808875

# convert data to a data frame
df_2a_Turbo <- data.frame(global = global, SR64_pulldn = SR64_pulldn)

# create the correlation plot
correlation_plot_2a_Turbo <- ggplot(df_2a_Turbo, aes(x = global, y = SR64_pulldn)) +
  geom_point(color = "black", fill = "brown1", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "HEK.global log2(intensity)",
       y = "HEK.SR64_pulldn log2(intensity)",
       title = "Linear correlation for uncollapsed intensities_TurboNorm",
       caption = paste("Correlation coefficient:", round(correlation_2a_Turbo, 2),
                       "\nNumber of proteins:", nrow(df_2a_Turbo))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 1, vjust = 1))  # Align the caption to the bottom right
#) +
#scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
#scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

# Print the plot
print(correlation_plot_2a_Turbo)

################################################################################
# b. HEK.SR71.global vs HEK.global

# subset into groups of interest and take the row mean
global <- rowMeans(cleanDat[,c(1:3)])
length(global)
# [1] 7333

SR71_pulldn <- rowMeans(cleanDat[,c(7:9)])
length(SR71_pulldn)
# [1] 7333

# generate a linear correlation
correlation_2b_Turbo <- cor(global, SR71_pulldn)
# [1] 0.6817794

# convert data to a data frame
df_2b_Turbo <- data.frame(global = global, SR71_pulldn = SR71_pulldn)

# create the correlation plot
correlation_plot_2b_Turbo <- ggplot(df_2b_Turbo, aes(x = global, y = SR71_pulldn)) +
  geom_point(color = "black", fill = "brown3", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "HEK.global log2(intensity)",
       y = "HEK.SR71_pulldn log2(intensity)",
       title = "Linear correlation for uncollapsed intensities",
       caption = paste("Correlation coefficient:", round(correlation_2b_Turbo, 2),
                       "\nNumber of proteins:", nrow(df_2b_Turbo))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 1, vjust = 1))  # Align the caption to the bottom right
#) +
#scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
#scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

# Print the plot
print(correlation_plot_2b_Turbo)

################################################################################
# c. HEK.SR71.pulldn vs HEK.SR71.pulldn

# subset into groups of interest and take the row mean
SR64_pulldn <- rowMeans(cleanDat[,c(4:6)])
length(SR64_pulldn)
# [1] 7333

SR71_pulldn <- rowMeans(cleanDat[,c(7:9)])
length(SR71_pulldn)
# [1] 7333

# generate a linear correlation
correlation_2c_Turbo <- cor(SR64_pulldn, SR71_pulldn)
# [1] 0.9264664

# convert data to a data frame
df_2c_Turbo <- data.frame(SR64_pulldn = SR64_pulldn, SR71_pulldn = SR71_pulldn)

# create the correlation plot
correlation_plot_2c_Turbo <- ggplot(df_2c_Turbo, aes(x = SR64_pulldn, y = SR71_pulldn)) +
  geom_point(color = "black", fill = "brown", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "HEK.SR64_pulldn log2(intensity)",
       y = "HEK.SR71_pulldn log2(intensity)",
       title = "Linear correlation for uncollapsed intensities_TurboNorm",
       caption = paste("Correlation coefficient:", round(correlation_2c_Turbo, 2),
                       "\nNumber of proteins:", nrow(df_2c_Turbo))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 1, vjust = 1))  # Align the caption to the bottom right
#) +
#scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
#scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

# Print the plot
print(correlation_plot_2c_Turbo)

dev.off()

#################################################################################
#################################################################################
################################################################################
## Part 2: Data analysis and visualization 
# Step 3. Volcanoes of collapsed intensity df

## 3. ANOVA / DiffEx -- 3 core user functions:
##
##   parANOVA.dex      -- create ANOVAout dataframe of tests for differential expression/abundance
##   plotVolc          -- create PDF and HTML Volcano Plots, output volcano settings to variables used later
##   DEXpercentStacked -- create PDF 
##
## By Eric Dammer, Duc Duong, and Qiudong Deng
################################################################################
#################################################################################
#################################################################################
## NOTES
##
## - Output ANOVA+Tukey pairwise stats for volcano and downstream analyses
## - if only 2 comparison groups are present, ANOVA overall p value is equivalent to T test;
##   FDR correction for all proteinwide comparisons provided in that special case (default twoGroupCorrMethod="BH")
## - can be repeat for different subgroup comparisons, if necessary
## - parallelized function. Requires R packages doParallel, parallel, and dependencies
## - writes DEX table typically to pipeline variable ANOVAout and .csv table
## - 0 Tukey values are inaccurate when <1e-9 or -10; Estimates of very small values become 0 in R base Tukey post hoc calculations.
## - So, we implement an option to fallback from Tukey p values less than 1e-10 to Bonferroni-corrected T test (unequal variance)
## - avoid dashes ("-") in strings representing comparison groups (Grouping vector); they will be substituted with '.'
##
##
################################################################################
## Required Loaded Data and Parameters ##
################################################################################
source("./parANOVA.dex.R")

#######################################################################################
## Volcano Parameters
FCmin=1                    # 0.25 for 25%, 0 for no threshold (vertical minimum FC threshold dashed lines)
selectComps="ALL"          # "ALL" for volcano output(s) on all pairwise comparisons in ANOVAout
flip=c(0)                  # ANOVAout column index numbers for p values in which to swap denominator of pair for x axis range (gene products high in denominator, will be on left)
# As a general rule, the group with less severe effects is usually set to be the denominator (represented by what is right of '-' in ANOVAout column names)
signifP=0.05               # p value threshold for counting Differential Expression points
useNETcolors=FALSE         # use module colors saved to ANOVAout, if available; otherwise when FALSE, specify downColor upColor, and NCcolor (must be valid R color specifications in quotes)
downColor="royalblue"      # significant points above/beyond thresholds on the upper left are this color if useNETcolors=FALSE
upColor="red"              # significant points above/beyond thresholds on the upper right are this color if useNETcolors=FALSE
NCcolor="grey"             # points not significant are this color if useNETcolors=FALSE
splitColors=FALSE          # create a separate volcano plot(s) for each color in an outputfigs/splitVolcanoes subfolder (folder created if it does not exist)
#highlightGeneProducts=c("Aqp4","Gfap","Aldh1l1","Camk2a","Ndrg2","S100b")  a list of uniqueID rownames to highlight as larger gold points. If symbolsOnly=TRUE, this can be a list of symbols, like c("APP","SMOC1","MAPT")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
symbolsOnly=TRUE           # for mouse-over HTML plots and the above highlight callouts, consider only displaying and using official gene symbol from first part of UniqueID rownames of ANOVAout.
labelTop=5                 # maximum p below which to label all points in the PDF output; OR an integer number of top ranked most significant points to label
labelSize=4.5              # text label font size, if any labels are found (when labelHighlighted=TRUE or labelTop>0)
sameScale=FALSE            # When multiple plots are drawn, should they all be the same scale with min and max x and y ranges?
HTMLout=TRUE               # output interactive HTML copies that can be opened in browser. Requires plotly package.
#outFilePrefix="5a"         # typically the step # in the pipeline being run
#outFileSuffix="AstroLPSvsNoLPS"
# A description of the project, used as a filename suffix
outputfigs=getwd()         # Location to save figure file output(s)

################################################################################
# Step 3: Diffabund volcanoes of intensity df

cleanDat <- read.csv("4c.cleanDat.bgSubtr.0centered.filtered.imputed-7333x9.csv_Perseus-Imputation_08302025.csv", header = TRUE, row.names = 1)
dim(cleanDat)
# [1] 7333    9

#cleanDat<-cleanDat.bgSubtr.filtered.imputed<-cleanDat.bgSubtr.meanMedian.unfiltered.imputed[rownames.toKeep,]

colnames(cleanDat)
# [1] "Input_1" "Input_2" "Input_3" "IP_4"    "IP_5"    "IP_6"    "IP_7"    "IP_8"    "IP_9"  

# Comparisons
  # a. HEK.SR64.pulldn vs HEK.global
  # b. HEK.SR71.global vs HEK.global
  # c. HEK.SR64.pulldn vs HEK.SR71.pulldn

numericMeta_nonegs <- numericMeta[c(1:3, 7:12),]
numericMeta_nonegs
#         DIA.MS_ID Sample.ID..           Group Replicate Group.simple
# Input_1   Input_1         CR1      HEK.global         1       Global
# Input_2   Input_2         CR2 HEK.SR64.global         2       Global
# Input_3   Input_3         CR3 HEK.SR71.global         3       Global
# IP_4         IP_4         CR7 HEK.SR64.pulldn         2  SR64_pulldn
# IP_5         IP_5         CR8 HEK.SR64.pulldn         3  SR64_pulldn
# IP_6         IP_6         CR9 HEK.SR64.pulldn         1  SR64_pulldn
# IP_7         IP_7        CR10 HEK.SR71.pulldn         2  SR71_pulldn
# IP_8         IP_8        CR11 HEK.SR71.pulldn         3  SR71_pulldn
# IP_9         IP_9        CR12 HEK.SR71.pulldn         1  SR71_pulldn

# a. HEK.SR64.pulldn vs HEK.global
SR64_pulldnVsglobal <- cleanDat[,which(numericMeta_nonegs$Group.simple=="SR64_pulldn" | numericMeta_nonegs$Group.simple=="Global")]
dim(SR64_pulldnVsglobal)
# [1] 7333    6
Grouping=numericMeta_nonegs$Group.simple[which(numericMeta_nonegs$Group.simple=="SR64_pulldn" | numericMeta_nonegs$Group.simple=="Global")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

head(SR64_pulldnVsglobal)
#                        Input_1    Input_2    Input_3       IP_4      IP_5       IP_6
# TAF11L11|A0A1W2PQ09 -2.2083268 -1.7193634 -2.3181059 -2.8600277 -3.834677 -3.8424965
# MCTS2|A0A3B3IRV3     0.3403200  0.3294251  0.5380254 -0.8692367 -1.406839 -0.2274323

cleanDat <- SR64_pulldnVsglobal
head(cleanDat)
#                        Input_1    Input_2    Input_3       IP_4      IP_5       IP_6
# TAF11L11|A0A1W2PQ09 -2.2083268 -1.7193634 -2.3181059 -2.8600277 -3.834677 -3.8424965

outFilePrefix="3a"
outFileSuffix= "SR64_pulldnVsglobal"
ANOVAout.SR64_pulldnVsglobal <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0) #3rd column of ANOVAout will be flipped (changes numerator vs denominator for volcano)
plotVolc()  #plots volcano using ANOVAout

################################################################################
# b. HEK.SR71.pulldn vs HEK.global

SR71_pulldnVsglobal <- cleanDat[,which(numericMeta_nonegs$Group.simple=="SR71_pulldn" | numericMeta_nonegs$Group.simple=="Global")]
dim(SR71_pulldnVsglobal)
# [1] 7333    6
Grouping=numericMeta_nonegs$Group.simple[which(numericMeta_nonegs$Group.simple=="SR71_pulldn" | numericMeta_nonegs$Group.simple=="Global")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

cleanDat <- SR71_pulldnVsglobal
head(cleanDat)
#                        Input_1    Input_2    Input_3       IP_7       IP_8       IP_9
# TAF11L11|A0A1W2PQ09 -2.2083268 -1.7193634 -2.3181059 -0.9435292 -1.7597154 -0.8521431

outFilePrefix="3b"
outFileSuffix= "SR71_pulldnVsglobal"
ANOVAout.SR71_pulldnVsglobal <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0) #3rd column of ANOVAout will be flipped (changes numerator vs denominator for volcano)
plotVolc()  #plots volcano using ANOVAout

################################################################################
# c. HEK.SR64.pulldn vs HEK.SR71.pulldn

SR64_pulldnVsSR71_pulldn <- cleanDat[,which(numericMeta_nonegs$Group.simple=="SR64_pulldn" | numericMeta_nonegs$Group.simple=="SR71_pulldn")]
dim(SR64_pulldnVsSR71_pulldn)
# [1] 7333    6
Grouping=numericMeta_nonegs$Group.simple[which(numericMeta_nonegs$Group.simple=="SR64_pulldn" | numericMeta_nonegs$Group.simple=="SR71_pulldn")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

cleanDat <- SR64_pulldnVsSR71_pulldn 
head(cleanDat)
#                           IP_4      IP_5       IP_6       IP_7       IP_8       IP_9
# TAF11L11|A0A1W2PQ09 -2.8600277 -3.834677 -3.8424965 -0.9435292 -1.7597154 -0.8521431

outFilePrefix="3c"
outFileSuffix= "SR64_pulldnVsSR71_pulldn"
ANOVAout.SR64_pulldnVsSR71_pulldn <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0) #3rd column of ANOVAout will be flipped (changes numerator vs denominator for volcano)
plotVolc()  #plots volcano using ANOVAout

################################################################################
################################################################################
## Part 2: Data analysis and visualization
# 4. GOparallel GSEA on volcano-defined DAP lists from uncollapsed intensity df
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

################################################################################
# a. HEK.SR64.pulldn vs HEK.global

outFilename <- "4a.SR64_pulldnVsglobal_0centered_uncollapsed.DAPs-GO"
GMTdatabaseFile="Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt"
GOparallel(ANOVAout.SR64_pulldnVsglobal)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all totals available.

################################################################################
# b. HEK.SR71.pulldn vs HEK.global

outFilename <- "4b.SR71_pulldnVsglobal_0centered_uncollapsed.DAPs-GO"
GMTdatabaseFile="Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt"
GOparallel(ANOVAout.SR71_pulldnVsglobal)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all totals available.

################################################################################
# c. HEK.SR64.pulldn vs HEK.SR71.pulldn

outFilename <- "4c.SR64_pulldnVsSR71_pulldn_0centered_uncollapsed.DAPs-GO"
GMTdatabaseFile="Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt"
GOparallel(ANOVAout.SR64_pulldnVsSR71_pulldn)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all totals available.


# end of data analysis and visualization
