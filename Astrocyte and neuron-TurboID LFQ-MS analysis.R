###################################################################################################
# Pipeline-Load+Norm+QC+FIlter+Impute+Volc...R
# Christina Ramelow, MS and Eric Dammer, PhD
# March 4th, 2024
###################################################################################################
# Based off Seyfried Lab Dammer "Loader-MaxQuant SummedIntensity (proteinGroups.txt).R" script
###################################################################################################
# ## Goals ###
# 1. Read in LFQ dataset from Exp. 9-2 and 10-1 new search w/0 Exp. 7-1 in vitro SPARO
# 2. Perform QC on log transformed dataset 
# 3. Determine which rows to keep to control missing values -- do not remove them yet
     # (not "<=" to 50%, but rather a minimum nonmissingness per group, and at least one group meeting its criteria)
# 4. Impute missing values so are able to do statistical analyses on the FULL dataset
     # THEN throw out rows that do not meet criteria from step #3.
# 5. Conduct t.test to compare samples + plot volcano
     # includes outputs to ANOVAout data frames, csv files, for subsetting diff abund proteins for further analysis.
###################################################################################################
options(stringsAsFactors=FALSE)

rootdir = "/Users/christina/Desktop/TurboID dual-omics/Exp. 9-2 and 10-1/2. LFQ-MS"
datadir=rootdir
setwd(rootdir)
###################################################################################################
## Part 1: Data processing and clean up
# steps:
# 1. Load and clean up proteinGroups.txt intensity matrix and load traits (numericMeta)
  # a. Remove categorical columns, contaminants and Camk2a missing mouse gene symbols
# 2. Normalization
  # a. pos.pulldown - neg.pulldown
  # b. pos.pulldown ribosomal protein normalization
# 3. Filter out missing values
  # Venns: Keep a protein if present in 2/3 samples/group
  # PCAs, correlations, diffabund: keep a protein if present in 2/3 pos.pulldown samples/group
# 4. Impute missing values
  # a. Perseus-style based on normal distribution
# 5. Collapse duplicate gene IDs 
  # a. keep a protein isoform with the most variance

## Part 2: Data analysis and visualization
  # 1. PCAs
  # 2. log2(intensity) correlations
  # 3. Volcanoes
  # 4. GOPar
####################################################################################################
####################################################################################################
## Part 1: Data processing and clean up
# steps:
# 1. Clean up proteinGroups.txt intensity matrix
  # a. Remove categorical columns, contaminants and Camk2a missing mouse gene symbols
####################################################################################################
#start fresh loading the case-sample raw summed intensity data
proteinGroups<-read.delim(file=paste0("Exp. 9-2 and 10-1 new search_proteinGroups.txt"),sep="\t",header=TRUE) #,row.names=1 #decoys can have same ID
dim(proteinGroups) 
# [1] 3179  239
decoyCount=length(which(proteinGroups$Reverse=="+"))
originalRows=nrow(proteinGroups)
decoyIndices=which(proteinGroups$Reverse=="+")
FDR=paste0(round(length(decoyIndices)/nrow(proteinGroups)*100,2),"% FDR")
LFQindices=which(grepl("Intensity.", colnames(proteinGroups)))
cat(paste0("Imported data has ",originalRows,"x",ncol(proteinGroups)," rows x columns; ",decoyCount," reverse hits, for a net ",FDR,".\n","Summed intensity is available for ",length(LFQindices)," experiment samples.\n"))

# Imported data has 3179x239 rows x columns; 52 reverse hits, for a net 1.64% FDR.
# Summed intensity is available for 22 experiment samples.

#make a matrix of just LFQ values with rownames as UniqueID (decoys removed, so they are truly unique now)
exprMat0<-proteinGroups[-decoyIndices,c(1,LFQindices)]
colnames(exprMat0)[1]<-"Protein.IDs"  #***

#Split Out Preferred Uniprot ID preference:(sp| > tr| > other)
IDdf<-as.data.frame(do.call(rbind,strsplit(as.character(exprMat0$Protein.IDs),";")))
IDdf2<-apply(IDdf,1,function(x) if(sum(as.numeric(grepl("sp\\|",x)))>0) { x[which(grepl("sp\\|",x))[1]] } else { if(sum(as.numeric(grepl("tr\\|",x)))>0) { x[which(grepl("tr\\|",x))[1]] } else { x[1] } } )
IDdf2<-apply(data.frame(col1=IDdf2),1,function(x) gsub("sp\\|","",gsub("tr\\|","",x)))
uniprotIDs.preferred<-as.data.frame(do.call(rbind,strsplit(IDdf2,"[|]")))[,1]

#Lookup HGNC Gene Symbols
symbolLookup<-readxl::read_excel("2020-symbolLookupTableForR_mouse.xlsx") #column names are "Uniprot.ID" and "Symbol"
symbolLookup<-as.data.frame(symbolLookup)
colnames(symbolLookup)<-c("Uniprot.ID","Symbol")

#lookup Uniprot symbol
symbols.lookup<-symbolLookup[match(uniprotIDs.preferred,symbolLookup[,1]),2]
symbols.lookup[is.na(symbols.lookup)]<-""

#make concatenated UniqueID for row.names of exprMat0
uniqueIDs<-paste0(symbols.lookup,"|",uniprotIDs.preferred)

#Set rownames to UniqueIDs and remove original column1
exprMat0<-exprMat0[,!colnames(exprMat0)=="Protein.IDs"]

rownames(exprMat0)<-uniqueIDs
exprMat0<-as.matrix(exprMat0)
dim(exprMat0)
# [1] 3127   22

colnames(exprMat0)<-gsub("Intensity.","Intensity ",colnames(exprMat0))

#clean up data structures we will not use again
rm(proteinGroups)
#rm(uniqueIDs)
rm(IDdf)
rm(IDdf2)

# Let's remove human contaminants, Trypsin, and LysC, but keep (probable) mouse gene products
conLookup.df<-read.csv(file="customMQcontamininantLookup.csv",header=T,row.names=1)
conLookup.rows.idx<-which(uniprotIDs.preferred %in% rownames(conLookup.df))
conLookups<-gsub("\\|","",rownames(exprMat0)[conLookup.rows.idx])
conLookups
# [1] "CON__O43790"           "CON__O76013"           "CON__O95678"           "CON__P00711"          
# [5] "CON__P00761"           "CON__P02533"           "CON__P02538"           "CON__P02662"          
# [9] "CON__P02754"           "CON__P02768-1"         "CON__P04259"           "CON__P04264"          
# [13] "CON__P05787"           "CON__P08727"           "CON__P08779"           "CON__P12035"          
# [17] "CON__P12763"           "CON__P13645"           "CON__P13646-1"         "CON__P13647"          
# [21] "CON__P15636"           "CON__P19013"           "CON__P35527"           "CON__P35908"          
# [25] "CON__P35908v2"         "CON__P48668"           "CON__P78386"           "CON__Q01546"          
# [29] "CON__Q04695"           "CON__Q14525"           "CON__Q9UE12"           "CON__Q3SY84"          
# [33] "CON__Q3SZH5"           "CON__Q5XQN5"           "CON__Q9NSB2"           "CON__Q6KB66-1"        
# [37] "CON__Q6NT21"           "CON__Q8N1N4-2"         "CON__Q7Z3Y8"           "CON__Q7Z3Z0"          
# [41] "CON__Q7Z794"           "CON__Q86YZ3"           "CON__Q9BYR9"           "CON__Q9NSB4"          
# [45] "CON__REFSEQ:XP_986630" "CON__Streptavidin"     "Q8BNF3"                "Q91WD4"   
for (uniprotID in conLookups) {
  if(conLookup.df$Keep.mouse[match(uniprotID,rownames(conLookup.df))]) {
    rownames(exprMat0)[which(uniprotIDs.preferred==uniprotID)] <- conLookup.df$UniqueID[which(rownames(conLookup.df)==uniprotID)]
  }
}
rownames(exprMat0)[conLookup.rows.idx]
# human contaminants removed later.
# Same as above       

# Lookup missing mouse symbols.
mmLookup.df<-read.csv(file="customMouseLookup.csv",header=T,row.names=1)
mmLookup.rows.idx<-which(uniprotIDs.preferred %in% rownames(mmLookup.df))
mmLookups<-gsub("\\|","",rownames(exprMat0)[mmLookup.rows.idx])
mmLookups
#  [1] "Q8C5J1"      "CON__Q3TTY5" "CON__A2A5Y0" "CON__P01966" "CON__P02769" "CON__P05784" "CON__P20930" "CON__Q9Z2K1"
# [9] "CON__Q5D862" "CON__Q61726" "CON__Q6IME9" "CON__Q9R0H5" "Q6Y642"      "Q3THL1"      "Q61674"      "Q3U4D1"     
# [17] "Q3UJM3"      "Q8BT09"      "Q3U7K1"      "Q8VER7"      "Q3TQQ9-3"    "Q3UJP5"      "Q61177"      "Q61934"     
# [25] "Q62013"      "Q8BR90"      "Q8C1L7"      "Q8C3W1"      "Q8CE18"      "Q8K1U6"      "Q9CY10"      "Q922R1"     
# [33] "Q9D454"      "Q9D9D9"      "V9GZG5"  

for (uniprotID in mmLookups) rownames(exprMat0)[which(uniprotIDs.preferred==uniprotID)] <- mmLookup.df$UniqueID[which(rownames(mmLookup.df)==uniprotID)]
rownames(exprMat0)[mmLookup.rows.idx]
#  [1] "Aldh3a2|Q8C5J1"       "Krt2|Q3TTY5"          "Krt31|A2A5Y0"         "Hba|P01966"           "Alb|P02769"          
# [6] "Krt18|P05784"         "Flg|P20930"           "Krt16|Q9Z2K1"         "Flg2|Q5D862"          "Krt83|Q61726"        
# [11] "Krt72|Q6IME9"         "Krt71|Q9R0H5"         "Gm20498|Q6Y642"       "Tmed2|Q3THL1"         "Hhex|Q61674"         
# [16] "Ahcy|Q3U4D1"          "Cox17|Q3UJM3"         "Rps6|Q8BT09"          "Ppp1cc|Q3U7K1"        "Olfr1308|Q8VER7"     
# [21] "Firrm|Q3TQQ9-3"       "2610301B20Rik|Q3UJP5" "Csnk2a1|Q61177"       "Myh4|Q61934"          "Luzp4|Q62013"        
# [26] "Rimoc1|Q8BR90"        "Rps21|Q8C1L7"         "2310022B05Rik|Q8C3W1" "Cdh13|Q8CE18"         "Hpx|Q8K1U6"          
# [31] "Hba|Q9CY10"           "Phaf1|Q922R1"         "4933412E24Rik|Q9D454" "Cfap141|Q9D9D9"       "Glul|V9GZG5"  

## Now remove human contaminants
dim(exprMat0)
# [1] 3127   22

exprMat0<-exprMat0[which(!rownames(exprMat0) %in% paste0("|",conLookups)), ]

dim(exprMat0)
# [1] 3079   22

##############################################################################
## Load traits (numericMeta)
numericMeta<-read.csv("Exp. 9-2 and 10-1_new search traits.csv", header = T)
rownames(numericMeta) <- numericMeta[,1]

dat.sumInt<-exprMat0[,match(rownames(numericMeta),colnames(exprMat0))]

dat.sumInt.log<-log2(dat.sumInt)
dat.sumInt.log[!is.finite(dat.sumInt.log)]<-NA

dim(dat.sumInt.log)
# [1] 3079   22

## Finalize Grouping of Samples for t.test
Grouping<-numericMeta$group.simple
Grouping
# [1] "neg_pulldn"  "neg_pulldn"  "neg_pulldn"  "neg_pulldn"  "neg_pulldn"  "Aldh1l1"     "Camk2a"      "Aldh1l1"    
# [9] "Aldh1l1.LPS" "Camk2a"      "Aldh1l1"     "Aldh1l1.LPS" "Camk2a"      "Aldh1l1.LPS" "Aldh1l1.LPS" "Aldh1l1.LPS"
# [17] "Camk2a"       "Camk2a"       "Camk2a"       "Camk2a"       "Camk2a"       "Camk2a"   

######################################################################################################
## Part 1: Data processing and clean up
# steps:
# 2. Normalization
  # a. pos.pulldown - neg.pulldown
  # b. pos.pulldown ribosomal protein normalization
######################################################################################################
# d. negative pulldown background subtraction 
## Background (negative pulldown) subtraction
dat.sumInt<-2^dat.sumInt.log

table(Grouping)
#Grouping
# Grouping
# Aldh1l1 Aldh1l1.LPS      Camk2a  neg_pulldn       bulk cortex 
# 3           5           3           5           6 

## From the above, select your group names for negative (.neg) and positive (.pos) in paired order to be used in the for loop within the apply below. Camk2a samples can be left out.
myGroups.neg=c("neg_pulldn")
myGroups.pos=c("Aldh1l1","Aldh1l1.LPS", "Camk2a")

dat.sumInt.bgSubtr <- t(apply(dat.sumInt,1,function(x) {
  for (i in 1:length(myGroups.pos)) {
    bkgr<-mean(x[which(Grouping==myGroups.neg[i])], na.rm=T)
    #bkgr[!is.finite(bkgr)]<- 0  # when no value in the background was present, we can't subtract anything (without first imputing) -- saves 31000 values from NaN status!
    #unfortunately, we cannot subtract background for any value in a row if one group is missed! -- keep NaN values in and fix by replacing the whole row back below.
    newRow=x
    newRow[which(Grouping==myGroups.pos[i])] <- x[which(Grouping==myGroups.pos[i])] - bkgr
  }
  newRow }))


## how many background subtracted values went <0 ? (should be a very small fraction of the large number)
length(which(dat.sumInt.bgSubtr<0))
#[1] 0

# we can set these to 0, for imputation later  (note: log2(0)=-Inf)
dat.sumInt.bgSubtr[dat.sumInt.bgSubtr<0]<- 0
#dat.sumInt.bgSubtr[!is.finite(dat.sumInt.bgSubtr)]<- NA


## find rows for which no background(s) were available -- they contain NaN -- and put back the unsubtracted values for the whole row in each case
length(which(is.nan(dat.sumInt.bgSubtr)))
# [1] 5627 # individual values in data

sum(unlist(apply(dat.sumInt.bgSubtr,1,function(x) if(length(which(is.nan(x)))>0) 1)))
# [1] 2149 # rows with at least 1 NaN

# perform cleanup, returning the unsubtracted matrix rows where there is at least 1 NaN
dat.sumInt.bgSubtr2<-t(sapply(rownames(dat.sumInt.bgSubtr),function(x) {
  if(length(which(is.nan(dat.sumInt.bgSubtr[x,])))>0) { 
    #cat(paste0(x,"\n"))
    dat.sumInt[x,]
  } else {
    dat.sumInt.bgSubtr[x,]
  }
}))

length(which(is.nan(dat.sumInt.bgSubtr2)))
#[1] 0  # should always be 0

dim(dat.sumInt.bgSubtr2)
# [1] 3079   22

dat.sumInt.bgSubtr.log<-log2(dat.sumInt.bgSubtr2)
dat.sumInt.bgSubtr.log[!is.finite(dat.sumInt.bgSubtr.log)]<-NA  # handles log2(0)

# subset to remove neg.pulldowns
#dat.sumInt.bgSubtr.nonegs.log <- dat.sumInt.bgSubtr.log[,c(1:6, 11:16)]
dat.sumInt.bgSubtr.nonegs.log <- dat.sumInt.bgSubtr.log[,which(!Grouping %in% myGroups.neg)]

dim(dat.sumInt.bgSubtr.nonegs.log)
# [1] 3079   17

######################################################################################################
## Ribosomal protein(nonmissing across all pulldown samples) for sum normalization after background subtraction

# Find nonmissing (in pos pulldowns) Rpl* and Rps* proteins
all.ribo<-which(grepl("^Rp[l|s]\\d+",rownames(dat.sumInt.bgSubtr2)))
rownames(dat.sumInt.bgSubtr2)[all.ribo]
# [1] "Rpl31|Q9CY93"          "Rps27|A0A0G2JDW7"      "Rpl9-ps6|A0A140T8T4"   "Rpl18|A0A1B0GQU8"     
# [5] "Rpl18a|A0A1D5RME4"     "Rps25|A0A1L1SQA8"      "Rps12|Q6ZWZ6"          "Rps24|A0A286YEB7"     
# [9] "Rpl36a-ps1|A0A2I3BPG9" "Rpl30|A0A2I3BQF4"      "Rpl10a|A0A3B2WBL1"     "Rpl19|A2A547"         
# [13] "Rpl26|B1ARA3"          "Rps15|D3YTQ9"          "Rps5|D3YYM6"           "Rps6ka1|P18653"       
# [17] "Rps28|G3UYV7"          "Rpl10|I7HLV2"          "Rpl21|Q9CQM8"          "Rpl35a|Q9DC85"        
# [21] "Rpl7a|Q6P1A9"          "Rpl27a|P14115"         "Rps16|Q5CZY9"          "Rpl7|Q3UBI6"          
# [25] "Rpl13a|P19253"         "Rps2|Q3TXS9"           "Rpl3|Q3UB90"           "Rpl12|Q3TIQ2"         
# [29] "Rpl28|Q5M9N5"          "Rpl6|P47911"           "Rpl5|Q3U850"           "Rpl13|P47963"         
# [33] "Rps20|P60867"          "Rpl27|P61358"          "Rpl37a|Q5M9N6"         "Rps7|P62082"          
# [37] "Rps8|P62242"           "Rps15a|P62245"         "Rps23|Q9CWI9"          "Rps18|P62270"         
# [41] "Rps11|Q9DB79"          "Rps13|P62301"          "Rps4x|Q545F8"          "Rpl23a|Q91YK6"        
# [45] "Rps6|Q8BT09"           "Rpl23|P62830"          "Rps26|P62855"          "Rpl39|P62892"         
# [49] "Rps3|P62908"           "Rpl32|P62911"          "Rpl8|Q3UJS0"           "Rps27a|P62983"        
# [53] "Rps17|P63276"          "Rps10|Q3U9P0"          "Rpl22|P67984"          "Rps3a1|Q3UAC2"        
# [57] "Rpl15|Q9CZM2"          "Rpl24|Q3UW40"          "Rpl17|Q6ZWZ7"          "Rpl36|Q5M9L1"         
# [61] "Rps19|Q5M9P3"          "Rps9|Q6ZWN5"           "Rpl35|Q6ZWV7"          "Rps27l|Q6ZWY3"        
# [65] "Rps6ka5|Q8C050-2"      "Rps21|Q8C1L7"          "Rpl11|Q8VC94"          "Rpl14|Q9CWK0"         
# [69] "Rpl34|Q9D1R9"          "Rpl4|Q9D8E6"           "Rpl38|Q9JJI8"  

all.ribo.posCols<-dat.sumInt.bgSubtr2[all.ribo, which(Grouping=="Aldh1l1" | Grouping=="Aldh1l1.LPS" | Grouping=="Camk2a")]
ribo.noMissingPos <- unlist( sapply(rownames(all.ribo.posCols),function(x) if(length(which(is.na(all.ribo.posCols[x,])))==0) x) )
names(ribo.noMissingPos)<-NULL

ribo.unmissed<-ribo.noMissingPos

#sanity check - no missing values here?
all.ribo.posCols[ribo.unmissed,]
# correct, no missing values

ribo.signals<-log2(colSums(dat.sumInt.bgSubtr2[which(rownames(dat.sumInt.bgSubtr2) %in% ribo.unmissed),]))
ribo.signals
# Intensity 21ip_01 Intensity 21ip_02 Intensity 21ip_04 Intensity 21ip_05 Intensity 21ip_07 Intensity 21ip_08 Intensity 21ip_09 
# NA                NA                NA                NA                NA          27.76628          27.10436 
# Intensity 21ip_11 Intensity 21ip_12 Intensity 21ip_13 Intensity 21ip_15 Intensity 21ip_16 Intensity 21ip_17 Intensity 21ip_18 
# 28.05976          27.62901          27.39127          27.06463          26.29016          27.40317          27.46719 
# Intensity 21ip_20 Intensity 21ip_21 Intensity 8tcl_01 Intensity 8tcl_02 Intensity 8tcl_04 Intensity 8tcl_05 Intensity 8tcl_07 
# 26.99784          25.55919          33.17859          33.51324          33.44218          33.62508          33.40082 
# Intensity 8tcl_08 
# 33.47098 

ribo.signals.posPulldown.mean<-mean(ribo.signals[6:16])  #mean of the 11 pos pulldowns
ribo.signals.posPulldown.mean
# [1] 27.15753
ribo.signals.log2.shifts = ribo.signals.posPulldown.mean - ribo.signals
# keep all samples that are not positive pulldowns anadjusted (0 shift)
ribo.signals.log2.shifts[c(1:5,17:22)] = rep(0, 11)  
ribo.signals.log2.shifts
# Intensity 21ip_01 Intensity 21ip_02 Intensity 21ip_04 Intensity 21ip_05 Intensity 21ip_07 Intensity 21ip_08 Intensity 21ip_09 
# 0.00000000        0.00000000        0.00000000        0.00000000        0.00000000       -0.60874894        0.05317646 
# Intensity 21ip_11 Intensity 21ip_12 Intensity 21ip_13 Intensity 21ip_15 Intensity 21ip_16 Intensity 21ip_17 Intensity 21ip_18 
# -0.90222539       -0.47147607       -0.23373711        0.09290736        0.86737116       -0.24563830       -0.30966115 
# Intensity 21ip_20 Intensity 21ip_21 Intensity 8tcl_01 Intensity 8tcl_02 Intensity 8tcl_04 Intensity 8tcl_05 Intensity 8tcl_07 
# 0.15969304        1.59833895        0.00000000        0.00000000        0.00000000        0.00000000        0.00000000 
# Intensity 8tcl_08 
# 0.00000000 

# apply shifts to the samples for normalized background subtracted matrix
dat.sumInt.bgSubtr.riboNorm.log2 <- t(t(dat.sumInt.bgSubtr.log) + ribo.signals.log2.shifts)
colMeans(dat.sumInt.bgSubtr.riboNorm.log2,na.rm=T) - colMeans(dat.sumInt.bgSubtr.log,na.rm=T)
# values of shift in column means should match the shifts above; CR: they match
# Intensity 21ip_01 Intensity 21ip_02 Intensity 21ip_04 Intensity 21ip_05 Intensity 21ip_07 Intensity 21ip_08 
# 0.00000000        0.00000000        0.00000000        0.00000000        0.00000000       -0.60874894 
# Intensity 21ip_09 Intensity 21ip_11 Intensity 21ip_12 Intensity 21ip_13 Intensity 21ip_15 Intensity 21ip_16 
# 0.05317646       -0.90222539       -0.47147607       -0.23373711        0.09290736        0.86737116 
# Intensity 21ip_17 Intensity 21ip_18 Intensity 21ip_20 Intensity 21ip_21 Intensity 8tcl_01 Intensity 8tcl_02 
# -0.24563830       -0.30966115        0.15969304        1.59833895        0.00000000        0.00000000 
# Intensity 8tcl_04 Intensity 8tcl_05 Intensity 8tcl_07 Intensity 8tcl_08 
# 0.00000000        0.00000000        0.00000000        0.00000000 


################################################################################
################################################################################
## QC Plot - step 2 - Unnorm Sum Intensity vs. Enforced Equal Loading Assumption WITHIN TREATMENT GROUP (samplewise summed intensity normalization) -- not MaxLFQ norm
pdf(file="2. Exp. 9-2 and 10-1 QC-runorder-SummedIntensity.Vs.ribosomal protein norm_03042024.pdf",width=21,height=14)
par(mar=c(12,5,3,3))
par(cex.axis=1.5)
par(cex.lab=1.5)
par(cex.main=3)

layout(matrix(c(1,2, 3,4, 5,6), nrow = 3, ncol = 2, byrow=TRUE),
       heights = c(0.95,0.95,0.95), # Heights of the rows
       widths = c(4.4,1)) # Widths of the columns

# colors for samples by group.simple
colvec=WGCNA::labels2colors(as.numeric(factor(Grouping)))
colvec.legend=WGCNA::labels2colors(as.numeric(factor(levels(factor(Grouping)))))

# boxplot 1: a. no normalization
boxplot(dat.sumInt.log, ylab=bquote("log"[2]~"MQ Summed Intensity"), main="a. Summed Intensity Protein Data (NO Normalization)", col=colvec, las=2)
par(mar=c(12,0,3,0))
plot.new()
legend("topleft", legend=levels(factor(Grouping)), fill=colvec.legend, cex=2.2)
par(mar=c(12,5,3,3))

# boxplot 2: b. neg.pulldown background subtraction
boxplot(dat.sumInt.bgSubtr.log[,c(6:16)], ylab=bquote("log"[2]~"MQ Summed Intensity"), main="Background Subtracted (Pulldowns Only)", col=colvec[c(6:16)], las=2)
plot.new()

# boxplot 3. c. ribosomal protein normalization
boxplot(dat.sumInt.bgSubtr.riboNorm.log2[,c(6:16)], ylab=bquote("log"[2]~"MQ Summed Intensity"), main="Ribosome Summed Intensity Samplewise Norm after Bkgr Subtr (Pulldowns Only)", col=colvec[c(6:16)], las=2)
plot.new()

dev.off()

dat.sumInt.log.norm <- dat.sumInt.bgSubtr.riboNorm.log2

###############################################################################
###############################################################################
###############################################################################
# 3. Determine which rows to keep to control missing values
#    (Enforce missingness across samples with logical criteria)
###############################################################################
###############################################################################
###############################################################################
rownames.toKeep<-unlist(sapply(rownames(dat.sumInt.bgSubtr.riboNorm.log2),function(x) {
  if(length(which(!is.na(dat.sumInt.bgSubtr.riboNorm.log2[x,which(Grouping=="Aldh1l1")]))) >= 2 | # at least 2 samples in each row kept are not NA
     length(which(!is.na(dat.sumInt.bgSubtr.riboNorm.log2[x,which(Grouping=="Camk2a")]))) >= 2 | # at least 2 samples in each row kept are not NA
     length(which(!is.na(dat.sumInt.bgSubtr.riboNorm.log2[x,which(Grouping=="Aldh1l1.LPS")]))) >= 3 ) x })) # at least 3 samples in each row kept are not NA

length(rownames.toKeep)
# [1] 2110


###############################################################################
# missingness threshold for Venn diagrams: keep a protein if present in 2/3 samples/group
# the input df does not have duplicate gene IDs collapsed
# Aldh1l1 filter
rownames.toKeep_Aldh1l1<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping=="Aldh1l1")]))) >= 2 ) x }))

length(rownames.toKeep_Aldh1l1)
# [1] 1967

cleanDat_Aldh1l1_filtered<-dat.sumInt.log.norm[rownames.toKeep_Aldh1l1,]

write.csv(cleanDat_Aldh1l1_filtered, file=paste0("Exp. 9-2 and 10-1_LFQ-MS_Aldh1l1_Venn filtered-",dim(cleanDat_Aldh1l1_filtered)[1],"x",dim(cleanDat_Aldh1l1_filtered)[2],".csv_04192024.csv"))

# Aldh1l1.LPS filter
rownames.toKeep_Aldh1l1.LPS<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping=="Aldh1l1.LPS")]))) >= 3 ) x }))

length(rownames.toKeep_Aldh1l1.LPS)
# [1] 1775

cleanDat_Aldh1l1.LPS_filtered<-dat.sumInt.log.norm[rownames.toKeep_Aldh1l1.LPS,]

write.csv(cleanDat_Aldh1l1.LPS_filtered, file=paste0("Exp. 9-2 and 10-1_LFQ-MS_Aldh1l1.LPS_Venn filtered-",dim(cleanDat_Aldh1l1.LPS_filtered)[1],"x",dim(cleanDat_Aldh1l1.LPS_filtered)[2],".csv_04192024.csv"))

# Camk2a filter
rownames.toKeep_Camk2a<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping=="Camk2a")]))) >= 2 ) x }))

length(rownames.toKeep_Camk2a)
# [1] 1881

cleanDat_input_filtered<-dat.sumInt.log.norm[rownames.toKeep_Camk2a,]

write.csv(cleanDat_input_filtered, file=paste0("Exp. 9-2 and 10-1_LFQ-MS_Camk2a_Venn filtered-",dim(cleanDat_input_filtered)[1],"x",dim(cleanDat_input_filtered)[2],".csv_04192024.csv"))

# bulk cortex filter
rownames.toKeep_bulk cortex<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping=="bulk cortex")]))) >= 2 ) x }))

length(rownames.toKeep_bulk cortex)
# [1] 2912

cleanDat_bulk cortex_filtered<-dat.sumInt.log.norm[rownames.toKeep_bulk cortex,]

write.csv(cleanDat_bulk cortex_filtered, file=paste0("Exp. 9-2 and 10-1_LFQ-MS_bulk cortex_Venn filtered-",dim(cleanDat_bulk cortex_filtered)[1],"x",dim(cleanDat_bulk cortex_filtered)[2],".csv_04192024.csv"))


# file names:
# Exp. 9-2 and 10-1_LFQ-MS_Aldh1l1_Venn filtered-1967x22.csv_04192024.csv
# Exp. 9-2 and 10-1_LFQ-MS_Aldh1l1.LPS_Venn filtered-1775x22.csv_04192024.csv
# Exp. 9-2 and 10-1_LFQ-MS_Camk2a_Venn filtered-1881x22.csv_04192024.csv
# Exp. 9-2 and 10-1_LFQ-MS_bulk cortex_Venn filtered-2912x22.csv_04192024.csv

#Load intensity files
Camk2a_protein<-read.csv("Exp. 9-2 and 10-1_LFQ-MS_Camk2a_Venn filtered-1881x22.csv_04192024.csv", header = TRUE)
bulk cortex_protein<-read.csv("Exp. 9-2 and 10-1_LFQ-MS_bulk cortex_Venn filtered-2912x22.csv_04192024.csv", header = TRUE)
Aldh1l1_protein<-read.csv("Exp. 9-2 and 10-1_LFQ-MS_Aldh1l1_Venn filtered-1967x22.csv_04192024.csv", header = TRUE)
Aldh1l1.LPS_protein<-read.csv("Exp. 9-2 and 10-1_LFQ-MS_Aldh1l1.LPS_Venn filtered-1775x22.csv_04192024.csv", header = TRUE)

## Subset the intensity matrices to extract only the protein ID columns to create Venn diagrams
Camk2a_protein_IDs <- data.frame(Camk2a_protein[,1])
colnames(Camk2a_protein_IDs)
colnames(Camk2a_protein_IDs) <- c("gene_protein_IDs")
library(tidyr)
Camk2a_protein_IDs_sep <- separate(Camk2a_protein_IDs, gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

bulk cortex_protein_IDs <- data.frame(bulk cortex_protein[,1])
colnames(bulk cortex_protein_IDs)
colnames(bulk cortex_protein_IDs) <- c("gene_protein_IDs")
bulk cortex_protein_IDs_sep <- separate(bulk cortex_protein_IDs, gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

Aldh1l1_protein_IDs <-data.frame(Aldh1l1_protein[,1])
colnames(Aldh1l1_protein_IDs)
colnames(Aldh1l1_protein_IDs) <- c("gene_protein_IDs")
Aldh1l1_protein_IDs_sep <- separate(Aldh1l1_protein_IDs, gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

Aldh1l1.LPS_protein_IDs <- data.frame(Aldh1l1.LPS_protein[,1])
colnames(Aldh1l1.LPS_protein_IDs)
colnames(Aldh1l1.LPS_protein_IDs) <- c("gene_protein_IDs")
Aldh1l1.LPS_protein_IDs_sep <- separate(Aldh1l1.LPS_protein_IDs, gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

Camk2a_protein_list <- Camk2a_protein_IDs_sep[,2]
bulk cortex_protein_list <- bulk cortex_protein_IDs_sep[,2]
Aldh1l1_protein_list <-  Aldh1l1_protein_IDs_sep[,2]
Aldh1l1.LPS_protein_list <- Aldh1l1.LPS_protein_IDs_sep[,2]
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
pdf(file = "Exp. 9-2 and 10-1_LFQ-MS_protein_Venns_CR_04192024.pdf", height = 11, width = 8.5, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# bulk cortex (x) vs Camk2a (y)
biovenn1 <- draw.venn(bulk cortex_protein_list, Camk2a_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="bulk cortex vs Camk2a protein", 
                      xtitle = "bulk cortex", 
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
                      x_c = "darkgrey",
                      y_c = "dodgerblue",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x bulk cortex: 2912"
# [1] "y bulk cortex: 1881"
# [1] "z bulk cortex: 0"
# [1] "x only: 1091"
# [1] "y only: 60"
# [1] "z only: 0"
# [1] "x-y bulk cortex overlap: 1821"

# bulk cortex (x) vs Aldh1l1 (y)
biovenn5 <- draw.venn(bulk cortex_protein_list, Aldh1l1_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="bulk cortex vs Aldh1l1 protein", 
                      xtitle = "bulk cortex", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "Aldh1l1",
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
                      y_c = "purple",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)
# [1] "x bulk cortex: 2912"
# [1] "y bulk cortex: 1967"
# [1] "z bulk cortex: 0"
# [1] "x only: 1008"
# [1] "y only: 63"
# [1] "z only: 0"
# [1] "x-y bulk cortex overlap: 1904"


# bulk cortex (x) vs Aldh1l1.LPS(y)
biovenn4 <- draw.venn(bulk cortex_protein_list, Aldh1l1.LPS_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="bulk cortex vs Aldh1l1.LPS protein", 
                      xtitle = "bulk cortex", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "Aldh1l1.LPS",
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
                      y_c = "darkturquoise",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x bulk cortex: 2912"
# [1] "y bulk cortex: 1775"
# [1] "z bulk cortex: 0"
# [1] "x only: 1192"
# [1] "y only: 55"
# [1] "z only: 0"
# [1] "x-y bulk cortex overlap: 1720"



# Aldh1l1 (x) vs Aldh1l1.LPS (y)
biovenn3 <- draw.venn(Aldh1l1_protein_list, Aldh1l1.LPS_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Aldh1l1 vs Aldh1l1.LPS protein", 
                      xtitle = "Aldh1l1", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "Aldh1l1.LPS",
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
                      y_c = "darkturquoise",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x bulk cortex: 1967"
# [1] "y bulk cortex: 1775"
# [1] "z bulk cortex: 0"
# [1] "x only: 244"
# [1] "y only: 52"
# [1] "z only: 0"
# [1] "x-y bulk cortex overlap: 1723"

# Camk2a (x) vs Aldh1l1 (y)
biovenn2 <- draw.venn(Camk2a_protein_list, Aldh1l1_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="Camk2a vs Aldh1l1 protein", 
                      xtitle = "Camk2a", 
                      xt_f = "Arial",
                      xt_fb = 2,
                      xt_s = 1,
                      xt_c = "black",
                      ytitle = "Aldh1l1",
                      yt_f = "Arial",
                      yt_fb = 2,
                      yt_s = 1,
                      yt_c = "black",
                      nrtype = "abs",
                      nr_f = "Arial",
                      nr_fb = 2,
                      nr_s = 1,
                      nr_c = "black",
                      x_c = "dodgerblue",
                      y_c = "purple",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)
# [1] "x bulk cortex: 1881"
# [1] "y bulk cortex: 1967"
# [1] "z bulk cortex: 0"
# [1] "x only: 116"
# [1] "y only: 202"
# [1] "z only: 0"
# [1] "x-y bulk cortex overlap: 1765"


dev.off()


# Create .csv files for common proteins between comparisons
# load library
library(dplyr)

# Venn overlap values
# bulk cortex (x) vs Camk2a (y)
# [1] "x-y bulk cortex overlap: 1821"

# bulk cortex (x) vs Aldh1l1 (y)
# [1] "x-y bulk cortex overlap: 1904"

# bulk cortex (x) vs Aldh1l1.LPS(y)
# [1] "x-y bulk cortex overlap: 1720"

# Aldh1l1 (x) vs Aldh1l1.LPS (y)
# [1] "x-y bulk cortex overlap: 1723"

# Camk2a (x) vs Aldh1l1 (y)
# [1] "x-y bulk cortex overlap: 1765"


# bulk cortex and Camk2a common proteins
common_proteins<-inner_join(bulk cortex_protein, Camk2a_protein)
dim(common_proteins)
# [1] 1821   23 CR: this matches the overlap number for the Venn
write.csv(common_proteins, file ="Exp. 9-2 and 10-1_LFQ-MS_bulk cortex and Camk2a_common proteins_CR_04192024.csv")

# bulk cortex & Aldh1l1 common proteins
common_proteins_4<-inner_join(bulk cortex_protein, Aldh1l1_protein)
dim(common_proteins_4)
# [1] 1904   13 CR: this matches the overlap number for the Venn
write.csv(common_proteins_4, file ="Exp. 9-2 and 10-1_LFQ-MS_bulk cortex and Aldh1l1_common proteins_CR_04192024.csv")

# bulk cortex & Aldh1l1.LPS common proteins
common_proteins_4<-inner_join(bulk cortex_protein, Aldh1l1.LPS_protein)
dim(common_proteins_4)
# [1] 1720   13 CR: this matches the overlap number for the Venn
write.csv(common_proteins_4, file ="Exp. 9-2 and 10-1_LFQ-MS_bulk cortex and Aldh1l1.LPS_common proteins_CR_04192024.csv")

# Aldh1l1 & Aldh1l1.LPS common proteins
common_proteins_2<-inner_join(Aldh1l1_protein, Aldh1l1.LPS_protein)
dim(common_proteins_2)
# [1] 1723  13 CR: this matches the overlap number for the Venn
write.csv(common_proteins_2, file ="Exp. 9-2 and 10-1_LFQ-MS_Aldh1l1 and Aldh1l1.LPS.pull_common proteins_CR_04192024.csv")

# Camk2a & Aldh1l1 common proteins
common_proteins_3<-inner_join(Camk2a_protein, Aldh1l1_protein)
dim(common_proteins_3)
# [1] 1765   17 CR: this matches the overlap number for the Venn
write.csv(common_proteins_3, file ="Exp. 9-2 and 10-1_LFQ-MS_Camk2a and Aldh1l1_common proteins_CR_04192024.csv")


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

# d. negative pulldown background subtraction
cleanDat.bgSubtr.unfiltered.imputed<-imputePerseusStyle(dat.sumInt.bgSubtr.log) #full data with imputation Perseus Style Column-wise Independent - all columns
head(cleanDat.bgSubtr.unfiltered.imputed)

dim(cleanDat.bgSubtr.unfiltered.imputed)
# [1] 3079   22

# write output with imputed values (log2-transformed), cleanDat.unfiltered.imputed
write.csv(cleanDat.bgSubtr.unfiltered.imputed, file=paste0("4d.cleanDat.bgSubtr.unfiltered.imputed-",dim(cleanDat.bgSubtr.unfiltered.imputed)[1],"x",dim(cleanDat.bgSubtr.unfiltered.imputed)[2],".csv_Perseus-Imputation_04192024.csv"))

# now includes imputed data points -- and we have just applied the filtering to keep only rows meeting criteria from step 3
cleanDat<-cleanDat.bgSubtr.filtered.imputed<-cleanDat.bgSubtr.unfiltered.imputed[rownames.toKeep,]
write.csv(cleanDat.bgSubtr.filtered.imputed, file=paste0("4d.cleanDat.bgSubtr.filtered.imputed-",dim(cleanDat.bgSubtr.filtered.imputed)[1],"x",dim(cleanDat.bgSubtr.filtered.imputed)[2],".csv_Perseus-Imputation_04192024.csv"))

dim(cleanDat)
# [1] 2110   22

# e. ribosomal protein normalization
cleanDat.bgSubtrRiboNorm.unfiltered.imputed<-imputePerseusStyle(dat.sumInt.bgSubtr.riboNorm.log2) #full data with imputation Perseus Style Column-wise Independent - all columns
head(cleanDat.bgSubtrRiboNorm.unfiltered.imputed)

dim(cleanDat.bgSubtrRiboNorm.unfiltered.imputed)
# [1] 3079   22

# write output with imputed values (log2-transformed), cleanDat.unfiltered.imputed
write.csv(cleanDat.bgSubtrRiboNorm.unfiltered.imputed, file=paste0("4e.cleanDat.bgSubtrRiboNorm.unfiltered.imputed-",dim(cleanDat.bgSubtrRiboNorm.unfiltered.imputed)[1],"x",dim(cleanDat.bgSubtrRiboNorm.unfiltered.imputed)[2],".csv_Perseus-Imputation_04192024.csv"))

# now includes imputed data points -- and we have just applied the filtering to keep only rows meeting criteria from step 3
cleanDat<-cleanDat.bgSubtrRiboNorm.filtered.imputed<-cleanDat.bgSubtrRiboNorm.unfiltered.imputed[rownames.toKeep,]
write.csv(cleanDat.bgSubtrRiboNorm.filtered.imputed, file=paste0("4e.cleanDat.bgSubtrRiboNorm.filtered.imputed-",dim(cleanDat.bgSubtrRiboNorm.filtered.imputed)[1],"x",dim(cleanDat.bgSubtrRiboNorm.filtered.imputed)[2],".csv_Perseus-Imputation_04192024.csv"))

dim(cleanDat)
#[1] 2110 22


###############################################################################
## Part 1: Data processing and clean up
# steps:
# 3. Collapse duplicate gene IDs 
  # a. Keep a protein isoform with the most variance
###############################################################################
###############################################################################
# normalized, non-filtered and non-imputed pulldowns intensity matrix as input
# rename normalized intensity matrix for simplicity
dat.sumInt.log.norm <- dat.sumInt.bgSubtr.riboNorm.log2
dim(dat.sumInt.log.norm)
# [1] 3079   22

# check to see if any gene symbols are missing, should be 0
badRows<-which(grepl("^\\|",rownames(dat.sumInt.log.norm)))
rownames(dat.sumInt.log.norm)[badRows]
# [1] "|Q91V76"

rownames(dat.sumInt.log.norm)[rownames(dat.sumInt.log.norm) == "|Q91V76"] <- "C11orf54|Q91V76"

badRows<-which(grepl("^\\|",rownames(dat.sumInt.log.norm)))
rownames(dat.sumInt.log.norm)[badRows]
# character(0)

myGroups.pos
# [1] "Aldh1l1"     "Aldh1l1.LPS" "Camk2a"

# subset to only include pos.pulldowns
#cleanDat <- dat.sumInt.log.norm[,c(6:16)]
cleanDat <- dat.sumInt.log.norm[,Grouping %in% myGroups.pos] # more efficient, less error-prone approach
dim(cleanDat)
# [1] 3079   11

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
# [1] 2214   11

###############################################################################
# a. normalized, filtered and imputed intensity matrix with pulldowns as input
cleanDat_pos <- cleanDat.bgSubtrRiboNorm.filtered.imputed[,Grouping %in% myGroups.pos] # more efficient, less error-prone approach
dim(cleanDat_pos)
# [1] 2110   11

# The norm intensity matrix row names are labeled gene symbol|protein ID, so we
# need to extract the first part of each row name to obtain gene symbols and put into a df
symbols <- as.data.frame(do.call(rbind, strsplit(rownames(cleanDat_pos), "[|]")))[, 1]
symbols <- as.data.frame(do.call(rbind, strsplit(symbols, "[;]")))[, 1]

# Collapse rows based on variance
collapsed_data_3 <- collapseRows(cleanDat_pos, symbols, rownames(cleanDat_pos), method = "MaxMean")

cleanDat.collapsed_3 <- collapsed_data_3$datETcollapsed
dim(cleanDat.collapsed_3)
# [1] 2060   11

write.csv(cleanDat.collapsed_3, file=paste0("5a. cleanDat.bgSubtrRiboNorm.filtered.imputed.collapsed-",dim(cleanDat.collapsed_3)[1],"x",dim(cleanDat.collapsed_3)[2],".csv_Perseus-Imputation_04192024.csv"))

# how many rows were collapsed?
2110 - 2060
# [1] 50

# what are the symbols|protein IDs for the collapsed rows?
which(!collapsed_data_3$selectedRow)
# Relch|A0A087WSS1 Camk2d|A0A0G2JGS4  Tanc2|A0A140LJD2 Camk2g|A0A286YCW8 Fam49b|A0A2I3BQK1 
# 26                62               126               188               198 
# Aak1|A0A571BDM4     Sptan1|A3KGU5      Vsnl1|B2L107        Hba|P01966    Epb41l3|D0VYV6 
# 244               328               367               414               421 
# Mapre3|D3Z6G3      Tubb3|D5MR34      Dnm1l|E9PUD2       Map4|E9PZ43       Myo6|E9Q3L1 
# 477               485               500               518               524 
# Actn4|E9Q2W9       Actb|E9Q5F4      Ncam1|E9QB01      Rims1|F6TZK4       Map2|G3UZJ2 
# 531               542               555               579               612 
# Stxbp1|O08599-2      Ap1b1|Q5SVG5      Mbp|P04370-6       Mapt|Q3UH19     Mapt|P10637-5 
# 663               688               753               781               782 
# Ncam1|P13595      Gnao1|P18872     Dnm1|P39053-4     Dnm1|P39053-6      Pkm|P52480-2 
# 802               832               912               913               977 
# Actg1|P63260   Ppp3ca|P63328-2      Rab5c|Q3TCT9      Vps35|Q3TJ43       Palm|Q3TRX4 
# 1107              1112              1218              1256              1272 
# Vim|Q3TWV0       Actb|Q3UBP6        Ina|Q3UMG4     Ogdh|Q60597-4     Syn2|Q64332-2 
# 1283              1322              1363              1441              1494 
# Cyfip2|Q6PGK0       Rtn1|Q7M6W1      Actn1|Q7TPR4      Cadps|Q80TJ1    Cadps|Q80TJ1-2 
# 1545              1570              1578              1590              1591 
# Ndrg4|Q8BTG7     Ank2|Q8C8R3-4        Dst|S4R1P5     Rtn4|Q99P72-1     Rtn3|Q9ES97-3 
# 1676              1715              1835              1882              2013 

# end of data processing and clean up pipeline

################################################################################
################################################################################
################################################################################
## Part 2: Data analysis and visualization
# 1. PCAs
# 2. log2(intensity) correlations
# 3. Volcanoes
# 4. GOparallel GSEA on volcano-defined DEP lists 
################################################################################
################################################################################
################################################################################
## Part 2: Data analysis and visualization
# 1. PCAs of uncollapsed data
# The cleanDat should be the normalized, filtered and imputed intensity matrix
# 4e.cleanDat.bgSubtrRiboNorm.filtered.imputed-2110x22.csv_Perseus-Imputation_04192024.csv
cleanDat <- read.csv("4e.cleanDat.bgSubtrRiboNorm.filtered.imputed-2110x22.csv_Perseus-Imputation_04192024.csv", header = TRUE, row.names = 1)
dim(cleanDat)
# [1] 2110   22

Exp9_2_and_10_1_intensities_nonegs <- cleanDat[,c(6:22)]

# load library
library(limma)

# EBD pdf formatting
#layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
#       heights = c(1.3,0.7,1.3), # Heights of the rows
#      widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

# CR pdf formatting
pdf(file = "1a. Exp. 9-2 and 10-1_LFQ-MS_all groups_nonegs_uncollapsed_Limma_PCA_CR_04182024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns


Grouping_vector_1 <- read.csv("Exp. 9-2 and 10-1_Grouping vector.csv", header = TRUE)
dim(Grouping_vector_1)
# [1] 22  1

Grouping_vector <- Grouping_vector_1[c(6:22),]


pch_values <- ifelse(Grouping_vector == "Camk2a", 16,
                     ifelse(Grouping_vector == "bulk cortex", 16,
                            ifelse(Grouping_vector == "Aldh1l1", 16, 
                                   ifelse(Grouping_vector == "Aldh1l1.LPS", 16,
                                          ifelse(Grouping_vector == "bulk.cortex.LPS ", 16, NA)))))

pt_colors <- ifelse(Grouping_vector == "Camk2a", "dodgerblue",
                    ifelse(Grouping_vector == "bulk cortex", "grey",
                           ifelse(Grouping_vector == "Aldh1l1", "purple3", 
                                  ifelse(Grouping_vector == "Aldh1l1.LPS", "darkcyan", 
                                         ifelse(Grouping_vector == "bulk.cortex.LPS ", "black", NA)))))

legend_groups <- c("bulk cortex", "bulk.cortex.LPS ", "Camk2a", "Aldh1l1", "Aldh1l1.LPS")

plotMDS_9_2_allgroups_nonegs <- plotMDS((Exp9_2_and_10_1_intensities_nonegs), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, ylim = c(-2, 3), xlim = c(-4, 7), plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for uncollapsed log2(intensity)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("grey", "black", "dodgerblue", "purple3", "darkcyan"), 
       pch = 16, title = "Proteome groups",cex=1.4)

dev.off()

###############################################################################
# Aldh1l1 vs Camk2a
Exp9_2_and_10_1_intensities_Aldh1l1VsCamk2a <- cleanDat[,c(6:8,10,11,13)]
dim(Exp9_2_and_10_1_intensities_Aldh1l1VsCamk2a)
# [1] 2110    6

colnames(Exp9_2_and_10_1_intensities_Aldh1l1VsCamk2a)
# [1] "Intensity 21ip_08" "Intensity 21ip_09" "Intensity 21ip_11" "Intensity 21ip_13" "Intensity 21ip_15" "Intensity 21ip_17"

#pdf(file ="1b. Exp. 9-2 and 10-1_LFQ-MS_Aldh vs Camk2a_uncollapsed_Limma_PCA_CR_04182024.pdf", width = 4, height = 4)

pdf(file = "1b. Exp. 9-2 and 10-1_LFQ-MS_Aldh vs Camk2a_uncollapsed_Limma_PCA_CR_04182024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

Grouping_vector <- c("Aldh1l1", "Camk2a","Aldh1l1","Camk2a","Aldh1l1","Camk2a")

pch_values <- ifelse(Grouping_vector== "Aldh1l1", 16,
                     ifelse(Grouping_vector == "Camk2a", 16, NA))

pt_colors <- ifelse(Grouping_vector== "Aldh1l1", "purple3",
                    ifelse(Grouping_vector == "Camk2a", "dodgerblue",NA))

legend_groups <- c("Aldh1l1", "Camk2a")


Exp9_2_and_10_1_symbols<- plotMDS(Exp9_2_and_10_1_intensities_Aldh1l1VsCamk2a, #top = 500, 
                                  labels = NULL, pch = pch_values,
                                  col = pt_colors, cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                  xlab = NULL, ylab = NULL,  ylim = c(-2, 2), xlim = c(-1.5, 1.5), plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for uncollapsed log2(intensity)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("purple3", "dodgerblue"), 
       pch = 16, title = "Proteome groups",cex=1.4)


dev.off()

################################################################################
# Aldh1l1 vs Aldh1l1.LPS
Exp9_2_and_10_1_intensities_Aldh1l1VsAldh1l1.LPS <- cleanDat[,c(6,8,9,11,12,14:16)]
dim(Exp9_2_and_10_1_intensities_Aldh1l1VsAldh1l1.LPS)
# [1] 2110    8

colnames(Exp9_2_and_10_1_intensities_Aldh1l1VsAldh1l1.LPS)
# [1] "Intensity 21ip_08" "Intensity 21ip_11" "Intensity 21ip_12" "Intensity 21ip_15" "Intensity 21ip_16" "Intensity 21ip_18"
# [7] "Intensity 21ip_20" "Intensity 21ip_21"

#pdf(file ="1b. Exp. 9-2 and 10-1_LFQ-MS_Aldh vs Camk2a_uncollapsed_Limma_PCA_CR_04182024.pdf", width = 4, height = 4)

pdf(file = "1c. Exp. 9-2 and 10-1_LFQ-MS_Aldh vs Aldh.LPS_uncollapsed_Limma_PCA_CR_04182024.pdf", width = 8, height = 8)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1,1,1), # Equal heights for all rows
       widths = c(1,1,1)) # Equal widths for all columns

Grouping_vector <- c("Aldh1l1", "Aldh1l1", "Aldh1l1.LPS", "Aldh1l1", "Aldh1l1.LPS", "Aldh1l1.LPS", "Aldh1l1.LPS", "Aldh1l1.LPS")

pch_values <- ifelse(Grouping_vector== "Aldh1l1", 16,
                     ifelse(Grouping_vector == "Aldh1l1.LPS", 16, NA))

pt_colors <- ifelse(Grouping_vector== "Aldh1l1", "purple3",
                    ifelse(Grouping_vector == "Aldh1l1.LPS", "darkcyan",NA))

legend_groups <- c("Aldh1l1", "Aldh1l1.LPS")


Exp9_2_and_10_1_symbols <- plotMDS(Exp9_2_and_10_1_intensities_Aldh1l1VsAldh1l1.LPS, #top = 500, 
                                  labels = NULL, pch = pch_values,
                                  col = pt_colors, cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                  xlab = NULL, ylab = NULL,  ylim = c(-2, 2), xlim = c(-2, 4), plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for uncollapsed log2(intensity)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("purple3", "darkcyan"), 
       pch = 16, title = "Proteome groups",cex=1.4)


dev.off()

################################################################################
## Part 2: Data analysis and visualization
# 2. log2(intensity) correlations of uncollapsed df

# Step 2. Correlations of mean intensities
  # a. Camk2a (neuron-TurboID) vs bulk cortex
  # b. Aldh1l1 (astrocyte-TurboID) vs bulk cortex
  # c. Aldh1l1.LPS vs bulk.cortex.LPS 
  # d. Aldh1l1 (astrocyte-TurboID) vs Camk2a (neuron-TurboID)
  # e. Aldh1l1.LPS (astrocyte-TurboID + LPS) vs Aldh1l1 (astrocyte-TurboID)

# load library to generate plot
library(ggplot2)

################################################################################
dim(cleanDat)
# [1] 2110   22

# 2a. Camk2a vs bulk cortex
pdf(file = "2a-e. Exp. 9-2 and 10-1_LFQ-MS_intensity correlations_uncollapsed_CR_04192024.pdf", width = 6, height = 6, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# subset into groups of interest and take the row mean
Camk2a <- cleanDat[,c(7,10,13)]
dim(Camk2a)
# [1] 2110    3
colnames(Camk2a)
# [1] "Intensity 21ip_09" "Intensity 21ip_13" "Intensity 21ip_17"
Camk2a_rowmeans <- rowMeans(Camk2a)

# Add row means as a new column to Camk2a data frame
Camk2a <- cbind(Camk2a, Camk2a_rowmeans)

bulk_cortex <- cleanDat[,c(17:20)]
dim(bulk_cortex)
# [1] 2110    4
colnames(bulk_cortex)
# [1] "Intensity 8tcl_01" "Intensity 8tcl_02" "Intensity 8tcl_04" "Intensity 8tcl_05"
bulk_cortex_rowmeans <- rowMeans(bulk_cortex)

# Add row means as a new column to Camk2a data frame
bulk_cortex <- cbind(bulk_cortex, bulk_cortex_rowmeans)

# generate a linear correlation
correlation_2a <- cor(Camk2a_rowmeans, bulk_cortex_rowmeans)
correlation_2a
# [1,] 0.4582397

# convert data to a data frame
df_2a <- data.frame(Camk2a_rowmeans = Camk2a_rowmeans, bulk_cortex_rowmeans = bulk_cortex_rowmeans)

# create the correlation plot
correlation_plot_2a <- ggplot(df_2a, aes(x = Camk2a_rowmeans, y = bulk_cortex_rowmeans)) +
  geom_point(color = "black", fill = "dodgerblue3", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Log2(mean intensity) Camk2a pulldown",
       y = "Log2(mean intensity) bulk cortex brain",
       title = "Linear correlation for uncollapsed intensities: 2a. Camk2a vs bulk cortex",
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
    plot.caption = element_text(hjust = 1, vjust = 1)  # Align the caption to the bottom right
  ) +
  scale_x_continuous(breaks = seq(18, 36, by = 9), limits = c(17, 36)) +
  scale_y_continuous(breaks = seq(18, 36, by = 9), limits = c(17, 36))

# Print the plot
print(correlation_plot_2a)


################################################################################
# b. Aldh1l1 vs bulk cortex

# subset into groups of interest and take the row mean
Aldh1l1 <- cleanDat[,c(6,8,11)]
dim(Aldh1l1)
# [1] 2110    3
colnames(Aldh1l1)
# [1] "Intensity 21ip_09" "Intensity 21ip_13" "Intensity 21ip_17"
Aldh1l1_rowmeans <- rowMeans(Aldh1l1)

# Add row means as a new column to Aldh1l1 data frame
Aldh1l1 <- cbind(Aldh1l1, Aldh1l1_rowmeans)

# bulk cortex group subsetted and rowmeaned in 2a above

# generate a linear correlation
correlation_2b <- cor(Aldh1l1_rowmeans, bulk_cortex_rowmeans)
correlation_2b
# [1,] 0.4791407

# convert data to a data frame
df_2b <- data.frame(Aldh1l1_rowmeans = Aldh1l1_rowmeans, bulk_cortex_rowmeans = bulk_cortex_rowmeans)

# create the correlation plot
correlation_plot_2b <- ggplot(df_2b, aes(x = Aldh1l1_rowmeans, y = bulk_cortex_rowmeans)) +
  geom_point(color = "black", fill = "purple4", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Log2(mean intensity) Aldh1l1 pulldown",
       y = "Log2(mean intensity) bulk cortex brain",
       title = "2b. Aldh1l1 vs bulk cortex",
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
    plot.caption = element_text(hjust = 1, vjust = 1)  # Align the caption to the bottom right
  ) +
  scale_x_continuous(breaks = seq(18, 36, by = 9), limits = c(17, 36)) +
  scale_y_continuous(breaks = seq(18, 36, by = 9), limits = c(17, 36))

# Print the plot
print(correlation_plot_2b)

################################################################################
# c. Aldh1l1.LPS vs bulk.cortex.LPS 

# subset into groups of interest and take the row mean
Aldh1l1.LPS <- cleanDat[,c(9,12,14:16)]
dim(Aldh1l1.LPS)
# [1] 2110    5
colnames(Aldh1l1.LPS)
# [1] "Intensity 21ip_12" "Intensity 21ip_16" "Intensity 21ip_18" "Intensity 21ip_20" "Intensity 21ip_21"
Aldh1l1.LPS_rowmeans <- rowMeans(Aldh1l1.LPS)

# Add row means as a new column to Aldh1l1.LPS data frame
Aldh1l1.LPS <- cbind(Aldh1l1.LPS, Aldh1l1.LPS_rowmeans)

bulk.cortex.LPS  <- cleanDat[,c(21:22)]
dim(bulk.cortex.LPS )
# [1] 2110    2
colnames(bulk.cortex.LPS )
# [1] "Intensity 8tcl_07" "Intensity 8tcl_08"
bulk.cortex.LPS_rowmeans <- rowMeans(bulk.cortex.LPS )

# Add row means as a new column to Aldh1l1.LPS data frame
bulk.cortex.LPS  <- cbind(bulk.cortex.LPS , bulk.cortex.LPS_rowmeans)

# generate a linear correlation
correlation_2c <- cor(Aldh1l1.LPS_rowmeans, bulk.cortex.LPS_rowmeans)
correlation_2c
# [1,] 0.4792233

# convert data to a data frame
df_2c <- data.frame(Aldh1l1.LPS_rowmeans = Aldh1l1.LPS_rowmeans, bulk.cortex.LPS_rowmeans = bulk.cortex.LPS_rowmeans)

# create the correlation plot
correlation_plot_2c <- ggplot(df_2c, aes(x = Aldh1l1.LPS_rowmeans, y = bulk.cortex.LPS_rowmeans)) +
  geom_point(color = "black", fill = "cyan4", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Log2(mean intensity) Aldh1l1.LPS pulldown",
       y = "Log2(mean intensity) bulk.cortex.LPS  brain",
       title = "2c. Aldh1l1.LPS vs bulk.cortex.LPS ",
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
    plot.caption = element_text(hjust = 1, vjust = 1)  # Align the caption to the bottom right
  ) +
  scale_x_continuous(breaks = seq(18, 36, by = 9), limits = c(17, 36)) +
  scale_y_continuous(breaks = seq(18, 36, by = 9), limits = c(17, 36))

# Print the plot
print(correlation_plot_2c)

################################################################################
# d. Aldh1l1 vs Camk2a

# Aldh1l1 and Camk2a groups subsetted and rowmean above

# generate a linear correlation
correlation_2d <- cor(Aldh1l1_rowmeans, Camk2a_rowmeans)
correlation_2d
# [1,] 0.9108584

# convert data to a data frame
df_2d <- data.frame(Aldh1l1_rowmeans = Aldh1l1_rowmeans, Camk2a_rowmeans = Camk2a_rowmeans)

# create the correlation plot
correlation_plot_2d <- ggplot(df_2d, aes(x = Aldh1l1_rowmeans, y = Camk2a_rowmeans)) +
  geom_point(color = "black", fill = "gray", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Log2(mean intensity) Aldh1l1 pulldown",
       y = "Log2(mean intensity) Camk2a pulldown",
       title = "2d. Aldh1l1 vs Camk2a",
       caption = paste("Correlation coefficient:", round(correlation_2d, 2),
                       "\nNumber of proteins:", nrow(df_2d))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 1, vjust = 1)  # Align the caption to the bottom right
  ) +
  scale_x_continuous(breaks = seq(15, 35, by = 10), limits = c(15, 35)) +
  scale_y_continuous(breaks = seq(15, 35, by = 10), limits = c(15, 35))

# Print the plot
print(correlation_plot_2d)

################################################################################
# e. Aldh1l1.LPS vs Aldh1l1

# Aldh1l1 and Aldh1l1.LPS groups subsetted and rowmeaned above

# generate a linear correlation
correlation_2e <- cor(Aldh1l1_rowmeans, Aldh1l1.LPS_rowmeans)
correlation_2e
# [1,] 0.9645978

# convert data to a data frame
df_2e <- data.frame(Aldh1l1_rowmeans = Aldh1l1_rowmeans, Aldh1l1.LPS_rowmeans = Aldh1l1.LPS_rowmeans)

# create the correlation plot
correlation_plot_2e <- ggplot(df_2e, aes(x = Aldh1l1_rowmeans, y = Aldh1l1.LPS_rowmeans)) +
  geom_point(color = "black", fill = "cyan2", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "Log2(mean intensity) Aldh1l1.LPS pulldown",
       y = "Log2(mean intensity) Aldh1l1 pulldown",
       title = "2e. Aldh1l1.LPS vs Aldh1l1",
       caption = paste("Correlation coefficient:", round(correlation_2e, 2),
                       "\nNumber of proteins:", nrow(df_2e))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 1, vjust = 1)  # Align the caption to the bottom right
  ) +
  scale_x_continuous(breaks = seq(15, 35, by = 10), limits = c(15, 35)) +
  scale_y_continuous(breaks = seq(15, 35, by = 10), limits = c(15, 35))

# Print the plot
print(correlation_plot_2e)

dev.off()

################################################################################
# histogram of Rpl22 intensity
library(ggplot2)
library(dplyr)

# Read in the data
protein_traits <- read.csv("Exp. 9-2 and 10-1_new search traits.csv", header = TRUE)
norm_intensities <- read.csv("4e.cleanDat.nobgSubtr.RiboNorm.filtered.imputed-2110x22.csv_Perseus-Imputation_03192024.csv", header = TRUE, row.names = 1)

match <- traits[match(rownames(norm_intensities), traits$group.simple), ]

if (!all(match$sample_id == rownames(counts))) {
  stop("Sample IDs in traits and counts files do not match.")
}

# Specify the protein of interest
protein_of_interest <- "Rpl22|P67984"

# Extract the intensities for the specific protein
Rpl22_intensities <- as.numeric(norm_intensities[protein_of_interest, ])
Rpl22_intensities 
# [1] 17.18389 15.40123 19.56493 16.08611 15.86526 18.07589
# [7] 17.48227 20.90022 20.92885 17.72620 17.49742 17.89642
# [13] 18.60430 19.42798 18.42505 16.57640 26.41774 26.66370
# [19] 26.79996 26.89862 26.43122 26.85007

# Combine protein intensities with experimental group information
data <- data.frame(sample_id = colnames(norm_intensities),
                   Rpl22_intensities = Rpl22_intensities,
                   group = protein_traits$group.simple)

# Check for missing values in the group column
if (any(is.na(data$group))) {
  stop("There are missing values in the group column.")
}

# Ensure the group column is a factor
data$group <- as.factor(data$group)

# Create histograms separated by experimental group
Rpl22_botplot <- ggplot(data, aes(x = group, y = Rpl22_intensities, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = paste("Box Plot of", protein_of_interest, "Protein abundance"),
       x = "Experimental Group",
       y = "Normalized Intensities") +
  theme_minimal()

print(Rpl22_botplot)

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
# Step 3: Diffabund volcanoes of collapsed intensity df
  # comparisons:
  # d. Aldh1l1 (astrocyte-TurboID) vs Camk2a (neuron-TurboID)
  # e. Aldh1l1.LPS (astrocyte-TurboID + LPS) vs Aldh1l1

numericMeta_parent <-read.csv("Exp. 9-2 and 10-1_new search traits.csv", header = T)

# subset the numericMeta to only include the positive pulldowns
numericMeta_nonegs <- numericMeta_parent[c(6:22),]
dim(numericMeta_nonegs)
# [1] 17  4

numericMeta <- numericMeta_nonegs[c(1:11),]
dim(numericMeta)
# [1] 11  4

# collapsed df
dim(cleanDat.collapsed_3)
# [1] 2060   11

## 3d Aldh1l1 vs Camk2a - ribo norm - collapsed
cleanDat<-cleanDat.collapsed_3[,which(numericMeta$group.simple=="Aldh1l1" | numericMeta$group.simple=="Camk2a")]
dim(cleanDat)
# [1] 2060    6
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="Aldh1l1" | numericMeta$group.simple=="Camk2a")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="3d_collapsed"
outFileSuffix=" Aldh1l1vsCamk2a.riboNorm"
ANOVAout.Aldh1l1vsCamk2a.riboNorm <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required Camk2a.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(3) #3rd column of ANOVAout will be flipped (changes numerator vs denominator for volcano)
plotVolc()  #plots volcano using ANOVAout

## 3e repeat for Aldh1l1.LPS vs Aldh1l1 - ribo norm - collapsed
cleanDat<-cleanDat.collapsed_3[,which(numericMeta$group.simple=="Aldh1l1.LPS" | numericMeta$group.simple=="Aldh1l1")]
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="Aldh1l1.LPS" | numericMeta$group.simple=="Aldh1l1")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="3e_collapsed"
outFileSuffix=" Aldh1l1.LPSvsAldh1l1.riboNorm_collapsed"
ANOVAout.Aldh1l1.LPSvsAldh1l1.riboNorm <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required Camk2a.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0)
plotVolc()  #plots volcano using ANOVAout

################################################################################
## Part 2: Data analysis and visualization
# 4. GOparallel GSEA on volcano-defined DEP lists from collapsed intensity df
################################################################################
source("GOparallel-FET.R")

## Parameters saved to variables in memory.

######################## EDIT THESE VARIABLES (USER PARAMETERS SET IN GLOBAL ENVIRONMENT) ############################################
#input File <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#input File <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
            # Input CSV FILE - in the filePath folder.
            #Can be formatted as Kme table from WGCNA pipeline, or
            #can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
            #in this case, the longest list is used as background or the "universe" for the FET contingencies
            #  For simple columnwise list Input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.

filePath <- getwd() 
            #Folder that (may) contain the Input file specified above, and which will contain the outFilename project Folder.
modulesInMemory=FALSE
ANOVAgroups=TRUE  #change to FALSE if using Input File
            #if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
            #input file will be ignored
maxBarsPerOntology=6
panelDimensions=c(2,2)    #dimensions of the individual parblots within a page of the main barplot PDF output
pageDimensions=c(8.5,11)  #main barplot PDF output page dimensions, in inches

############ MUST HAVE AT LEAST 2 THREADS ENABLED TO RUN ############################################################################
parallelThreads=7
removeRedundantGOterms=TRUE
            #if true, the 3 GO ontology types are collapased into a minimal set of less redundant terms using the below OBO file
cocluster=FALSE
            #If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)

######################## END OF PARAMETER VARIABLES ###################################################################################
#to lookup rows with missing genes symbols that halt GOPar from running
badRows<-which(grepl("^\\|",rownames(cleanDat)))
rownames(cleanDat)[badRows]
#character(0), so we don't have any missing gene symbols :)

rm(ANOVAout)
# rm(dexComps)
# Warning message:
# In rm(dexComps) : object 'dexComps' not found
rm(comparisonIDs)
rm(testIndexMasterList)


###############################################################################
## 3d: Aldh1l1 vs Camk2a DEP list GSEA - ribo norm - collapsed
flip=c(3) #3rd column of ANOVAout will be flipped
outFilename <- "3d_collapsed_Aldh1l1VsCamk2a_riboNorm_collapsed_DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.Aldh1l1vsCamk2a.riboNorm)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all Camk2as available.

## 3e: Aldh1l1.LPS vs Aldh1l1 DEP list GSEA - ribo norm - collapsed
flip=c(0)
outFilename <- "3e_collapsed_Aldh1l1.LPSVsAldh1l1_riboNorm_DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.Aldh1l1.LPSvsAldh1l1.riboNorm)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all Camk2as available.


################################################################################
################################################################################
# Step 3: Diffabund analysis and GOPar of uncollapsed intensity df

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
# Step 3: Diffabund analysis and GOPar of uncollapsed intensity df
# comparisons:
  # a: Camk2a vs bulk cortex
  #    Camk2avsbulk cortex
  # b. Aldh1l1 vs bulk cortex
  #    Aldh1l1vsbulk cortex
  # c. Aldh1l1.LPS vs bulk.cortex.LPS 
  #    Aldh1l1.LPSvsbulk.cortex.LPS 
  # d. Aldh1l1 vs Camk2a
  #    Aldh1l1vsCamk2a
  # e. Aldh1l1.LPS vs Aldh1l1
  #    Aldh1l1.LPSvsAldh1l1`

# update numericMeta
numericMeta <- numericMeta_parent
dim(numericMeta)
# [1] 22  4

numericMeta
#               label                                  group replicate group.simple
# 1  Intensity 21ip_01              Aldh1l1_negative_pulldown         1   neg_pulldn
# 2  Intensity 21ip_02               Camk2a_negative_pulldown         1   neg_pulldn
# 3  Intensity 21ip_04          Aldh1l1.LPS_negative_pulldown         1   neg_pulldn
# 4  Intensity 21ip_05               Camk2a_negative_pulldown         2   neg_pulldn
# 5  Intensity 21ip_07          Aldh1l1.LPS_negative pulldown         2   neg_pulldn
# 6  Intensity 21ip_08              Aldh1l1_positive_pulldown         1      Aldh1l1
# 7  Intensity 21ip_09               Camk2a_positive_pulldown         1       Camk2a
# ...

# update cleanDat
uncollapsed.cleanDat <- read.csv("4e.cleanDat.bgSubtrRiboNorm.filtered.imputed-2110x22.csv_Perseus-Imputation_04192024.csv", header = TRUE, row.names = 1)
dim(uncollapsed.cleanDat)
# [1] 2110   22

# 3d. Aldh1l1 vs Camk2a
#    Aldh1l1vsCamk2a
cleanDat<-uncollapsed.cleanDat[,which(numericMeta$group.simple=="Aldh1l1" | numericMeta$group.simple=="Camk2a")]
dim(cleanDat)
# [1] 2110    6

Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="Aldh1l1" | numericMeta$group.simple=="Camk2a")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.
Grouping
# [1] "Aldh1l1" "Camk2a"  "Aldh1l1" "Camk2a"  "Aldh1l1" "Camk2a" 

outFilePrefix="3d"
outFileSuffix="Aldh1l1vsCamk2a.riboNorm.uncollapsed"
ANOVAout.Aldh1l1vsCamk2a.riboNorm.uncollapsed <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required Camk2a.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(3)
plotVolc()  #plots volcano using ANOVAout

# 3e. Aldh1l1.LPS vs Aldh1l1
#    Aldh1l1.LPSvsAldh1l1`
cleanDat<-uncollapsed.cleanDat[,which(numericMeta$group.simple=="Aldh1l1" | numericMeta$group.simple=="Aldh1l1.LPS")]
dim(cleanDat)
# [1] 2110    8

Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="Aldh1l1" | numericMeta$group.simple=="Aldh1l1.LPS")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.
Grouping
# [1] "Aldh1l1"     "Aldh1l1"     "Aldh1l1.LPS" "Aldh1l1"     "Aldh1l1.LPS" "Aldh1l1.LPS" "Aldh1l1.LPS"
# [8] "Aldh1l1.LPS"

outFilePrefix="3e"
outFileSuffix="Aldh1l1.LPSvsAldh1l1.riboNorm.uncollapsed"
ANOVAout.Aldh1l1.LPSvsAldh1l1.riboNorm.uncollapsed <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required Camk2a.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0)
plotVolc()  #plots volcano using ANOVAout

################################################################################
## Part 2: Data analysis and visualization
# Step 4. GOparallel GSEA on volcano-defined DEP lists from uncollapsed intensity df
################################################################################
source("GOparallel-FET.R")

## Parameters saved to variables in memory.

######################## EDIT THESE VARIABLES (USER PARAMETERS SET IN GLOBAL ENVIRONMENT) ############################################
#Input File <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#InputFile <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
#InputCSV FILE - in the filePath folder.
#Can be formatted as Kme table from WGCNA pipeline, or
#can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
#in this case, the longest list is used as background or the "universe" for the FET contingencies
#  For simple columnwise list Input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.

filePath <- getwd() 
#Folder that (may) contain the Input file specified above, and which will contain the outFilename project Folder.
modulesInMemory=FALSE
ANOVAgroups=TRUE  #change to FALSE if using Camk2aFile
#if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
#Camk2aFile will be ignored
maxBarsPerOntology=6
panelDimensions=c(2,2)    #dimensions of the individual parblots within a page of the main barplot PDF output
pageDimensions=c(8.5,11)  #main barplot PDF output page dimensions, in inches

############ MUST HAVE AT LEAST 2 THREADS ENABLED TO RUN ############################################################################
parallelThreads=7
removeRedundantGOterms=TRUE
#if true, the 3 GO ontology types are collapased into a minimal set of less redundant terms using the below OBO file
cocluster=FALSE
#If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)

######################## END OF PARAMETER VARIABLES ###################################################################################
#to lookup rows with missing genes symbols that halt GOPar from running
badRows<-which(grepl("^\\|",rownames(cleanDat)))
rownames(cleanDat)[badRows]
#character(0), so we don't have any missing gene symbols :)

rm(ANOVAout)

###############################################################################
# Step 4. GOparallel GSEA on volcano-defined DEP lists from uncollapsed intensity df
# comparisons:
# 4d. Aldh1l1 (astrocyte-TurboID) vs Camk2a (neuron-TurboID)
# 4e. Aldh1l1.LPS (astrocyte-TurboID + LPS) vs Aldh1l1 (astrocyte-TurboID)


# 4d. Aldh1l1 vs Camk2a
#    Aldh1l1vsCamk2a
flip=c(3)
outFilename <- "4d_Aldh1l1vsCamk2a.riboNorm_uncollapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.Aldh1l1vsCamk2a.riboNorm.uncollapsed)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all Camk2as available.

#################################################################################
#################################################################################
# 4e. Aldh1l1.LPS vs Aldh1l1
#    Aldh1l1.LPSvsAldh1l1`
flip=c(0)
outFilename <- "4e_Aldh1l1.LPSvsAldh1l1.riboNorm_uncollapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.Aldh1l1.LPSvsAldh1l1.riboNorm.uncollapsed)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all Camk2as available.


#################################################################################
# end of data analysis and visualization

