###################################################################################################
# Pipeline-Load+Norm+QC+FIlter+Impute+Volc...R
# Christina Ramelow, MS and Eric Dammer, PhD
# February 13th, 2024
###################################################################################################
# Based off Seyfried Lab Dammer "Loader-MaxQuant SummedIntensity (proteinGroups.txt).R" script
###################################################################################################
# ## Goals ###
# 1. Read in LFQ dataset from Exp. 1-4 new search w/0 Exp. 7-1 in vitro CRAPP
# 2. Perform QC on log transformed dataset 
# 3. Determine which rows to keep to control missing values -- do not remove them yet
     # (not "<=" to 50%, but rather a minimum nonmissingness per group, and at least one group meeting its criteria)
# 4. Impute missing values so are able to do statistical analyses on the FULL dataset
     # THEN throw out rows that do not meet criteria from step #3.
# 5. Conduct t.test to compare samples + plot volcano
     # includes outputs to ANOVAout data frames, csv files, for subsetting diff abund proteins for further analysis.
###################################################################################################
options(stringsAsFactors=FALSE)

rootdir = "/Users/christina/Desktop/TurboID dual-omics/Exp. 1-4/2. LFQ-MS"
datadir=rootdir
setwd(rootdir)
###################################################################################################
## Part 1: Data processing and clean up
# steps:
# 1. Load and clean up proteinGroups.txt intensity matrix and load traits (numericMeta)
  # a. Remove categorical columns, contaminants and global missing mouse gene symbols
# 2. Normalization
  # a. pos.pulldown - neg.pulldown
  # b. globals and (pos.pulldown - neg.pulldown) median zero-centered
# 3. Filter out missing values
  # Venns: Keep a protein if present in 2/3 samples/ group
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
  # a. Remove categorical columns, contaminants and global missing mouse gene symbols
####################################################################################################
#start fresh loading the case-sample raw summed intensity data
proteinGroups<-read.delim(file=paste0("Exp. 1-4 new search_no Exp. 7-1_proteinGroups.txt"),sep="\t",header=TRUE) #,row.names=1 #decoys can have same ID
dim(proteinGroups) 
# [1] 2873  185
decoyCount=length(which(proteinGroups$Reverse=="+"))
originalRows=nrow(proteinGroups)
decoyIndices=which(proteinGroups$Reverse=="+")
FDR=paste0(round(length(decoyIndices)/nrow(proteinGroups)*100,2),"% FDR")
LFQindices=which(grepl("Intensity.", colnames(proteinGroups)))
cat(paste0("Imported data has ",originalRows,"x",ncol(proteinGroups)," rows x columns; ",decoyCount," reverse hits, for a net ",FDR,".\n","Summed intensity is available for ",length(LFQindices)," experiment samples.\n"))

# Imported data has 2873x185 rows x columns; 36 reverse hits, for a net 1.25% FDR.
# Summed intensity is available for 16 experiment samples.

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
# [1] 2837   16

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
# [1] "CON__Q3TTY5"       "CON__O43790"       "CON__O76013"       "CON__O95678"       "CON__P00761"      
# [6] "CON__P01966"       "CON__P02533"       "CON__P02538"       "CON__P02662"       "CON__P02663"      
# [11] "CON__P02754"       "CON__P02768-1"     "CON__P02769"       "CON__P04259"       "CON__P04264"      
# [16] "CON__P05784"       "CON__P08727"       "CON__P08779"       "CON__P12035"       "CON__P12763"      
# [21] "CON__P13645"       "CON__P13646-1"     "CON__P13647"       "CON__P15636"       "CON__P19013"      
# [26] "CON__P20930"       "CON__P34955"       "CON__P35527"       "CON__P35908"       "CON__P35908v2"    
# [31] "CON__P48668"       "CON__P78386"       "CON__Q04695"       "CON__Q29443"       "CON__Q14525"      
# [36] "CON__Q9UE12"       "CON__Q3SY84"       "CON__Q9Z2K1"       "CON__Q3ZBS7"       "CON__Q5D862"      
# [41] "CON__Q5XQN5"       "CON__Q9NSB2"       "CON__Q6KB66-1"     "CON__Q8N1N4-2"     "CON__Q7Z3Y8"      
# [46] "CON__Q7Z3Z0"       "CON__Q7Z794"       "CON__Q86YZ3"       "CON__Q8IUT8"       "CON__Q9C075"      
# [51] "CON__Q9R0H5"       "CON__Streptavidin" "Q8BNF3"            "Q9D002"            "Q9D4W6"           
# [56] "Q8CBU3"   
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
# [1] "Q1KYM0"                "A1E5T5"                "B2CSK2"                "CON__Q61726"          
# [5] "CON__Q922U2"           "CON__REFSEQ:XP_986630" "Q3TBT0"                "Q8C1Z9"               
# [9] "Q8BZ07"                "Q9CY10"                "Q3U4D1"                "Q8BRA9"               
# [13] "Q3UJM3"                "Q3U7K1"                "Q3V0Y9"                "Q61177"               
# [17] "Q6YIY0"                "Q8BW35"                "Q8C1L7"                "Q8C5X4"               
# [21] "Q99M08"                "Q9D454"                "Q9QYS0"    

for (uniprotID in mmLookups) rownames(exprMat0)[which(uniprotIDs.preferred==uniprotID)] <- mmLookup.df$UniqueID[which(rownames(mmLookup.df)==uniprotID)]
rownames(exprMat0)[mmLookup.rows.idx]
# [1] "Env polyprotein|Q1KYM0"                            "Tctex-1|A1E5T5"                                   
# [3] "Hsp70 1|B2CSK2"                                    "5430421N21Rik|CON__Q61726"                        
# [5] "Krt5|Q922U2"                                       "Krt33b|CON__REFSEQ:XP_986630"                     
# [7] "VWFA domain-containing protein|Q3TBT0"             "Tmed2|Q8C1Z9"                                     
# [9] "Signal recognition particle 19 kDa protein|Q8BZ07" "Hba-a1|Q9CY10"                                    
# [11] "Achy|Q3U4D1"                                       "UV excision repair protein RAD23 homolog A|Q8BRA9"
# [13] "Cox17|Q3UJM3"                                      "Ppp1cc|Q3U7K1"                                    
# [15] "Adam1a|Q3V0Y9"                                     "Csnk2a1|Q61177"                                   
# [17] "Mela; Gag|Q6YIY0"                                  "Acox1|Q8BW35"                                     
# [19] "Rps21|Q8C1L7"                                      "2-Hacid_dh|Q8C5X4"                                
# [21] "Uncharacterized protein C4orf3 homolog|Q99M08"     "4933412E24Rik|Q9D454"                             
# [23] "ATP sulfurylase/APS kinase isoform SK2|Q9QYS0" 

## Now remove human contaminants
dim(exprMat0)
# [1] 2837   16

exprMat0<-exprMat0[which(!rownames(exprMat0) %in% paste0("|",conLookups)), ]

dim(exprMat0)
# [1] 2777   16

##############################################################################
## Load traits (numericMeta)
numericMeta<-read.csv("Exp. 1-4_proteome traits.csv", header = T)
rownames(numericMeta) <- numericMeta[,1]

dat.sumInt<-exprMat0[,match(rownames(numericMeta),colnames(exprMat0))]

dat.sumInt.log<-log2(dat.sumInt)
dat.sumInt.log[!is.finite(dat.sumInt.log)]<-NA

dim(dat.sumInt.log)
# [1] 2777   16

## Finalize Grouping of Samples for t.test
Grouping<-numericMeta$group.simple
Grouping
# [1] "global"                 "LPS.global"             "global"                 "LPS.global"            
# [5] "global"                 "LPS.global"             "neg.pulldown"          "LPS.neg.pulldown"     
# [9] "neg.pulldown"          "LPS.neg.pulldown"      "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown"
# [13] "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown" "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown"

######################################################################################################
## Part 1: Data processing and clean up
# steps:
# 2. Normalization
  # a. pos.pulldown - neg.pulldown
  # b. globals and (pos.pulldown - neg.pulldown) median zero-centered
######################################################################################################
# d. negative pulldown background subtraction 
## Background (negative pulldown) subtraction
dat.sumInt<-2^dat.sumInt.log

table(Grouping)
#BV2T.LPS.pos.pulldown     BV2T.pos.pulldown      LPS.neg.pulldown             LPS.global          neg.pulldown                 global 
#                    3                     3                     2                     3                     2                     3

## From the above, select your group names for negative (.neg) and positive (.pos) in paired order to be used in the for loop within the apply below. global samples can be left out.
myGroups.neg=c("neg.pulldown","LPS.neg.pulldown")
myGroups.pos=c("BV2T.pos.pulldown","BV2T.LPS.pos.pulldown")

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
#[1] 160

# we can set these to 0, for imputation later  (note: log2(0)=-Inf)
dat.sumInt.bgSubtr[dat.sumInt.bgSubtr<0]<- 0
#dat.sumInt.bgSubtr[!is.finite(dat.sumInt.bgSubtr)]<- NA


## find rows for which no background(s) were available -- they contain NaN -- and put back the unsubtracted values for the whole row in each case
length(which(is.nan(dat.sumInt.bgSubtr)))
#[1] 2422  # individual values in data

sum(unlist(apply(dat.sumInt.bgSubtr,1,function(x) if(length(which(is.nan(x)))>0) 1)))
#[1] 1010  # rows with at least 1 NaN

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
# [1] 2777   16

dat.sumInt.bgSubtr.log<-log2(dat.sumInt.bgSubtr2)
dat.sumInt.bgSubtr.log[!is.finite(dat.sumInt.bgSubtr.log)]<-NA  # handles log2(0)

# subset to remove neg.pulldowns
#dat.sumInt.bgSubtr.nonegs.log <- dat.sumInt.bgSubtr.log[,c(1:6, 11:16)]
dat.sumInt.bgSubtr.nonegs.log <- dat.sumInt.bgSubtr.log[,which(!Grouping %in% myGroups.neg)]

dim(dat.sumInt.bgSubtr.nonegs.log)
# [1] 2777   12

######################################################################################################
# c. median samplewise zero centered normalization
## Simple median samplewise subtraction will center the distributions all at 0.

# Calculate median sample-wise
sampleMedianLFQ <- apply(dat.sumInt.bgSubtr.nonegs.log, 2, function(x) median(x, na.rm = TRUE))

# Calculate the mean of sample medians
meanMedian <- mean(sampleMedianLFQ)

# Apply median sample-wise subtraction
dat.sumInt.log.meanMedian<-sapply(colnames(dat.sumInt.bgSubtr.nonegs.log),function(x) dat.sumInt.bgSubtr.nonegs.log[,x] -sampleMedianLFQ[x] )  # keeps values positive: -(sampleMedianLFQ[x]-meanMedian) )


## QC Plot - step 2 - Unnorm Sum Intensity vs. Enforced Equal Loading Assumption WITHIN TREATMENT GROUP (samplewise summed intensity normalization) -- not MaxLFQ norm
pdf(file="2. Exp. 1-4 QC-runorder-SummedIntensity.Vs.meanMedian zero-centered norm_02122024.pdf",width=21,height=14)
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
boxplot(dat.sumInt.bgSubtr.log[,which(Grouping %in% myGroups.pos)], ylab=bquote("log"[2]~"MQ Summed Intensity"), main="b. Background Subtracted (Pulldowns Only)", col=colvec[which(Grouping %in% myGroups.pos)], las=2)
plot.new()

# adjust "Grouping" to remove neg.pulldowns
traits_nonegs <- numericMeta[which(!Grouping %in% myGroups.neg),]

Grouping_2<-traits_nonegs$group.simple
Grouping_2
# [1] "global"                 "LPS.global"             "global"                
# [4] "LPS.global"             "global"                 "LPS.global"            
# [7] "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown" "BV2T.pos.pulldown"    
# [10] "BV2T.LPS.pos.pulldown" "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown"

# colors for samples by group.simple
#colvec_2=WGCNA::labels2colors(as.numeric(factor(Grouping_2)))
#colvec.legend=WGCNA::labels2colors(as.numeric(factor(levels(factor(Grouping_2)))))
colvec_2=colvec[which(!Grouping %in% myGroups.neg)]

# boxplot 3. c. median samplewise zero centered
boxplot(dat.sumInt.log.meanMedian, ylab=bquote("log"[2]~"MQ Summed Intensity"), main="c. Median Samplewise (Zero-)Centered Protein Data", col=colvec_2, las=2)
plot.new()

dev.off()

dat.sumInt.log.norm <- dat.sumInt.log.meanMedian
dim(dat.sumInt.log.norm)
# [1] 2777   12

#EBD - not edited below

###############################################################################
###############################################################################
###############################################################################
# 3. Determine which rows to keep to control missing values
#    (Enforce missingness across samples with logical criteria)
###############################################################################
###############################################################################
###############################################################################
Grouping_2 # need the traits file w/o neg.pulldowns
# [1] "global"                 "LPS.global"             "global"                 "LPS.global"             "global"                
# [6] "LPS.global"             "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown" "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown"
# [11] "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown"

# missingness threshold for diffabund analysis: keep a protein if present in 2/3 samples/pos.pulldown group
# the input df does not have duplicate gene IDs collapsed
rownames.toKeep<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="BV2T.pos.pulldown")]))) >= 2 |
     length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="BV2T.LPS.pos.pulldown")]))) >= 2) x }))

# Note from Eric: at least 2 BV2T.pos.pulldown samples in each row kept are not NA

length(rownames.toKeep)
# [1] 1455

rownames.notKept_presentglobal<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(is.na(dat.sumInt.log.norm[x,which(Grouping_2=="BV2T.pos.pulldown")]) & 
                  !is.na(dat.sumInt.log.norm[x,which(grepl("global",Grouping_2))]))) >= 1 &
     length(which(is.na(dat.sumInt.log.norm[x,which(Grouping_2=="BV2T.LPS.pos.pulldown")]))) >= 2) x }))


length(rownames.notKept_presentglobal)
# [1] 1310

dat.sumInt.log.norm[rownames.notKept_presentglobal[1],] # should have values across globals
# Intensity 10tcl_01 Intensity 10tcl_02 Intensity 10tcl_03 Intensity 10tcl_04 Intensity 10tcl_05 Intensity 10tcl_06 
# 0.30589531         0.61564452         0.07279568         0.41894738        -2.49427235         0.29624871 
# Intensity 21ip_07  Intensity 21ip_08  Intensity 21ip_09  Intensity 21ip_10  Intensity 21ip_11  Intensity 21ip_12 
# NA                 NA                 NA                 NA                 NA                 NA 



rownames.toKeep_2<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="BV2T.pos.pulldown")]))) >= 1 |
     length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="BV2T.LPS.pos.pulldown")]))) >= 1) x }))

# Note from Eric: at least 2 BV2T.pos.pulldown samples in each row kept are not NA

length(rownames.toKeep_2)
# [1] 1834

###############################################################################
# missingness threshold for Venn diagrams: keep a protein if present in 2/3 samples/group
# the input df does not have duplicate gene IDs collapsed
# BV2T.pos.pulldown filter
rownames.toKeep_BV2T.pos.pulldown<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="BV2T.pos.pulldown")]))) >= 2 ) x }))

length(rownames.toKeep_BV2T.pos.pulldown)
# [1] 1140

cleanDat_BV2T_filtered<-dat.sumInt.log.norm[rownames.toKeep_BV2T.pos.pulldown,]

write.csv(cleanDat_BV2T_filtered, file=paste0("Exp. 1-4_LFQ-MS_BV2T.pos.pulldown_Venn filtered-",dim(cleanDat_BV2T_filtered)[1],"x",dim(cleanDat_BV2T_filtered)[2],".csv_02152024.csv"))

# BV2T.LPS.pos.pulldown filter
rownames.toKeep_BV2T.pos.pulldown.LPS.pos.pulldown<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="BV2T.LPS.pos.pulldown")]))) >= 2 ) x }))

length(rownames.toKeep_BV2T.pos.pulldown.LPS.pos.pulldown)
# [1] 1333

cleanDat_BV2T.LPS.pos.pulldown_filtered<-dat.sumInt.log.norm[rownames.toKeep_BV2T.pos.pulldown.LPS.pos.pulldown,]

write.csv(cleanDat_BV2T.LPS.pos.pulldown_filtered, file=paste0("Exp. 1-4_LFQ-MS_BV2T.LPS.pos.pulldown_Venn filtered-",dim(cleanDat_BV2T.LPS.pos.pulldown_filtered)[1],"x",dim(cleanDat_BV2T.LPS.pos.pulldown_filtered)[2],".csv_02152024.csv"))

# global filter
rownames.toKeep_global<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="global")]))) >= 2 ) x }))

length(rownames.toKeep_global)
# [1] 2493

cleanDat_input_filtered<-dat.sumInt.log.norm[rownames.toKeep_global,]

write.csv(cleanDat_input_filtered, file=paste0("Exp. 1-4_LFQ-MS_global_Venn filtered-",dim(cleanDat_input_filtered)[1],"x",dim(cleanDat_input_filtered)[2],".csv_02152024.csv"))

# LPS.global filter
rownames.toKeep_LPS.global<-unlist(sapply(rownames(dat.sumInt.log.norm),function(x) {
  if(length(which(!is.na(dat.sumInt.log.norm[x,which(Grouping_2=="LPS.global")]))) >= 2 ) x }))

length(rownames.toKeep_LPS.global)
# [1] 2521

cleanDat_LPS.global_filtered<-dat.sumInt.log.norm[rownames.toKeep_LPS.global,]

write.csv(cleanDat_LPS.global_filtered, file=paste0("Exp. 1-4_LFQ-MS_LPS.global_Venn filtered-",dim(cleanDat_LPS.global_filtered)[1],"x",dim(cleanDat_LPS.global_filtered)[2],".csv_02152024.csv"))


# file names:
# Exp. 1-4_LFQ-MS_global_Venn filtered-2493x12.csv_02152024.csv
# Exp. 1-4_LFQ-MS_LPS.global_Venn filtered-2521x12.csv_02152024.csv
# Exp. 1-4_LFQ-MS_BV2T.pos.pulldown_Venn filtered-1140x12.csv_02152024.csv
# Exp. 1-4_LFQ-MS_BV2T.LPS.pos.pulldown_Venn filtered-1333x12.csv_02152024.csv

#Load intensity files
global_protein<-read.csv("Exp. 1-4_LFQ-MS_global_Venn filtered-2493x12.csv_02152024.csv", header = TRUE)
LPS.global_protein<-read.csv("Exp. 1-4_LFQ-MS_LPS.global_Venn filtered-2521x12.csv_02152024.csv", header = TRUE)
BV2T.pos.pulldown_protein<-read.csv("Exp. 1-4_LFQ-MS_BV2T.pos.pulldown_Venn filtered-1140x12.csv_02152024.csv", header = TRUE)
BV2T.LPS.pos.pulldown_protein<-read.csv("Exp. 1-4_LFQ-MS_BV2T.LPS.pos.pulldown_Venn filtered-1333x12.csv_02152024.csv", header = TRUE)

## Subset the intensity matrices to extract only the protein ID columns to create Venn diagrams
global_protein_IDs <- data.frame(global_protein[,1])
colnames(global_protein_IDs)
colnames(global_protein_IDs) <- c("gene_protein_IDs")
library(tidyr)
global_protein_IDs_sep <- separate(global_protein_IDs, gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

LPS.global_protein_IDs <- data.frame(LPS.global_protein[,1])
colnames(LPS.global_protein_IDs)
colnames(LPS.global_protein_IDs) <- c("gene_protein_IDs")
LPS.global_protein_IDs_sep <- separate(LPS.global_protein_IDs, gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

BV2T.pos.pulldown_protein_IDs <-data.frame(BV2T.pos.pulldown_protein[,1])
colnames(BV2T.pos.pulldown_protein_IDs)
colnames(BV2T.pos.pulldown_protein_IDs) <- c("gene_protein_IDs")
BV2T.pos.pulldown_protein_IDs_sep <- separate(BV2T.pos.pulldown_protein_IDs, gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

BV2T.LPS.pos.pulldown_protein_IDs <- data.frame(BV2T.LPS.pos.pulldown_protein[,1])
colnames(BV2T.LPS.pos.pulldown_protein_IDs)
colnames(BV2T.LPS.pos.pulldown_protein_IDs) <- c("gene_protein_IDs")
BV2T.LPS.pos.pulldown_protein_IDs_sep <- separate(BV2T.LPS.pos.pulldown_protein_IDs, gene_protein_IDs, into = c("col1", "col2"), sep ="\\|")

global_protein_list <- global_protein_IDs_sep[,2]
LPS.global_protein_list <- LPS.global_protein_IDs_sep[,2]
BV2T.pos.pulldown_protein_list <-  BV2T.pos.pulldown_protein_IDs_sep[,2]
BV2T.LPS.pos.pulldown_protein_list <- BV2T.LPS.pos.pulldown_protein_IDs_sep[,2]
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
pdf(file = "Exp. 1-4_LFQ-MS_protein_Venns_CR_02152024.pdf", height = 11, width = 8.5, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# global (x) vs LPS.global (y)
biovenn1 <- draw.venn(global_protein_list, LPS.global_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="BV2T global vs BV2T LPS global protein", 
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

# [1] "x global: 2493"
# [1] "y global: 2521"
# [1] "z global: 0"
# [1] "x only: 54"
# [1] "y only: 82"
# [1] "z only: 0"
# [1] "x-y global overlap: 2439"


# global (x) vs BV2T pulldown (y)
biovenn2 <- draw.venn(global_protein_list, BV2T.pos.pulldown_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="BV2T global vs BV2T pulldown protein", 
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
                      x_c = "darkgrey",
                      y_c = "azure",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)
# [1] "x global: 2493"
# [1] "y global: 1140"
# [1] "z global: 0"
# [1] "x only: 1455"
# [1] "y only: 102"
# [1] "z only: 0"
# [1] "x-y global overlap: 1038"

# BV2T.pos.pulldown (x) vs BV2T.LPS.pos.pulldown (y)
biovenn3 <- draw.venn(BV2T.pos.pulldown_protein_list, BV2T.LPS.pos.pulldown_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="BV2T vs BV2T LPS pulldown protein", 
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
                      x_c = "azure3",
                      y_c = "darkorange",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x global: 1140"
# [1] "y global: 1333"
# [1] "z global: 0"
# [1] "x only: 122"
# [1] "y only: 315"
# [1] "z only: 0"
# [1] "x-y global overlap: 1018"

# LPS.global (x) vs BV2T.LPS pulldown (y)
biovenn4 <- draw.venn(LPS.global_protein_list, BV2T.LPS.pos.pulldown_protein_list, list_z,  
                      t_fb = 2,
                      t_s = 1.5,
                      t_c = "black",
                      t_f = "Arial",
                      title="BV2T LPS global vs BV2T LPS pulldown protein", 
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
                      x_c = "darkgrey",
                      y_c = "darkorange",
                      # z_c = "blue",
                      bg_c = "white",
                      width = 1000,
                      height = 1000,
                      #output = "screen",
                      filename = NULL,
                      map2ens = FALSE
)

# [1] "x global: 2521"
# [1] "y global: 1333"
# [1] "z global: 0"
# [1] "x only: 1289"
# [1] "y only: 101"
# [1] "z only: 0"
# [1] "x-y global overlap: 1232"

dev.off()


# Create .csv files for common proteins between comparisons
# load library
library(dplyr)

# Venn overlap values
# global (x) vs LPS.global (y)
  # [1] "x-y global overlap: 2439"

# global (x) vs BV2T pulldown (y)
# [1] "x-y global overlap: 1038"

# BV2T.pos.pulldown (x) vs BV2T.LPS.pos.pulldown (y)
# [1] "x-y global overlap: 1018"

# LPS.global (x) vs BV2T.LPS pulldown (y)
  # [1] "x-y global overlap: 1232"


# global & LPS.global common proteins
common_proteins<-inner_join(global_protein,LPS.global_protein)
dim(common_proteins)
# [1] 2439   13 CR: this matches the overlap number for the Venn
write.csv(common_proteins, file ="Exp. 1-4_LFQ-MS_global & global.LPS_common proteins_CR_02152024.csv")

# global & BV2T.pos.pulldown common proteins
common_proteins_3<-inner_join(global_protein,BV2T.pos.pulldown_protein)
dim(common_proteins_3)
# [1] 1038   17 CR: this matches the overlap number for the Venn
write.csv(common_proteins_3, file ="Exp. 1-4_LFQ-MS_global & BV2T.pull_common proteins_CR_02152024.csv")

# BV2T.pos.pulldown & BV2T.LPS.pos.pulldown common proteins
common_proteins_2<-inner_join(BV2T.pos.pulldown_protein,BV2T.LPS.pos.pulldown_protein)
dim(common_proteins_2)
# [1] 1018   13 CR: this matches the overlap number for the Venn
write.csv(common_proteins_2, file ="Exp. 1-4_LFQ-MS_BV2T.pull & BV2T.LPS.pull_common proteins_CR_02152024.csv")

# LPS.global & BV2T.LPS.pos.pulldown common proteins
common_proteins_4<-inner_join(LPS.global_protein,BV2T.LPS.pos.pulldown_protein)
dim(common_proteins_4)
# [1] 1232   13 CR: this matches the overlap number for the Venn
write.csv(common_proteins_4, file ="Exp. 1-4_LFQ-MS_global.LPS & BV2T.LPS.pull_common proteins_CR_02152024.csv")


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


# c. median samplewise zero centered
cleanDat.0centered.unfiltered.imputed<-imputePerseusStyle(dat.sumInt.log.meanMedian) #full data with imputation Perseus Style Column-wise Independent - all columns
head(cleanDat.0centered.unfiltered.imputed)

# write output with imputed values (log2-transformed), cleanDat.unfiltered.imputed
write.csv(cleanDat.0centered.unfiltered.imputed, file=paste0("4c.cleanDat.0centered.unfiltered.imputed-",dim(cleanDat.0centered.unfiltered.imputed)[1],"x",dim(cleanDat.0centered.unfiltered.imputed)[2],".csv_Perseus-Imputation_02132024.csv"))

dim(cleanDat)
# [1] 2777   12

# now includes imputed data points -- and we have just applied the filtering to keep only rows meeting criteria from step 3
cleanDat<-cleanDat.0centered.filtered.imputed<-cleanDat.0centered.unfiltered.imputed[rownames.toKeep,]
write.csv(cleanDat.0centered.filtered.imputed, file=paste0("4c.cleanDat.0centered.filtered.imputed-",dim(cleanDat.0centered.filtered.imputed)[1],"x",dim(cleanDat.0centered.filtered.imputed)[2],".csv_Perseus-Imputation_02132024.csv"))

dim(cleanDat)
# [1] 1455   12

# d. negative pulldown background subtraction
cleanDat.bgSubtr.unfiltered.imputed<-imputePerseusStyle(dat.sumInt.bgSubtr.log) #full data with imputation Perseus Style Column-wise Independent - all columns
head(cleanDat.bgSubtr.unfiltered.imputed)

# write output with imputed values (log2-transformed), cleanDat.unfiltered.imputed
write.csv(cleanDat.bgSubtr.unfiltered.imputed, file=paste0("4d.cleanDat.bgSubtr.unfiltered.imputed-",dim(cleanDat.bgSubtr.unfiltered.imputed)[1],"x",dim(cleanDat.bgSubtr.unfiltered.imputed)[2],".csv_Perseus-Imputation_02132024.csv"))

dim(cleanDat)
# [1] 1455   12

# now includes imputed data points -- and we have just applied the filtering to keep only rows meeting criteria from step 3
cleanDat<-cleanDat.bgSubtr.filtered.imputed<-cleanDat.bgSubtr.unfiltered.imputed[rownames.toKeep,]
write.csv(cleanDat.bgSubtr.filtered.imputed, file=paste0("4d.cleanDat.bgSubtr.filtered.imputed-",dim(cleanDat.bgSubtr.filtered.imputed)[1],"x",dim(cleanDat.bgSubtr.filtered.imputed)[2],".csv_Perseus-Imputation_02132024.csv"))

dim(cleanDat)
# [1] 1455   16

###############################################################################
## Part 1: Data processing and clean up
# steps:
# 3. Collapse duplicate gene IDs 
  # a. Keep a protein isoform with the most variance
###############################################################################
###############################################################################
# rename normalized intensity matrix for simplicity and subset to only include pos.pulldowns
# normalized, non-filtered and non-imputed intensity matrix as input
cleanDat <- dat.sumInt.log.meanMedian[,c(7:12)]
dim(cleanDat)
# [1] 2777   6

# check to see if any gene symbols are missing, should be 0
badRows<-which(grepl("^\\|",rownames(cleanDat)))
rownames(cleanDat)[badRows]
# character(0)

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
# [1] 1807    6

Grouping
# [1] "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown" "BV2T.pos.pulldown"     "BV2T.LPS.pos.pulldown" "BV2T.pos.pulldown"    
# [6] "BV2T.LPS.pos.pulldown"

rownames.toKeep<-unlist(sapply(rownames(cleanDat.collapsed),function(x) {
  if(length(which(!is.na(cleanDat.collapsed[x,which(Grouping=="BV2T.pos.pulldown")]))) >= 2 |
     length(which(!is.na(cleanDat.collapsed[x,which(Grouping=="BV2T.LPS.pos.pulldown")]))) >= 2) x }))

# Note from Eric: at least 2 BV2T.pos.pulldown samples in each row kept are not NA

length(rownames.toKeep)
# [1] 1436



###############################################################################
# normalized, filtered and imputed intensity matrix with pos.pulldowns as input
cleanDat <- cleanDat.0centered.filtered.imputed[,c(7:12)]
dim(cleanDat)
# [1] 1455    6

# The norm intensity matrix row names are labeled gene symbol|protein ID, so we
# need to extract the first part of each row name to obtain gene symbols and put into a df
symbols <- as.data.frame(do.call(rbind, strsplit(rownames(cleanDat), "[|]")))[, 1]
symbols <- as.data.frame(do.call(rbind, strsplit(symbols, "[;]")))[, 1]

# Collapse rows based on variance
collapsed_data_2 <- collapseRows(cleanDat, symbols, rownames(cleanDat), method = "MaxMean")

cleanDat.collapsed_2 <- collapsed_data_2$datETcollapsed
dim(cleanDat.collapsed_2)
# [1] 1436    6

###############################################################################
# a. normalized, filtered and imputed intensity matrix with global and pos.pulldowns as input
cleanDat <- cleanDat.0centered.filtered.imputed
dim(cleanDat)
# [1] 1455    12

# The norm intensity matrix row names are labeled gene symbol|protein ID, so we
# need to extract the first part of each row name to obtain gene symbols and put into a df
symbols <- as.data.frame(do.call(rbind, strsplit(rownames(cleanDat), "[|]")))[, 1]
symbols <- as.data.frame(do.call(rbind, strsplit(symbols, "[;]")))[, 1]

# Collapse rows based on variance
collapsed_data_3 <- collapseRows(cleanDat, symbols, rownames(cleanDat), method = "MaxMean")

cleanDat.collapsed_3 <- collapsed_data_3$datETcollapsed
dim(cleanDat.collapsed_3)
# [1] 1436    12

write.csv(cleanDat.collapsed_3, file=paste0("5a. cleanDat.bgSubtr.filtered.imputed.collapsed-",dim(cleanDat.collapsed_3)[1],"x",dim(cleanDat.collapsed_3)[2],".csv_Perseus-Imputation_02132024.csv"))

cleanDat <- cleanDat.0centered.filtered.imputed[,c(1:6)]
dim(cleanDat)
# [1] 1455    6

# The norm intensity matrix row names are labeled gene symbol|protein ID, so we
# need to extract the first part of each row name to obtain gene symbols and put into a df
symbols <- as.data.frame(do.call(rbind, strsplit(rownames(cleanDat), "[|]")))[, 1]
symbols <- as.data.frame(do.call(rbind, strsplit(symbols, "[;]")))[, 1]

# Collapse rows based on variance
collapsed_data_globals <- collapseRows(cleanDat, symbols, rownames(cleanDat), method = "MaxMean")

cleanDat.collapsed_globals <- collapsed_data_globals$datETcollapsed
dim(cleanDat.collapsed_globals)
# [1] 1436    6

intersect(which(collapsed_data_globals$selectedRow), which(!collapsed_data_2$selectedRow))
# [1] 946

intersect(which(!collapsed_data_globals$selectedRow), which(collapsed_data_2$selectedRow))
# [1] 946 row index ID for the gene retained in one df over the other

which(!collapsed_data_globals$selectedRow)
# Map4|A0A0G2JFH2 Map4|A0A140T8T5 Aak1|A0A571BEI2  Hnrnpa3|Q6P6I7   Ahnak2|E9PYB0    Capzb|F6YHZ8 
# 33              75             162             197             337             392 
# Ahnak2|F7CVJ5    Birc6|S4R1L5      Vim|Q3TFD9    Coro7|Q3TU50      Vim|Q3TWV0     Actb|Q3UBP6 
# 399             502             849             886             896             928 
# Ina|Q3UMG4    Haus6|Q6P252     Plec|Q6S393 Ubap2l|Q80X50-5     Cast|Q8CE04     Git2|Q9JLQ2 
# 969            1081            1098            1139            1194            1407 
# Plec|Q9QXS1-3 
# 1416 

which(!collapsed_data_2$selectedRow)
# Map4|A0A0G2JFH2 Map4|A0A140T8T5  Hnrnpa3|Q6P6I7   Ahnak2|E9PYB0    Capzb|F6YHZ8   Ahnak2|F7CVJ5 
# 33              75             197             337             392             399 
# Birc6|S4R1L5      Vim|Q3TFD9    Coro7|Q3TU50      Vim|Q3TWV0     Actb|Q3UBP6   Aak1|Q3UHJ0-2 
# 502             849             886             896             928             946 
# Ina|Q3UMG4    Haus6|Q6P252     Plec|Q6S393 Ubap2l|Q80X50-5     Cast|Q8CE04     Git2|Q9JLQ2 
# 969            1081            1098            1139            1194            1407 
# Plec|Q9QXS1-3 
# 1416 


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
# change the cleanDat to be the normalized intensity matrix
cleanDat <- dat.sumInt.log.norm
dim(cleanDat)
# [1] 2777   12

# load library
library(limma)

pdf(file = "1a. Exp. 1-4_LFQ-MS_all groups_nonegs_uncollapsed_Limma_PCA_CR_02132024.pdf", width = 8, height = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- read.csv("Exp. 1-4_Grouping vector.csv", header = TRUE)

pch_values <- ifelse(Grouping_vector == "global", 16,
                     ifelse(Grouping_vector == "LPS.global", 16,
                            ifelse(Grouping_vector == "BV2T.pos.pulldown", 16, 
                                   ifelse(Grouping_vector == "BV2T.LPS.pos.pulldown", 16, NA))))


pt_colors <- ifelse(Grouping_vector == "global", "orange",
                    ifelse(Grouping_vector == "LPS.global", "purple",
                           ifelse(Grouping_vector == "BV2T.pos.pulldown", "dodgerblue", 
                                  ifelse(Grouping_vector == "BV2T.LPS.pos.pulldown",  "magenta", NA))))

legend_groups <- c("global", "LPS.global", "BV2T.pos.pulldown","BV2T.LPS.pos.pulldown")

plotMDS_1_4_allgroups_nonegs <- plotMDS((cleanDat), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for uncollapsed log2(intensity)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("orange", "purple", "dodgerblue", "magenta"), 
       pch = 16, title = "Proteome groups",cex=1.4)

dev.off()


###############################################################################
# 1. PCAs of collapsed data

# load cleaned, filter, imputed and collapsed intensity matrix
cleanDat <- read.csv("5a. cleanDat.bgSubtr.filtered.imputed.collapsed-1436x12.csv_Perseus-Imputation_02132024.csv", header = TRUE, row.names = 1)
dim(cleanDat)
# [1] 1436   12 should match cleanDat.collapsed_3 CR: it does

pdf(file = "1b. Exp. 1-4_LFQ-MS_all groups_nonegs_collapsed_Limma_PCA_CR_02132024.pdf", width = 8, height = 10.5)
layout(matrix(c(1,1,2, 1,1,3, 4,5,6), nrow = 3, ncol = 3, byrow=TRUE),
       heights = c(1.3,0.7,1.3), # Heights of the rows
       widths = c(0.65,0.65,0.7)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

Grouping_vector <- read.csv("Exp. 1-4_Grouping vector.csv", header = TRUE)

pch_values <- ifelse(Grouping_vector == "global", 16,
                     ifelse(Grouping_vector == "LPS.global", 16,
                            ifelse(Grouping_vector == "BV2T.pos.pulldown", 16, 
                                   ifelse(Grouping_vector == "BV2T.LPS.pos.pulldown", 16, NA))))


pt_colors <- ifelse(Grouping_vector == "global", "orange",
                    ifelse(Grouping_vector == "LPS.global", "purple",
                           ifelse(Grouping_vector == "BV2T.pos.pulldown", "dodgerblue", 
                                  ifelse(Grouping_vector == "BV2T.LPS.pos.pulldown",  "magenta", NA))))

legend_groups <- c("global", "LPS.global", "BV2T.pos.pulldown","BV2T.LPS.pos.pulldown")

plotMDS_1_4_allgroups_nonegs <- plotMDS((cleanDat), #top = 500, 
                                        labels = NULL, pch = pch_values, col = pt_colors, 
                                        cex = 3, dim.plot = c(1,2), gene.selection = "pairwise",
                                        xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)

mtext(side=3, text="MDS Plot for collapsed log2(intensity)", cex=1, line=2)

plot.new()
legend("topright", legend = legend_groups, col = c("orange", "purple", "dodgerblue", "magenta"), 
       pch = 16, title = "Proteome groups",cex=1.4)

dev.off()



################################################################################
## Part 2: Data analysis and visualization
# 2. log2(intensity) correlations of uncollapsed df

# Step 2. Correlations of mRNA counts
# a. BV2T.global vs BV2T.LPS.global
# b. BV2T.global vs BV2T.pos.pulldown
# c. BV2T.LPS.global vs BV2T.LPS.pos.pulldown
# d. BV2T.LPS.pos.pulldown vs BV2T.pos.pulldown

# load library to generate plot
library(ggplot2)

# change the cleanDat to be the normalized, filtered and imputed intensity matrix
cleanDat <- read.csv("4c.cleanDat.0centered.filtered.imputed-1455x12.csv_Perseus-Imputation_02132024.csv", header = TRUE, row.names = 1)
dim(cleanDat)
# [1] 1455   12

################################################################################
# 2a. BV2T.global vs BV2T.LPS.global

pdf(file = "2a-d. Exp. 1-4_LFQ-MS_intensity correlations_uncollapsed_CR_02132024.pdf", width = 6, height = 6, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# subset into groups of interest and take the row mean
global <- rowMeans(cleanDat[,c(1,3,5)])
length(global)
# [1] 1455

LPS.global <- rowMeans(cleanDat[,c(2,4,6)])
length(LPS.global)
# [1] 1455

# generate a linear correlation
correlation_2a <- cor(global, LPS.global)

# convert data to a data frame
df_2a <- data.frame(global = global, LPS.global = LPS.global)

# create the correlation plot
correlation_plot_2a <- ggplot(df_2a, aes(x = global, y = LPS.global)) +
  geom_point(color = "black", fill = "darkorange", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "global log2(intensity)",
       y = "LPS global log2(intensity)",
       title = "Linear correlation for uncollapsed intensities",
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
  scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

# Print the plot
print(correlation_plot_2a)


################################################################################
# 2b.  BV2T.global vs BV2T.pos.pulldown
# subset into groups of interest and take the row mean
global <- rowMeans(cleanDat[,c(1,3,5)])
BV2T.pos.pulldown <- rowMeans(cleanDat[,c(7,9,11)])

# convert data to a data frame
df_2b <- data.frame(global = global, BV2T.pos.pulldown = BV2T.pos.pulldown)

# calculate correlation coefficient
correlation_2b <- cor(global, BV2T.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others

correlation_plot_2b <- ggplot(df_2b, aes(x = global, y = BV2T.pos.pulldown)) +
  geom_point(color = "black", fill = "darkorange", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "global log2(intensity)",
       y = "BV2T pulldown log2(intensity)",
       title = "Linear correlation for uncollapsed intensities",
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
  scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

print(correlation_plot_2b)


################################################################################
# 2c. BV2T.LPS.global vs BV2T.LPS.pos.pulldown
# subset into groups of interest and take the row mean
LPS.global <- rowMeans(cleanDat[,c(2,4,6)])
BV2T.LPS.pos.pulldown <- rowMeans(cleanDat[,c(8,10,12)])

# convert data to a data frame
df_2c <- data.frame(LPS.global = LPS.global, BV2T.LPS.pos.pulldown = BV2T.LPS.pos.pulldown)

# calculate correlation coefficient
correlation_2c <- cor(LPS.global, BV2T.LPS.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others

correlation_plot_2c <- ggplot(df_2c, aes(x = LPS.global, y = BV2T.LPS.pos.pulldown)) +
  geom_point(color = "black", fill = "darkorange", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "LPS global log2(intensity)",
       y = "BV2T LPS pulldown log2(intensity)",
       title = "Linear correlation for uncollapsed intensities",
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
  scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

print(correlation_plot_2c)


################################################################################
# 2d. BV2T.LPS.pos.pulldown vs BV2T.pos.pulldown
# subset into groups of interest and take the row mean
BV2T.pos.pulldown <- rowMeans(cleanDat[,c(7,9,11)])
BV2T.LPS.pos.pulldown <- rowMeans(cleanDat[,c(8,10,12)])

# convert data to a data frame
df_2d <- data.frame(BV2T.pos.pulldown = BV2T.pos.pulldown, BV2T.LPS.pos.pulldown  = BV2T.LPS.pos.pulldown)

# calculate correlation coefficient
correlation_2d <- cor(BV2T.pos.pulldown, BV2T.LPS.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others

correlation_plot_2d <- ggplot(df_2d, aes(x = BV2T.pos.pulldown, y = BV2T.LPS.pos.pulldown)) +
  geom_point(color = "black", fill = "darkorange", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "BV2T pulldown log2(intensity)",
       y = "BV2T LPS pulldown log2(intensity)",
       title = "Linear correlation for uncollapsed intensities",
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
  scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

print(correlation_plot_2d)

dev.off()

###############################################################################
# log2(intensity) correlations of collapsed intensity df
# change the cleanDat
cleanDat <- cleanDat.collapsed_3
dim(cleanDat)
# [1] 1436   12

pdf(file = "2a-d. Exp. 1-4__LFQ-MS_intensity correlations_collapsed_CR_02132024.pdf", width = 6, height = 6, family = "Arial")
par(mfrow= c(2,1)) # centers plot on single page (I think)

# subset into groups of interest and take the row mean
global <- rowMeans(cleanDat[,c(1,3,5)])
length(global)
# [1] 1436

LPS.global <- rowMeans(cleanDat[,c(2,4,6)])
length(LPS.global)
# [1] 1436

# generate a linear correlation
correlation_2a <- cor(global, LPS.global)

# convert data to a data frame
df_2a <- data.frame(global = global, LPS.global = LPS.global)

# create the correlation plot
correlation_plot_2a <- ggplot(df_2a, aes(x = global, y = LPS.global)) +
  geom_point(color = "black", fill = "darkorange", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "global log2(intensity)",
       y = "LPS global log2(intensity)",
       title = "Linear correlation for collapsed intensities",
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
  scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

# Print the plot
print(correlation_plot_2a)


################################################################################
# 2b.  BV2T.global vs BV2T.pos.pulldown
# subset into groups of interest and take the row mean
global <- rowMeans(cleanDat[,c(1,3,5)])
BV2T.pos.pulldown <- rowMeans(cleanDat[,c(7,9,11)])

# convert data to a data frame
df_2b <- data.frame(global = global, BV2T.pos.pulldown = BV2T.pos.pulldown)

# calculate correlation coefficient
correlation_2b <- cor(global, BV2T.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others

correlation_plot_2b <- ggplot(df_2b, aes(x = global, y = BV2T.pos.pulldown)) +
  geom_point(color = "black", fill = "darkorange", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "global log2(intensity)",
       y = "BV2T pulldown log2(intensity)",
       title = "Linear correlation for collapsed intensities",
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
  scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

print(correlation_plot_2b)


################################################################################
# 2c. BV2T.LPS.global vs BV2T.LPS.pos.pulldown
# subset into groups of interest and take the row mean
LPS.global <- rowMeans(cleanDat[,c(2,4,6)])
BV2T.LPS.pos.pulldown <- rowMeans(cleanDat[,c(8,10,12)])

# convert data to a data frame
df_2c <- data.frame(LPS.global = LPS.global, BV2T.LPS.pos.pulldown = BV2T.LPS.pos.pulldown)

# calculate correlation coefficient
correlation_2c <- cor(LPS.global, BV2T.LPS.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others

correlation_plot_2c <- ggplot(df_2c, aes(x = LPS.global, y = BV2T.LPS.pos.pulldown)) +
  geom_point(color = "black", fill = "darkorange", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "LPS global log2(intensity)",
       y = "BV2T LPS pulldown log2(intensity)",
       title = "Linear correlation for collapsed intensities",
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
  scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

print(correlation_plot_2c)

dev.off()

################################################################################
# 2d. BV2T.LPS.pos.pulldown vs BV2T.pos.pulldown
# subset into groups of interest and take the row mean
BV2T.pos.pulldown <- rowMeans(cleanDat[,c(7,9,11)])
BV2T.LPS.pos.pulldown <- rowMeans(cleanDat[,c(8,10,12)])

# convert data to a data frame
df_2d <- data.frame(BV2T.pos.pulldown = BV2T.pos.pulldown, BV2T.LPS.pos.pulldown  = BV2T.LPS.pos.pulldown )

# calculate correlation coefficient
correlation_2d <- cor(BV2T.pos.pulldown, BV2T.LPS.pos.pulldown) # cor() is based on Pearson's correlation by default, can change if needed by adding method = "spearman" or others

correlation_plot_2d <- ggplot(df_2d, aes(x = BV2T.pos.pulldown, y = BV2T.LPS.pos.pulldown)) +
  geom_point(color = "black", fill = "darkorange", pch = 21, alpha = 0.5, size = 3) +
  labs(x = "BV2T pulldown log2(intensity)",
       y = "BV2T LPS pulldown log2(intensity)",
       title = "Linear correlation for collapsed intensities",
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
  scale_x_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 12, by = 6), limits = c(0, 12))

print(correlation_plot_2d)

dev.off()

################################################################################
## Part 2: Data analysis and visualization 
# 3. Volcanoes of collapsed intensity df

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
# diffabund for collapsed intensity df
# update numericMeta
numericMeta <- traits_nonegs
dim(numericMeta)
# [1] 12  4

# comparisons:
# 1: BV2T LPS global vs global
#    BV2TLPSglobalvsglobal
# 2: BV2T pulldown vs global
#    BV2Tpulldownvsglobal
# 3: BV2T LPS pulldown vs LPS global
#    BV2TLPSpulldownvsLPSglobal
# 4: BV2T LPS pulldown vs non-LPS pulldown
#    BV2TLPSpulldownvsnonLPSpulldown

## c. median samplewise zero centered
## 5c1  BV2T LPS global vs global - 0 centered norm
cleanDat<-cleanDat.collapsed_3[,which(numericMeta$group.simple=="LPS.global" | numericMeta$group.simple=="global")]
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="LPS.global" | numericMeta$group.simple=="global")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="5c1"
outFileSuffix=" BV2TLPSglobalvsglobal.0centered_collapsed"
ANOVAout.BV2TLPSglobalvsglobal.0centered <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required global.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0)
plotVolc()  #plots volcano using ANOVAout

## 5c2 repeat for BV2T pulldown vs global- 0 centered norm
cleanDat<-cleanDat.collapsed_3[,which(numericMeta$group.simple=="BV2T.pos.pulldown" | numericMeta$group.simple=="global")]
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="BV2T.pos.pulldown" | numericMeta$group.simple=="global")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="5c2"
outFileSuffix=" BV2Tpulldownvsglobal.0centered_collapsed"
ANOVAout.BV2Tpulldownvsglobal.0centered <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required global.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0)
plotVolc()  #plots volcano using ANOVAout

## 5c3 repeat for BV2T LPS pulldown vs LPS global - 0 centered norm
cleanDat<-cleanDat.collapsed_3[,which(numericMeta$group.simple=="BV2T.LPS.pos.pulldown" | numericMeta$group.simple=="LPS.global")]
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="BV2T.LPS.pos.pulldown" | numericMeta$group.simple=="LPS.global")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="5c3"
outFileSuffix=" BV2TLPSpulldownvsLPSglobal.0centered_collapsed"
ANOVAout.BV2TLPSpulldownvsLPSglobal.0centered <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required global.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0)
plotVolc()  #plots volcano using ANOVAout

## 5c4 repeat for BV2T LPS pulldown vs non-LPS pulldown - 0 centered norm
cleanDat<-cleanDat.collapsed_3[,which(numericMeta$group.simple=="BV2T.LPS.pos.pulldown" | numericMeta$group.simple=="BV2T.pos.pulldown")]
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="BV2T.LPS.pos.pulldown" | numericMeta$group.simple=="BV2T.pos.pulldown")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="5c4"
outFileSuffix=" BV2TLPSpulldownvsnonLPSpulldown.0centered_collapsed"
ANOVAout.BV2TLPSpulldownvsnonLPSpulldown.0centered <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required global.
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
#globalFile <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#globalFile <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
            #global CSV FILE - in the filePath folder.
            #Can be formatted as Kme table from WGCNA pipeline, or
            #can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
            #in this case, the longest list is used as background or the "universe" for the FET contingencies
            #  For simple columnwise list global, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.

filePath <- getwd() 
            #Folder that (may) contain the global file specified above, and which will contain the outFilename project Folder.
modulesInMemory=FALSE
ANOVAgroups=TRUE  #change to FALSE if using globalFile
            #if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
            #globalFile will be ignored
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
# c. median samplewise zero centered

## 4c1 (BV2-TurboID LPS vs no LPS)pulldown DEP list GSEA
flip=c(0)
outFilename <- "4c1.BV2TLPSglobalvsglobal.0centered_collapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.BV2TLPSglobalvsglobal.0centered)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all globals available.

## 4c2 BV2-TurboID pulldown vs global DEP list GSEA
flip=c(3)
outFilename <- "4c2.BV2Tpulldownvsglobal.0centered_collapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.BV2Tpulldownvsglobal.0centered)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all globals available.

## 4c3 BV2-TurboID LPS pulldown vs LPS global DEP list GSEA
flip=c(3)
outFilename <- "4c3.BV2TLPSpulldownvsLPSglobal.0centered_collapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.BV2TLPSpulldownvsLPSglobal.0centered)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all globals available.

## 4c4 BV2-TurboID LPS pulldown vs LPS global DEP list GSEA
flip=c(3)
outFilename <- "4c4.BV2TLPSpulldownvsnonLPSpulldown.0centered_collapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.BV2TLPSpulldownvsnonLPSpulldown.0centered)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all globals available.


################################################################################
################################################################################
# diffabundance analysis and GOPar of uncollapsed intensity matrix

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
# diffabund for collapsed intensity df
# update numericMeta
numericMeta <- traits_nonegs
dim(numericMeta)
# [1] 12  4

# update cleanDat
uncollapsed.cleanDat <- read.csv("4c.cleanDat.0centered.filtered.imputed-1455x12.csv_Perseus-Imputation_02132024.csv", header = TRUE, row.names = 1)
dim(cleanDat)
# [1] 1455   12

# comparisons:
# 1: BV2T LPS global vs global
#    BV2TLPSglobalvsglobal
# 2: BV2T pulldown vs global
#    BV2Tpulldownvsglobal
# 3: BV2T LPS pulldown vs LPS global
#    BV2TLPSpulldownvsLPSglobal
# 4: BV2T LPS pulldown vs non-LPS pulldown
#    BV2TLPSpulldownvsnonLPSpulldown

## c. median samplewise zero centered
## 5c1  BV2T LPS global vs global - 0 centered norm
cleanDat<-uncollapsed.cleanDat[,which(numericMeta$group.simple=="LPS.global" | numericMeta$group.simple=="global")]
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="LPS.global" | numericMeta$group.simple=="global")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="5c1"
outFileSuffix=" BV2TLPSglobalvsglobal.0centered_uncollapsed"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
ANOVAout.BV2TLPSglobalvsglobal.0centered <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required global.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0)
plotVolc()  #plots volcano using ANOVAout

# kept incase I need to run down the line
# library(permFDP)
# myDesign <- numericMeta$group.simple[match(colnames(cleanDat), rownames(numericMeta))]
# myDesign[myDesign=="global"]<-1
# myDesign[myDesign=="LPS.global"]<-2
# c1permAdjp <- permFDP.adjust.threshold(ANOVAout.BV2TLPSglobalvsglobal.0centered[,3], threshold = 0.05, 
#                                       as.numeric(myDesign), intOnly = cleanDat, nPerms = 1000)
# signifP=c1permAdjp
# flip=c(0)
# plotVolc()  #plots volcano using ANOVAout
# c1permAdjp
#run GoPar for each comparison after the volcano since the adj p will be different for each comparison

## 5c2 repeat for BV2T pulldown vs global- 0 centered norm
cleanDat<-uncollapsed.cleanDat[,which(numericMeta$group.simple=="BV2T.pos.pulldown" | numericMeta$group.simple=="global")]
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="BV2T.pos.pulldown" | numericMeta$group.simple=="global")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="5c2"
outFileSuffix=" BV2Tpulldownvsglobal.0centered_uncollapsed"
ANOVAout.BV2Tpulldownvsglobal.0centered <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required global.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0)
plotVolc()  #plots volcano using ANOVAout

## 5c3 repeat for BV2T LPS pulldown vs LPS global - 0 centered norm
cleanDat<-uncollapsed.cleanDat[,which(numericMeta$group.simple=="BV2T.LPS.pos.pulldown" | numericMeta$group.simple=="LPS.global")]
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="BV2T.LPS.pos.pulldown" | numericMeta$group.simple=="LPS.global")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="5c3"
outFileSuffix=" BV2TLPSpulldownvsLPSglobal.0centered_uncollapsed"
ANOVAout.BV2TLPSpulldownvsLPSglobal.0centered <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required global.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0)
plotVolc()  #plots volcano using ANOVAout

## 5c4 repeat for BV2T LPS pulldown vs non-LPS pulldown - 0 centered norm
cleanDat<-uncollapsed.cleanDat[,which(numericMeta$group.simple=="BV2T.LPS.pos.pulldown" | numericMeta$group.simple=="BV2T.pos.pulldown")]
Grouping=numericMeta$group.simple[which(numericMeta$group.simple=="BV2T.LPS.pos.pulldown" | numericMeta$group.simple=="BV2T.pos.pulldown")]                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.

outFilePrefix="5c4"
outFileSuffix=" BV2TLPSpulldownvsnonLPSpulldown.0centered_uncollapsed"
ANOVAout.BV2TLPSpulldownvsnonLPSpulldown.0centered <- ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required global.
#...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: 0 [0%]

flip=c(0)
plotVolc()  #plots volcano using ANOVAout

################################################################################
## Part 2: Data analysis and visualization
# 4. GOparallel GSEA on volcano-defined DEP lists from uncollapsed intensity df
################################################################################
source("GOparallel-FET.R")

## Parameters saved to variables in memory.

######################## EDIT THESE VARIABLES (USER PARAMETERS SET IN GLOBAL ENVIRONMENT) ############################################
#globalFile <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#globalFile <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
#global CSV FILE - in the filePath folder.
#Can be formatted as Kme table from WGCNA pipeline, or
#can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
#in this case, the longest list is used as background or the "universe" for the FET contingencies
#  For simple columnwise list global, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.

filePath <- getwd() 
#Folder that (may) contain the global file specified above, and which will contain the outFilename project Folder.
modulesInMemory=FALSE
ANOVAgroups=TRUE  #change to FALSE if using globalFile
#if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
#globalFile will be ignored
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
# c. median samplewise zero centered

## 4c1 (BV2-TurboID LPS vs no LPS)pulldown DEP list GSEA
flip=c(0)
outFilename <- "4c1.BV2TLPSglobalvsglobal.0centered_uncollapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.BV2TLPSglobalvsglobal.0centered)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all globals available.

## 4c2 BV2-TurboID pulldown vs global DEP list GSEA
flip=c(3)
outFilename <- "4c2.BV2Tpulldownvsglobal.0centered_uncollapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.BV2Tpulldownvsglobal.0centered)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all globals available.

## 4c3 BV2-TurboID LPS pulldown vs LPS global DEP list GSEA
flip=c(3)
outFilename <- "4c3.BV2TLPSpulldownvsLPSglobal.0centered_uncollapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.BV2TLPSpulldownvsLPSglobal.0centered)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all globals available.

## 4c4 BV2-TurboID LPS pulldown vs LPS global DEP list GSEA
flip=c(3)
outFilename <- "4c4.BV2TLPSpulldownvsnonLPSpulldown.0centered_uncollapsed.DEPs-GO"
GMTdatabaseFile="Mouse_GO_AllPathways_with_GO_iea_November_07_2023_symbol.gmt"
GOparallel(ANOVAout.BV2TLPSpulldownvsnonLPSpulldown.0centered)  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all globals available.


# end of data analysis and visualization

