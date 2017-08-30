####################### Average LD Decay ###############################
# February 2017
# author: Annabel Beichman
# contact: annabel.beichman@gmail.com
# purpose: This script will average r2 values within bins of physical distance 
# by dividing the sum of all r2 values in that bin by the sum of SNP pairs in that bin

#################### GUTENKUNST ########################
ceu_gut_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/gutenkunst.ceu.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(ceu_gut_in) # 2000000 4
colnames(ceu_gut_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
ceu_gut_avg <- aggregate(ceu_gut_in[,3:4],by=list(factor(ceu_gut_in$binStop)),sum)

# now get average:
ceu_gut_avg$avgR2 <- ceu_gut_avg$sumR2/ceu_gut_avg$sumPairs
colnames(ceu_gut_avg) <- c("binStop","sumR2","sumPairs","avgR2")


###################### MSMC: ###############################################
##################### CEU 2 ##############################################
ceu_msmc_2_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/msmc.ceu_2.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(ceu_msmc_2_in) # 2000000 4
head(ceu_msmc_2_in) # checked against hoffman. CORRECT.
colnames(ceu_msmc_2_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
ceu_msmc_2_avg <- aggregate(ceu_msmc_2_in[,3:4],by=list(factor(ceu_msmc_2_in$binStop)),sum)
head(ceu_msmc_2_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
ceu_msmc_2_avg$avgR2 <- ceu_msmc_2_avg$sumR2/ceu_msmc_2_avg$sumPairs
colnames(ceu_msmc_2_avg) <- c("binStop","sumR2","sumPairs","avgR2")

##################### CEU 4 ##############################################
ceu_msmc_4_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/msmc.ceu_4.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(ceu_msmc_4_in) # 500000 4
head(ceu_msmc_4_in) # checked against hoffman. CORRECT.
colnames(ceu_msmc_4_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
ceu_msmc_4_avg <- aggregate(ceu_msmc_4_in[,3:4],by=list(factor(ceu_msmc_4_in$binStop)),sum)
head(ceu_msmc_4_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
ceu_msmc_4_avg$avgR2 <- ceu_msmc_4_avg$sumR2/ceu_msmc_4_avg$sumPairs
colnames(ceu_msmc_4_avg) <- c("binStop","sumR2","sumPairs","avgR2")


########## RECALL THAT THIS IS ONLY 5000 blocks not 20000

##################### CEU 8 ##############################################
ceu_msmc_8_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/msmc.ceu_8.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(ceu_msmc_8_in) # 2000000 4
head(ceu_msmc_8_in) # checked against hoffman. CORRECT.
colnames(ceu_msmc_8_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
ceu_msmc_8_avg <- aggregate(ceu_msmc_8_in[,3:4],by=list(factor(ceu_msmc_8_in$binStop)),sum)
head(ceu_msmc_8_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
ceu_msmc_8_avg$avgR2 <- ceu_msmc_8_avg$sumR2/ceu_msmc_8_avg$sumPairs
colnames(ceu_msmc_8_avg) <- c("binStop","sumR2","sumPairs","avgR2")

#################################### SMC++ ###################################
#ceu_smcpp_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/smcpp.ceu.ALL.LD.BINNED.OUTPUT.20170211.txt")
####**NEW** smcpp with time points: 20170214
ceu_smcpp_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/smcpp.ceu.ALL.LD.BINNED.OUTPUT.20170215.txt")

dim(ceu_smcpp_in) # 2000000 4
head(ceu_smcpp_in) # checked against hoffman. CORRECT.
colnames(ceu_smcpp_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
ceu_smcpp_avg <- aggregate(ceu_smcpp_in[,3:4],by=list(factor(ceu_smcpp_in$binStop)),sum)
head(ceu_smcpp_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
ceu_smcpp_avg$avgR2 <- ceu_smcpp_avg$sumR2/ceu_smcpp_avg$sumPairs
colnames(ceu_smcpp_avg) <- c("binStop","sumR2","sumPairs","avgR2")




##################################################################################################
##################################### CHB #######################################################
##################################################################################################
####################### Plottind LD Decay From Scratch ###############################

#################### GUTENKUNST ########################
chb_gut_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/gutenkunst.chb.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(chb_gut_in) # 2000000 4
head(chb_gut_in) # checked against hoffman. CORRECT.
colnames(chb_gut_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
chb_gut_avg <- aggregate(chb_gut_in[,3:4],by=list(factor(chb_gut_in$binStop)),sum)
head(chb_gut_avg) ### CHECK THIS SUPER CAREFULLY 

#   0 3113717  5402412 *see below * CORRECT
# test:
sum(chb_gut_in[chb_gut_in$binStop==1000,]$sumR2) # 3113717 CORRECT test what first bin of R2 aggregate should be
sum(chb_gut_in[chb_gut_in$binStop==1000,]$sumPairs) # 5402412 CORRECT test what first bin of Pairs aggregate should be
tail(chb_gut_avg)
#  95000 58343.464   203060
sum(chb_gut_in[chb_gut_in$binStop==96000,]$sumR2) #  58343.46 CORRECT test what first bin of R2 aggregate should be
sum(chb_gut_in[chb_gut_in$binStop==96000,]$sumPairs) # 203060 CORRECT test what first bin of Pairs aggregate should be
# SO AGGREGATE WORKS THE WAY I THINK IT DOES 

# now get average:
chb_gut_avg$avgR2 <- chb_gut_avg$sumR2/chb_gut_avg$sumPairs
colnames(chb_gut_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(chb_gut_avg)
chb_gut_forPlot <- chb_gut_avg[,c(1,4)]
head(chb_gut_forPlot)
# add label:
chb_gut_forPlot$label <- "2. Gutenkunst "
chb_gut_forPlot$shape <- "average"

###################### MSMC: ###############################################
##################### CHB 2 ##############################################
chb_msmc_2_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/msmc.chb_2.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(chb_msmc_2_in) # 2000000 4
head(chb_msmc_2_in) # checked against hoffman. CORRECT.
colnames(chb_msmc_2_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
chb_msmc_2_avg <- aggregate(chb_msmc_2_in[,3:4],by=list(factor(chb_msmc_2_in$binStop)),sum)
head(chb_msmc_2_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
chb_msmc_2_avg$avgR2 <- chb_msmc_2_avg$sumR2/chb_msmc_2_avg$sumPairs
colnames(chb_msmc_2_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(chb_msmc_2_avg)
chb_msmc_2_forPlot <- chb_msmc_2_avg[,c(1,4)]
head(chb_msmc_2_forPlot)
# add label:
chb_msmc_2_forPlot$label <- "3. MSMC 2 Hap "
chb_msmc_2_forPlot$shape <- "average"
##################### CHB 4 ##############################################
chb_msmc_4_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/msmc.chb_4.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(chb_msmc_4_in) # 500000 4
head(chb_msmc_4_in) # checked against hoffman. CORRECT.
colnames(chb_msmc_4_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
chb_msmc_4_avg <- aggregate(chb_msmc_4_in[,3:4],by=list(factor(chb_msmc_4_in$binStop)),sum)
head(chb_msmc_4_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
chb_msmc_4_avg$avgR2 <- chb_msmc_4_avg$sumR2/chb_msmc_4_avg$sumPairs
colnames(chb_msmc_4_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(chb_msmc_4_avg)
chb_msmc_4_forPlot <- chb_msmc_4_avg[,c(1,4)]
head(chb_msmc_4_forPlot)
# add label:
chb_msmc_4_forPlot$label <- "4. MSMC 4 Hap "
chb_msmc_4_forPlot$shape <- "average"
########## RECALL THAT THIS IS ONLY 5000 blocks not 20000

##################### CHB 8 ##############################################
chb_msmc_8_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/msmc.chb_8.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(chb_msmc_8_in) # 2000000 4
head(chb_msmc_8_in) # checked against hoffman. CORRECT.
colnames(chb_msmc_8_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
chb_msmc_8_avg <- aggregate(chb_msmc_8_in[,3:4],by=list(factor(chb_msmc_8_in$binStop)),sum)
head(chb_msmc_8_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
chb_msmc_8_avg$avgR2 <- chb_msmc_8_avg$sumR2/chb_msmc_8_avg$sumPairs
colnames(chb_msmc_8_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(chb_msmc_8_avg)
chb_msmc_8_forPlot <- chb_msmc_8_avg[,c(1,4)]
head(chb_msmc_8_forPlot)
# add label:
chb_msmc_8_forPlot$label <- "5. MSMC 8 Hap "
chb_msmc_8_forPlot$shape <- "average"
#################################### SMC++ ###################################
### NEW smcpp: 20170214
chb_smcpp_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/smcpp.chb.ALL.LD.BINNED.OUTPUT.20170215.txt")
dim(chb_smcpp_in) # 2000000 4
head(chb_smcpp_in) # checked against hoffman. CORRECT.
colnames(chb_smcpp_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
chb_smcpp_avg <- aggregate(chb_smcpp_in[,3:4],by=list(factor(chb_smcpp_in$binStop)),sum)
head(chb_smcpp_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
chb_smcpp_avg$avgR2 <- chb_smcpp_avg$sumR2/chb_smcpp_avg$sumPairs
colnames(chb_smcpp_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(chb_smcpp_avg)
chb_smcpp_forPlot <- chb_smcpp_avg[,c(1,4)]
head(chb_smcpp_forPlot)
# add label:
chb_smcpp_forPlot$label <- "6. SMC++ "
chb_smcpp_forPlot$shape <- "average"

##################################################################################################
##################################### YRI #######################################################
##################################################################################################
####################### Plottind LD Decay From Scratch ###############################


#################### GUTENKUNST ########################
yri_gut_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/gutenkunst.yri.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(yri_gut_in) # 2000000 4
head(yri_gut_in) # checked against hoffman. CORRECT.
colnames(yri_gut_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
yri_gut_avg <- aggregate(yri_gut_in[,3:4],by=list(factor(yri_gut_in$binStop)),sum)
head(yri_gut_avg) ### CHECK THIS SUPER CAREFULLY 

#   0 3113717  5402412 *see below * CORRECT
# test:
sum(yri_gut_in[yri_gut_in$binStop==1000,]$sumR2) # 3113717 CORRECT test what first bin of R2 aggregate should be
sum(yri_gut_in[yri_gut_in$binStop==1000,]$sumPairs) # 5402412 CORRECT test what first bin of Pairs aggregate should be
tail(yri_gut_avg)
#  95000 58343.464   203060
sum(yri_gut_in[yri_gut_in$binStop==96000,]$sumR2) #  58343.46 CORRECT test what first bin of R2 aggregate should be
sum(yri_gut_in[yri_gut_in$binStop==96000,]$sumPairs) # 203060 CORRECT test what first bin of Pairs aggregate should be
# SO AGGREGATE WORKS THE WAY I THINK IT DOES 

# now get average:
yri_gut_avg$avgR2 <- yri_gut_avg$sumR2/yri_gut_avg$sumPairs
colnames(yri_gut_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(yri_gut_avg)
yri_gut_forPlot <- yri_gut_avg[,c(1,4)]
head(yri_gut_forPlot)
# add label:
yri_gut_forPlot$label <- "2. Gutenkunst "
yri_gut_forPlot$shape <- "average"

###################### MSMC: ###############################################
##################### YRI 2 ##############################################
yri_msmc_2_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/msmc.yri_2.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(yri_msmc_2_in) # 2000000 4
head(yri_msmc_2_in) # checked against hoffman. CORRECT.
colnames(yri_msmc_2_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
yri_msmc_2_avg <- aggregate(yri_msmc_2_in[,3:4],by=list(factor(yri_msmc_2_in$binStop)),sum)
head(yri_msmc_2_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
yri_msmc_2_avg$avgR2 <- yri_msmc_2_avg$sumR2/yri_msmc_2_avg$sumPairs
colnames(yri_msmc_2_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(yri_msmc_2_avg)
yri_msmc_2_forPlot <- yri_msmc_2_avg[,c(1,4)]
head(yri_msmc_2_forPlot)
# add label:
yri_msmc_2_forPlot$label <- "3. MSMC 2 Hap "
yri_msmc_2_forPlot$shape <- "average"
##################### YRI 4 ##############################################
yri_msmc_4_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/msmc.yri_4.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(yri_msmc_4_in) # 500000 4
head(yri_msmc_4_in) # checked against hoffman. CORRECT.
colnames(yri_msmc_4_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
yri_msmc_4_avg <- aggregate(yri_msmc_4_in[,3:4],by=list(factor(yri_msmc_4_in$binStop)),sum)
head(yri_msmc_4_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
yri_msmc_4_avg$avgR2 <- yri_msmc_4_avg$sumR2/yri_msmc_4_avg$sumPairs
colnames(yri_msmc_4_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(yri_msmc_4_avg)
yri_msmc_4_forPlot <- yri_msmc_4_avg[,c(1,4)]
head(yri_msmc_4_forPlot)
# add label:
yri_msmc_4_forPlot$label <- "4. MSMC 4 Hap "
yri_msmc_4_forPlot$shape <- "average"
########## RECALL THAT THIS IS ONLY 5000 blocks not 20000

##################### YRI 8 ##############################################
yri_msmc_8_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/msmc.yri_8.ALL.LD.BINNED.OUTPUT.20170211.txt")
dim(yri_msmc_8_in) # 500000 4
head(yri_msmc_8_in) # checked against hoffman. CORRECT.
colnames(yri_msmc_8_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
yri_msmc_8_avg <- aggregate(yri_msmc_8_in[,3:4],by=list(factor(yri_msmc_8_in$binStop)),sum)
head(yri_msmc_8_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
yri_msmc_8_avg$avgR2 <- yri_msmc_8_avg$sumR2/yri_msmc_8_avg$sumPairs
colnames(yri_msmc_8_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(yri_msmc_8_avg)
yri_msmc_8_forPlot <- yri_msmc_8_avg[,c(1,4)]
head(yri_msmc_8_forPlot)
# add label:
yri_msmc_8_forPlot$label <- "5. MSMC 8 Hap "
yri_msmc_8_forPlot$shape <- "average"
#################################### SMC++ ###################################
###**NEW** smc++ with time points 
yri_smcpp_in <- read.table("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/1000GenomesData/LDDecay_NEWgenor2Results_round2_20170211/NewSimResults_round2_genor2_mac2_spotchecked/smcpp.yri.ALL.LD.BINNED.OUTPUT.20170215.txt")
dim(yri_smcpp_in) # 2000000 4
head(yri_smcpp_in) # checked against hoffman. CORRECT.
colnames(yri_smcpp_in) <- c("binStart","binStop","sumR2","sumPairs")

# use binStop as the bin ID # : 
yri_smcpp_avg <- aggregate(yri_smcpp_in[,3:4],by=list(factor(yri_smcpp_in$binStop)),sum)
head(yri_smcpp_avg) ### CHECK THIS SUPER CAREFULLY 


# now get average:
yri_smcpp_avg$avgR2 <- yri_smcpp_avg$sumR2/yri_smcpp_avg$sumPairs
colnames(yri_smcpp_avg) <- c("binStop","sumR2","sumPairs","avgR2")

# want just binStop,avgR2 and label for plotting:
head(yri_smcpp_avg)
yri_smcpp_forPlot <- yri_smcpp_avg[,c(1,4)]
head(yri_smcpp_forPlot)
# add label:
yri_smcpp_forPlot$label <- "6. SMC++ "
yri_smcpp_forPlot$shape <- "average"
