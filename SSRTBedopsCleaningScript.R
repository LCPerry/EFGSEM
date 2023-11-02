library(data.table)
library(stringr)
library(caroline)
library(dplyr)

setwd("~/SummaryStatistics/")

Arnat_SSRT <- fread('C:/Users/lucas/Documents/SummaryStatistics/ArnatSumstats/EUR_results_SSRT_210706.txt.gz')

#Make bedops input file

Arnat_SSRTChrom <- data.frame(strsplit(as.character(Arnat_SSRT$MarkerName),':',fixed=TRUE))
rownames(Arnat_SSRTChrom) <- c("chromosome","position","A1","A2")
Arnat_SSRTChromt <- transpose(Arnat_SSRTChrom)
Arnat_SSRTChromt <- Arnat_SSRTChromt[,1:2]
Arnat_SSRTChromt[,3] <- Arnat_SSRTChromt[,2]
write.delim(Arnat_SSRTChromt, file = 'positionsSSRT.txt', row.names=FALSE, col.names=FALSE, quote = FALSE)

#Run Bedops, use positionsSSRT.txt as input.txt
#Combine rsid results

SSRTanswer <- fread('C:/Users/lucas/Documents/SummaryStatistics/BEDOPSanswers/Arnat_SSRT_answer.txt')
colnames(SSRTanswer) <- c('Chromosome','Start','End','rsid')
SSRTanswer$crpos <- paste(SSRTanswer$Chromosome, SSRTanswer$End, sep=":")
SSRTanswer <-SSRTanswer[,4:5]
SSRTanswer2 <- as.data.frame(gsub(".*;", "", SSRTanswer$rsid))
colnames(SSRTanswer2) <- c('rsid')
SSRTanswer3 <- as.data.frame(gsub(".*;", "", SSRTanswer2$rsid))
colnames(SSRTanswer3) <- c('rsid')
SSRTanswer4 <- as.data.frame(gsub(".*;", "", SSRTanswer3$rsid))
colnames(SSRTanswer4) <- c('rsid')
SSRTanswer5 <- as.data.frame(gsub(".*;", "", SSRTanswer4$rsid))
colnames(SSRTanswer5) <- c('rsid')
SSRTanswer2$crpos <- SSRTanswer$crpos
SSRTanswer3$crpos <- SSRTanswer$crpos
SSRTanswer4$crpos <- SSRTanswer$crpos
SSRTanswer5$crpos <- SSRTanswer$crpos
SSRTanswer$rsid <- gsub("\\;.*", "", SSRTanswer$rsid)
SSRTanswer2$rsid <- gsub("\\;.*", "", SSRTanswer2$rsid)
SSRTanswer3$rsid <- gsub("\\;.*", "", SSRTanswer3$rsid)
SSRTanswer4$rsid <- gsub("\\;.*", "", SSRTanswer4$rsid)
SSRTanswer5$rsid <- gsub("\\;.*", "", SSRTanswer5$rsid)
SSRTanswerfull<- rbind(SSRTanswer, SSRTanswer2)
SSRTanswerfull <- distinct(SSRTanswerfull)
SSRTanswerfull<- rbind(SSRTanswerfull, SSRTanswer3)
SSRTanswerfull <- distinct(SSRTanswerfull)
SSRTanswerfull<- rbind(SSRTanswerfull, SSRTanswer4)
SSRTanswerfull <- distinct(SSRTanswerfull)
SSRTanswerfull<- rbind(SSRTanswerfull, SSRTanswer5)
SSRTanswerfull <- distinct(SSRTanswerfull)

Arnat_SSRTChromt<- fread('C:/Users/lucas/Documents/SummaryStatistics/positionsSSRT.txt')
Arnat_SSRT$crpos <- paste(Arnat_SSRTChromt$V1, Arnat_SSRTChromt$V2, sep=":")
Arnat_SSRTFull <- merge(Arnat_SSRT, SSRTanswerfull, by='crpos')
Arnat_SSRTFull[Arnat_SSRTFull == ''] <- NA
Arnat_SSRTFull <- na.omit(Arnat_SSRTFull)
Arnat_SSRTFull <- Arnat_SSRTFull[,-1]
colnames(Arnat_SSRTFull) <- c("MarkerName","A1","A2","effect","StdErr","P","Direction","N","MAX_N","SNP")
Arnat_SSRTFull$A1 <- str_sub(Arnat_SSRTFull$A1, start = -1)
Arnat_SSRTFull$A2 <- str_sub(Arnat_SSRTFull$A2, start = -1)
Arnat_SSRTFull$A1 <- toupper(Arnat_SSRTFull$A1)
Arnat_SSRTFull$A2 <- toupper(Arnat_SSRTFull$A2)
Arnat_SSRTFull$effect <- Arnat_SSRTFull$effect*(-1)
write.table(Arnat_SSRTFull, file = 'Arnat_SSRT.txt')
