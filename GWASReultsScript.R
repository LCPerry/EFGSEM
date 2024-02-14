library(data.table)
library(GenomicSEM)

setwd("~/SummaryStatistics/")
load("~/SummaryStatistics/CommonFactor/CommonFactor.RData")
load("~/SummaryStatistics/EFGSEMResults/CorrelatedFactors.RData")

EFResults <- CommonFactor
SSResults <- CorrelatedFactors[[1]]
WMResults <- CorrelatedFactors[[2]]

EFResultsClean <- as.data.frame(filter(EFResults, EFResults$Z_smooth < 1))
SSResultsClean <- as.data.frame(filter(SSResults, SSResults$Z_smooth < 1))  
SSResultsClean <- SSResultsClean[- grep('lavaan WARNING:*',SSResultsClean$warning),]
WMResultsClean <- as.data.frame(filter(WMResults, WMResults$Z_smooth < 1))
WMResultsClean <- WMResultsClean[- grep('lavaan WARNING:*',WMResultsClean$warning),]


CorrelatedFactorsClean <- list(SSResultsClean, WMResultsClean)
save(CorrelatedFactorsClean,file='CorrelatedFactorsClean.RData')

load("~/SummaryStatistics/CorrelatedFactorsClean.RData")

##Calculate Effective Sample Size for Factor 1
#restrict to MAF of 40% and 10%

EFResultsClean<-subset(EFResultsClean, EFResultsClean$MAF <= .4 & EFResultsClean$MAF >= .1)
CorrelatedFactorsSub1<-subset(SSResultsClean, SSResultsClean$MAF <= .4 & SSResultsClean$MAF >= .1)
CorrelatedFactorsSub2<-subset(WMResultsClean, WMResultsClean$MAF <= .4 & WMResultsClean$MAF >= .1)

EFResultsClean <- na.omit(EFResultsClean)
N_hat_EF<-mean(1/((2*EFResultsClean$MAF*(1-EFResultsClean$MAF))*EFResultsClean$se_c^2))

CorrelatedFactorsSub1 <- na.omit(CorrelatedFactorsSub1)
N_hat_SS<-mean(1/((2*CorrelatedFactorsSub1$MAF*(1-CorrelatedFactorsSub1$MAF))*CorrelatedFactorsSub1$SE^2))

CorrelatedFactorsSub2 <- na.omit(CorrelatedFactorsSub2)
N_hat_WM<-mean(1/((2*CorrelatedFactorsSub2$MAF*(1-CorrelatedFactorsSub2$MAF))*CorrelatedFactorsSub2$SE^2))

colnames(EFResultsClean)[14] <- 'P'
colnames(CorrelatedFactorsSub1)[15] <- 'P'
colnames(CorrelatedFactorsSub2)[15] <- 'P'

write.table(EFResultsClean, file = 'CommonFactor.txt', quote = FALSE, row.names = FALSE)
write.table(CorrelatedFactorsSub1, file = 'CorrelatedFactorsSub1.txt')
write.table(CorrelatedFactorsSub2, file = 'CorrelatedFactorsSub2.txt')


# Munge
files<-c('CommonFactor.txt','CorrelatedFactorsSub1.txt','CorrelatedFactorsSub2.txt')
hm3<-"C:/Users/lucas/Documents/SummaryStatistics/reference.1000G.maf.0.005.txt"
trait.names<-c('Perry_EF','Perry_SS','Perry_WM')
N=c(4409915,77687,43156)
info.filter=0.9
maf.filter=0.01

munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

Haroum_Full <- fread('C:/Users/lucas/Documents/SummaryStatistics/HaroumSumstats/GCST90162547/GCST90162547_buildGRCh37.tsv')
Haroum_Reaction <- fread('C:/Users/lucas/Documents/SummaryStatistics/HaroumSumstats/GCST90179115/GCST90179115_buildGRCh37.tsv')
Okbay_EA <- fread('C:/Users/lucas/Documents/SummaryStatistics/Other Trait Sumstats/EA4_additive_excl_23andMe.txt.gz')
Demange_CP <- fread('C:/Users/lucas/Documents/SummaryStatistics/Other Trait Sumstats/GWAS_CP_all.txt')
PGC_ADHD <- fread('C:/Users/lucas/Documents/SummaryStatistics/ADHDSumstats/ADHD2022_iPSYCH_deCODE_PGC.meta.gz')
PGC_Autism <- fread('C:/Users/lucas/Documents/SummaryStatistics/AutismSumstats/iPSYCH-PGC_ASD_Nov2017.gz')
PGC_SCZ <- fread('C:/Users/lucas/Documents/SummaryStatistics/SCZSumstats/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz')
PGC_BiP <- fread('C:/Users/lucas/Documents/SummaryStatistics/BiPolar/daner_bip_pgc3_nm_noukbiobank.gz',fill=TRUE)
PGC_Anx <- fread('C:/Users/lucas/Documents/SummaryStatistics/ANX/pgc-panic2019.vcf.tsv.gz')
PGC_AN <- fread('C:/Users/lucas/Documents/SummaryStatistics/Anorexia/pgcAN2.2019-07.vcf.tsv.gz')
PGC_MDD <- fread('C:/Users/lucas/Documents/SummaryStatistics/MDD/PGC_UKB_depression_genome-wide.txt')
PGC_OCD <- fread('C:/Users/lucas/Documents/SummaryStatistics/OCD/ocd_aug2017.gz')
PGC_TS <- fread('C:/Users/lucas/Documents/SummaryStatistics/Tourette/TS_Oct2018.gz')
PGC_AUD <- fread('C:/Users/lucas/Documents/SummaryStatistics/BiPolar/daner_bip_pgc3_nm_noukbiobank.gz')
PGC_CUD <- fread('C:/Users/lucas/Documents/SummaryStatistics/BiPolar/daner_bip_pgc3_nm_noukbiobank.gz')
PGC_OUD <- fread('C:/Users/lucas/Documents/SummaryStatistics/BiPolar/daner_bip_pgc3_nm_noukbiobank.gz')

write.table(Haroum_Full, file = 'Haroum_Full.txt')
write.table(Haroum_Reaction, file = 'Haroum_Reaction.txt')
write.table(Okbay_EA, file = 'Okbay_EA.txt')
write.table(Demange_CP, file = 'Demange_CP.txt')
write.table(PGC_ADHD, file = 'PGC_ADHD.txt')
write.table(PGC_Autism, file = 'PGC_Autism.txt')
write.table(PGC_SCZ, file = 'PGC_SCZ.txt')
write.table(PGC_BiP, file = 'PGC_BiP.txt')
write.table(PGC_Anx, file = 'PGC_Anx.txt')
write.table(PGC_AN, file = 'PGC_AN.txt')
write.table(PGC_MDD, file = 'PGC_MDD.txt')
write.table(PGC_OCD, file = 'PGC_OCD.txt')
write.table(PGC_TS, file = 'PGC_TS.txt')

files<-c('Haroum_Full.txt', 'Haroum_Reaction.txt','Demange_CP.txt','Okbay_EA.txt','PGC_ADHD.txt','PGC_Autism.txt','PGC_SCZ.txt','PGC_BiP.txt','PGC_Anx.txt','PGC_AN.txt','PGC_MDD.txt','PGC_OCD.txt','PGC_TS.txt')
hm3<-"reference.1000G.maf.0.005.txt"
trait.names<-c('Haroum_Full','Haroum_Reaction','Demange_CP','Okbay_EA','PGC_ADHD','PGC_Autism','PGC_SCZ','PGC_BiP','PGC_Anx','PGC_AN','PGC_MDD','PGC_OCD','PGC_TS')
N=c(427037,432297,11351,765283,51568,41457,58746,101574,46322,97171,430424,7281,12140)
info.filter=0.9
maf.filter=0.01

munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)


setwd("~/SummaryStatistics/")

traits<-c('Perry_EF.sumstats.gz','Perry_SS.sumstats.gz','Perry_WM.sumstats.gz','Haroum_Full.sumstats.gz','Haroum_Reaction.sumstats.gz','Okbay_EA.sumstats.gz','Demange_CP.sumstats.gz','PGC_ADHD.sumstats.gz','PGC_Autism.sumstats.gz','PGC_SCZ.sumstats.gz','PGC_BiP.sumstats.gz','PGC_Anx.sumstats.gz','PGC_AN.sumstats.gz','PGC_MDD.sumstats.gz','PGC_OCD.sumstats.gz','PGC_TS.sumstats.gz')
sample.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
population.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA,0.05,0.012,0.01,0.02,0.1,0.01,0.302,0.025,0.008)
ld<-"eur_w_ld_chr/"
wld<-"eur_w_ld_chr/"
trait.names<-c('Perry_EF','Perry_SS','Perry_WM','Haroum_Full','Haroum_Reaction','Okbay_EA','Demange_CP','PGC_ADHD','PGC_Autism','PGC_SCZ','PGC_BiP','PGC_Anx','PGC_AN','PGC_MDD','PGC_OCD','PGC_TS')

LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)
save(LDSCoutput,file="PGCLDSCoutput.RData")
load("C:/Users/lucas/Documents/SummaryStatistics/PGCLDSCoutput.RData")



SSmodel<-'Okbay_EA ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results

SSmodel<-'PGC_ADHD ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results

SSmodel<-'PGC_Autism ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results

SSmodel<-'PGC_SCZ ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results


SSmodel<-'PGC_BiP ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results

SSmodel<-'PGC_Anx ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results

SSmodel<-'PGC_AN ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results

SSmodel<-'PGC_MDD ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results

SSmodel<-'PGC_OCD ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results

SSmodel<-'PGC_TS ~ Perry_SS + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_SS
         Haroum_Reaction ~ Perry_SS
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=SSmodel)
output$results



WMmodel<-'Okbay_EA ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results

WMmodel<-'PGC_ADHD ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results

WMmodel<-'PGC_Autism ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results

WMmodel<-'PGC_SCZ ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results

WMmodel<-'PGC_BiP ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results

WMmodel<-'PGC_Anx ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results

WMmodel<-'PGC_AN ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results

WMmodel<-'PGC_MDD ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results

WMmodel<-'PGC_OCD ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results

WMmodel<-'PGC_TS ~ Perry_WM + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_WM
         Haroum_Reaction ~ Perry_WM
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=WMmodel)
output$results


EFmodel<-'Okbay_EA ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results

EFmodel<-'PGC_ADHD ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results

EFmodel<-'PGC_Autism ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results

EFmodel<-'PGC_SCZ ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results

EFmodel<-'PGC_BiP ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results

EFmodel<-'PGC_Anx ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results

EFmodel<-'PGC_AN ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results

EFmodel<-'PGC_MDD ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results

EFmodel<-'PGC_OCD ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results

EFmodel<-'PGC_TS ~ Perry_EF + Demange_CP + Haroum_Reaction
         Demange_CP ~ Perry_EF
         Haroum_Reaction ~ Perry_EF
         '

output<-usermodel(LDSCoutput,estimation="DWLS",model=EFmodel)
output$results