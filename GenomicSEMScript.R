library(GenomicSEM)
library(data.table)

setwd("~/SummaryStatistics/")

Ibrahim_Stroop3 <- fread('~/SummaryStatistics/IbrahimSumstats/ExecutiveGWAS_Hapmap/Stroop3_45agesex_IV.csv.gz')
Ibrahim_DSST <- fread('~/SummaryStatistics/IbrahimSumstats/ExecutiveGWAS_Hapmap/DSST_45agesex_effN.txt.gz')
Hatoum_Trails <- fread('~/SummaryStatistics/HatoumSumstats/GCT90179118/GCST90179118_buildGRCh37.tsv')
Hatoum_Pairs <- fread('~/SummaryStatistics/HatoumSumstats/GCST90179117/GCST90179117_buildGRCh37.tsv')
Hatoum_SDST <- fread('~/SummaryStatistics/HatoumSumstats/GCST90179119/GCST90179119_buildGRCh37.tsv')
Hatoum_Digit <- fread('~/SummaryStatistics/HatoumSumstats/GCST90179120/GCST90179120_buildGRCh37.tsv')
Hatoum_Memory <- fread('~/SummaryStatistics/HatoumSumstats/GCST90179116/GCST90179116_buildGRCh37.tsv')
Donati_WM <- fread('~/SummaryStatistics/DonatiSumstats/Working_Memory.txt')
Rahman_G6 <- fread('~/SummaryStatistics/RahmanSumstats/g6_GWAS_Sumstats_Cleaned.txt')
Rahman_G4 <- fread('~/SummaryStatistics/RahmanSumstats/g4_GWAS_Sumstats_Cleaned.txt')

colnames(Hatoum_Trails)[6] <- 'SE'
colnames(Hatoum_Pairs)[6] <- 'SE'
colnames(Hatoum_SDST)[6] <- 'SE'
colnames(Hatoum_Digit)[6] <- 'SE'
colnames(Hatoum_Memory)[6] <- 'SE'
colnames(Donati_WM) <- c( "alternate_ids","rsid","chromosome","position","A2","A1","info","all_maf","missing_data_proportion","pvalue","add_info","effect","SE")
Ibrahim_Stroop3$Zscore <- (Ibrahim_Stroop3$Zscore)*(-1)
Hatoum_Pairs$beta <- (Hatoum_Pairs$beta)*(-1)
Rahman_G6$BETA <- (Rahman_G6$BETA)*(-1)
Rahman_G4$BETA <- (Rahman_G4$BETA)*(-1)

write.table(Hatoum_Trails, file = 'Hatoum_Trails.txt')
write.table(Hatoum_Pairs, file = 'Hatoum_Pairs.txt')
write.table(Hatoum_SDST, file = 'Hatoum_SDST.txt')
write.table(Hatoum_Digit, file = 'Hatoum_Digit.txt')
write.table(Hatoum_Memory, file = 'Hatoum_Memory.txt')
write.table(Ibrahim_DSST, file = 'Ibrahim_DSST.txt')
write.table(Ibrahim_Stroop3, file = 'Ibrahim_Stroop3.txt')
write.table(Donati_WM, file = 'Donati_WM.txt')
write.table(Rahman_G6, file = 'Rahman_G6.txt', quote = FALSE, row.names = FALSE)
write.table(Rahman_G4, file = 'Rahman_G4.txt', quote = FALSE, row.names = FALSE)

#Use SSRTBedopsCleaningScript to obtain the file for the Stop Signal task

# Munge
files<-c('Hatoum_Trails.txt', 'Hatoum_Pairs.txt', 'Hatoum_SDST.txt', 'Hatoum_Digit.txt','Hatoum_Memory.txt', 'Ibrahim_Stroop3.txt', 'Ibrahim_DSST.txt', 'Donati_WM.txt','Arnat_SSRT.txt','Rahman_G6.txt','Rahman_G4.txt')
hm3<-"reference.1000G.maf.0.005.txt"
trait.names<-c('Hatoum_Trails','Hatoum_Pairs','Hatoum_SDST','Hatoum_Digit','Hatoum_Memory','Ibrahim_Stroop3','Ibrahim_DSST' ,'Donati_WM','Arnat_SSRT','Rahman_G6','Rahman_G4')
N=c(93024,81701,84125,81701,162335,12866,32070,4611,11715,9879,9879)
info.filter=0.9
maf.filter=0.01

munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

#Executive Function LDSC

traits<-c('Hatoum_Trails.sumstats.gz', 'Hatoum_Pairs.sumstats.gz', 'Hatoum_SDST.sumstats.gz', 'Hatoum_Digit.sumstats.gz','Hatoum_Memory.sumstats.gz', "Ibrahim_Stroop3.sumstats.gz", 'Ibrahim_DSST.sumstats.gz','Donati_WM.sumstats.gz','Arnat_SSRT.sumstats.gz','Rahman_G6.sumstats.gz','Rahman_G4.sumstats.gz')
sample.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
population.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
ld<-"eur_w_ld_chr/"
wld<-"eur_w_ld_chr/"
trait.names<-c('Hatoum_Trails','Hatoum_Pairs','Hatoum_SDST','Hatoum_Digit','Hatoum_Memory', "Ibrahim_Stroop3", 'Ibrahim_DSST','Donati_WM','Arnat_SSRT','Rahman_G6','Rahman_G4')

LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)
save(LDSCoutput,file="EFLDSCoutput.RData")
load("~/SummaryStatistics/EFLDSCoutput.RData")

#LDSC remove G6 and Stop Signal
traits<-c('Hatoum_Trails.sumstats.gz', 'Hatoum_Pairs.sumstats.gz', 'Hatoum_SDST.sumstats.gz', 'Hatoum_Digit.sumstats.gz','Hatoum_Memory.sumstats.gz', "Ibrahim_Stroop3.sumstats.gz", 'Ibrahim_DSST.sumstats.gz','Donati_WM.sumstats.gz','Rahman_G4.sumstats.gz')
sample.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
population.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
ld<-"eur_w_ld_chr/"
wld<-"eur_w_ld_chr/"
trait.names<-c('Hatoum_Trails','Hatoum_Pairs','Hatoum_SDST','Hatoum_Digit','Hatoum_Memory', "Ibrahim_Stroop3", 'Ibrahim_DSST','Donati_WM','Rahman_G4')


LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)
save(LDSCoutput,file="CLNLDSCoutput.RData")
load("~/SummaryStatistics/CLNLDSCoutput.RData")


# Common factor
CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
CommonFactor_DWLS$modelfit
CommonFactor_DWLS$results

# CFA
# Accepted Model, Bifactor with WM
cEF <- 'cEF =~ NA*Hatoum_Trails + Hatoum_SDST + Hatoum_Pairs + Ibrahim_Stroop3 + Hatoum_Memory + Ibrahim_DSST + Donati_WM + Hatoum_Digit + Rahman_G4
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs 
        RT =~ NA*Ibrahim_DSST + Hatoum_SDST
        Donati_WM ~~ NA*Hatoum_Memory
        cEF ~~ 0*RT
        cEF ~~ 0*WM
        Hatoum_SDST ~~ a*Hatoum_SDST 
        a > .001
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Three Factor
cEF <- 'SS =~ NA*Hatoum_Trails + Hatoum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory 
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs 
        RT =~ NA*Ibrahim_DSST + Hatoum_SDST
        Donati_WM ~~ NA*Hatoum_Memory
        SS ~~ 0*RT
        Ibrahim_DSST ~~ 0*Ibrahim_DSST
        Hatoum_Trails ~~ b*Hatoum_Trails
        b > .01
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Three Factor Pairs is SS
cEF <- 'SS =~ NA*Hatoum_Trails + Hatoum_SDST + Ibrahim_DSST + Hatoum_Pairs
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory 
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 
        RT =~ NA*Ibrahim_DSST + Hatoum_SDST
        Donati_WM ~~ NA*Hatoum_Memory
        SS ~~ 0*RT
        Ibrahim_DSST ~~ a*Ibrahim_DSST
        a > .001
        Hatoum_Trails ~~ b*Hatoum_Trails
        b > .01
        Hatoum_SDST ~~ c*Hatoum_SDST
        c > .001
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor With WM
cEF <- 'cEF =~ Hatoum_SDST + Hatoum_Trails + Hatoum_Pairs + Ibrahim_Stroop3 + Hatoum_Memory + Ibrahim_DSST + Donati_WM + Hatoum_Digit + Rahman_G4
        WM =~ Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs
        Donati_WM ~~ Hatoum_Memory
        cEF ~~ 0*WM
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor with WM, Substitution is Shifting
cEF <- 'cEF =~ NA*Hatoum_SDST + Hatoum_Trails + Hatoum_Pairs + Ibrahim_Stroop3 + Ibrahim_DSST + Hatoum_Memory + Donati_WM + Hatoum_Digit + Rahman_G4
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs
        SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST
        Donati_WM ~~ Hatoum_Memory
        cEF ~~ 0*WM
        cEF ~~ 0*SS
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor No Pairs on WM
cEF <- 'cEF =~ NA*Hatoum_Trails + Hatoum_SDST + Hatoum_Pairs + Ibrahim_Stroop3 + Hatoum_Memory + Ibrahim_DSST + Donati_WM + Hatoum_Digit + Rahman_G4
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4
        RT =~ NA*Ibrahim_DSST + Hatoum_SDST
        Donati_WM ~~ NA*Hatoum_Memory
        cEF ~~ 0*RT
        cEF ~~ 0*WM
        Hatoum_SDST ~~ a*Hatoum_SDST 
        a > .001
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor Prospective loads onto WM
cEF <- 'cEF =~ NA*Hatoum_Trails + Hatoum_SDST + Hatoum_Pairs + Ibrahim_Stroop3 + Hatoum_Memory + Ibrahim_DSST + Donati_WM + Hatoum_Digit + Rahman_G4
        WM =~ NA*Donati_WM + Hatoum_Digit + Hatoum_Pairs + Hatoum_Memory
        RT =~ NA*Ibrahim_DSST + Hatoum_SDST
        cEF ~~ 0*RT
        cEF ~~ 0*WM
        Hatoum_SDST ~~ a*Hatoum_SDST 
        a > .001
        Hatoum_Digit ~~ b*Hatoum_Digit
        b > .001
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Substitution Correlated Residual
cEF <- 'cEF =~ NA*Hatoum_Trails + Hatoum_SDST + Hatoum_Pairs + Ibrahim_Stroop3 + Hatoum_Memory + Ibrahim_DSST + Donati_WM + Hatoum_Digit + Rahman_G4
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs 
        Ibrahim_DSST ~~ Hatoum_SDST
        Donati_WM ~~ NA*Hatoum_Memory
        cEF ~~ 0*WM
        Hatoum_SDST ~~ a*Hatoum_SDST 
        a > .001
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results


#SSRT Models from supplemental methods
#Three Factor
cEF <- 'SS =~ NA*Hatoum_Trails + Hatoum_SDST + Ibrahim_DSST + Hatoum_Pairs 
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory + Arnat_SSRT
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Three Factor, Donati crossloading
cEF <- 'SS =~ Hatoum_Trails + Hatoum_SDST + Ibrahim_DSST
        IC =~ Ibrahim_Stroop3 + Hatoum_Memory + Arnat_SSRT + Donati_WM 
        WM =~ Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs
        
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Three Factor with RT
cEF <- 'SS =~ NA*Hatoum_Trails + Hatoum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory + Arnat_SSRT
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs 
        RT =~ NA*Ibrahim_DSST + Hatoum_SDST
        SS ~~ 0*RT
        Ibrahim_DSST ~~ a*Ibrahim_DSST
        a > .001
        Hatoum_Trails ~~ b*Hatoum_Trails
        b > .01
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Three Factor with RT, Donati Crossloading
cEF <- 'SS =~ NA*Hatoum_Trails + Hatoum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs 
        RT =~ NA*Ibrahim_DSST + Hatoum_SDST
        SS ~~ 0*RT
        Hatoum_Trails ~~ b*Hatoum_Trails
        b > .005
        Ibrahim_DSST ~~ 0*Ibrahim_DSST
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor
cEF <- 'cEF =~ NA*Hatoum_SDST + Hatoum_Trails + Hatoum_Pairs + Ibrahim_Stroop3 + Ibrahim_DSST + Hatoum_Memory + Donati_WM + Hatoum_Digit + Rahman_G4 + Arnat_SSRT
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs
        RT =~ NA*Hatoum_SDST+ Ibrahim_DSST
        cEF ~~ 0*WM
        cEF ~~ 0*RT
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        Arnat_SSRT ~~ b*Arnat_SSRT
        b > .001

' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor with Inhibition, Donati Crossloading
cEF <- 'cEF =~ NA*Hatoum_SDST + Hatoum_Trails + Hatoum_Pairs + Ibrahim_Stroop3 + Ibrahim_DSST + Hatoum_Memory + Donati_WM + Hatoum_Digit + Rahman_G4 + Arnat_SSRT
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs
        RT =~ NA*Hatoum_SDST+ Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory + Arnat_SSRT + Donati_WM
        cEF ~~ 0*WM
        cEF ~~ 0*RT
        cEF ~~ 0*IC
        WM ~~ 0*IC
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        Arnat_SSRT ~~ b*Arnat_SSRT
        b > .001

' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor with Inhibition, drop non-sig
cEF <- 'cEF =~ NA*Hatoum_SDST + Hatoum_Trails + Hatoum_Pairs + Ibrahim_Stroop3 + Ibrahim_DSST + Hatoum_Memory + Donati_WM + Hatoum_Digit + Rahman_G4 + Arnat_SSRT
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs
        RT =~ NA*Hatoum_SDST+ Ibrahim_DSST
        IC =~ NA*Donati_WM + Arnat_SSRT
        cEF ~~ 0*WM
        cEF ~~ 0*RT
        cEF ~~ 0*IC
        WM ~~ 0*IC
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        Arnat_SSRT ~~ b*Arnat_SSRT
        b > .001

' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Four Factor
cEF <- 'SS =~ NA*Hatoum_Trails + Hatoum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory
        BI =~ NA*Hatoum_Memory + Arnat_SSRT + Donati_WM 
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs
        RT =~ NA*Ibrahim_DSST + Hatoum_SDST
        SS ~~ 0*RT
        Ibrahim_Stroop3 ~~ a*Ibrahim_Stroop3
        a > .001
        Hatoum_Trails ~~ b*Hatoum_Trails
        b > .005
        Ibrahim_DSST ~~ 0*Ibrahim_DSST
        
' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor Two Inhibitions
cEF <- 'cEF =~ NA*Hatoum_SDST + Hatoum_Trails + Hatoum_Pairs + Ibrahim_Stroop3 + Ibrahim_DSST + Hatoum_Memory + Donati_WM + Hatoum_Digit + Rahman_G4 + Arnat_SSRT
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs
        RT =~ NA*Hatoum_SDST+ Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory
        RI =~ NA*Hatoum_Memory + Arnat_SSRT + Donati_WM
        cEF ~~ 0*WM
        cEF ~~ 0*RT
        cEF ~~ 0*IC
        cEF ~~ 0*RI
        WM ~~ 0*IC
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        Arnat_SSRT ~~ b*Arnat_SSRT
        b > .001
        Donati_WM ~~ c*Donati_WM
        c > .001

' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor just RI
cEF <- 'cEF =~ NA*Hatoum_SDST + Hatoum_Trails + Hatoum_Pairs + Ibrahim_Stroop3 + Ibrahim_DSST + Hatoum_Memory + Donati_WM + Hatoum_Digit + Rahman_G4 + Arnat_SSRT
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs
        RT =~ NA*Hatoum_SDST+ Ibrahim_DSST
        RI =~ NA*Hatoum_Memory + Arnat_SSRT + Donati_WM
        cEF ~~ 0*WM
        cEF ~~ 0*RT
        cEF ~~ 0*RI
        WM ~~ 0*RI
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        Arnat_SSRT ~~ b*Arnat_SSRT
        b > .001

' 

EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results


# GWAS
N<-c(93024,81701,84125,81701,162335,12866,32070,4611,9879)
se.logit <-c(F,F,F,F,F,F,F,F,F)
Hail<-c(F,F,F,F,F,T,T,F,F)
files <- c('Hatoum_Trails.txt','Hatoum_Pairs.txt','Hatoum_SDST.txt','Hatoum_Digit.txt','Hatoum_Memory.txt','Ibrahim_Stroop3.txt','Ibrahim_DSST.txt','Donati_WM.txt','Rahman_G4.txt')
ref <- "reference.1000G.maf.0.005.txt"
trait.names<-c('Hatoum_Trails','Hatoum_Pairs','Hatoum_SDST','Hatoum_Digit','Hatoum_Memory','Ibrahim_Stroop3','Ibrahim_DSST','Donati_WM','Rahman_G4')


EF_sumstats <- sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS = NULL,linprob=Hail,N=N,betas=NULL, info.filter=.6, maf.filter=0.01,keep.indel=FALSE,parallel=FALSE,cores=NULL)
save(EF_sumstats,file="EF_sumstats.RData")
load("~/SummaryStatistics/EF_sumstats.RData")

library(doParallel)

cl <- makePSOCKcluster(12, outfile='')
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("C:/Program Files/R/R-4.2.3/library"))
clusterEvalQ(cl, library(doParallel))

CommonFactors<-commonfactorGWAS(covstruc = LDSCoutput, SNPs = EF_sumstats, parallel = TRUE, cores=12, toler=1e-60, smooth_check=TRUE)
save(CorrelatedFactors,file="CorrelatedFactors.RData")

EF <- 'cEF =~ NA*Hatoum_Trails + Hatoum_SDST + Hatoum_Pairs + Ibrahim_Stroop3 + Hatoum_Memory + Ibrahim_DSST + Donati_WM + Hatoum_Digit + Rahman_G4
        WM =~ NA*Donati_WM + Hatoum_Digit + Rahman_G4 + Hatoum_Pairs 
        RT =~ NA*Ibrahim_DSST + Hatoum_SDST
        Donati_WM ~~ NA*Hatoum_Memory
        cEF ~~ 0*RT
        cEF ~~ 0*WM
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        Hatoum_Trails ~~ b*Hatoum_Trails
        b > .001
        Ibrahim_DSST ~~ c*Ibrahim_DSST
        c > .001
        Hatoum_Pairs ~~ d*Hatoum_Pairs
        d > .001
        Ibrahim_Stroop3 ~~ e*Ibrahim_Stroop3
        e > .001
        Hatoum_Memory ~~ f*Hatoum_Memory
        f > .001
        Donati_WM ~~ g*Donati_WM
        g > .001
        Hatoum_Digit ~~ h*Hatoum_Digit
        h > .001
        Rahman_G4 ~~ i*Rahman_G4
        i > .001
        cEF ~ SNP
        WM ~ SNP
        RT ~ SNP
        
' 

#run the multivariate GWAS
CorrelatedFactors<-userGWAS(covstruc = LDSCoutput, SNPs = EF_sumstatsSub, model = EF, parallel = TRUE, cores = 16, sub = c('cEF ~ SNP','WM ~ SNP','RT ~ SNP'),toler=1e-60,std.lv=TRUE,smooth_check=TRUE,Q_SNP=TRUE)
save(CorrelatedFactors,file="CorrelatedFactors.RData")
