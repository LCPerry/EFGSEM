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

colnames(Hatoum_Trails)[6] <- 'SE'
colnames(Hatoum_Pairs)[6] <- 'SE'
colnames(Hatoum_SDST)[6] <- 'SE'
colnames(Hatoum_Digit)[6] <- 'SE'
colnames(Hatoum_Memory)[6] <- 'SE'
colnames(Donati_WM) <- c( "alternate_ids","rsid","chromosome","position","A2","A1","info","all_maf","missing_data_proportion","pvalue","add_info","effect","SE")
Ibrahim_Stroop3$Zscore <- (Ibrahim_Stroop3$Zscore)*(-1)
Hatoum_Pairs$beta <- (Hatoum_Pairs$beta)*(-1)

write.table(Hatoum_Trails, file = 'Hatoum_Trails.txt')
write.table(Hatoum_Pairs, file = 'Hatoum_Pairs.txt')
write.table(Hatoum_SDST, file = 'Hatoum_SDST.txt')
write.table(Hatoum_Digit, file = 'Hatoum_Digit.txt')
write.table(Hatoum_Memory, file = 'Hatoum_Memory.txt')
write.table(Ibrahim_DSST, file = 'Ibrahim_DSST.txt')
write.table(Ibrahim_Stroop3, file = 'Ibrahim_Stroop3.txt')
write.table(Donati_WM, file = 'Donati_WM.txt')

#Use SSRTBedopsCleaningScript to obtain the file for the Stop Signal task

# Munge
files<-c('Hatoum_Trails.txt', 'Hatoum_Pairs.txt', 'Hatoum_SDST.txt', 'Hatoum_Digit.txt','Hatoum_Memory.txt', 'Ibrahim_Stroop3.txt', 'Ibrahim_DSST.txt', 'Donati_WM.txt','Arnat_SSRT.txt')
hm3<-"reference.1000G.maf.0.005.txt"
trait.names<-c('Hatoum_Trails','Hatoum_Pairs','Hatoum_SDST','Hatoum_Digit','Hatoum_Memory','Ibrahim_Stroop3','Ibrahim_DSST' ,'Donati_WM','Arnat_SSRT')
N=c(93024,81701,84125,81701,162335,12866,32070,4611,11715)
info.filter=0.9
maf.filter=0.01

munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

#Executive Function LDSC

traits<-c('Hatoum_Trails.sumstats.gz', 'Hatoum_Pairs.sumstats.gz', 'Hatoum_SDST.sumstats.gz', 'Hatoum_Digit.sumstats.gz','Hatoum_Memory.sumstats.gz', "Ibrahim_Stroop3.sumstats.gz", 'Ibrahim_DSST.sumstats.gz','Donati_WM.sumstats.gz','Arnat_SSRT.sumstats.gz')
sample.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA)
population.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA)
ld<-"eur_w_ld_chr/"
wld<-"eur_w_ld_chr/"
trait.names<-c('Hatoum_Trails','Hatoum_Pairs','Hatoum_SDST','Hatoum_Digit','Hatoum_Memory', "Ibrahim_Stroop3", 'Ibrahim_DSST','Donati_WM','Arnat_SSRT')

LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)
save(LDSCoutput,file="EFLDSCoutput.RData")
load("~/SummaryStatistics/EFLDSCoutput.RData")


# Common factor and EFA
CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
CommonFactor_DWLS$modelfit
CommonFactor_DWLS$results

# CFA
# Accepted Model, SS=IC

cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs + Ibrahim_Stroop3 + Hatoum_Memory
        WM =~ NA*Donati_WM + Hatoum_Digit
        RT =~ NA*Hatoum_SDST + Ibrahim_DSST
        Donati_WM ~~ Hatoum_Memory
        SS ~~ 0*RT
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results


# Three Factor
cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs
        WM =~ NA*Donati_WM + Hatoum_Digit
        RT =~ NA*Hatoum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory
        Donati_WM ~~ Hatoum_Memory
        SS ~~ 0*RT
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# NO RT
cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs
        WM =~ NA*Donati_WM + Hatoum_Digit
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory
        Donati_WM ~~ Hatoum_Memory
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# WM=SS

cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs + Donati_WM + Hatoum_Digit
        RT =~ NA*Hatoum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory
        Donati_WM ~~ Hatoum_Memory
        SS ~~ 0*RT
        Hatoum_SDST ~~ a*Hatoum_SDST
        Ibrahim_DSST ~~ a*Ibrahim_DSST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# WM=IC

cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs
        WM =~ NA*Donati_WM + Hatoum_Digit + Ibrahim_Stroop3 + Hatoum_Memory
        RT =~ NA*Hatoum_SDST + Ibrahim_DSST
        Donati_WM ~~ Hatoum_Memory
        SS ~~ 0*RT
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# Common factor with RT

cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs + Ibrahim_Stroop3 + Hatoum_Memory + Donati_WM + Hatoum_Digit
        RT =~ NA*Hatoum_SDST + Ibrahim_DSST
        SS ~~ 0*RT
        Hatoum_SDST ~~ a*Hatoum_SDST
        Ibrahim_DSST ~~ a*Ibrahim_DSST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# Bifactor
cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs
        WM =~ NA*Donati_WM + Hatoum_Digit
        RT =~ NA*Hatoum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Hatoum_Memory
        cEF =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs + Donati_WM + Hatoum_Digit + Ibrahim_Stroop3 + Hatoum_Memory
        Donati_WM ~~ Hatoum_Memory
        SS ~~ 0*RT
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        cEF ~~ 0*SS
        cEF ~~ 0*WM
        cEF ~~ 0*IC
        WM ~~ 0*SS
        WM ~~ 0*IC
        SS ~~ 0*IC
        RT ~~ 0*cEF
        RT ~~ 0*IC
        RT ~~ 0*WM
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# Bifactor no Inhibition
cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs
        WM =~ NA*Donati_WM + Hatoum_Digit
        RT =~ NA*Hatoum_SDST + Ibrahim_DSST
        cEF =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs + Donati_WM + Hatoum_Digit + Ibrahim_Stroop3 + Hatoum_Memory
        Donati_WM ~~ Hatoum_Memory
        SS ~~ 0*RT
        Hatoum_SDST ~~ a*Hatoum_SDST
        Ibrahim_DSST ~~ a*Ibrahim_DSST
        a > .001
        cEF ~~ 0*SS
        cEF ~~ 0*WM
        WM ~~ 0*SS
        RT ~~ 0*cEF
        RT ~~ 0*WM
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# Bifactor, No Inhibition or Shifting

cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs + Donati_WM + Hatoum_Digit + Hatoum_Memory
        WM =~ NA*Donati_WM + Hatoum_Digit
        RT =~ NA*Hatoum_SDST + Ibrahim_DSST
        Donati_WM ~~ Hatoum_Memory
        SS ~~ 0*RT
        Hatoum_SDST ~~ a*Hatoum_SDST
        Hatoum_Digit ~~ a*Hatoum_Digit
        a > .001
        WM ~~ 0*SS
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor Without Digit-Specific

cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs
        WM =~ NA*Donati_WM + Hatoum_Digit
        cEF =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs + Donati_WM + Hatoum_Digit + Ibrahim_Stroop3 + Hatoum_Memory
        Donati_WM ~~ Hatoum_Memory
        Hatoum_SDST ~~ a*Hatoum_SDST
        a > .001
        cEF ~~ 0*SS
        cEF ~~ 0*WM
        WM ~~ 0*SS
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results


# GWAS
N<-c(93024,81701,84125,81701,162335,12866,32070,4611)
se.logit <-c(F,F,F,F,F,F,F,F)
Hail<-c(F,F,F,F,F,T,T,F)
files <- c('Hatoum_Trails.txt','Hatoum_Pairs.txt','Hatoum_SDST.txt','Hatoum_Digit.txt','Hatoum_Memory.txt','Ibrahim_Stroop3.txt','Ibrahim_DSST.txt','Donati_WM.txt')
ref <- "reference.1000G.maf.0.005.txt"
trait.names<-c('Hatoum_Trails','Hatoum_Pairs','Hatoum_SDST','Hatoum_Digit','Hatoum_Memory','Ibrahim_Stroop3','Ibrahim_DSST','Donati_WM')

EF_sumstats <- sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS = NULL,linprob=Hail,N=N,betas=NULL, info.filter=.6, maf.filter=0.01,keep.indel=FALSE,parallel=FALSE,cores=NULL)
save(EF_sumstats,file="EF_sumstats.RData")
load("~/SummaryStatistics/EF_sumstats.RData")

library(doParallel)

cl <- makePSOCKcluster(12, outfile='')
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("C:/Program Files/R/R-4.2.3/library"))
clusterEvalQ(cl, library(doParallel))


cEF <- 'SS =~ NA*Hatoum_SDST + Hatoum_Trails + Ibrahim_DSST + Hatoum_Pairs + Ibrahim_Stroop3 + Hatoum_Memory
        WM =~ NA*Donati_WM + Hatoum_Digit
        RT =~ NA*Hatoum_SDST + Ibrahim_DSST
        Donati_WM ~~ Hatoum_Memory
        SS ~~ 0*RT
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
        SS ~~ i*SS
        i > .001
        WM ~~ j*WM
        j > .001
        RT ~~ k*RT
        k > .001
        SS ~ SNP
        WM ~ SNP

' 

#run the multivariate GWAS
CorrelatedFactors<-userGWAS(covstruc = LDSCoutput, SNPs = EF_sumstats, model = cEF, parallel = TRUE, cores=12, sub = c('SS ~ SNP','WM ~ SNP'),toler=1e-60,std.lv=TRUE,)
save(CorrelatedFactors,file="CorrelatedFactors.RData")
