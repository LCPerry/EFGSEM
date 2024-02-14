library(GenomicSEM)
library(data.table)

setwd("~/SummaryStatistics/")

Ibrahim_Stroop3 <- fread('~/SummaryStatistics/IbrahimSumstats/ExecutiveGWAS_Hapmap/Stroop3_45agesex_IV.csv.gz')
Ibrahim_DSST <- fread('~/SummaryStatistics/IbrahimSumstats/ExecutiveGWAS_Hapmap/DSST_45agesex_effN.txt.gz')
Haroum_Trails <- fread('~/SummaryStatistics/HaroumSumstats/GCT90179118/GCST90179118_buildGRCh37.tsv')
Haroum_Pairs <- fread('~/SummaryStatistics/HaroumSumstats/GCST90179117/GCST90179117_buildGRCh37.tsv')
Haroum_SDST <- fread('~/SummaryStatistics/HaroumSumstats/GCST90179119/GCST90179119_buildGRCh37.tsv')
Haroum_Digit <- fread('~/SummaryStatistics/HaroumSumstats/GCST90179120/GCST90179120_buildGRCh37.tsv')
Haroum_Memory <- fread('~/SummaryStatistics/HaroumSumstats/GCST90179116/GCST90179116_buildGRCh37.tsv')
Donati_WM <- fread('~/SummaryStatistics/DonatiSumstats/Working_Memory.txt')

colnames(Haroum_Trails)[6] <- 'SE'
colnames(Haroum_Pairs)[6] <- 'SE'
colnames(Haroum_SDST)[6] <- 'SE'
colnames(Haroum_Digit)[6] <- 'SE'
colnames(Haroum_Memory)[6] <- 'SE'
colnames(Donati_WM) <- c( "alternate_ids","rsid","chromosome","position","A2","A1","info","all_maf","missing_data_proportion","pvalue","add_info","effect","SE")
Ibrahim_Stroop3$Zscore <- (Ibrahim_Stroop3$Zscore)*(-1)
Haroum_Pairs$beta <- (Haroum_Pairs$beta)*(-1)

write.table(Haroum_Trails, file = 'Haroum_Trails.txt')
write.table(Haroum_Pairs, file = 'Haroum_Pairs.txt')
write.table(Haroum_SDST, file = 'Haroum_SDST.txt')
write.table(Haroum_Digit, file = 'Haroum_Digit.txt')
write.table(Haroum_Memory, file = 'Haroum_Memory.txt')
write.table(Ibrahim_DSST, file = 'Ibrahim_DSST.txt')
write.table(Ibrahim_Stroop3, file = 'Ibrahim_Stroop3.txt')
write.table(Donati_WM, file = 'Donati_WM.txt')

#Use SSRTBedopsCleaningScript to obtain the file for the Stop Signal task

# Munge
files<-c('Haroum_Trails.txt', 'Haroum_Pairs.txt', 'Haroum_SDST.txt', 'Haroum_Digit.txt','Haroum_Memory.txt', 'Ibrahim_Stroop3.txt', 'Ibrahim_DSST.txt', 'Donati_WM.txt','Arnat_SSRT.txt')
hm3<-"reference.1000G.maf.0.005.txt"
trait.names<-c('Haroum_Trails','Haroum_Pairs','Haroum_SDST','Haroum_Digit','Haroum_Memory','Ibrahim_Stroop3','Ibrahim_DSST' ,'Donati_WM','Arnat_SSRT')
N=c(93024,81701,84125,81701,162335,12866,32070,4611,11715)
info.filter=0.9
maf.filter=0.01

munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

#Executive Function LDSC

traits<-c('Haroum_Trails.sumstats.gz', 'Haroum_Pairs.sumstats.gz', 'Haroum_SDST.sumstats.gz', 'Haroum_Digit.sumstats.gz','Haroum_Memory.sumstats.gz', "Ibrahim_Stroop3.sumstats.gz", 'Ibrahim_DSST.sumstats.gz','Donati_WM.sumstats.gz','Arnat_SSRT.sumstats.gz')
sample.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA)
population.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA)
ld<-"eur_w_ld_chr/"
wld<-"eur_w_ld_chr/"
trait.names<-c('Haroum_Trails','Haroum_Pairs','Haroum_SDST','Haroum_Digit','Haroum_Memory', "Ibrahim_Stroop3", 'Ibrahim_DSST','Donati_WM','Arnat_SSRT')

LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)
save(LDSCoutput,file="EFLDSCoutput.RData")
load("~/SummaryStatistics/EFLDSCoutput.RData")


# Common factor and EFA
CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
CommonFactor_DWLS$modelfit
CommonFactor_DWLS$results

# CFA
# Accepted Model, SS=IC

cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs + Ibrahim_Stroop3 + Haroum_Memory
        WM =~ NA*Donati_WM + Haroum_Digit
        RT =~ NA*Haroum_SDST + Ibrahim_DSST
        Donati_WM ~~ Haroum_Memory
        SS ~~ 0*RT
        Haroum_SDST ~~ a*Haroum_SDST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results


# Three Factor
cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs
        WM =~ NA*Donati_WM + Haroum_Digit
        RT =~ NA*Haroum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Haroum_Memory
        Donati_WM ~~ Haroum_Memory
        SS ~~ 0*RT
        Haroum_SDST ~~ a*Haroum_SDST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# NO RT
cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs
        WM =~ NA*Donati_WM + Haroum_Digit
        IC =~ NA*Ibrahim_Stroop3 + Haroum_Memory
        Donati_WM ~~ Haroum_Memory
        Haroum_SDST ~~ a*Haroum_SDST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# WM=SS

cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs + Donati_WM + Haroum_Digit
        RT =~ NA*Haroum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Haroum_Memory
        Donati_WM ~~ Haroum_Memory
        SS ~~ 0*RT
        Haroum_SDST ~~ a*Haroum_SDST
        Ibrahim_DSST ~~ a*Ibrahim_DSST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# WM=IC

cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs
        WM =~ NA*Donati_WM + Haroum_Digit + Ibrahim_Stroop3 + Haroum_Memory
        RT =~ NA*Haroum_SDST + Ibrahim_DSST
        Donati_WM ~~ Haroum_Memory
        SS ~~ 0*RT
        Haroum_SDST ~~ a*Haroum_SDST
        a > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# Common factor with RT

cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs + Ibrahim_Stroop3 + Haroum_Memory + Donati_WM + Haroum_Digit
        RT =~ NA*Haroum_SDST + Ibrahim_DSST
        SS ~~ 0*RT
        Haroum_SDST ~~ a*Haroum_SDST
        Ibrahim_DSST ~~ b*Ibrahim_DSST
        a > .001
        b > .001
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# Bifactor
cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs
        WM =~ NA*Donati_WM + Haroum_Digit
        RT =~ NA*Haroum_SDST + Ibrahim_DSST
        IC =~ NA*Ibrahim_Stroop3 + Haroum_Memory
        cEF =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs + Donati_WM + Haroum_Digit + Ibrahim_Stroop3 + Haroum_Memory
        Donati_WM ~~ Haroum_Memory
        SS ~~ 0*RT
        Haroum_SDST ~~ a*Haroum_SDST
        a > .001
        Haroum_Trails ~~ b*Haroum_Trails
        b > .001
        Ibrahim_DSST ~~ c*Ibrahim_DSST
        c > .001
        Haroum_Pairs ~~ d*Haroum_Pairs
        d > .001
        Ibrahim_Stroop3 ~~ e*Ibrahim_Stroop3
        e > .001
        Haroum_Memory ~~ f*Haroum_Memory
        f > .001
        Donati_WM ~~ g*Donati_WM
        g > .001
        Haroum_Digit ~~ h*Haroum_Digit
        h > .001
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
cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs
        WM =~ NA*Donati_WM + Haroum_Digit
        RT =~ NA*Haroum_SDST + Ibrahim_DSST
        cEF =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs + Donati_WM + Haroum_Digit + Ibrahim_Stroop3 + Haroum_Memory
        Donati_WM ~~ Haroum_Memory
        SS ~~ 0*RT
       Haroum_SDST ~~ a*Haroum_SDST
        a > .001
        Haroum_Trails ~~ b*Haroum_Trails
        b > .001
        Ibrahim_DSST ~~ c*Ibrahim_DSST
        c > .001
        Haroum_Pairs ~~ d*Haroum_Pairs
        d > .001
        Ibrahim_Stroop3 ~~ e*Ibrahim_Stroop3
        e > .001
        Haroum_Memory ~~ f*Haroum_Memory
        f > .001
        Donati_WM ~~ g*Donati_WM
        g > .001
        Haroum_Digit ~~ h*Haroum_Digit
        h > .001
        cEF ~~ 0*SS
        cEF ~~ 0*WM
        WM ~~ 0*SS
        RT ~~ 0*cEF
        
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

# Bifactor, No Inhibition or Shifting
cEF <- 'WM =~ NA*Donati_WM + Haroum_Digit
        RT =~ NA*Haroum_SDST + Ibrahim_DSST
        cEF =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs + Donati_WM + Haroum_Digit + Ibrahim_Stroop3 + Haroum_Memory
        Donati_WM ~~ Haroum_Memory
        Haroum_SDST ~~ a*Haroum_SDST
        a > .001
        Haroum_Trails ~~ b*Haroum_Trails
        b > .001
        Ibrahim_DSST ~~ c*Ibrahim_DSST
        c > .001
        Haroum_Pairs ~~ d*Haroum_Pairs
        d > .001
        Ibrahim_Stroop3 ~~ e*Ibrahim_Stroop3
        e > .001
        Haroum_Memory ~~ f*Haroum_Memory
        f > .001
        Donati_WM ~~ g*Donati_WM
        g > .001
        Haroum_Digit ~~ h*Haroum_Digit
        h > .001
        cEF ~~ 0*WM
        RT ~~ 0*cEF
        
        
' 
EF<-usermodel(LDSCoutput, estimation = "DWLS", model = cEF, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
EF$modelfit
EF$results

#Bifactor Without Digit-Specific

cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs
        WM =~ NA*Donati_WM + Haroum_Digit
        cEF =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs + Donati_WM + Haroum_Digit + Ibrahim_Stroop3 + Haroum_Memory
        Donati_WM ~~ Haroum_Memory
        Haroum_SDST ~~ a*Haroum_SDST
        a > .001
        Haroum_Trails ~~ b*Haroum_Trails
        b > .001
        Ibrahim_DSST ~~ c*Ibrahim_DSST
        c > .001
        Haroum_Pairs ~~ d*Haroum_Pairs
        d > .001
        Ibrahim_Stroop3 ~~ e*Ibrahim_Stroop3
        e > .001
        Haroum_Memory ~~ f*Haroum_Memory
        f > .001
        Donati_WM ~~ g*Donati_WM
        g > .001
        Haroum_Digit ~~ h*Haroum_Digit
        h > .001
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
files <- c('Haroum_Trails.txt','Haroum_Pairs.txt','Haroum_SDST.txt','Haroum_Digit.txt','Haroum_Memory.txt','Ibrahim_Stroop3.txt','Ibrahim_DSST.txt','Donati_WM.txt')
ref <- "reference.1000G.maf.0.005.txt"
trait.names<-c('Haroum_Trails','Haroum_Pairs','Haroum_SDST','Haroum_Digit','Haroum_Memory','Ibrahim_Stroop3','Ibrahim_DSST','Donati_WM')

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


cEF <- 'SS =~ NA*Haroum_SDST + Haroum_Trails + Ibrahim_DSST + Haroum_Pairs + Ibrahim_Stroop3 + Haroum_Memory
        WM =~ NA*Donati_WM + Haroum_Digit
        RT =~ NA*Haroum_SDST + Ibrahim_DSST
        Donati_WM ~~ Haroum_Memory
        SS ~~ 0*RT
        Haroum_SDST ~~ a*Haroum_SDST
        a > .001
        Haroum_Trails ~~ b*Haroum_Trails
        b > .001
        Ibrahim_DSST ~~ c*Ibrahim_DSST
        c > .001
        Haroum_Pairs ~~ d*Haroum_Pairs
        d > .001
        Ibrahim_Stroop3 ~~ e*Ibrahim_Stroop3
        e > .001
        Haroum_Memory ~~ f*Haroum_Memory
        f > .001
        Donati_WM ~~ g*Donati_WM
        g > .001
        Haroum_Digit ~~ h*Haroum_Digit
        h > .001
        SNP ~ SS
        SNP ~ WM
        SNP ~ RT

' 

#run the multivariate GWAS
CorrelatedFactors<-userGWAS(covstruc = LDSCoutput, SNPs = EF_sumstats, model = cEF, parallel = TRUE, cores=12, sub = c('SS ~ SNP','WM ~ SNP'),toler=1e-60,std.lv=TRUE,)
save(CorrelatedFactors,file="CorrelatedFactors.RData")
