#Script authors: Carolina Reyes-Puig and Antigoni Kaliontzopoulou 
#This script is part of the paper Rensch’s Rule: linking intraspecific to evolutionary allometry
#Reyes-Puig, C., Adams, D. C., Enriquez-Urzelai, U. &  Kaliontzopoulou, A. (2023). Evolution. 

rm(list = ls())

library(ape)
library(RRPP)
library(nlme)
require(gridExtra)
library(tidyverse)


# Read species data and separate phenotypes and tree
load("SPdata.bin")
phy <- sp.dt$phy 
phy
sp.means.bs <- sp.dt$sp.means[1] ## body size means 

#### trait species means, which are already standardized (trait/SVL) at individual level##
bd.trait <- sp.dt$dt.trait

lin.sex <- as.factor(paste(bd.trait$sp_lin, bd.trait$SEX, sep = "_"))
trait.means <- lapply(3:12, function(x){
  stats <- matrix(tapply(bd.trait[,x], lin.sex, mean)[c(1:11, 14, 12, 13, 15:38)], 
                  ncol = 2, byrow = T)
  colnames(stats) <- levels(bd.trait$SEX)
  rownames(stats) <- levels(bd.trait$sp_lin)
  return(stats)
})
names(trait.means) <- colnames(bd.trait)[3:12]

##SD RATIO for body size as log(Y/X)####

SD.r.bs <- log(sp.means.bs$SVL[,2]/sp.means.bs$SVL[,1])


# Calculate species mean size as log(mean species size) ####
sp.sz <- log(apply(sp.means.bs$SVL, 1, mean))


### Calculate SD RATIO for traits ####

SD.ratio.t <- data.frame(lapply(trait.means, function (x) {log(x[,2]/ x[,1])}) [-1])

SD.ratio.all <- cbind(SD.r.bs, SD.ratio.t) 
colnames(SD.ratio.all)[1] <- "SVL"

##include slopes difference ###### 
slopesTRL <- read.table("slopeTRL_sex.txt", sep = ",")
slopesHS <- read.table("slopeHS_sex.txt", sep= ",", row.names = NULL)[-1]
slopesPL <- read.table("slopePL_sex.txt", sep= ",", row.names = NULL)[-1]
slopesMO <- read.table("slopeMO_sex.txt", sep= ",", row.names = NULL)[-1]
slopesHLL<-read.table("slopeHLL_sex.txt", sep= ",", row.names = NULL)[-1]
slopesFLL<-read.table("slopeFLL_sex.txt", sep= ",", row.names = NULL)[-1]

## include intercepts difference##
interTRL <- read.table("interceptTRL_sex.txt", sep = ",", row.names = NULL)[-1]
interHS <- read.table("intercepHS_sex.txt", sep= ",",row.names = NULL)[-1]
interPL <- read.table("intercepPL_sex.txt", sep= ",",row.names = NULL)[-1]
interMO <- read.table("intercepMO_sex.txt", sep= ",",row.names = NULL)[-1]
interHLL<-read.table("intercepHLL_sex.txt", sep= ",",row.names = NULL)[-1]
interFLL<-read.table("intercepFLL_sex.txt", sep= ",",row.names = NULL)[-1]


SD.ratio.all <- cbind(SD.ratio.all, slopesTRL[1], slopesHS[1], slopesPL[1], slopesMO[1], 
                  slopesHLL[1], slopesFLL[1], interTRL[1], interHS[1], interPL[1], 
                  interMO[1], interHLL [1], interFLL[1], sp.sz)

colnames(SD.ratio.all)[11] <- "slopesTRL"
colnames(SD.ratio.all)[17] <- "interTRL"

#write.table(SD.ratio.all,"SD.ratio.all.txt",sep=",")

SD.ratio.all <- cbind(SD.ratio.all, rownames(SD.ratio.all))
names(SD.ratio.all)[names(SD.ratio.all) == "rownames(SD.ratio.all)"] <- "tip.label"

# Test for RR with SD RATIO####
# SIZE 
fit.SizeSD <- lm.rrpp(SVL ~ sp.sz, Cov = vcv.phylo(phy), data = SD.ratio.all)
anova(fit.SizeSD)
summary(fit.SizeSD)
#write.table if necessary 


par(mar = c(4, 4, 1, 1), mgp = c(2, 0.5, 0))
plot(sp.sz, SD.ratio.all$SVL, pch=21, cex=2, bg="black", 
     xlab = "log(Species size)", ylab = "SD ratio Size", cex.lab = 2.5, cex.axis = 2, ylim = c(-0.06, 0.2)) 
abline(h = 0, lwd = 2.5, lty = 2)
abline(a = coef(fit.SizeSD)[1], 
       b = coef(fit.SizeSD)[2], col = "red", lwd = 3.5)
# HS 

fit.HS <- lm.rrpp(HS ~ sp.sz * slopesHS * interHS, Cov = vcv.phylo(phy), data = SD.ratio.all)
anova(fit.HS)
summary(fit.HS)

# TRL
fit.TRL <- lm.rrpp(TRL ~ sp.sz * slopesTRL * interTRL, Cov = vcv.phylo(phy), data = SD.ratio.all)
anova(fit.TRL)
summary(fit.TRL)

# HLL
fit.HLL <- lm.rrpp(HLL ~ sp.sz * slopesHLL * interHLL, Cov = vcv.phylo(phy), data = SD.ratio.all)
anova(fit.HLL)
summary(fit.HLL)

# FLL
fit.FLL <- lm.rrpp(FLL ~ sp.sz * slopesFLL *interFLL, Cov = vcv.phylo(phy), data = SD.ratio.all)
anova(fit.FLL)
summary(fit.FLL)
