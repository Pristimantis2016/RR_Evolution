#Authors: Carolina Reyes-Puig and Antigoni Kaliontzopoulou 
#This script is part of the paper Rensch’s Rule: linking intraspecific to evolutionary allometry
#Reyes-Puig, C., Adams, C., Enriquez-Urzelai, U. &  Kaliontzopoulou, A. (2023). Evolution. 

rm(list = ls())
library("RRPP")

dt <- read.csv("dt.15logmax.csv", sep = ",", header = T)
str(dt)
dt$SEX <- as.factor(dt$SEX)
dt$linAK <- as.factor(dt$linAK)

# 1) TEST WHETHER MALES AND FEMALES DIFER IN THEIR ALLOMETRIC SLOPE WITHIN EACH LINEAGE ####
#__________________________________________________________________________________________#
#TRL
anova.statsTRL <- NULL
summary.statsTRL <- NULL
for (l in levels(dt$linAK)){
  dt.sub <- subset(dt, dt$linAK == l)
  lm.lTRL <- lm.rrpp (dt.sub$logTRL ~ dt.sub$logSVL*dt.sub$SEX)
  stats.sumTRL <- cbind(coef(lm.lTRL, test = T), l)
  stats.subTRL <- cbind(anova(lm.lTRL)[[1]], l)
  summary.statsTRL <- rbind(summary.statsTRL, stats.subTRL)
  anova.statsTRL <- rbind(anova.statsTRL, stats.subTRL)
}
##write table if necessary

#HS
anova.statsHS <- NULL
for (l in levels(dt$linAK)){
  dt.sub <- subset(dt, dt$linAK == l)
  lm.lHS <- lm.rrpp(dt.sub$logHS ~ dt.sub$logSVL*dt.sub$SEX)
  stats.subHS <- cbind(anova(lm.lHS)[[1]], l)
  anova.statsHS <- rbind(anova.statsHS, stats.subHS)
}
##write table if necessary

#HLL
anova.statsHLL <- NULL
for (l in levels(dt$linAK)){
  dt.sub <- subset(dt, dt$linAK == l)
  lm.lHLL <- lm.rrpp(dt.sub$logHLL ~ dt.sub$logSVL*dt.sub$SEX)
  stats.subHLL <- cbind(anova(lm.lHLL)[[1]], l)
  anova.statsHLL <- rbind(anova.statsHLL, stats.subHLL)
}


#FLL
anova.statsFLL <- NULL
for (l in levels(dt$linAK)){
  dt.sub <- subset(dt, dt$linAK == l)
  lm.lFLL <- lm.rrpp(dt.sub$logFLL ~ dt.sub$logSVL*dt.sub$SEX)
  stats.subFLL <- cbind(anova(lm.lFLL)[[1]], l)
  anova.statsFLL <- rbind(anova.statsFLL, stats.subFLL)
}


# 2) TEST WHETHER THE DEGREE OF DIMORPHISM IN ALLOMETRIC SLOPE VARIES ACROSS LINEAGES ####
#________________________________________________________________________________________#
logTRL <- dt$logTRL
logHS <- dt$logHS
logHLL <-dt$logHLL
logFLL <-dt$logFLL
logSVL <-dt$logSVL
logMO <- dt$logMO
logPL <- dt$logPL
sex <- dt$SEX
lineage <- dt$linAK

#TRL
lm.globalTRL <- lm.rrpp(logTRL ~ logSVL*sex*lineage)
anovTRL<-anova(lm.globalTRL)
names(anovTRL)
anovTRL$table
##write table if necessary

#HS
lm.globalHS <- lm.rrpp(logHS ~ logSVL*sex*lineage)
anovHS<-anova(lm.globalHS)
anovHS$table

#HLL
lm.globalHLL <- lm.rrpp(logHLL ~ logSVL*sex*lineage)
anovHLL<-anova(lm.globalHLL)
anovHLL$table
write.table(anovHLL$table,"Allometricslopes_linagesHLL.txt",sep = ",")
summary(lm.globalHLL)

#FLL
lm.globalFLL <- lm.rrpp(logFLL ~ logSVL*sex*lineage)
anovFLL<-anova(lm.globalFLL)
anovFLL$table
write.table(anovFLL$table,"Allometricslopes_linagesFLL.txt",sep = ",")
summary(lm.globalFLL)
