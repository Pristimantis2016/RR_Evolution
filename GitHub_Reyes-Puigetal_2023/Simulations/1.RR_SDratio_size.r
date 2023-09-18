#Authors: Carolina Reyes-Puig and Antigoni Kaliontzopoulou 
#This script is part of the paper Rensch’s Rule: linking intraspecific to evolutionary allometry
#Reyes-Puig, C., Adams, C., Enriquez-Urzelai, U. &  Kaliontzopoulou, A. (2023). Evolution.

library(ape)
library(phytools)
library(nlme)
library(geiger)
library(geomorph)
library(RRPP)
library(tidyverse)
rm(list = ls())

ntips <- 50
my_tree <- pbtree(n=ntips)
nsim <- 1000

## Fast BM of female means body size ##

spF <- lapply(1:nsim, function(x) fastBM(tree = my_tree,mu = 4,a = 0))

## rnorm for individuals per species ##
n <- 50 
sizeF <- list()
for (i in seq_along(spF)){
  phen <- NULL
  for (j in names(spF[[i]])) {
    dt <- rnorm(n, mean = spF[[i]][j], sd = 0.1)
    dtf <- data.frame(sp = rep(j, n))
    dtf$SVL <- dt
    phen <- rbind(phen, dtf)
  }
  sizeF [[i]] <- phen
}


###variation in SSD### 
ssd <- seq(0.95, 1.05, by = 0.001)

###Male data, males will be created from Huxley equation##  
bd.sd <- list()
for (d in ssd) {
  sizeM <- NULL
  for (i in seq_along(spF)){
    phen <- NULL 
    for (j in names(spF[[i]])) {
      dt.F <- sizeF[[i]]
      dt <- rnorm(n, mean = (spF[[i]][j]^d), sd = 0.1)
      dtm <- data.frame(sp = rep(j, n))
      dtm$SVL <- dt
      phen <- rbind(phen, dtm)
    }
    bd <- rbind(dt.F, phen)
    bd$sex <- as.factor(rep(c("F", "M"), each = nrow(phen)))
    sizeM[[i]] <- bd 
  }
  list.name <- paste("ssd = ", d, sep="")
  bd.sd[[list.name]] <- sizeM
}

##RR rule###
#RRsize.SM <- list()
pval.sim.SM <- NULL
prob.T <- NULL
slopes.sim.sd <- NULL
for (d in names(bd.sd)) {
  #RRsize <- NULL
  pval.sim <- NULL
  slopes.sim <- NULL
  for (i in seq_along(bd.sd[[d]])) {
    bd <- bd.sd[[d]][[i]]
    svl <- bd$SVL
    sp <- bd$sp
    sex <- bd$sex
    ###Means per lineage, per sex###
    lin.sex <- as.factor(paste(sp, sex, sep = "_"))
    sp_means <- matrix(tapply(bd[,2], lin.sex, mean), 
                      ncol = 2, byrow = T)
    colnames(sp_means) <- levels(sex)
    rownames(sp_means) <- levels(as.factor(sp))
    sp_means <- data.frame(sp_means)
    sp_means <- cbind(sp_means, rownames(sp_means))
    names(sp_means[,2]) <- "tip.label"
    names(sp_means)[names(sp_means) == "rownames(sp_means)"] <- "tip.label"
    sp_means$mean.sp <- apply(sp_means[,c(1:2)], 1, mean)
    sp_means$ratio <- apply(sp_means[,c(1:2)], 1, function (x) {x[2]/x[1]})
    lm.SD <- gls(log(ratio) ~ log((mean.sp)), correlation = corBrownian(phy=my_tree, form = ~ tip.label), data = sp_means)
    sum <- summary(lm.SD)
    T.values <- sum$tTable[,4][2]
    prob.T [[i]] <- unlist(T.values)
    pval.sim  <- c(pval.sim, prob.T[[i]])
    slopes.sim <- c(slopes.sim, coef(lm.SM) [2])
  }
  pval.sim.SM[[d]] <- pval.sim
  slopes.sim.sd[[d]] <- slopes.sim 
}

pval.means <- lapply(pval.sim.SM, mean)
pval.means.un <- unlist(pval.means)
pval.means.slop <- unlist(lapply(slopes.sim.sd, mean))
pval.means.slop_sd <- unlist(lapply(slopes.sim.sd, sd))
pval.sd <- unlist(lapply(pval.sim.SM, sd))
CI.sim <- lapply(pval.sim.SM, function (x) { (1.96 * sd(x))/ sqrt(nsim) })
CI.sim.un <- unlist(CI.sim)

