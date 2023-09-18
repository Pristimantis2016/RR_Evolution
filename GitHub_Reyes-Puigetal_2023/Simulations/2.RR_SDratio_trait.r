#Authors: Carolina Reyes-Puig and Antigoni Kaliontzopoulou 
#This script is part of the paper Rensch’s Rule: linking intraspecific to evolutionary allometry
#Reyes-Puig, C., Adams, C., Enriquez-Urzelai, U. &  Kaliontzopoulou, A. (2023). Evolution.

library(ape)
library(phytools)
library(nlme)
library(geiger)
library(geomorph)
library(RRPP)
library(dplyr)
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
    dt <- rnorm(n, mean = spF[[i]][[j]], sd = 0.1) 
    dtf <- data.frame(sp = rep(j, n))
    dtf$SVL <- dt
    phen <- rbind(phen, dtf)
  }
  sizeF [[i]] <- phen
}

###ssd = 1 isometry, < 1 FM, > 1 MB###
ssd <- 1 

aM <- 0.583011 # Intercept for males - constant for generating the trait
bM <- 0.25843 # Slope for males - constant for generating the trait
interceptD <- seq(0.8, 1.2, by = 0.004)*(aM) # variation in intercept

###Male data##  
bd.sd <- list()
for (d in ssd) {
  sizeM <- NULL
  for (i in seq_along(spF)){
    phen <- NULL 
    for (j in names(spF[[i]])) {
      dt.F <- sizeF[[i]]
      dt <- rnorm(n, mean = (spF[[i]][[j]] ^ d), sd = 0.1) + rnorm(n = n, 0, 0.1)
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

bd.b0 <- list()
for (a in interceptD) {  
  tempo <- NULL
  for (i in seq_along(bd.sd[[1]])) { 
    bd <- bd.sd[[1]][[i]]
    ma <- bd[which(bd$sex == "M"), ]
    fe <- bd[which(bd$sex == "F"), ]
    traitM <- a * ma$SVL^bM + rnorm(n = n, 0, 0.01)
    dt.M <- cbind(ma, traitM, traitM/ma$SVL) 
    colnames(dt.M) <- c("sp", "SVL","sex", "trait", "std.trait")
    traitF <- aM * fe$SVL^bM + rnorm(n = n, 0, 0.01)
    dt.F <- cbind(fe, traitF, traitF/fe$SVL) 
    colnames(dt.F) <- c("sp", "SVL","sex", "trait", "std.trait")
    temp <- rbind(dt.F,dt.M)
    list.name <- paste("sim = ", i, sep=" ")
    tempo[[list.name]] <- temp
  }
  list.name <- paste("intercept = ", a, sep=" ") 
  bd.b0[[list.name]] <- tempo
}

##RR rule###
#RR.trait <- list()
pval.sim <- NULL
slopes.sim.sd <- NULL
for (a in names(bd.b0)) {
  RRtraiti <- NULL
  prob.T <- NULL
  pval <- NULL
  slopes.sim <- NULL
  for (i in 1:length(bd.b0[[a]])) {
    bd <- bd.b0[[a]][[i]]
    trait <- bd$std.trait
    svl <- bd$SVL
    sp <- bd$sp
    sex <- bd$sex
    #######Calculate means per lineage and sex#####
    lin.sex <- as.factor(paste(sp, sex, sep = "_"))
    sp_means <- lapply(c(c(2,5)), function(x){ 
      stats <- matrix(tapply(bd[,x], lin.sex, mean), 
                      ncol = 2, byrow = T)
      colnames(stats) <- levels(sex)
      rownames(stats) <- levels(as.factor(sp))
      return(stats)
    })
    names(sp_means) <- colnames(bd)[c(2,5)]
    sp_means <- data.frame(sp_means)
    sp_means <- cbind(sp_means, rownames(sp_means))
    names(sp_means)[names(sp_means) == "rownames(sp_means)"] <- "tip.label"
    sp_means$mean.sp <- apply(sp_means[,c(1:2)], 1, mean)
    sp_means$mean.rtr <- apply(sp_means[,c(3:4)], 1, function (x) {x[2]/x[1]})
    lm.SM <- gls(log(std.trait.M/std.trait.F) ~ log(mean.sp) , correlation = corBrownian(phy=my_tree, form = ~ tip.label), data = sp_means)
    sum <- summary(lm.SM)
    T.values <- sum$tTable[,4][2]
    prob.T [[i]] <- unlist(T.values)
    pval  <- c(pval, prob.T[[i]])
    slopes.sim <- c(slopes.sim, coef(lm.SM)[2])
  }
  pval.sim[[a]] <- pval
  slopes.sim.sd[[a]] <- slopes.sim
}

pval.means <- lapply(pval.sim, mean)
pval.means.un <- unlist(pval.means)
pval.means.slop <- unlist(lapply(slopes.sim.sd, mean))
pval.means.slop_sd <- unlist(lapply(slopes.sim.sd, sd))
pval.sd <- unlist(lapply(pval.sim, sd))
CI.sim <- lapply(pval.sim, function (x)  1.96 * (sd(x)/ sqrt(nsim)) )
CI.sim.un <- unlist(CI.sim)

