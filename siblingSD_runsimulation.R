


############################### SETUP: load packages and functions #####################

# load relevant packages
library(dplyr)
library(abind)
library(reshape2)
#library(caret)
library(stats)
library(stargazer)
#library(viridis)
library(stringr)
library(lme4)
library(ade4)
library(plm)
library(scales)
library(ggplot2)
library(data.table)

# load functions used at various points

##load other functions
##simMD function used to create family
## genotypes
simMD <- function (N, P, L, p = NULL, c.vec1, c.vec2 = 1, ac = 2, beta = 1) 
{
  rdirichlet <- function(n, a) {
    l <- length(a)
    x <- matrix(rgamma(l * n, a), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
  }
  rmultinom <- function(p) length(p) - sum(runif(1) < cumsum(p)) + 1
  if (is.null(p)) {
    p <- matrix(0, ac, L)
    for (l in 1:L) p[, l] <- rdirichlet(1, rep(beta, ac))
  }
  X <- array(0, dim = c(N * P * length(c.vec2), 2, L))
  for (nn in 1:P) {
    for (l in 1:L) {
      temp1 <- rdirichlet(1, p[, l] * (1 - c.vec1[nn])/c.vec1[nn])
      for (j in 1:length(c.vec2)) {
        if (length(c.vec2 == 1)) {
          temp2 = temp1
        }
        else {
          temp2 = rdirichlet(1, temp1 * (1 - c.vec2[j])/c.vec2[j])
        }
        for (i in 1:N) 
          for (a in 1:2) {
            X[i + (nn - 1) * N * length(c.vec2) + (j - 1) * N, a, l] <- rmultinom(temp2)
          }
      }
    }
  }
  return(X)
}

##function to convert 1, 2, 3 into factor variable
##for what they correspond to, and relevel
##so that reference category is major homozygotes
recode_factor_allele <- function(variable){
  return(factor(variable,
                levels = c(1, 2, 3),
                labels = c("min_h", "hz",
                           "maj_h")))
}

##function to generalize summing minor alleles
##across snps
minorallele_sum <- function(snpname){
  return(paste("2*", paste(snpname, "min_h_offspring1", sep = "_"),
               "+", "2*", paste(snpname, "min_h_offspring2", sep = "_"),
               "+", paste(snpname, "hz_offspring1", sep = "_"),
               "+", paste(snpname, "hz_offspring2", sep = "_"),
               sep = ""))
}

## script to check location of current script
## and then source resulting scripts from there

LocationOfThisScript <- function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file)) return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}
current.dir <- LocationOfThisScript()

########################## STEP ONE: GENERATE GENOTYPES AND OUTCOMES ##############
########################## WITH GIVEN CORRELATION STRUCTURE AND SIZE ##############

## source the script that runs the 
## code to generate genotypes of a given sample size

source("vGWAS_1000replicates_createdf.R") 


######################## STEP TWO: RUN DIFFERENT TYPES OF REGRESSIONS #############
######################### ON MEAN OUTCOMES TO SHOW BIASED ESTIMATES AT NON-ZERO COR ###

## run mean regs with correlation 0
source("meanregs_cor0.R")

## run mean regs with correlation 0.01
source("meanregs_cor001.R")

## run mean regs with correlation 0.05
source("meanregs_cor005.R")

## run mean regs with correlation 0.1
source("meanregs_cor01.R")

## creates files with:
## 1. coefficient estimates, and 
## 2. hausman test results

## can then create plot shown in paper
## using estimates and hausman test results 



#################### STEP THREE: RUN SQUARED Z-SCORE METHOD ##########
#################### ON DEMEANED AND NON-DEMEANED DATA ###############



## run squared Z-score with two datasets:
## 1. data set with no confounding (cor0)
## 2. dataset with medium confounding (cor 0.05)

## for each dataset, run
## 1. with and without controls for ancestry (non-demeaned)
## 2. without controls for ancestry (demeaned)


source("zscore_createdf_runregs.R") 

## creates files with results of squared Z-score 
## regressions on non-demeaned and demeaned 
## data for supplementary tables


############################ STEP FOUR: RUN DGLM METHOD #####################
############################ ON NON-DEMEANED AND DEMEANED DATA ###############


## run DGLM models

source("dglm_runregs.R")


########################## STEP FIVE: RUN SIBLING SD METHOD ############
###########################################################################


## run script to create sibling-format data
## with sibling mean, SD, minor allele count,
## and parent mean, SD, minor allele count
source("siblingsd_createdf.R")

## run script to estimate models on the
## sibling-format data
source("siblingsd_runregs.R")


## See other code file: siblingSD_finalfigures_tables
## that loads the output of above
## into R and creates the figures













