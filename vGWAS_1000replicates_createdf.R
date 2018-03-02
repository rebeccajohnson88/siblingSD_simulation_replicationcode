
##eventually add seed and round again
##add seed here: (randomly drew from 1:10000000)
set.seed(7276720)

## initialize parameters
## for now, just did one causal snp
###number of snp's (just did two total, since 
###Cao method for generating phenotype 
### from genotype just uses all snps as 
### causal but generates different outcome variables)
Lu <- 0
Lc <- 1
Lcv <- 1

###number of families per subpopulation
## did 250 fams/subpop x 4 subpops x 2 = 2000 total
## 1000 sibling pairs
N <- 250
s <- 4
o <- 2

###range for uniform distribution to generate
###probability of allele A v. G, etc. at locus
###took from parameter values set in other simulation
pmin=0.1
pmax=0.9

### c parameter of subpopulations differing in variance 
c=0.01


##


##create function
generate_genphen <- function(Lu, Lc, Lcv, f, 
                             N, s, pmin,
                             pmax, c,
                             intercept_cor){
  # STEP 1: Generate genotypes of parent generation
  L <- Lu + Lc + Lcv                 # total number of SNPs (unassociated, causal
  # for mean value, causal for variance)
  f <- N*s                           # total number of families
  
  # allele frequencies in ancestral population (drawn from uniform[pmin,pmax] distribution 
  # as per Price et al. 2006 and Wu et al. 2011 who used [0.1,0.9])
  p <- runif(L,min=pmin,max=pmax)    
  p2 <- rbind(p,1-p)                 # need to specify frequency of each allele
  
  # generate parent genotypes (ac=2 for 2 alleles/locus, beta parameter not given
  # because allele frequencies given)
  parents <- simMD(N*2,s,L,p=p2,c.vec1=rep(c,s),ac=2) - 1 
  
  
  # returns f*2 X 2 (diploid) X L matrix of parent genotypes [SNPs 1 to Lu = 
  # unassociated, Lu+1 to Lu+Lc = causal for mean, Lu+Lc+1 to Lu+Lc+Lcv = causal for 
  # variance]
  # original output is 1 or 2 for allele names; changed to 0 or 1 by subtracting 1
  
  # rearrange to 2-D matrix where columns are loci and rows go: parent 1 allele 1,
  # parent 1 allele 2, parent 2 allele 1, etc [parents 1+2 same family, 
  # 3+4 same family, etc]
  parents2 <- matrix(aperm(parents,c(2,1,3)),nrow=2*2*f)  
  
  
  
  # parent genotypes of causal SNPs as 1=AA, 2=AB, 3=BB (used for determining phenotypes)   
  parentsC <- parents2[seq(1,2*f*2,by=2),(Lu+1):L] + 
    parents2[seq(2,2*f*2,by=2),(Lu+1):L] + 1 
  
  
  # STEP 2: Generate genotypes of offspring assuming random mating and independent segregation
  # f*o*2 X L matrix of offspring genotypes [same row structure as parents2 but instead of just 2 parents per family have o offspring per family]
  offspring <- matrix(parents2[cbind((runif(2*f*o*L)<0.5) + rep(rep(seq(1,(f-1)*2*2+1,by=2*2),each=2*o),L) + rep(c(0,2),f*o*L),rep(seq(1:L),each=2*f*o))],nrow=2*f*o)  
  
  # offspring genotypes of causal SNPs as 1=AA, 2=AB, 3=BB (used for determining phenotypes)
  offspringC <- offspring[seq(1,2*f*o,by=2),(Lu+1):L] + 
    offspring[seq(2,2*f*o,by=2),(Lu+1):L] + 1 
  
  # turn into df
  offspring_df_4_Cao <- as.data.frame(offspringC)
  
  # name columns with snps
  colnames(offspring_df_4_Cao) <- paste("snp", seq(1, Lc+ Lcv, by = 1), sep = "")
  
  # add individual ID's and family ID
  offspring_df_4_Cao <- offspring_df_4_Cao %>%
    mutate(IID = seq(from = 1, to = N*s*o,
                     by = 1),
           FID = rep(seq(from = 1, to = N*s,
                         by = 1), 
                     each = o),
           subpop = rep(seq(from = 1, to = s,
                            by = 1), 
                        each = N*o)) 
  
  # recode snps into factor variables
  offspring_factor_alleles <- apply(offspring_df_4_Cao[, grepl("snp",
                                                               colnames(offspring_df_4_Cao))],
                                    2, recode_factor_allele) 
     
  # use dummyvars function within caret package
  # to create dummy indicators for each of the 
  # types of snps 
  ##use function for ade4
  indic_snps <- acm.disjonctif(offspring_factor_alleles)
  colnames(indic_snps) <- paste(rep(grep("snp", 
                                         colnames(offspring_df_4_Cao), 
                                         value = TRUE),
                                    each = 3),
                                gsub("^\\.", "", colnames(indic_snps)),
                                sep = "_")
  indic_snps_df <- cbind.data.frame(offspring_df_4_Cao[,
                                  c("FID", "IID", "subpop")],
                                    indic_snps)
  
  # for now, generate simulated sex and age
  # since used as covariates in Cao phenotype generation
  # step
  gendf <- indic_snps_df %>%
    mutate(sex = sample(c(rep(0, (N*s*o)/2),
                          rep(1, (N*s*o)/2)),
                        N*s*o, replace = FALSE),
           age = rnorm(N*s*o, mean = 50, sd = 10))
  
  # merge in parent count of minor alleles at snp
  parent_df <- as.data.frame(parentsC)
  
  # name columns with snps
  colnames(parent_df) <- paste("snp", seq(1, Lc+ Lcv, by = 1), 
                               sep = "")
  
  # add individual ID's and family ID
  parent_df <- parent_df %>%
    mutate(IID = seq(from = 1, to = N*s*o,
                     by = 1),
           FID = rep(seq(from = 1, to = N*s,
                         by = 1), 
                     each = 2),
           subpop = rep(seq(from = 1, to = s,
                            by = 1), 
                        each = N*o))
  
  
  # obtain minor allele count of parents for each family ID
  
  # right now, parentsC coded as 1 = AA, 2 = AB, 3 = BB
  # recode so that 2 = AA, 1 = AB, and 0 = BB to make
  # summing across both parents easier
  parent_df_minor_recode <- parent_df %>%
    mutate_if(grepl("snp", names(.)),
              funs (ifelse(. == 1, 
                           2, 
                           ifelse(. == 2, 
                                  1, 0))))
  
  print("made past mutate_if")
  
  # generate count of minor alleles for each parent
  # sums across any snp generated 
  parent_sum_minor <- parent_df_minor_recode %>%
    group_by(FID) %>%
    summarize_at(vars(contains("snp")),
                 sum) 
  
  print("made past summarize_at")
  
  
  colnames(parent_sum_minor) <- ifelse(grepl("snp",
                                             colnames(parent_sum_minor)),
                                       paste(colnames(parent_sum_minor),
                                             "minorcount_parents",
                                             sep = "_"),
                                       colnames(parent_sum_minor)) 
  
  # merge with offspring genotype- each offspring
  # in a family will have same count of parental minor alleles
  gendf_withparents <- merge(gendf, parent_sum_minor,
                             by = "FID")
  
  # create count of sibling minor alleles at snp1
  snp1_siblingcount <- gendf %>%
    group_by(factor(FID)) %>%
    summarise(snpcount = sum(2 * snp1_min_h +
                               snp1_hz))
  snp1_siblingcount_vec <- as.vector(snp1_siblingcount$snpcount)
  
  # generate family-level intercept correlated with snp1 count
  # use general cholesky decomposition method 
  # used in samfim where we first generate
  # a non-correalted intercept with mean = 0, sd = 1
  # (startingintercept) and then scale
  # by cholesky decomp of correlation matrix
  # to get new intercept correlated with snp1
  # minor allele count
  cor_mat <- matrix(intercept_cor, nrow = 2, ncol = 2)
  diag(cor_mat) <- 1
  cholesky_cor_mat <- chol(cor_mat)
  startingintercept <- rnorm(f)
  startingintercept_sibcount <- cbind(snp1_siblingcount_vec,
                                      startingintercept)
  # induce correlation using the cholesky decomp
  correlatedintercept_sibcount <- startingintercept_sibcount %*% 
    cholesky_cor_mat
  
  
  # add to data.frame
  gendf_withparents <- gendf_withparents %>%
    mutate(family_intercept = 
             rep(correlatedintercept_sibcount[, 2],
                 each = o))
  
  
  # generate non-confounded and confounded phenotypes
  # according to method outlined in Cao (with modification)
  # for confounding
  gendf_withparents <- gendf_withparents %>%
    mutate(base_outcome = 0.5 * sex + 0.05 * age,
           
           ##neither mean nor variance effects
           normal_neithereffect_eithersnp = base_outcome +
             rnorm(N*s*o, mean = 0, sd = 1),
           
           ##mean effects but not variance effects
           normal_meaneffect_snp1 = base_outcome + 
             0.35*snp1_min_h + 0.15*snp1_hz +
             rnorm(N*s*o, mean = 0, sd = 1),
           
           ##variance effects but not mean effects
           normal_vareffect_snp1 = base_outcome +
             ifelse(snp1_min_h == 1,
                    rnorm(nrow(gendf[gendf$snp1_min_h == 1, ]), 
                          mean = 0, sd = 1.4),
                    ifelse(snp1_hz == 1,
                           rnorm(nrow(gendf[gendf$snp1_hz == 1, ]), mean = 0, 
                                 sd = 1.15),
                           rnorm(nrow(gendf[gendf$snp1_maj_h == 1, ]), mean = 0, 
                                 sd = 1))),
           
           ##mean and variance effects
           normal_meanvareffect_snp1 = base_outcome +
             0.35*snp1_min_h + 0.15*snp1_hz + 
             ifelse(snp1_min_h == 1,
                    rnorm(nrow(gendf[gendf$snp1_min_h == 1, ]), 
                          mean = 0, sd = 1.4),
                    ifelse(snp1_hz == 1,
                           rnorm(nrow(gendf[gendf$snp1_hz == 1, ]), mean = 0, 
                                 sd = 1.15),
                           rnorm(nrow(gendf[gendf$snp1_maj_h == 1, ]), mean = 0, 
                                 sd = 1))),
           
           ##repeating with family-level confounder
           base_outcome_withconfound = 0.5 * sex + 0.05 * age +
             family_intercept,
           
           ##neither mean nor variance effects
           normal_neithereffect_eithersnp_confounded = base_outcome_withconfound +
             rnorm(N*s*o, mean = 0, sd = 1),
           
           ##mean effects but not variance effects
           normal_meaneffect_snp1_confounded = base_outcome_withconfound + 
             0.35*snp1_min_h + 0.15*snp1_hz +
             rnorm(N*s*o, mean = 0, sd = 1),
           
           ##variance effects but not mean effects
           normal_vareffect_snp1_confounded = base_outcome_withconfound +
             ifelse(snp1_min_h == 1,
                    rnorm(nrow(gendf[gendf$snp1_min_h == 1, ]), 
                          mean = 0, sd = 1.4),
                    ifelse(snp1_hz == 1,
                           rnorm(nrow(gendf[gendf$snp1_hz == 1, ]), mean = 0, 
                                 sd = 1.15),
                           rnorm(nrow(gendf[gendf$snp1_maj_h == 1, ]), mean = 0, 
                                 sd = 1))),
           
           ##mean and variance effects
           normal_meanvareffect_snp1_confounded = base_outcome_withconfound +
             0.35*snp1_min_h + 0.15*snp1_hz + 
             ifelse(snp1_min_h == 1,
                    rnorm(nrow(gendf[gendf$snp1_min_h == 1, ]), 
                          mean = 0, sd = 1.4),
                    ifelse(snp1_hz == 1,
                           rnorm(nrow(gendf[gendf$snp1_hz == 1, ]), mean = 0, 
                                 sd = 1.15),
                           rnorm(nrow(gendf[gendf$snp1_maj_h == 1, ]), mean = 0, 
                                 sd = 1)))) 
  
  return(gendf_withparents)  
}


##run for correlation = 0.1
nsims <- 1000
simulated_df <- list(length = nsims)
for(i in 1:nsims){
  simulated_df[[i]] <- generate_genphen(Lu = Lu, Lc = Lc, Lcv = Lcv, 
                                        f = f, 
                                        N = N, s = s, 
                                        pmin = pmin,
                                        pmax = pmax, c = c,
                                        intercept_cor = 0.1)
  print(paste("Generated ",
              round(i/nsims, 3), " of iterations",
              sep = ""))
}

##generate additive coding for snp
simulated_df_withadd_0.1 <- lapply(simulated_df, function(x){
  return(updated_df <- x %>%
           mutate(snp1_additive = 2 * snp1_min_h +
                    snp1_hz)) 
})


##save resulting data
saveRDS(simulated_df_withadd_0.1,
        file = "simulated_1000rep_interceptcor0.1.RDS")


##run for correlation = 0.05
simulated_df <- list(length = nsims)
for(i in 1:nsims){
  simulated_df[[i]] <- generate_genphen(Lu = Lu, Lc = Lc, Lcv = Lcv, 
                                        f = f, 
                                        N = N, s = s, 
                                        pmin = pmin,
                                        pmax = pmax, c = c,
                                        intercept_cor = 0.05)
  print(paste("Generated ",
              round(i/nsims, 3), " of iterations",
              sep = ""))
}

##generate additive coding for snp
simulated_df_withadd_0.05 <- lapply(simulated_df, function(x){
  return(updated_df <- x %>%
           mutate(snp1_additive = 2 * snp1_min_h +
                    snp1_hz)) 
})


##save resulting data
saveRDS(simulated_df_withadd_0.05,
        file = "simulated_1000rep_interceptcor0.05.RDS")

print("saved second RData object")

## run for correlation 0.01
##run for correlation = 0.1
nsims <- 1000
simulated_df <- list(length = nsims)
for(i in 1:nsims){
  simulated_df[[i]] <- generate_genphen(Lu = Lu, Lc = Lc, Lcv = Lcv, 
                                        f = f, 
                                        N = N, s = s, 
                                        pmin = pmin,
                                        pmax = pmax, c = c,
                                        intercept_cor = 0.1)
  print(paste("Generated ",
              round(i/nsims, 3), " of iterations",
              sep = ""))
}

##generate additive coding for snp
simulated_df_withadd_0.1 <- lapply(simulated_df, function(x){
  return(updated_df <- x %>%
           mutate(snp1_additive = 2 * snp1_min_h +
                    snp1_hz)) 
})


##save resulting data
saveRDS(simulated_df_withadd_0.1,
        file = "simulated_1000rep_interceptcor0.1.RDS")


##run for correlation = 0.01
simulated_df <- list(length = nsims)
for(i in 1:nsims){
  simulated_df[[i]] <- generate_genphen(Lu = Lu, Lc = Lc, Lcv = Lcv, 
                                        f = f, 
                                        N = N, s = s, 
                                        pmin = pmin,
                                        pmax = pmax, c = c,
                                        intercept_cor = 0.01)
  print(paste("Generated ",
              round(i/nsims, 3), " of iterations",
              sep = ""))
}

##generate additive coding for snp
simulated_df_withadd_0.01 <- lapply(simulated_df, function(x){
  return(updated_df <- x %>%
           mutate(snp1_additive = 2 * snp1_min_h +
                    snp1_hz)) 
})


##save resulting data
saveRDS(simulated_df_withadd_0.01,
        file = "simulated_1000rep_interceptcor0.01.RDS")

print("saved third RData object")


##run for correlation = 0
simulated_df <- list(length = nsims)
for(i in 1:nsims){
  simulated_df[[i]] <- generate_genphen(Lu = Lu, Lc = Lc, Lcv = Lcv, 
                                        f = f, 
                                        N = N, s = s, 
                                        pmin = pmin,
                                        pmax = pmax, c = c,
                                        intercept_cor = 0)
  print(paste("Generated ",
              round(i/nsims, 3), " of iterations",
              sep = ""))
}


simulated_df_withadd_0 <- lapply(simulated_df, function(x){
  return(updated_df <- x %>%
           mutate(snp1_additive = 2 * snp1_min_h +
                    snp1_hz)) 
})

saveRDS(simulated_df_withadd_0,
        file = "simulated_1000rep_interceptcor0.RDS")

print("saved fourth RData object")


