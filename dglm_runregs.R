

##load packages
library(dplyr)
library(dglm)
library(data.table)

##load data
simulated_df_withadd_confound <- readRDS("simulated_1000rep_interceptcor0.1.RDS")
simulated_df_withadd_noconfound <- readRDS("simulated_1000rep_interceptcor0.RDS")

## run DGLM on:
## 1. NON-confounded null, mean, variance
## 2. WITHOUT demeaning
nulloutcome_dglm_noconfound <- lapply(simulated_df_withadd_noconfound,
                                      function(x) {
                                        dglm(normal_neithereffect_eithersnp_confounded ~ 
                                               snp1_additive + sex + age, 
                                             ~ snp1_additive + sex + age,
                                             data=x,
                                             family= gaussian(link = identity),
                                             method = "REML") })
meanoutcome_dglm_noconfound <- lapply(simulated_df_withadd_noconfound,
                                    function(x) {
                                      dglm(normal_meaneffect_snp1_confounded ~ 
                                      snp1_additive + sex + age, 
                                    ~ snp1_additive + sex + age,
                                    data=x,
                                    family= gaussian(link = identity),
                                    method = "REML") })
varoutcome_dglm_noconfound <- lapply(simulated_df_withadd_noconfound,
                                      function(x) {
                                        dglm(normal_vareffect_snp1_confounded ~ 
                                               snp1_additive + sex + age, 
                                             ~ snp1_additive + sex + age,
                                             data=x,
                                             family= gaussian(link = identity),
                                             method = "REML") })


## run DGLM on:
## 1. CONFOUNDED null, mean, and variance 
## 2. without demeaning
nulloutcome_dglm_confound <- lapply(simulated_df_withadd_confound,
                                      function(x) {
                                        dglm(normal_neithereffect_eithersnp_confounded ~ 
                                               snp1_additive + sex + age, 
                                             ~ snp1_additive + sex + age,
                                             data=x,
                                             family= gaussian(link = identity),
                                             method = "REML") })
meanoutcome_dglm_confound <- lapply(simulated_df_withadd_confound,
                                      function(x) {
                                        dglm(normal_meaneffect_snp1_confounded ~ 
                                               snp1_additive + sex + age, 
                                             ~ snp1_additive + sex + age,
                                             data=x,
                                             family= gaussian(link = identity),
                                             method = "REML") })
varoutcome_dglm_confound <- lapply(simulated_df_withadd_confound,
                                     function(x) {
                                       dglm(normal_vareffect_snp1_confounded  ~ 
                                              snp1_additive + sex + age, 
                                            ~ snp1_additive + sex + age,
                                            data=x,
                                            family= gaussian(link = identity),
                                            method = "REML") })


## Extract coefs for mean and variance parameters
## for each of the simulated DVs and save
## 1. NON-CONFOUNDED

### extract coefs for null DV
nulloutcome_dglm_noconfound_meancoef <- do.call("rbind.data.frame",
                                                lapply(nulloutcome_dglm_noconfound,
                                                       function(x) {
                                                         summary(x)$coefficients["snp1_additive",
                                                                                 c(1, 4)]
                                                       })) %>%
  mutate(coeftype = "mean",
         outcome = "null",
         confounding = "no",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))
nulloutcome_dglm_noconfound_varcoef <- do.call("rbind.data.frame",
                                               lapply(nulloutcome_dglm_noconfound,
                                                      function(x) {
                                                        summary(x$dispersion.fit)$coefficients["snp1_additive",
                                                                                               c(1, 4)]
                                                      })) %>%
  mutate(coeftype = "var",
         outcome = "null",
         confounding = "no",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))


### extract coefs for mean DV
meanoutcome_dglm_noconfound_meancoef <- do.call("rbind.data.frame",
                                      lapply(meanoutcome_dglm_noconfound,
                                      function(x) {
                                        summary(x)$coefficients["snp1_additive",
                                            c(1, 4)]
                                      })) %>%
                                      mutate(coeftype = "mean",
                                             outcome = "mean",
                                             confounding = "no",
                                             data = "normal",
                                             sims = 1:length(simulated_df_withadd_noconfound))
meanoutcome_dglm_noconfound_varcoef <- do.call("rbind.data.frame",
                                      lapply(meanoutcome_dglm_noconfound,
                                        function(x) {
                                          summary(x$dispersion.fit)$coefficients["snp1_additive",
                                             c(1, 4)]
                                               })) %>%
                                        mutate(coeftype = "var",
                                        outcome = "mean",
                                        confounding = "no",
                                        data = "normal",
                                        sims = 1:length(simulated_df_withadd_noconfound))


### extract coefs for variance DV
varoutcome_dglm_noconfound_meancoef <- do.call("rbind.data.frame",
                       lapply(varoutcome_dglm_noconfound,
                        function(x) {
                      summary(x)$coefficients["snp1_additive",
                           c(1, 4)] })) %>%
  mutate(coeftype = "mean",
         outcome = "var",
         confounding = "no",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))
varoutcome_dglm_noconfound_varcoef <- do.call("rbind.data.frame",
          lapply(varoutcome_dglm_noconfound,
           function(x) {
            summary(x$dispersion.fit)$coefficients["snp1_additive",
           c(1, 4)]
            })) %>%
  mutate(coeftype = "var",
         outcome = "var",
         confounding = "no",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))


## Extract coefs for mean and variance parameters
## for each of the simulated DVs and save
## 2. CONFOUNDED

### null outcome
nulloutcome_dglm_confound_meancoef <- do.call("rbind.data.frame",
                                                lapply(nulloutcome_dglm_confound,
                                                       function(x) {
                                                         summary(x)$coefficients["snp1_additive",
                                                                                 c(1, 4)]
                                                       })) %>%
  mutate(coeftype = "mean",
         outcome = "null",
         confounding = "yes",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))
nulloutcome_dglm_confound_varcoef <- do.call("rbind.data.frame",
                                               lapply(nulloutcome_dglm_confound,
                                                      function(x) {
                                                        summary(x$dispersion.fit)$coefficients["snp1_additive",
                                                                                               c(1, 4)]
                                                      })) %>%
  mutate(coeftype = "var",
         outcome = "null",
         confounding = "yes",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))


meanoutcome_dglm_confound_meancoef <- do.call("rbind.data.frame",
                           lapply(meanoutcome_dglm_confound,
                          function(x) {
                           summary(x)$coefficients["snp1_additive",
                            c(1, 4)]
                             })) %>%
  mutate(coeftype = "mean",
         outcome = "mean",
         confounding = "yes",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))
meanoutcome_dglm_confound_varcoef <- do.call("rbind.data.frame",
                lapply(meanoutcome_dglm_confound,
                 function(x) {
                 summary(x$dispersion.fit)$coefficients["snp1_additive",
                 c(1, 4)]
                })) %>%
  mutate(coeftype = "var",
         outcome = "mean",
         confounding = "yes",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))
varoutcome_dglm_confound_meancoef <- do.call("rbind.data.frame",
              lapply(varoutcome_dglm_confound,
              function(x) {
              summary(x)$coefficients["snp1_additive",
              c(1, 4)] })) %>%
  mutate(coeftype = "mean",
         outcome = "var",
         confounding = "yes",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))
varoutcome_dglm_confound_varcoef <- do.call("rbind.data.frame",
          lapply(varoutcome_dglm_confound,
          function(x) {
          summary(x$dispersion.fit)$coefficients["snp1_additive",
           c(1, 4)] })) %>%
  mutate(coeftype = "var",
         outcome = "var",
         confounding = "yes",
         data = "normal",
         sims = 1:length(simulated_df_withadd_noconfound))


var_reg_results_nondemean_list <- list(nulloutcome_dglm_noconfound_meancoef,
                                      nulloutcome_dglm_noconfound_varcoef,
                                      meanoutcome_dglm_noconfound_meancoef,
                                       meanoutcome_dglm_noconfound_varcoef,
                                      varoutcome_dglm_noconfound_meancoef,
                                       varoutcome_dglm_noconfound_varcoef,
                                      nulloutcome_dglm_confound_meancoef,
                                      nulloutcome_dglm_confound_varcoef,
                                      meanoutcome_dglm_confound_meancoef,
                                      meanoutcome_dglm_confound_varcoef,
                                       varoutcome_dglm_confound_meancoef,
                                       varoutcome_dglm_confound_varcoef)

var_reg_results_nondemean_list_rename <- lapply(var_reg_results_nondemean_list,
                                                setNames,
                                                c("beta", "p",
                                                  "coeftype",
                                                  "outcome",
                                                  "confounding",
                                                  "data",
                                                  "sims")) 

var_reg_results_nondemean <- do.call("rbind.data.frame", 
                                     var_reg_results_nondemean_list_rename)

print("non-demeaned regressions finished")
saveRDS(var_reg_results_nondemean, 
        "var_reg_results_nondemean_20180214_updated.RDS")


############################# REPEAT ABOVE WITH DEMEANED DATA #########################

##repeat all above with demeaned data
##create vector of variables to demean
to_demean <- c("snp1_additive", "normal_meaneffect_snp1_confounded",
               "normal_vareffect_snp1_confounded",
               "normal_neithereffect_eithersnp_confounded",
               "sex", "age")


## demean each of the two datasets
df_noconfound_demean <- lapply(simulated_df_withadd_noconfound,
                                      function(x) {
                                        x[, 
                                          c('FID', to_demean)] %>%
                                          group_by(factor(FID)) %>%
                                          mutate_at(.cols = to_demean,
                                                    .funs = funs(. - mean(.))) })  

df_confound_demean <- lapply(simulated_df_withadd_confound,
                               function(x) {
                                 x[, 
                                   c('FID', to_demean)] %>%
                                   group_by(factor(FID)) %>%
                                   mutate_at(.cols = to_demean,
                                             .funs = funs(. - mean(.))) })  


## use those demeaned df in subsequent regressions 
## run DGLM on:
## 1. NON-confounded null, mean, variance
## 2. WITHOUT demeaning
nulloutcome_dglm_noconfound_demean <- lapply(df_noconfound_demean,
                                             function(x) {
                                               dglm(normal_neithereffect_eithersnp_confounded ~ 
                                                      snp1_additive + sex + age, 
                                                    ~ snp1_additive + sex + age,
                                                    data=x,
                                                    family= gaussian(link = identity),
                                                    method = "REML") })
meanoutcome_dglm_noconfound_demean <- lapply(df_noconfound_demean,
                                             function(x) {
                                               dglm(normal_meaneffect_snp1_confounded ~ 
                                                      snp1_additive + sex + age, 
                                                    ~ snp1_additive + sex + age,
                                                    data=x,
                                                    family= gaussian(link = identity),
                                                    method = "REML") })
varoutcome_dglm_noconfound_demean <- lapply(df_noconfound_demean,
                                            function(x) {
                                              dglm(normal_vareffect_snp1_confounded ~ 
                                                     snp1_additive + sex + age, 
                                                   ~ snp1_additive + sex + age,
                                                   data=x,
                                                   family= gaussian(link = identity),
                                                   method = "REML") })


## run DGLM on:
## 1. CONFOUNDED null, mean, and variance 
## 2. without demeaning
nulloutcome_dglm_confound_demean <- lapply(df_confound_demean,
                                           function(x) {
                                             dglm(normal_neithereffect_eithersnp_confounded ~ 
                                                    snp1_additive + sex + age, 
                                                  ~ snp1_additive + sex + age,
                                                  data=x,
                                                  family= gaussian(link = identity),
                                                  method = "REML") })
meanoutcome_dglm_confound_demean <- lapply(df_confound_demean,
                                           function(x) {
                                             dglm(normal_meaneffect_snp1_confounded ~ 
                                                    snp1_additive + sex + age, 
                                                  ~ snp1_additive + sex + age,
                                                  data=x,
                                                  family= gaussian(link = identity),
                                                  method = "REML") })
varoutcome_dglm_confound_demean <- lapply(df_confound_demean,
                                          function(x) {
                                            dglm(normal_vareffect_snp1_confounded  ~ 
                                                   snp1_additive + sex + age, 
                                                 ~ snp1_additive + sex + age,
                                                 data=x,
                                                 family= gaussian(link = identity),
                                                 method = "REML") })


## Extract coefs for mean and variance parameters
## for each of the simulated DVs and save
## 1. NON-CONFOUNDED

### extract coefs for null DV
nulloutcome_dglm_noconfound_meancoef_demean <- do.call("rbind.data.frame",
                                                       lapply(nulloutcome_dglm_noconfound_demean,
                                                              function(x) {
                                                                summary(x)$coefficients["snp1_additive",
                                                                                        c(1, 4)]
                                                              })) %>%
  mutate(coeftype = "mean",
         outcome = "null",
         confounding = "no",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))
nulloutcome_dglm_noconfound_varcoef_demean <- do.call("rbind.data.frame",
                                                      lapply(nulloutcome_dglm_noconfound_demean,
                                                             function(x) {
                                                               summary(x$dispersion.fit)$coefficients["snp1_additive",
                                                                                                      c(1, 4)]
                                                             })) %>%
  mutate(coeftype = "var",
         outcome = "null",
         confounding = "no",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))


### extract coefs for mean DV
meanoutcome_dglm_noconfound_meancoef_demean <- do.call("rbind.data.frame",
                                                       lapply(meanoutcome_dglm_noconfound_demean,
                                                              function(x) {
                                                                summary(x)$coefficients["snp1_additive",
                                                                                        c(1, 4)]
                                                              })) %>%
  mutate(coeftype = "mean",
         outcome = "mean",
         confounding = "no",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))
meanoutcome_dglm_noconfound_varcoef_demean <- do.call("rbind.data.frame",
                                                      lapply(meanoutcome_dglm_noconfound_demean,
                                                             function(x) {
                                                               summary(x$dispersion.fit)$coefficients["snp1_additive",
                                                                                                      c(1, 4)]
                                                             })) %>%
  mutate(coeftype = "var",
         outcome = "mean",
         confounding = "no",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))


### extract coefs for variance DV
varoutcome_dglm_noconfound_meancoef_demean <- do.call("rbind.data.frame",
                                                      lapply(varoutcome_dglm_noconfound_demean,
                                                             function(x) {
                                                               summary(x)$coefficients["snp1_additive",
                                                                                       c(1, 4)] })) %>%
  mutate(coeftype = "mean",
         outcome = "var",
         confounding = "no",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))
varoutcome_dglm_noconfound_varcoef_demean <- do.call("rbind.data.frame",
                                                     lapply(varoutcome_dglm_noconfound_demean,
                                                            function(x) {
                                                              summary(x$dispersion.fit)$coefficients["snp1_additive",
                                                                                                     c(1, 4)]
                                                            })) %>%
  mutate(coeftype = "var",
         outcome = "var",
         confounding = "no",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))


## Extract coefs for mean and variance parameters
## for each of the simulated DVs and save
## 2. CONFOUNDED

### null outcome
nulloutcome_dglm_confound_meancoef_demean <- do.call("rbind.data.frame",
                                                     lapply(nulloutcome_dglm_confound_demean,
                                                            function(x) {
                                                              summary(x)$coefficients["snp1_additive",
                                                                                      c(1, 4)]
                                                            })) %>%
  mutate(coeftype = "mean",
         outcome = "null",
         confounding = "yes",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))
nulloutcome_dglm_confound_varcoef_demean <- do.call("rbind.data.frame",
                                                    lapply(nulloutcome_dglm_confound_demean,
                                                           function(x) {
                                                             summary(x$dispersion.fit)$coefficients["snp1_additive",
                                                                                                    c(1, 4)]
                                                           })) %>%
  mutate(coeftype = "var",
         outcome = "null",
         confounding = "yes",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))


meanoutcome_dglm_confound_meancoef_demean <- do.call("rbind.data.frame",
                                                     lapply(meanoutcome_dglm_confound_demean,
                                                            function(x) {
                                                              summary(x)$coefficients["snp1_additive",
                                                                                      c(1, 4)]
                                                            })) %>%
  mutate(coeftype = "mean",
         outcome = "mean",
         confounding = "yes",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))
meanoutcome_dglm_confound_varcoef_demean <- do.call("rbind.data.frame",
                                                    lapply(meanoutcome_dglm_confound_demean,
                                                           function(x) {
                                                             summary(x$dispersion.fit)$coefficients["snp1_additive",
                                                                                                    c(1, 4)]
                                                           })) %>%
  mutate(coeftype = "var",
         outcome = "mean",
         confounding = "yes",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))
varoutcome_dglm_confound_meancoef_demean <- do.call("rbind.data.frame",
                                                    lapply(varoutcome_dglm_confound_demean,
                                                           function(x) {
                                                             summary(x)$coefficients["snp1_additive",
                                                                                     c(1, 4)] })) %>%
  mutate(coeftype = "mean",
         outcome = "var",
         confounding = "yes",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))
varoutcome_dglm_confound_varcoef_demean <- do.call("rbind.data.frame",
                                                   lapply(varoutcome_dglm_confound_demean,
                                                          function(x) {
                                                            summary(x$dispersion.fit)$coefficients["snp1_additive",
                                                                                                   c(1, 4)] })) %>%
  mutate(coeftype = "var",
         outcome = "var",
         confounding = "yes",
         data = "demeaned",
         sims = 1:length(df_noconfound_demean))


var_reg_results_demean_list <- list(nulloutcome_dglm_noconfound_meancoef_demean,
                                    nulloutcome_dglm_noconfound_varcoef_demean,
                                    meanoutcome_dglm_noconfound_meancoef_demean,
                                    meanoutcome_dglm_noconfound_varcoef_demean,
                                    varoutcome_dglm_noconfound_meancoef_demean,
                                    varoutcome_dglm_noconfound_varcoef_demean,
                                    nulloutcome_dglm_confound_meancoef_demean,
                                    nulloutcome_dglm_confound_varcoef_demean,
                                    meanoutcome_dglm_confound_meancoef_demean,
                                    meanoutcome_dglm_confound_varcoef_demean,
                                    varoutcome_dglm_confound_meancoef_demean,
                                    varoutcome_dglm_confound_varcoef_demean)

var_reg_results_demean_list_rename <- lapply(var_reg_results_demean_list,
                                             setNames,
                                             c("beta", "p",
                                               "coeftype",
                                               "outcome",
                                               "confounding",
                                               "data",
                                               "sims")) 

var_reg_results_demean <- do.call("rbind.data.frame", 
                                  var_reg_results_demean_list_rename)

print("demeaned regressions finished")
saveRDS(var_reg_results_demean, 
        "var_reg_results_demean_20180215_updated.RDS")

