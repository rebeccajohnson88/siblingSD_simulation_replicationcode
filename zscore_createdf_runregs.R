

################################ STEP ZERO: SETUP ################################

## load packages
library(dplyr)
library(plm)
library(data.table)

## function that 
## takes in data.frame and names
## of outcome variables and creates 
## squared Z score by sex
trait_Zbysex_func <- function(data, traitnames){
  df_trait_z <- as.data.frame(data %>%
                                group_by(factor(sex)) %>%
                                mutate_at(.cols = traitnames,
     .funs = funs(squared_zscore = ((. - mean(.))/sd(.))^2)) %>%
                                ungroup())
  return(df_trait_z)
  
}

##load data with medium level of confounding
simulated_df_withadd_confound <- readRDS("simulated_1000rep_interceptcor0.1.RDS")

print("read in data")


##create vector of trait names
trait_names <- grep("effect", colnames(simulated_df_withadd_confound[[1]]), 
                    value = TRUE)

## apply squared Z function to subset 
simulated_1000rep_interceptcor0.1_withsqZ <- lapply(simulated_df_withadd_confound,
                     function(x) 
                      trait_Zbysex_func(data = x,
                  traitnames = trait_names)) 

print("calculated squared Z for confound df")

## save result
saveRDS(simulated_1000rep_interceptcor0.1_withsqZ,
        "simulated_1000rep_interceptcor0.1_withsqZ.RDS")



## repeat for the data with zero correlation
##load data with medium level of confounding
simulated_df_withadd_noconfound <- readRDS("simulated_1000rep_interceptcor0.RDS")

print("read in no confound data")


trait_names_noconfound <- grep("effect", colnames(simulated_df_withadd_noconfound[[1]]), 
                    value = TRUE)

## apply squared Z function to subset 
simulated_1000rep_interceptcor0_withsqZ <- lapply(simulated_df_withadd_noconfound,
                             function(x)  trait_Zbysex_func(data = x,
                               traitnames = trait_names)) 

print("calculated squared Z for non-confound df")

## save result
saveRDS(simulated_1000rep_interceptcor0_withsqZ,
        "simulated_1000rep_interceptcor0_withsqZ.RDS")


### testing code: commented out before server
## read in the two datasets
## z_noconf <- readRDS("simulated_1000rep_interceptcor0_withsqZ.RDS")
## z_conf <- readRDS("simulated_1000rep_interceptcor0.05_withsqZ.RDS")


dfwithZ_confound <- simulated_1000rep_interceptcor0.1_withsqZ
dfwithZ_noconfound <- simulated_1000rep_interceptcor0_withsqZ



############################## STEP ONE: RUN REGRESSIONS WITH NON-DEMEANED DATA #########

#### neither mean nor variance effects
squaredZreg_null_noconfound_nondemeaned_cont <- do.call("rbind.data.frame", 
                 lapply(dfwithZ_noconfound,
                 function(x){
                  summary(lm(normal_neithereffect_eithersnp_confounded_squared_zscore ~
                        snp1_additive + sex + age + factor(subpop),
                    data = x))$coefficients["snp1_additive", c(1, 4)]
                                                          })) %>%
  mutate(outcome = "neither",
         confounding = "no",
         data = "normal",
         popstratcontrol = "yes",
         sims = 1:length(dfwithZ_noconfound)) 

squaredZreg_null_noconfound_nondemeaned_nocont <- do.call("rbind.data.frame", 
               lapply(dfwithZ_noconfound,
              function(x){
              summary(lm(normal_neithereffect_eithersnp_confounded_squared_zscore ~
              snp1_additive + sex + age,
          data = x))$coefficients["snp1_additive", c(1, 4)]
                                                                 })) %>%
  mutate(outcome = "neither",
         confounding = "no",
         data = "normal",
         popstratcontrol = "no",
         sims = 1:length(dfwithZ_noconfound)) 


squaredZreg_null_confound_nondemeaned_cont <- do.call("rbind.data.frame", 
             lapply(dfwithZ_confound, function(x){ 
              summary(lm(normal_neithereffect_eithersnp_confounded_squared_zscore ~
                  snp1_additive + sex + age +  factor(subpop),
           data = x))$coefficients["snp1_additive", c(1, 4)]  })) %>%
  mutate(outcome = "neither",
         confounding = "yes",
         data = "normal",
         popstratcontrol = "yes",
         sims = 1:length(dfwithZ_noconfound)) 

squaredZreg_null_confound_nondemeaned_nocont <- do.call("rbind.data.frame", 
                        lapply(dfwithZ_confound,
                  function(x){  summary(lm(normal_neithereffect_eithersnp_confounded_squared_zscore ~
                         snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                                             })) %>%
  mutate(outcome = "neither",
         confounding = "yes",
         data = "normal",
         popstratcontrol = "no",
         sims = 1:length(dfwithZ_noconfound)) 

## mean effects only
squaredZreg_mean_noconfound_nondemeaned_cont <- do.call("rbind.data.frame", 
                      lapply(dfwithZ_noconfound, function(x){
                       summary(lm(normal_meaneffect_snp1_confounded_squared_zscore ~
                          snp1_additive + sex + age + factor(subpop),
                       data = x))$coefficients["snp1_additive", c(1, 4)]
                                                          })) %>%
  mutate(outcome = "mean",
         confounding = "no",
         data = "normal",
         popstratcontrol = "yes",
         sims = 1:length(dfwithZ_noconfound)) 

squaredZreg_mean_noconfound_nondemeaned_nocont <- do.call("rbind.data.frame", 
              lapply(dfwithZ_noconfound, 
            function(x){ 
        summary(lm(normal_meaneffect_snp1_confounded_squared_zscore ~
             snp1_additive + sex + age,   data = x))$coefficients["snp1_additive", c(1, 4)]
                                                                 })) %>%
  mutate(outcome = "mean",
         confounding = "no",
         data = "normal",
         popstratcontrol = "no",
         sims = 1:length(dfwithZ_noconfound)) 


squaredZreg_mean_confound_nondemeaned_cont <- do.call("rbind.data.frame", 
            lapply(dfwithZ_confound, function(x){ 
              summary(lm(normal_meaneffect_snp1_confounded_squared_zscore ~
                    snp1_additive + sex + age +  factor(subpop),
                    data = x))$coefficients["snp1_additive", c(1, 4)]  })) %>%
  mutate(outcome = "mean",
         confounding = "yes",
         data = "normal",
         popstratcontrol = "yes",
         sims = 1:length(dfwithZ_noconfound)) 

squaredZreg_mean_confound_nondemeaned_nocont <- do.call("rbind.data.frame", 
      lapply(dfwithZ_confound,   
          function(x){  
        summary(lm(normal_meaneffect_snp1_confounded_squared_zscore ~
        snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                                               })) %>%
  mutate(outcome = "mean",
         confounding = "yes",
         data = "normal",
         popstratcontrol = "no",
         sims = 1:length(dfwithZ_noconfound)) 

## var effects only
squaredZreg_var_noconfound_nondemeaned_cont <- do.call("rbind.data.frame", 
                      lapply(dfwithZ_noconfound, function(x){
             summary(lm(normal_vareffect_snp1_confounded_squared_zscore ~
               snp1_additive + sex + age + factor(subpop),
               data = x))$coefficients["snp1_additive", c(1, 4)]
                                                   })) %>%
  mutate(outcome = "var",
         confounding = "no",
         data = "normal",
         popstratcontrol = "yes",
         sims = 1:length(dfwithZ_noconfound)) 

squaredZreg_var_noconfound_nondemeaned_nocont <- do.call("rbind.data.frame", 
                 lapply(dfwithZ_noconfound,   
                        function(x){ 
                  summary(lm(normal_vareffect_snp1_confounded_squared_zscore ~
                snp1_additive + sex + age,   data = x))$coefficients["snp1_additive", c(1, 4)]
                                                                 })) %>%
  mutate(outcome = "var",
         confounding = "no",
         data = "normal",
         popstratcontrol = "no",
         sims = 1:length(dfwithZ_noconfound)) 


squaredZreg_var_confound_nondemeaned_cont <- do.call("rbind.data.frame", 
           lapply(dfwithZ_confound, function(x) {
             summary(lm(normal_vareffect_snp1_confounded_squared_zscore ~
         snp1_additive + sex + age +  factor(subpop), data = x))$coefficients["snp1_additive", c(1, 4)]  })) %>%
  mutate(outcome = "var",
         confounding = "yes",
         data = "normal",
         popstratcontrol = "yes",
         sims = 1:length(dfwithZ_noconfound)) 

squaredZreg_var_confound_nondemeaned_nocont <- do.call("rbind.data.frame", 
        lapply(dfwithZ_confound,  function(x){
          summary(lm(normal_vareffect_snp1_confounded_squared_zscore ~
              snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                                               })) %>%
  mutate(outcome = "var",
         confounding = "yes",
         data = "normal",
         popstratcontrol = "no",
         sims = 1:length(dfwithZ_noconfound)) 

## mean and var. effects 
squaredZreg_meanvar_noconfound_nondemeaned_cont <- do.call("rbind.data.frame", 
           lapply(dfwithZ_noconfound, function(x){
            summary(lm(normal_meanvareffect_snp1_confounded_squared_zscore ~
           snp1_additive + sex + age + factor(subpop),
           data = x))$coefficients["snp1_additive", c(1, 4)]
                                                  })) %>%
  mutate(outcome = "meanvar",
         confounding = "no",
         data = "normal",
         popstratcontrol = "yes",
         sims = 1:length(dfwithZ_noconfound)) 

squaredZreg_meanvar_noconfound_nondemeaned_nocont <- do.call("rbind.data.frame", 
       lapply(dfwithZ_noconfound,  function(x){ 
         summary(lm(normal_meanvareffect_snp1_confounded_squared_zscore ~
        snp1_additive + sex + age,   data = x))$coefficients["snp1_additive", c(1, 4)]
                                                                })) %>%
  mutate(outcome = "meanvar",
         confounding = "no",
         data = "normal",
         popstratcontrol = "no",
         sims = 1:length(dfwithZ_noconfound)) 


squaredZreg_meanvar_confound_nondemeaned_cont <- do.call("rbind.data.frame", 
              lapply(dfwithZ_confound, function(x) {
            summary(lm(normal_meanvareffect_snp1_confounded_squared_zscore ~
          snp1_additive + sex + age +  factor(subpop), 
          data = x))$coefficients["snp1_additive", c(1, 4)]  })) %>%
  mutate(outcome = "meanvar",
         confounding = "yes",
         data = "normal",
         popstratcontrol = "yes",
         sims = 1:length(dfwithZ_noconfound)) 

squaredZreg_meanvar_confound_nondemeaned_nocont <- do.call("rbind.data.frame", 
               lapply(dfwithZ_confound,  function(x){
                 summary(lm(normal_meanvareffect_snp1_confounded_squared_zscore ~
            snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                                       })) %>%
  mutate(outcome = "meanvar",
         confounding = "yes",
         data = "normal",
         popstratcontrol = "no",
         sims = 1:length(dfwithZ_noconfound)) 


## save non-demeaned results as a list before moving onto demeaning stage
squaredZ_allregslist_nondemean <- list(squaredZreg_null_noconfound_nondemeaned_cont,
                                    squaredZreg_null_noconfound_nondemeaned_nocont,
                                    squaredZreg_null_confound_nondemeaned_cont,
                                    squaredZreg_null_confound_nondemeaned_nocont,
                                    squaredZreg_mean_noconfound_nondemeaned_cont,
                                    squaredZreg_mean_noconfound_nondemeaned_nocont,
                                    squaredZreg_mean_confound_nondemeaned_cont,
                                    squaredZreg_mean_confound_nondemeaned_nocont,
                                    squaredZreg_var_noconfound_nondemeaned_cont,
                                    squaredZreg_var_noconfound_nondemeaned_nocont,
                                    squaredZreg_var_confound_nondemeaned_cont,
                                    squaredZreg_var_confound_nondemeaned_nocont,
                                    squaredZreg_meanvar_noconfound_nondemeaned_cont,
                                    squaredZreg_meanvar_noconfound_nondemeaned_nocont,
                                    squaredZreg_meanvar_confound_nondemeaned_cont,
                                    squaredZreg_meanvar_confound_nondemeaned_nocont)
                                    
## standardize variable names
squaredZ_allregs_nondemean_listrename <- lapply(squaredZ_allregslist_nondemean,
                                             setNames,
                                             c("beta", "p",
                                               "outcome",
                                               "confounding",
                                               "data",
                                               "sims")) 

## write to file
squaredZ_allregs_df_nondemean <- do.call("rbind.data.frame", 
                                      squaredZ_allregs_nondemean_listrename)

print("non-demeaned results combined")

## save results

saveRDS(squaredZ_allregs_df_nondemean, "squaredZ_allregs_non_demean_updated20180214.RDS")

print("non-demeaned results saved")

################################## STEP TWO: DEMEAN DATA AND RUN REGRESSIONS WITH THOSE ##########

## select variables to demean
## and demean those variables for
## the df with confounding

to_demean <- c("sex", "age", "snp1_additive", "FID", 
                                        grep("confounded\\_squared\\_zscore",
                                      colnames(dfwithZ_confound[[1]]), value = TRUE))


dfwithZ_confound_todemean <- lapply(dfwithZ_confound,
                            function(x) x[, to_demean])
dfwithZ_confound_demean <- lapply(dfwithZ_confound_todemean,
                          function(x) 
                        as.data.table(x)[, lapply(.SD, function(x) x - mean(x)), 
                                             by = "FID"])

head(dfwithZ_confound_demean[[1]])

dfwithZ_noconfound_todemean <- lapply(dfwithZ_noconfound,
                                    function(x) x[, to_demean])
dfwithZ_noconfound_demean <- lapply(dfwithZ_noconfound_todemean,
                                  function(x) 
                                    as.data.table(x)[, lapply(.SD, function(x) x - mean(x)), 
                                                     by = "FID"])

## now feed that data to all models

#### neither mean nor variance effects
squaredZreg_null_demeaned <- do.call("rbind.data.frame", 
            lapply(dfwithZ_confound_demean,
                 function(x){
                   summary(lm(normal_neithereffect_eithersnp_confounded_squared_zscore ~
                         snp1_additive + sex + age,
                 data = x))$coefficients["snp1_additive", c(1, 4)]
                                                            })) %>%
  mutate(outcome = "neither",
         data = "demean",
         confound = "yes",
         sims = 1:length(dfwithZ_confound_demean)) 

## mean effects
squaredZreg_mean_demeaned <- do.call("rbind.data.frame", 
                  lapply(dfwithZ_confound_demean, function(x){
             summary(lm(normal_meaneffect_snp1_confounded_squared_zscore ~
           snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                                     })) %>%
  mutate(outcome = "mean",
         data = "demean",
         confound = "yes",
         sims = 1:length(dfwithZ_confound_demean)) 

## variance effects
squaredZreg_var_demeaned <- do.call("rbind.data.frame", 
            lapply(dfwithZ_confound_demean, 
          function(x){summary(lm(normal_vareffect_snp1_confounded_squared_zscore ~
     snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                                    })) %>%
  mutate(outcome = "var",
         data = "demean",
         confound = "yes",
         sims = 1:length(dfwithZ_confound_demean)) 


## both
squaredZreg_meanvar_demeaned <- do.call("rbind.data.frame", 
          lapply(dfwithZ_confound_demean, 
        function(x){summary(lm(normal_meanvareffect_snp1_confounded_squared_zscore
            ~ snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                                        })) %>%
  mutate(outcome = "meanvar",
         data = "demean",
         confound = "yes",
         sims = 1:length(dfwithZ_confound_demean)) 


##### repeat for non-confounded data 

#### neither mean nor variance effects
squaredZreg_null_demeaned_nc <- do.call("rbind.data.frame", 
                                     lapply(dfwithZ_noconfound_demean,
                                            function(x){
                                              summary(lm(normal_neithereffect_eithersnp_confounded_squared_zscore ~
                                                           snp1_additive + sex + age,
                                                         data = x))$coefficients["snp1_additive", c(1, 4)]
                                            })) %>%
  mutate(outcome = "neither",
         data = "demean",
         confound = "no",
         sims = 1:length(dfwithZ_confound_demean)) 

## mean effects
squaredZreg_mean_demeaned_nc <- do.call("rbind.data.frame", 
                                     lapply(dfwithZ_noconfound_demean, function(x){
                                       summary(lm(normal_meaneffect_snp1_confounded_squared_zscore ~
                                                    snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                     })) %>%
  mutate(outcome = "mean",
         data = "demean",
         confound = "no",
         sims = 1:length(dfwithZ_confound_demean)) 


## variance effects
squaredZreg_var_demeaned_nc <- do.call("rbind.data.frame", 
                       lapply(dfwithZ_noconfound_demean, 
                function(x){summary(lm(normal_vareffect_snp1_confounded_squared_zscore ~
                 snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                           })) %>%
  mutate(outcome = "var",
         data = "demean",
         confound = "no",
         sims = 1:length(dfwithZ_confound_demean)) 


## both
squaredZreg_meanvar_demeaned_nc <- do.call("rbind.data.frame", 
                                        lapply(dfwithZ_noconfound_demean, 
                                               function(x){summary(lm(normal_meanvareffect_snp1_confounded_squared_zscore
                                                                      ~ snp1_additive + sex + age, data = x))$coefficients["snp1_additive", c(1, 4)]
                                               })) %>%
  mutate(outcome = "meanvar",
         data = "demean",
         confound = "no",
         sims = 1:length(dfwithZ_confound_demean)) 


########################## STEP THREE: SAVE RESULTS ################ 

## bind into a list and export

### store in a list
squaredZ_allregslist_demean <- list(squaredZreg_mean_demeaned,
                                    squaredZreg_null_demeaned,
                                    squaredZreg_var_demeaned,
                                    squaredZreg_meanvar_demeaned,
                                    squaredZreg_mean_demeaned_nc,
                                    squaredZreg_null_demeaned_nc,
                                    squaredZreg_var_demeaned_nc,
                                    squaredZreg_meanvar_demeaned_nc)

### name all list elements

squaredZ_allregs_demean_listrename <- lapply(squaredZ_allregslist_demean,
                                             setNames,
                                             c("beta", "p",
                                               "outcome",
                                               "data",
                                               "confound",
                                               "sims")) 



## bind into a single data.frame

squaredZ_allregs_df_demean <- do.call("rbind.data.frame", 
                                      squaredZ_allregs_demean_listrename)

print("demeaned results combined")

## save results

saveRDS(squaredZ_allregs_df_demean, "squaredZ_allregs_demean_updated20180214_newdemean.RDS")

print("demeaned results saved")
