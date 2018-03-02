
##load RDS for given correlation
simulated_df_withadd <- readRDS("simulated_1000rep_interceptcor0.1.RDS")


##subset to first five and make sure below code runs ok
##first check length to make sure created properly


###### Pooled regressions: no effect outcome 

####### first subset to one randomly selected sibling
####### and then run regressions on that subset
simulated_df_withadd_sampleonesib <- lapply(simulated_df_withadd,
                                            function(x) {
                                              shuffle_df <- as.data.table(x)[sample(dim(x)[1])]
                                              return(unique(shuffle_df,
                                                            by = "FID") %>%
                                                       arrange(FID))
                                            })


print("sampled one sib")

pool_lm_noeffect_confound_bothsibs_pvals <- lapply(simulated_df_withadd_sampleonesib,
                                                   function(x){
                                                     return(summary(lm(normal_neithereffect_eithersnp_confounded ~  
                                                                         snp1_additive + sex + age,
                                                                       data = x))$coefficients["snp1_additive",
                                                                                               c(1, 4)])    
                                                   }) 

pool_lm_noeffect_confound_stratcontrol_bothsibs_pvals <- lapply(simulated_df_withadd_sampleonesib,
                                                                function(x){
                                                                  return(summary(lm(normal_neithereffect_eithersnp_confounded ~  
                                                                                      snp1_additive + sex + age + factor(subpop),
                                                                                    data = x))$coefficients["snp1_additive",
                                                                                                            c(1, 4)])}) 

pool_lm_noeffect_confound_bothsibs_pvals_df <- do.call("rbind.data.frame",
                                                       pool_lm_noeffect_confound_bothsibs_pvals) 
colnames(pool_lm_noeffect_confound_bothsibs_pvals_df) <- c("beta", "pval")

pool_lm_noeffect_confound_stratcontrol_bothsibs_pvals_df <- do.call("rbind.data.frame",
                                                                    pool_lm_noeffect_confound_stratcontrol_bothsibs_pvals) 
colnames(pool_lm_noeffect_confound_stratcontrol_bothsibs_pvals_df) <- c("beta", "pval")

print("pooled finished")

###### write results to Rdata file
saveRDS(pool_lm_noeffect_confound_bothsibs_pvals_df, 
        file = "pool_lm_mean_simres_cor01.RDS") 
saveRDS(pool_lm_noeffect_confound_stratcontrol_bothsibs_pvals_df, 
        file = "pool_lm_mean_simres_stratcontrol_cor01.RDS") 


#### Fixed effects regressions: no effects outcome

fe_lm_noeffect_confound_bothsibs_fullreg <- lapply(simulated_df_withadd,
                                                   function(x){
                                                     return(plm(normal_neithereffect_eithersnp_confounded ~  
                                                                  snp1_additive + age  + sex,
                                                                data = x,
                                                                model = "within",
                                                                index = "FID"))})  

fe_lm_noeffect_confound_bothsibs_pvals <- lapply(simulated_df_withadd,
                                                 function(x){
                                                   return(summary(plm(normal_neithereffect_eithersnp_confounded ~  
                                                                        snp1_additive + age  + sex,
                                                                      data = x,
                                                                      model = "within",
                                                                      index = "FID"))$coefficients["snp1_additive",
                                                                                                   c(1, 4)])    
                                                 }) 
fe_lm_noeffect_confound_bothsibs_pvals_df <- do.call("rbind.data.frame",
                                                     fe_lm_noeffect_confound_bothsibs_pvals) 
colnames(fe_lm_noeffect_confound_bothsibs_pvals_df) <- c("beta", "pval")

print("fe finished")
saveRDS(fe_lm_noeffect_confound_bothsibs_pvals_df, file = "fe_lm_mean_simres_cor01.RDS")   


#### Random effects regressions: no effects outcome
re_lm_noeffect_confound_bothsibs_fullreg <-  lapply(simulated_df_withadd,
                                                    function(x){
                                                      return(plm(normal_neithereffect_eithersnp_confounded ~  
                                                                   snp1_additive + sex + age,
                                                                 index = "FID",
                                                                 model = "random",
                                                                 data = x))})
re_lm_noeffect_confound_stratcontrol_bothsibs_fullreg <-  lapply(simulated_df_withadd,
                                                                 function(x){
                                                                   return(plm(normal_neithereffect_eithersnp_confounded ~  
                                                                                snp1_additive + sex + age + factor(subpop),
                                                                              index = "FID", 
                                                                              model = "random",
                                                                              data = x))})

re_lm_noeffect_confound_bothsibs_pvals <-  lapply(simulated_df_withadd,
                                                  function(x){
                                                    return(summary(plm(normal_neithereffect_eithersnp_confounded ~  
                                                                         snp1_additive + sex + age, index = "FID",
                                                                       model = "random",
                                                                       data = x))$coefficients["snp1_additive",
                                                                                               c(1, 4)])})

re_lm_noeffect_confound_stratcontrol_bothsibs_pvals <-  lapply(simulated_df_withadd,
                                                               function(x){
                                                                 return(summary(plm(normal_neithereffect_eithersnp_confounded ~  
                                                                                      snp1_additive + sex + age +
                                                                                      factor(subpop),
                                                                                    index = "FID",
                                                                                    model = "random",
                                                                                    data = x))$coefficients["snp1_additive",
                                                                                                            c(1, 4)])})

re_lm_noeffect_confound_bothsibs_pvals_df <- do.call("rbind.data.frame",
                                                     re_lm_noeffect_confound_bothsibs_pvals) 
colnames(re_lm_noeffect_confound_bothsibs_pvals_df) <- c("beta", "pval")

re_lm_noeffect_confound_stratcontrol_bothsibs_pvals_df <- do.call("rbind.data.frame",
                                                                  re_lm_noeffect_confound_stratcontrol_bothsibs_pvals) 
colnames(re_lm_noeffect_confound_stratcontrol_bothsibs_pvals_df) <- c("beta", "pval")

print("re finished")

#### save the results
saveRDS(re_lm_noeffect_confound_bothsibs_pvals_df, 
        file = "re_lm_mean_simres_cor01.RDS")  
saveRDS(re_lm_noeffect_confound_stratcontrol_bothsibs_pvals_df, 
        file = "re_lm_mean_simres_stratcontrol_cor01.RDS")  

### Conduct a hausman test comparing the outcomes
haus_refe_nopopstratcontrol <- mapply(function(x, y) phtest(x, y),
                                      x = re_lm_noeffect_confound_bothsibs_fullreg,
                                      y = fe_lm_noeffect_confound_bothsibs_fullreg) 


##before running the hausman test with the popstrat
##control included in the RE, 
haus_refe_popstratcontrol <- mapply(function(x, y) phtest(x, y),
                                    x = re_lm_noeffect_confound_stratcontrol_bothsibs_fullreg,
                                    y = fe_lm_noeffect_confound_bothsibs_fullreg) 


### Save the results of the no effects regressions
saveRDS(haus_refe_nopopstratcontrol, file = "haus_refe_nopopstratcontrol_cor01.RDS")
saveRDS(haus_refe_nopopstratcontrol, file = "haus_refe_popstratcontrol_cor01.RDS")

