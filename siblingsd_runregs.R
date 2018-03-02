

## set seed to randomly draw one 
## replicate
set.seed(06141971)
random_replicate <- sample(1:1000,1, replace = FALSE)

##load packages
library(dplyr)

##load data
siblingsd_df_noconfound <- readRDS("simulated_1000rep_sibwide_rho0.RDS")
siblingsd_df_confound <- readRDS("simulated_1000rep_sibwide_rho0.1.RDS")

print("data loaded")

## save one replicate for Figure 1
alldf_onesim_noconfound <- siblingsd_df_noconfound[random_replicate, ]
saveRDS(alldf_onesim_noconfound, "alldf_onesim_noconfound.RDS") 
alldf_onesim_confound <- siblingsd_df_confound[random_replicate, ]
saveRDS(alldf_onesim_confound, "alldf_onesim_confound.RDS") 




##create base of formulas
reg_predictors <- paste(c("snp1_minorcount",
                          "sex_offspring1", "age_offspring1",
                          "sex_offspring2", "age_offspring2"), 
                        collapse = "+")
reg_predictors_parents <- paste(c("snp1_minorcount",
                                  "snp1_minorcount_parents",
                                  "sex_offspring1", "age_offspring1",
                                  "sex_offspring2", "age_offspring2"),
                                collapse = "+") 


##### regressions: no confounding data

no_effect_noconfound <- lapply(siblingsd_df_noconfound,
                               function(x){
                                 lm(paste(paste("normal_neithereffect_eithersnp_confounded_siblingsd",
                                                reg_predictors, sep = "~"),
                                          "normal_neithereffect_eithersnp_confounded_siblingmean", 
                                          sep = "+"),
                                    data = x)
                               })

mean_effect_noconfound <- lapply(siblingsd_df_noconfound,
                                 function(x){
                                   lm(paste(paste("normal_meaneffect_snp1_confounded_siblingsd",
                                                  reg_predictors, sep = "~"),
                                            "normal_meaneffect_snp1_confounded_siblingmean",
                                            sep = "+"),
                                      data = x) 
                                 }) 

var_effect_noconfound <- lapply(siblingsd_df_noconfound,
                                function(x) {
                                  lm(paste(paste("normal_vareffect_snp1_confounded_siblingsd",
                                                 reg_predictors, sep = "~"),
                                           "normal_vareffect_snp1_confounded_siblingmean",
                                           sep  = "+"),
                                     data = x)
                                })

meanandvar_effect_noconfound <- lapply(siblingsd_df_noconfound,
                                       function(x) {
                                         lm(paste(paste("normal_meanvareffect_snp1_confounded_siblingsd",
                                                        reg_predictors, sep = "~"),
                                                  "normal_meanvareffect_snp1_confounded_siblingmean", 
                                                  sep = "+"),
                                            data = x)
                                         
                                       })


print("no parent control regressions completed")

##repeat regressions with control for parental genotype
##repeat w. control for parent genotype
no_effect_parents_noconfound <- lapply(siblingsd_df_noconfound,
                                       function(x) {
                                         lm(paste(paste("normal_neithereffect_eithersnp_confounded_siblingsd",
                                                        reg_predictors_parents, sep = "~"),
                                                  "normal_neithereffect_eithersnp_confounded_siblingmean", sep = "+"),
                                            data = x)
                                       })

##regression for outcome with mean effects- should not be significant
mean_effect_parents_noconfound <- lapply(siblingsd_df_noconfound,
                                         function(x) {
                                           lm(paste(paste("normal_meaneffect_snp1_confounded_siblingsd",
                                                          reg_predictors_parents, sep = "~"),
                                                    "normal_meaneffect_snp1_confounded_siblingmean",
                                                    sep = "+"),
                                              data = x) 
                                         })

##regression for outcome with var effects- should be significant
var_effect_parents_noconfound <- lapply(siblingsd_df_noconfound,
                                        function(x){
                                          lm(paste(paste("normal_vareffect_snp1_confounded_siblingsd",
                                                         reg_predictors_parents, sep = "~"),
                                                   "normal_vareffect_snp1_confounded_siblingmean",
                                                   sep  = "+"),
                                             data = x)
                                        }) 

##regression for outcome with var and mean effects- should be significant
meanandvar_effect_parents_noconfound <- lapply(siblingsd_df_noconfound,
                                               function(x) {
                                                 lm(paste(paste("normal_meanvareffect_snp1_confounded_siblingsd",
                                                                reg_predictors_parents, sep = "~"),
                                                          "normal_meanvareffect_snp1_confounded_siblingmean", 
                                                          sep = "+"),
                                                    data = x)
                                               })

##print status
print("parent control regressions completed")

##save full results for one randomly drawn replicate

all_results_noconfound <- list(no_effect_noconfound, mean_effect_noconfound, 
                               var_effect_noconfound,
                               meanandvar_effect_noconfound, no_effect_parents_noconfound,
                               mean_effect_parents_noconfound,
                               var_effect_parents_noconfound,
                               meanandvar_effect_parents_noconfound)


all_results_onesim_noconfound <- sapply(all_results_noconfound, "[", 
                                        random_replicate)
print("full results from one simulation generated")
saveRDS(all_results_onesim_noconfound, "all_results_onesim_cor0.RDS")

# save data associated with that replicate



##for other simulations, extract coefficient and p value on minor allele count
minorallele_coef_allsims_noconfound <- sapply(all_results_noconfound,
                                              function(x)
                                                lapply(x,
                                                       function(y)
                                                         summary(y)$coefficients["snp1_minorcount",
                                                                                 1]))
minorallele_p_allsims_noconfound <- sapply(all_results_noconfound,
                                           function(x)
                                             lapply(x,
                                                    function(y)
                                                      summary(y)$coefficients["snp1_minorcount",
                                                                              4]))
colnames(minorallele_coef_allsims_noconfound) <- paste(c("no_effect", "mean_effect", 
                                                         "var_effect",
                                                         "meanandvar_effect", 
                                                         "no_effect_parents",
                                                         "mean_effect_parents",
                                                         "var_effect_parents",
                                                         "meanandvar_effect_parents"),
                                                       "beta", sep = "_") 
colnames(minorallele_p_allsims_noconfound) <- paste(c("no_effect", "mean_effect", 
                                                      "var_effect",
                                                      "meanandvar_effect", 
                                                      "no_effect_parents",
                                                      "mean_effect_parents",
                                                      "var_effect_parents",
                                                      "meanandvar_effect_parents"),
                                                    "pval", sep = "_") 
minorallele_coef_df_noconfound <- cbind.data.frame(minorallele_coef_allsims_noconfound,
                                                   minorallele_p_allsims_noconfound) %>%
  mutate(sims = 1:length(no_effect_noconfound))

print("coef and p from all sim generated")

saveRDS(minorallele_coef_df_noconfound, "minorallele_coefandp_cor0.RDS")


### repeat process for ones with confounding
no_effect_confound <- lapply(siblingsd_df_confound,
                             function(x){
                               lm(paste(paste("normal_neithereffect_eithersnp_confounded_siblingsd",
                                              reg_predictors, sep = "~"),
                                        "normal_neithereffect_eithersnp_confounded_siblingmean", 
                                        sep = "+"),
                                  data = x)
                             })

mean_effect_confound <- lapply(siblingsd_df_confound,
                               function(x){
                                 lm(paste(paste("normal_meaneffect_snp1_confounded_siblingsd",
                                                reg_predictors, sep = "~"),
                                          "normal_meaneffect_snp1_confounded_siblingmean",
                                          sep = "+"),
                                    data = x) 
                               }) 

var_effect_confound <- lapply(siblingsd_df_confound,
                              function(x) {
                                lm(paste(paste("normal_vareffect_snp1_confounded_siblingsd",
                                               reg_predictors, sep = "~"),
                                         "normal_vareffect_snp1_confounded_siblingmean",
                                         sep  = "+"),
                                   data = x)
                              })

meanandvar_effect_confound <- lapply(siblingsd_df_confound,
                                     function(x) {
                                       lm(paste(paste("normal_meanvareffect_snp1_confounded_siblingsd",
                                                      reg_predictors, sep = "~"),
                                                "normal_meanvareffect_snp1_confounded_siblingmean", 
                                                sep = "+"),
                                          data = x)
                                       
                                     })


print("no parent control regressions completed")

##repeat regressions with control for parental genotype
##repeat w. control for parent genotype
no_effect_parents_confound <- lapply(siblingsd_df_confound,
                                     function(x) {
                                       lm(paste(paste("normal_neithereffect_eithersnp_confounded_siblingsd",
                                                      reg_predictors_parents, sep = "~"),
                                                "normal_neithereffect_eithersnp_confounded_siblingmean", sep = "+"),
                                          data = x)
                                     })

##regression for outcome with mean effects- should not be significant
mean_effect_parents_confound <- lapply(siblingsd_df_confound,
                                       function(x) {
                                         lm(paste(paste("normal_meaneffect_snp1_confounded_siblingsd",
                                                        reg_predictors_parents, sep = "~"),
                                                  "normal_meaneffect_snp1_confounded_siblingmean",
                                                  sep = "+"),
                                            data = x) 
                                       })

##regression for outcome with var effects- should be significant
var_effect_parents_confound <- lapply(siblingsd_df_confound,
                                      function(x){
                                        lm(paste(paste("normal_vareffect_snp1_confounded_siblingsd",
                                                       reg_predictors_parents, sep = "~"),
                                                 "normal_vareffect_snp1_confounded_siblingmean",
                                                 sep  = "+"),
                                           data = x)
                                      }) 

##regression for outcome with var and mean effects- should be significant
meanandvar_effect_parents_confound <- lapply(siblingsd_df_confound,
                                             function(x) {
                                               lm(paste(paste("normal_meanvareffect_snp1_confounded_siblingsd",
                                                              reg_predictors_parents, sep = "~"),
                                                        "normal_meanvareffect_snp1_confounded_siblingmean", 
                                                        sep = "+"),
                                                  data = x)
                                             })

##print status
print("parent control regressions completed")

##save full results for one randomly drawn replicate
## (drawn randomly above)
all_results_confound <- list(no_effect_confound, mean_effect_confound, 
                             var_effect_confound,
                             meanandvar_effect_confound, no_effect_parents_confound,
                             mean_effect_parents_confound,
                             var_effect_parents_confound,
                             meanandvar_effect_parents_confound)

all_results_onesim_confound <- sapply(all_results_confound, "[", random_replicate)
print("full results from one simulation generated")
saveRDS(all_results_onesim_confound, "all_results_onesim_cor0.1.RDS")






##for other simulations, extract coefficient and p value on minor allele count
minorallele_coef_allsims_confound <- sapply(all_results_confound,
                                            function(x)
                                              lapply(x,
                                                     function(y)
                                                       summary(y)$coefficients["snp1_minorcount",
                                                                               1]))
minorallele_p_allsims_confound <- sapply(all_results_confound,
                                         function(x)
                                           lapply(x,
                                                  function(y)
                                                    summary(y)$coefficients["snp1_minorcount",
                                                                            4]))
colnames(minorallele_coef_allsims_confound) <- paste(c("no_effect", "mean_effect", 
                                                       "var_effect",
                                                       "meanandvar_effect", 
                                                       "no_effect_parents",
                                                       "mean_effect_parents",
                                                       "var_effect_parents",
                                                       "meanandvar_effect_parents"),
                                                     "beta", sep = "_") 
colnames(minorallele_p_allsims_confound) <- paste(c("no_effect", "mean_effect", 
                                                    "var_effect",
                                                    "meanandvar_effect", 
                                                    "no_effect_parents",
                                                    "mean_effect_parents",
                                                    "var_effect_parents",
                                                    "meanandvar_effect_parents"),
                                                  "pval", sep = "_") 
minorallele_coef_df_confound <- cbind.data.frame(minorallele_coef_allsims_confound,
                                                 minorallele_p_allsims_confound) %>%
  mutate(sims = 1:length(no_effect_confound))

print("coef and p from all sim generated")

saveRDS(minorallele_coef_df_confound, "minorallele_coefandp_cor0.1.RDS")





