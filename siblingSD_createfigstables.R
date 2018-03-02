


## load packages
library(ggplot2)
library(dplyr)
library(xtable)
library(reshape2)
library("pwr")

## load custom theme
theme_new <- function(base_size = 16, base_family = "Helvetica"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank(),   
      panel.border = element_rect(fill = NA, colour = "black", size=1),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = NA),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black")
    )
}

## In-text Figure 1

### load in data from 
## one randomly selected simulation
siblingsd_rawdata_onesim <- readRDS("alldf_onesim_confound.RDS")

### subset to minor allele count 

meanandci_simulatedDV <- siblingsd_rawdata_onesim %>%
  group_by(factor(snp1_minorcount)) %>%
  summarise(meansd_noeffects = mean(normal_neithereffect_eithersnp_siblingsd),
            meansd_varianceoutcome = mean(normal_vareffect_snp1_siblingsd),
            meansd_meanoutcome = mean(normal_meaneffect_snp1_siblingsd),
            meansd_varmeanoutcome = mean(normal_meanvareffect_snp1_siblingsd),
            se_noeffects = sd(normal_neithereffect_eithersnp_siblingsd)/
              sqrt(nrow(siblingsd_rawdata_onesim)),
            se_varianceoutcome = sd(normal_vareffect_snp1_siblingsd)/
              sqrt(nrow(siblingsd_rawdata_onesim)),
            se_meanoutcome = sd(normal_meaneffect_snp1_siblingsd)/
              sqrt(nrow(siblingsd_rawdata_onesim)),
            se_varmeanoutcome = sd(normal_meanvareffect_snp1_siblingsd)/
              sqrt(nrow(siblingsd_rawdata_onesim)))
colnames(meanandci_simulatedDV)[1] <- "snp1_minorcount"
meanandci_simulatedDV_forplot <- meanandci_simulatedDV %>%
  reshape2::melt(id.vars = "snp1_minorcount") %>%
  mutate(parameter = gsub("\\_.*", "", 
                          variable),
         simulatedDV = gsub(".*\\_",
                            "",
                            variable)) %>%
  dplyr::select(-variable) %>%
  reshape(, idvar = c("snp1_minorcount",
                      "simulatedDV"),
          timevar = "parameter",
          direction = "wide",
          sep = "_") %>%
  mutate(lowerci = value_meansd - 1.96*value_se,
         upperci = value_meansd + 1.96 *value_se)
ggplot(meanandci_simulatedDV_forplot, 
       aes(x = snp1_minorcount,
           y = value_meansd,
           color = factor(simulatedDV,
                          levels = c("noeffects", "meanoutcome",
                                     "varianceoutcome", "varmeanoutcome"),
                          labels = c("SNP has no effect",
                                     "SNP affects mean only",
                                     "SNP affects variance only",
                                     "SNP affects mean and variance")),
           group = factor(simulatedDV))) +
  geom_point(position = "dodge") +
  geom_line() +
  geom_errorbar(aes(ymin = lowerci, ymax = upperci, width = 0.1)) +
  xlab("Sibling count of minor alleles") +
  ylab("Sibling standard \n deviation of outcome \n (non-adjusted mean and 95% CI)") +
  labs(color = "Type of simulated outcome") +
  theme_bw(base_size = 0) +
  ylim(0.5, 1.5) +
  theme_new() +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_blank()) +
  scale_color_manual(values = c("firebrick4", 
                                 "firebrick1",
                                 "springgreen4",
                                 "chartreuse2"))

ggsave("../finalfigures/Fig1.eps",
       plot = last_plot(),
       device = "eps",
       dpi = 400)

ggsave("../finalfigures/Fig1.pdf",
       plot = last_plot(),
       device = "pdf")

## In-text Figure 2

### load results for the sibling
## SD method for no correlation
## and correlation = 0.1
siblingsd_minorallele_cor0 <- readRDS("minorallele_coefandp_cor0.RDS") %>%
  mutate(correlation = 0)
siblingsd_minorallele_cor01 <- readRDS("minorallele_coefandp_cor0.1.RDS") %>%
  mutate(correlation = 0.1) 


siblingsd_minorallele_allregs <- rbind.data.frame(siblingsd_minorallele_cor0,
                                                  siblingsd_minorallele_cor01)


### plot coefficients for mean and variance across all simulations
siblingsd_minorallele_4plot <- siblingsd_minorallele_allregs %>%
  dplyr::select(sims, correlation,
                mean_effect_beta,
                var_effect_beta) %>%
  mutate(mean_effect_beta = unlist(mean_effect_beta),
         var_effect_beta = unlist(var_effect_beta)) %>%
  reshape2::melt(id.vars = c("sims", "correlation"))


ggplot(siblingsd_minorallele_4plot, aes(x = value,
                                        fill = factor(variable,
                                                      levels = c("mean_effect_beta", 
                                                                 "var_effect_beta"),
                                                      labels = c("Mean \n effects only",
                                                                 "Variance \n effects only")))) +
  geom_density() +
  facet_wrap(~ factor(correlation,
                      levels = c(0, 0.1),
                      labels = c("No confounding",
                                 "Some confounding (rho = 0.1)")),
             ncol = 1) +
  theme_bw(base_size = 20) +
  theme_new() +
  theme(legend.position = c(0.65, 0.9),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  labs(fill = "Simulated DV") +
  xlab("Coefficient on sibling count of minor alleles") +
  scale_fill_manual(values = c("firebrick",
                               "seagreen1"))

ggsave("../finalfigures/Fig2.eps",
       plot = last_plot(),
       device = "eps",
       dpi = 400)

ggsave("../finalfigures/Fig2.pdf",
       plot = last_plot(),
       device = "pdf")




## S1 fig: results for mean effects regressions of estimate
## on minor allele regression on minor allele count
## have ones but with and without
## control for popstrat
## but since previous results
## showed that popstrat control doesn't make
## much diff, just present ones 
## with popstrat control

## load data 

### pooled regression results
pool_corhigh_control <- readRDS("pool_lm_mean_simres_stratcontrol_cor01.RDS") %>%
  mutate(method = "pooled",
         correlation = 0.1,
         controlsforpop = "yes",
         sims = 1:1000)
pool_cormedhigh_control <- readRDS("pool_lm_mean_simres_stratcontrol_cor005.RDS") %>%
  mutate(method = "pooled",
         correlation = 0.05,
         controlsforpop = "yes",
         sims = 1:1000)
pool_cormedlow_control <- readRDS("pool_lm_mean_simres_stratcontrol_cor001.RDS") %>%
  mutate(method = "pooled",
         correlation = 0.01,
         controlsforpop = "yes",
         sims = 1:1000)
pool_cornone_control <- readRDS("pool_lm_mean_simres_stratcontrol_cor0.RDS") %>%
  mutate(method = "pooled",
         correlation = 0.0,
         controlsforpop = "yes",
         sims = 1:1000)

### random effects regression results
re_corhigh_control <- readRDS("re_lm_mean_simres_stratcontrol_cor01.RDS") %>%
  mutate(method = "RE",
         correlation = 0.1,
         controlsforpop = "yes",
         sims = 1:1000)
re_cormedhigh_control <- readRDS("re_lm_mean_simres_stratcontrol_cor005.RDS") %>%
  mutate(method = "RE",
         correlation = 0.05,
         controlsforpop = "yes",
         sims = 1:1000)
re_cormedlow_control <- readRDS("re_lm_mean_simres_stratcontrol_cor001.RDS") %>%
  mutate(method = "RE",
         correlation = 0.01,
         controlsforpop = "yes",
         sims = 1:1000)
re_cornone_control <- readRDS("re_lm_mean_simres_stratcontrol_cor0.RDS") %>%
  mutate(method = "RE",
         correlation = 0.0,
         controlsforpop = "yes",
         sims = 1:1000)


### fixed effects
fe_corhigh_control <- readRDS("fe_lm_mean_simres_cor01.RDS") %>%
  mutate(method = "FE",
         correlation = 0.1,
         controlsforpop = "yes",
         sims = 1:1000)
fe_cormedhigh_control <- readRDS("fe_lm_mean_simres_cor005.RDS") %>%
  mutate(method = "FE",
         correlation = 0.05,
         controlsforpop = "yes",
         sims = 1:1000)
fe_cormedlow_control <- readRDS("fe_lm_mean_simres_cor001.RDS") %>%
  mutate(method = "FE",
         correlation = 0.01,
         controlsforpop = "yes",
         sims = 1:1000)
fe_cornone_control <- readRDS("fe_lm_mean_simres_cor0.RDS") %>%
  mutate(method = "FE",
         correlation = 0.0,
         controlsforpop = "yes",
         sims = 1:1000)

## bind into one data.frame
reg_results_obs <- list(fe_corhigh_control,
                        fe_cormedhigh_control,
                        fe_cormedlow_control,
                        fe_cornone_control,
                        re_corhigh_control,
                        re_cormedhigh_control,
                        re_cormedlow_control,
                        re_cornone_control,
                        pool_corhigh_control,
                        pool_cormedhigh_control,
                        pool_cormedlow_control,
                        pool_cornone_control)

##bind into one dataframe
meanregs_allmethods <- do.call("rbind.data.frame", reg_results_obs)


##regressions with no controls for po
meanregs_allmethods_popstrat <- meanregs_allmethods %>%
  filter(controlsforpop == "yes") %>%
  mutate(correlation_factor = 
           factor(correlation,
       levels = c(0, 0.01, 0.05, 0.1),
     labels = c("Correlation between family intercept and genotype: 0",
                "Correlation between family intercept and genotype: 0.01",
                "Correlation between family intercept and genotype: 0.05",
                "Correlation between family intercept and genotype: 0.1"))) 


##plot results
ggplot(meanregs_allmethods_popstrat,
       aes(x = beta, 
           fill = factor(method,
                         levels = c("FE", "pooled", "RE"),
                         labels = c("FE", "Pooled", "RE")))) +
  geom_density(alpha = 0.6) +
  facet_wrap(~correlation_factor, ncol = 1) +
  xlab("Beta on individual minor allele count") +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "red") +
  theme_bw(base_size = 16) +
  labs(fill = "Method") +
  theme(legend.position = c(0.75, 0.65),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  scale_fill_brewer(palette = "Dark2")


ggsave("meanregs_allmethods_densitycompare.pdf",
       plot = last_plot(),
       device = "pdf")


## S4 table: comparing results of Hausman test 
## for fe, re, and pooled models
## compare models without control for population
## stratification because otherwise not 
## nested

### load hausman test data
hausman_cor0 <- as.data.frame(t(readRDS("haus_refe_nopopstratcontrol_cor0.RDS"))) %>%
  dplyr::select(statistic, p.value) %>%
  mutate(statistic = unlist(statistic),
         p.value= unlist(p.value),
         correlation = 0) 

hausman_cor0.01 <- as.data.frame(t(readRDS("haus_refe_nopopstratcontrol_cor001.RDS"))) %>%
  dplyr::select(statistic, p.value) %>%
  mutate(statistic = unlist(statistic),
         p.value= unlist(p.value),
         correlation = 0.01) 

hausman_cor0.05 <- as.data.frame(t(readRDS("haus_refe_nopopstratcontrol_cor005.RDS"))) %>%
  dplyr::select(statistic, p.value) %>%
  mutate(statistic = unlist(statistic),
         p.value= unlist(p.value),
         correlation = 0.05) 

hausman_cor0.1 <- as.data.frame(t(readRDS("haus_refe_nopopstratcontrol_cor01.RDS"))) %>%
  dplyr::select(statistic, p.value) %>%
  mutate(statistic = unlist(statistic),
         p.value= unlist(p.value),
         correlation = 0.1) 

### bind into single df
hausman_all <- rbind.data.frame(hausman_cor0,
                                hausman_cor0.01,
                                hausman_cor0.05,
                                hausman_cor0.1)


### summarize percent p value < 0.05
hausman_rejectrate <- hausman_all %>%
  group_by(factor(correlation)) %>%
  summarise(percent_reject = sum(p.value < 0.05)/1000)

### print table
xtable(hausman_rejectrate, digits = c(0, 2, 4))


## S5 and S6 Table: results of squared Z-score in 
## non-demeaned and demeaned data

### read in results
squaredz_results_nondemean <- readRDS("squaredZ_allregs_non_demean_updated20180212.RDS")
squaredz_results_demean <- readRDS("squaredZ_allregs_demean_updated20180214_newdemean.RDS")
colnames(squaredz_results_nondemean)[6:7] <- c("popstratcontrol", "sims")


## S5 Table. print results for non-demeaned-- percent of sims with p < 0.05
## on minor allele count coefficient
squaredz_results_fortable_nondemean <- squaredz_results_nondemean %>%
  dplyr::select(-sims) %>%
  group_by(outcome, popstratcontrol, confounding) %>%
  summarise(percent_sig = (sum(p < 0.05)/n())* 100) %>%
  arrange(confounding) %>%
  dplyr::select(-confounding)

print(xtable(squaredz_results_fortable_nondemean, 
              digits = c(4, 4, 4, 4)),
     include.rownames = FALSE)

## S6 Table. print results for demeaned

squaredz_results_fortable_demean <- squaredz_results_demean %>%
  dplyr::select(-sims) %>%
  group_by(outcome, confound) %>%
  summarise(percent_sig = (sum(p < 0.05)/n())* 100) %>%
  arrange(confound)

print(xtable(squaredz_results_fortable_demean, 
             digits = c(4, 4, 4, 4)),
      include.rownames = FALSE)


## S7 table: results of DGLM on non-demeaned and demeaned data

### first load in the results for non-demeaned data: look
### at three outcomes: null DV (both coef are false positive),
### mean DV (mean coef is true positive; var coef is false positive)
### var DV (mean coef is false positive; var coef is true positive)


dglm_nondemean <- readRDS("var_reg_results_nondemean_20180214_updated.RDS")

dglm_results_fortable_nondemean <- dglm_nondemean %>%
  dplyr::select(-sims, -data) %>%
  group_by(outcome, confounding, coeftype) %>%
  summarise(percent_sig = (sum(p < 0.05)/n())* 100) %>%
  arrange(outcome, coeftype)


print(xtable(dglm_results_fortable_nondemean, 
             digits = c(4, 4, 4, 4, 4)),
      include.rownames = FALSE)


### read in demean results and do the same
dglm_demean <- readRDS("var_reg_results_demean_20180215_updated.RDS")

dglm_results_fortable_demean <- dglm_demean %>%
  dplyr::select(-sims, -data) %>%
  group_by(outcome, confounding, coeftype) %>%
  summarise(percent_sig = (sum(p < 0.05)/n())* 100) %>%
  arrange(outcome, confounding, coeftype)


print(xtable(dglm_results_fortable_demean, 
             digits = c(4, 4, 4, 4, 4)),
      include.rownames = FALSE)


## S10 Table. Results across
## replicates for sibling SD method

### get p values from regressions
siblingsd_pval <- grep("pval", colnames(siblingsd_minorallele_allregs), 
                       value = TRUE)
siblingsd_coef2summarize <- siblingsd_minorallele_allregs[, 
                    c(siblingsd_pval, "sims", "correlation")] %>%
  mutate_if(is.list, as.numeric) %>%
  melt(, id.vars =  c("sims", "correlation"))  

## group by variable and summarize percent of p under 0.05
siblingsd_error_rate <- siblingsd_coef2summarize %>%
  mutate(factor_correlation = factor(correlation)) %>%
  group_by(variable, factor_correlation) %>%
  summarise(posrate = sum(value < 0.05)/1000 * 100) %>%
  mutate(parent_control = ifelse(grepl("parents",
                          variable), "Yes", "No")) %>%
  arrange(factor_correlation, variable) %>%
  dplyr::select(variable, parent_control, posrate, factor_correlation)



print(xtable(siblingsd_error_rate, digits = c(2, 2, 2, 2, 2)), 
      include.rownames = FALSE)

## S11 Table: correlation between 
## pop indicators and family intercepts

##load datasets with 4 correlation levels
## load these 
df_rho0 <- readRDS("simulated_1000rep_interceptcor0.RDS") 
df_rho001 <- readRDS("simulated_1000rep_interceptcor0.01.RDS")
df_rho005 <- readRDS("simulated_1000rep_interceptcor0.05.RDS")
df_rho01 <- readRDS("simulated_1000rep_interceptcor0.1.RDS")



df_rho0_selectvars <- lapply(df_rho0,
                              function(x)
                                x <- x %>%
                                group_by(factor(FID)) %>%
                                mutate(snp1_minorcount = 
                                         sum(2 * snp1_min_h +
                                               snp1_hz)) %>%
                                ungroup() %>%
                                dplyr::select(FID,
                                              IID,
                                              subpop,
                                              family_intercept,
                                              snp1_minorcount))

df_rho001_selectvars <- lapply(df_rho001,
                              function(x)
                                x <- x %>%
                                group_by(factor(FID)) %>%
                                mutate(snp1_minorcount = sum(2 * snp1_min_h +
                                                               snp1_hz)) %>%
                                ungroup() %>%
                                dplyr::select(FID,
                                              IID,
                                              subpop,
                                              family_intercept,
                                              snp1_minorcount))

df_rho005_selectvars <- lapply(df_rho05,
                               function(x)
                                 x <- x %>%
                                 group_by(factor(FID)) %>%
                                 mutate(snp1_minorcount = sum(2 * snp1_min_h +
                                                                snp1_hz)) %>%
                                 ungroup() %>%
                                 dplyr::select(FID,
                                               IID,
                                               subpop,
                                               family_intercept,
                                               snp1_minorcount))

df_rho01_selectvars <- lapply(df_rho01,
                               function(x)
                                 x <- x %>%
                                 group_by(factor(FID)) %>%
                                 mutate(snp1_minorcount = sum(2 * snp1_min_h +
                                                                snp1_hz)) %>%
                                 ungroup() %>%
                                 dplyr::select(FID,
                                               IID,
                                               subpop,
                                               family_intercept,
                                               snp1_minorcount))




## regress family intercept on subpop indicator
intercept_pop_rho0 <- unlist(lapply(df_rho0_selectvars,
                                     function(x)
                                       summary(lm(family_intercept ~ subpop,
                                                  data = x))$coefficients[2, 4])) 
intercept_pop_rho001 <- unlist(lapply(df_rho001_selectvars,
                                      function(x)
                                        summary(lm(family_intercept ~ subpop,
                                                   data = x))$coefficients[2, 4])) 

intercept_pop_rho005 <- unlist(lapply(df_rho005_selectvars,
                                      function(x)
                                        summary(lm(family_intercept ~ subpop,
                                                   data = x))$coefficients[2, 4])) 


intercept_pop_rho01 <- unlist(lapply(df_rho01_selectvars,
                                      function(x)
                                        summary(lm(family_intercept ~ subpop,
                                                   data = x))$coefficients[2, 4])) 

##bind results into a data.frame
corr_intercept_subpopindicator <- data.frame(correlation_level = 
                                               c(0, 0.01, 0.05, 0.1),
                                             percentsig = c((sum(intercept_pop_rho0 < 0.05)/1000),
                                                            (sum(intercept_pop_rho001 < 0.05)/1000),
                                                            (sum(intercept_pop_rho005 < 0.05)/1000),
                                                            (sum(intercept_pop_rho01 < 0.05)/1000))) 

## print 
xtable(corr_intercept_subpopindicator,
       digits = c(4, 4, 4))


## S12 table: relationship between intercept
## and observed genotype 


## change this code
allele_pop_rho0 <- unlist(lapply(df_rho0_selectvars,
                                  function(x)
                                    summary(lm(snp1_minorcount ~ family_intercept,
                                               data = x))$coefficients[2, 1])) 

allele_pop_rho001 <- unlist(lapply(df_rho001_selectvars,
                                   function(x)
                                     summary(lm(snp1_minorcount ~ family_intercept,
                                                data = x))$coefficients[2, 1]))  

allele_pop_rho005 <- unlist(lapply(df_rho005_selectvars,
                                 function(x)
                                   summary(lm(snp1_minorcount ~ family_intercept,
                                              data = x))$coefficients[2, 1]))

allele_pop_rho01 <- unlist(lapply(df_rho01_selectvars,
                                   function(x)
                                     summary(lm(snp1_minorcount ~ family_intercept,
                                                data = x))$coefficients[2, 1]))




corr_minorcount_familyint <- data.frame(correlation_level = 
                                          c(0, 0.01, 0.05, 0.1),
                                        percentsig = c(mean(allele_pop_rho0),
                                                       mean(allele_pop_rho001),
                                                       mean(allele_pop_rho005),
                                                       mean(allele_pop_rho01))) 

xtable(corr_minorcount_familyint,
       digits = c(4, 4, 4))


## S4 Fig: power analysis



## first replicate power analysis run in paper

## iterated through r^2 of putative effect up to point 1
## for ukb, only looked at full siblings and got 
## N from folliwng link: https://www.biorxiv.org/content/biorxiv/early/2017/07/20/166298.full.pdf
## for addhealth, got from parenting memo in biosoc 
## doc
r2_pwranalysis <- seq(0.00000001, 0.01, by = 0.0001)
n_pairlevel_FHS <- 583
n_pairlevel_addhealth <- 852
n_pairlevel_ukb <- 22666 
n_total <- n_pairlevel * 2
stringent_sig <- 0.00001
looser_sig <- 0.05 


## turn R2 into f2
f2_pwranalysis <- r2_pwranalysis/(1 - r2_pwranalysis)

## apply to each combination of sample sizes
pwr_FHS_stringent <- lapply(f2_pwranalysis,
                            pwr.f2.test,
                            u = 2, v = n_pairlevel_FHS - 2 -1, 
                            sig.level = stringent_sig,
                            power = NULL)
pwr_addH_stringent <- lapply(f2_pwranalysis,
                             pwr.f2.test,
                             u = 2, v = n_pairlevel_addhealth - 2 -1, 
                             sig.level = stringent_sig,
                             power = NULL)
pwr_ukb_stringent <- lapply(f2_pwranalysis,
                            pwr.f2.test,
                            u = 2, v = n_pairlevel_ukb - 2 -1, 
                            sig.level = stringent_sig,
                            power = NULL)


## get pwr from all
pwr_FHS_stringent_df <- data.frame(r2 = r2_pwranalysis,
                                   power = unlist(lapply(pwr_FHS_stringent,
                                                         function(x) x$power)),
                                   sig = stringent_sig,
                                   siblingpairs = n_pairlevel_FHS)
pwr_addh_stringent_df <- data.frame(r2 = r2_pwranalysis,
                                    power = unlist(lapply(pwr_addH_stringent,
                                                          function(x) x$power)),
                                    sig = stringent_sig,
                                    siblingpairs = n_pairlevel_addhealth)
pwr_ukb_stringent_df <- data.frame(r2 = r2_pwranalysis,
                                   power = unlist(lapply(pwr_ukb_stringent,
                                                         function(x) x$power)),
                                   sig = stringent_sig,
                                   siblingpairs = n_pairlevel_ukb)


## repeat with less stringent p-value threshold
pwr_FHS_replication <- lapply(f2_pwranalysis,
                              pwr.f2.test,
                              u = 2, v = n_pairlevel_FHS - 2 -1, 
                              sig.level = looser_sig,
                              power = NULL)
pwr_addH_replication <- lapply(f2_pwranalysis,
                               pwr.f2.test,
                               u = 2, v = n_pairlevel_addhealth - 2 -1, 
                               sig.level = looser_sig,
                               power = NULL)
pwr_ukb_replication <- lapply(f2_pwranalysis,
                              pwr.f2.test,
                              u = 2, v = n_pairlevel_ukb - 2 -1, 
                              sig.level = looser_sig,
                              power = NULL)


## create df
pwr_FHS_replication_df <- data.frame(r2 = r2_pwranalysis,
                                     power = unlist(lapply(pwr_FHS_replication,
                                                           function(x) x$power)),
                                     sig = looser_sig,
                                     siblingpairs = n_pairlevel_FHS)
pwr_addh_replication_df <- data.frame(r2 = r2_pwranalysis,
                                      power = unlist(lapply(pwr_addH_replication,
                                                            function(x) x$power)),
                                      sig = looser_sig,
                                      siblingpairs = n_pairlevel_addhealth)
pwr_ukb_replication_df <- data.frame(r2 = r2_pwranalysis,
                                     power = unlist(lapply(pwr_ukb_replication,
                                                           function(x) x$power)),
                                     sig = looser_sig,
                                     siblingpairs = n_pairlevel_ukb)

## combine into separate df's for now
pwr_alldf <- rbind.data.frame(pwr_FHS_stringent_df,
                              pwr_addh_stringent_df,
                              pwr_ukb_stringent_df,
                              pwr_FHS_replication_df,
                              pwr_addh_replication_df,
                              pwr_ukb_replication_df)

## plot to see similarity with table
ggplot(pwr_alldf, 
       aes(x = r2,
           y = power,
           color = factor(siblingpairs,
                          levels = c(583,
                                     852, 22666),
                          labels = c("FHS: 583",
                                     "AddHealth: 852",
                                     "UK Biobank: 22,666")),
           linetype = factor(sig,
                             levels = c(0.05,
                                        0.00001),
                             labels = c("p < 0.05",
                                        "p < 10^-5"))))  +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 1, 
                                  by = 0.2)) +
  xlab("R-squared of Putative Effect\n(depicted on log scale)") +
  ylab("Power") +
  labs(color = "Sample size\n (# of sibling pairs)",
       linetype = "Sig. threshold") +
  theme_bw(base_size = 16) +
  theme(legend.position = "right",
        legend.background = element_blank(),
        axis.text.x = element_text(color= "black"),
        axis.text.y = element_text(color = "black"))+
  geom_hline(yintercept = 0.8,
             linetype = "dashed",
             color = "red") +
  scale_x_log10() +
  scale_color_brewer(palette = "Dark2")

ggsave("finalfigures/updated_power_analysis.pdf",
       plot = last_plot(),
       device = "pdf")




