## Ellen Brandell
## May 2021

## GOAL: initial modeling for FULL *genotype success model* in Bayesian framework
## * not final

# load libraries for analysis
library(brms) # Bayesian model
library(sjstats) # calculate intra-class correlation (ICC)
library(tidybayes) # for analysis of posterior draws of a Bayesian model
# load libraries for plotting
library(tidyverse)
library(viridis)
library(cowplot)
library(bayesplot) # for plotting posterior draws of a Bayesian model

# load data:
setwd("~/Desktop/scat_project/DATA/data_summaries")
# individual scat data
# geno <- read.csv("all_geno_sex_results.csv")
# environmental data
success <- read.csv("sent_successful.csv")
success <- success[,c(1,8)] # just need scatID and geno
meta <- read.csv("scat_metadata.csv")
location <- read.csv("locations.csv")
# combine
data <- merge(success, meta, all.x=T, by="scatID")
data <- merge(data, location, all.x=T, by="scatID")

# remove those without meta
data <- data[complete.cases(data),]

# load data from Maddy's csv
# data.all <- read.csv("success_data_all.csv")
# keep variables of interest
# vars <- c("scatID","geno","days_elapsed","cover","period","high_temp_mean","total_precip","site")
# data <- data.all[,colnames(data.all) %in% vars]
# make sure variables are correct class
data <- data %>% mutate_at(vars("scatID","cover","period","site"), as.factor)
# cover 0:open, 1:covered.
data <- data %>% mutate_at(vars("geno","days_elapsed","high_temp_mean","total_precip"), as.numeric)
# make sure there is not missing data
summary(is.na(data))
str(data)
summary(data)



################################### FULL MODEL - default priors and MCMC settings
mod0 <- brm(formula= geno ~ days_elapsed + cover + high_temp_mean + total_precip + period + (1|site),
            data=data, family=bernoulli(link="logit"))
summary(mod0)
posterior_summary(mod0, probs=c(.025, .25, .75, .975))
plot(mod0)

################################### ASSESS MODEL CONVERGENCE
mcmc_plot(mod0, type="neff_hist")    # effective sample size
# ratio: Neff/N -> closer to 1 is better (lighter); light: between 0.5 and 1 (high); mid: between 0.1 and 0.5 (good); dark: below 0.1 (low)

mcmc_plot(mod0, type="rhat")         # R hat (~1)
mcmc_plot(mod0, type="rhat_hist")    # R hat (~1)
# pull out specific rhat estimates of interest
rhat(mod0)

mcmc_plot(mod0, type="trace")        # trace plots
mcmc_plot(mod0, type="acf_bar")      # look at autocorrelation
mcmc_plot(mod0, type="dens_overlay") # density plots for each chain

################################### PARAMETER ESTIMATES
posterior_summary(mod0, probs=c(.025, .25, .75, .975))
posterior <- as.array(mod0)
dim(posterior)
dimnames(posterior)

color_scheme_set("red")
# line plots
mcmc_intervals(posterior, pars=c("b_Intercept","b_days_elapsed","b_cover1","b_high_temp_mean","b_total_precip","b_periodwinter","b_periodz_den","sd_site__Intercept"))
# density plots
mcmc_areas(
  posterior, 
  pars=c("b_Intercept","b_days_elapsed","b_cover1","b_high_temp_mean","b_total_precip","b_periodwinter","b_periodz_den","sd_site__Intercept"),
  prob=0.5, # 50% intervals
  prob_outer=0.95, # 95%  intervals
  point_est="mean"
)
# histograms
mcmc_hist(posterior, pars=c("b_Intercept","b_days_elapsed","b_cover1","b_high_temp_mean","b_total_precip","b_periodwinter","b_periodz_den","sd_site__Intercept"))

# manual plot of line estimate + density with shading
post <- posterior_samples(mod0, add_chain=T)
post %>% 
  ggplot(aes(x=b_Intercept, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Intercept",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_cover1, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Covered",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_high_temp_mean, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Mean high temperature (F)",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_days_elapsed, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Days elapsed",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_total_precip, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Total precipitation (inches)",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_periodwinter, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Winter collection",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_periodz_den, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Den collection",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()



################################### FULL MODEL - set priors and MCMC settings
get_prior(formula= geno ~ period + cover + high_temp_mean*days_elapsed + total_precip + (1|site),
          data=data, family=bernoulli(link="logit"))

## PRIORS
prior.list <- c(set_prior("normal(0,10)", class="b", coef="cover1"),
                set_prior("normal(0,10)", class="b", coef="high_temp_mean"),
                set_prior("normal(0,10)", class="b", coef="total_precip"),
                set_prior("normal(0,10)", class="b", coef="periodwinter"),
                set_prior("normal(0,10)", class="b", coef="periodz_den"),
                set_prior("normal(0,10)", class="b", coef="days_elapsed"),
                prior_string("uniform(0,20)", class="sd", coef="Intercept", group="site", ub=NA))

## MCMC setting
nburn <- 5000
niter <- 10000
nchains <- 3


###### RUN MODEL
full.mod <- brm(formula= geno ~ period + cover + high_temp_mean + days_elapsed + total_precip + (1|site),
                data=data, family=bernoulli,
                prior=prior.list,
                warmup=nburn, iter=niter, chains=nchains)
summary(full.mod)
posterior_summary(full.mod, probs=c(.025, .25, .75, .975))
plot(full.mod)


################################### ASSESS MODEL CONVERGENCE
mcmc_plot(full.mod, type="neff_hist")    # effective sample size
# ratio: Neff/N -> closer to 1 is better (lighter); light: between 0.5 and 1 (high); mid: between 0.1 and 0.5 (good); dark: below 0.1 (low)

mcmc_plot(full.mod, type="rhat")         # R hat (~1)
mcmc_plot(full.mod, type="rhat_hist")    # R hat (~1)
# pull out specific rhat estimates of interest
rhat(full.mod)

mcmc_plot(full.mod, type="trace")        # trace plots
mcmc_plot(full.mod, type="acf_bar")      # look at autocorrelation
mcmc_plot(full.mod, type="dens_overlay") # density plots for each chain

################################### PARAMETER ESTIMATES
posterior_summary(full.mod, probs=c(.025, .25, .75, .975))
posterior <- as.array(full.mod)
dim(posterior)
dimnames(posterior)

color_scheme_set("red")
# line plots
mcmc_intervals(posterior, pars=c("b_Intercept","b_days_elapsed","b_cover1","b_high_temp_mean","b_total_precip","b_periodwinter","b_periodz_den","sd_site__Intercept"))
# density plots
mcmc_areas(
  posterior, 
  pars=c("b_Intercept","b_days_elapsed","b_cover1","b_high_temp_mean","b_total_precip","b_periodwinter","b_periodz_den","sd_site__Intercept"),
  prob=0.5, # 50% intervals
  prob_outer=0.95, # 95%  intervals
  point_est="mean"
)
# histograms
mcmc_hist(posterior, pars=c("b_Intercept","b_days_elapsed","b_cover1","b_high_temp_mean","b_total_precip","b_periodwinter","b_periodz_den","sd_site__Intercept"))

# manual plot of line estimate + density with shading
post <- posterior_samples(full.mod, add_chain=T)
post %>% 
  ggplot(aes(x=b_Intercept, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Intercept",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_cover1, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Covered",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_high_temp_mean, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Mean high temperature (F)",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_days_elapsed, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Days elapsed",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_total_precip, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Total precipitation (inches)",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_periodwinter, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Winter collection",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()
post %>% 
  ggplot(aes(x=b_periodz_den, y=0)) +
  stat_halfeye(point_interval=mode_hdi, .width=c(.95, .5)) +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(title="Den collection",
       x    =expression(beta)) +
  theme(panel.grid=element_blank()) +
  vline_0()


############################################# EXPONENTIATING = ODDS
# https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
exp(fixef(full.mod))
mcmc_plot(full.mod, 
         type="areas",
         prob=0.95,
         transformations="exp") +
  geom_vline(xintercept=1, color="grey") +
  xlim(0,3)

##### MORE TO DO HERE





sessionInfo()
