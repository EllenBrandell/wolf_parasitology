## Ellen Brandell
## May 2021

## GOAL: initial modeling for *genotype success model* in Bayesian framework

# load libraries for analysis
library(brms) # Bayesian model
library(sjstats) # calculate intra-class correlation (ICC)
library(tidybayes) # for analysis of posterior draws of a Bayesian model
library(reshape2)
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


###### MODEL COMPONENTS
## the model is a Bayesian hierarchical model; response = genotype success, Binomial

######################## START WITH A BASIC MODEL WITH DEFAULT PARAMETERS ########################
mod0 <- brm(formula= geno_success ~ cover + high_temp_mean + (1|group),
            data=data, family=bernoulli(link="logit"), sample_prior=T)
summary(mod0)
posterior_summary(mod0, probs=c(.025, .25, .75, .975))
plot(mod0)

## SOURCES:
# https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html 
# https://bookdown.org/ajkurz/DBDA_recoded/jags-brms.html


################################### ASSESS MODEL CONVERGENCE
mcmc_plot(mod0, type="neff_hist")    # effective sample size
# ratio: Neff/N -> closer to 1 is better (lighter); light: between 0.5 and 1 (high); mid: between 0.1 and 0.5 (good); dark: below 0.1 (low)

mcmc_plot(mod0, type="rhat")         # R hat (~1)
mcmc_plot(mod0, type="rhat_hist")    # R hat (~1)
# pull out specific rhat estimates of interest
rhat(mod0)
rhat(mod0)["b_Intercept"]
rhat(mod0)["b_cover1"]

mcmc_plot(mod0, type="trace")        # trace plots
mcmc_plot(mod0, type="acf_bar")      # look at autocorrelation
mcmc_plot(mod0, type="dens_overlay") # density plots for each chain

# plot some diagnostics specific to the NUTS sampler - WHAT DO THESE MEAN
mcmc_plot(mod0, type="nuts_acceptance")
mcmc_plot(mod0, type="nuts_divergence")

## alternative way to look at model convergence by drawing from posterior
post <- posterior_samples(mod0, add_chain=T)
mcmc_dens_overlay(post, pars=c("b_Intercept")) +
  theme(panel.grid=element_blank())
mcmc_acf(post, pars="b_Intercept", lags=20)


################################### PARAMETER ESTIMATES
posterior_summary(mod0, probs=c(.025, .25, .75, .975))
posterior <- as.array(mod0)
dim(posterior)
dimnames(posterior)

color_scheme_set("red")
# line plots
mcmc_intervals(posterior, pars=c("b_Intercept","b_cover1","b_high_temp_mean","sd_group__Intercept"))
# density plots
mcmc_areas(
  posterior, 
  pars=c("b_Intercept","b_cover1","b_high_temp_mean","sd_group__Intercept"),
  prob=0.5, # 50% intervals
  prob_outer=0.95, # 95%  intervals
  point_est="mean"
)
# histograms
mcmc_hist(posterior, pars=c("b_Intercept","b_cover1","b_high_temp_mean","sd_group__Intercept"))

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


################################### SAMPLE FROM PRIOR DISTRIBUTION
prior <- prior_samples(mod0)
head(prior)

# doesn't make a ton of sense with the defaults ?










######################## ADDING SPECIFICATIONS TO MODEL ########################
## see which parameters can have priors
get_prior(formula= geno_success ~ cover + high_temp_mean + (1|group),
    data=data, family=bernoulli(link="logit"))


## PRIORS
# it was a bit tricky to get the syntax for the random effect correct here
prior.list <- c(prior(normal(0,10), class=Intercept),
                prior(normal(0,10), class=b, coef=cover1),
                prior(normal(0,10), class=b, coef=high_temp_mean),
                prior_string("uniform(0, 20)", class="sd", coef="Intercept", group="group", ub=NA) )

## MCMC setting
nburn <- 5000
niter <- 10000
nchains <- 3

###### RUN MODEL
mod1 <- brm(formula= geno_success ~ cover + high_temp_mean + (1|group),
                data=data, family=bernoulli(link="logit"),
                prior=prior.list, sample_prior=T,
                warmup=nburn, iter=niter, chains=nchains)

summary(mod1)
head(posterior_summary(mod1, probs=c(.025, .25, .75, .975)))
plot(mod1)


################################### PARAMETER ESTIMATES
posterior1 <- as.array(mod1)
head(posterior1)
post1 <- posterior_samples(mod1, add_chain=T)
head(post1)
summary(post1$chain)

################################### SAMPLE FROM PRIOR DISTRIBUTION
prior1 <- prior_samples(mod1)
head(prior1)

prior1 %>% 
  rename(Intercept        = Intercept,
         Covered          = b_cover1,
         High_Temperature = b_high_temp_mean,
         Site_Error       = sd_group__Intercept)  %>% 
  gather() %>% 
  mutate(key=factor(key, levels=c("Intercept","Covered","High_Temperature","Site_Error"))) %>%
  ggplot(aes(x=value, group=key)) +
  geom_histogram(color="grey92", fill="grey67", size=.2) +
  stat_pointinterval(aes(y=0), point_interval=mode_hdi, .width=c(.95, .50)) +
  scale_y_continuous(NULL, breaks=NULL) +
  theme_grey() +
  theme(panel.grid=element_blank()) +
  facet_wrap(~key, scales="free_x") 
# the uniform intervals look odd

## plot the PRIOR and the POSTERIOR by variable to compare
cover.pp <- data.frame(prior=prior1$b_cover1, posterior=post1$b_cover1, chain=post1$chain)
cover.pp2 <- melt(cover.pp)

# plot each chain separately
ggplot(data=cover.pp2, aes(x=value, fill=variable, color=variable)) +
  geom_density(size=1, alpha=0.5) +
  facet_wrap(~chain) + # , labeller=c("Chain 1","Chain 2","Chain 3")
  theme(panel.grid=element_blank()) +
  theme_light() +
  labs(title="Covered",
       x    =expression(beta)) +
  xlim(-20,20)

# plot all chains and 1 prior
cover.pp2.post <- cover.pp2[cover.pp2$variable=='posterior',]
cover.pp2.prior <- cover.pp2[cover.pp2$variable=='prior' & cover.pp2$chain==1 ,]

ggplot(data=cover.pp2.post, aes(x=value, fill=chain, color=chain)) +
  geom_density(size=0.7, alpha=0.1) +
  theme(panel.grid=element_blank()) +
  theme_light() +
  geom_density(data=cover.pp2.prior, aes(x=value), color='black', fill=NA, size=1) +
  xlim(-8,8)
  


######################################################################
sessionInfo()
