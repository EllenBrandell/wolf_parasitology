## Ellen Brandell
## May 2021


# GOAL: initial modeling for *genotype success model* in Bayesian framework

# load libraries for analysis
library(brms) # Bayesian model
library(sjstats) # calculate intra-class correlation (ICC)
library(tidybayes) # for analysis of posterior draws of a Bayesian model
# load libraries for plotting
library(tidyverse)
library(viridis)
library(cowplot)

# load data
data.all <- read.csv("success_data_all.csv")
# keep variables of interest
vars <- c("scat_id","geno_success","days_elapsed","cover","collection_period","high_temp_mean","total_precip","group")
data <- data.all[,colnames(data) %in% vars]
# make sure variables are correct class
data <- data %>% mutate_at(vars("scat_id","cover","collection_period","group"), as.factor)
            # cover 0:open, 1:covered.
data <- data %>% mutate_at(vars("geno_success","days_elapsed","high_temp_mean","total_precip"), as.numeric)
# make sure there is not missing data
summary(is.na(data))
str(data)
summary(data)


###### MODEL COMPONENTS
## the model is a Bayesian hierarchical model; response = genotype success, Binomial

################################### START WITH A BASIC MODEL WITH DEFAULT PARAMETERS
mod0 <- brm(formula= geno_success ~ cover + high_temp_mean + (1|group),
            data=data, family=bernoulli(link = "logit"))
summary(mod0)
plot(mod0)

# assess model convergence
mcmc_plot(mod0, type="neff_hist")    # effective sample size
# ratio: Neff/N -> closer to 1 is better (lighter); light: between 0.5 and 1 (high); mid: between 0.1 and 0.5 (good); dark: below 0.1 (low)
mcmc_plot(mod0, type="rhat")         # R hat (~1)
mcmc_plot(mod0, type="rhat_hist")    # R hat (~1)
mcmc_plot(mod0, type="trace")        # traceplots
mcmc_plot(mod0, type="acf_bar")      # look at autocorrelation
mcmc_plot(mod0, type="dens_overlay") # density plots for each chain

# plot some diagnostics specific to the NUTS sampler - WHAT DO THESE MEAN
mcmc_plot(mod0, type="nuts_acceptance")
mcmc_plot(mod0, type="nuts_divergence")












################################### ADDING SPECIFICATIONS TO MODEL

## PRIORS
prior.list <- c(prior_string("normal(0,10)", class="Intercept"),
                prior_string("normal(0,10)", class="b", coef="cover"),
                prior_string("normal(0,10)", class="b", coef="high_temp_mean"),
                prior_string("uniform(0,20)", class="sd", group="group", coef="group") )

## MCMC setting
nburn <- 5000
niter <- 10000
nchains <- 3

###### RUN MODEL
mod1 <- brm(formula= geno_success ~ cover + high_temp_mean + (1|group),
                data=data, family=bernoulli,
                prior=prior.list,
                warmup=nburn, iter=niter, chains=nchains)

















## PRIORS
prior.list <- c(set_prior("normal(0,10)", class="b", coef="cover"),
                set_prior("normal(0,10)", class="b", coef="high_temp_mean"),
                set_prior("lognormal(0,10)", class="b", coef="total_precip"),
                set_prior("normal(0,10)", class="b", coef="collection_period"),
                set_prior("normal(0,10)", class="b", coef="days_elapsed"),
                set_prior("normal(0,20)", class="b", coef="group"))

## MCMC setting
nburn <- 5000
niter <- 10000
nchains <- 3


###### RUN MODEL
full.mod <- brm(formula= geno_success ~ cover + high_temp_mean*days_elapsed + total_precip + (1|group),
            data=data, family=bernoulli,
            prior=prior.list,
            warmup=nburn, iter=niter, chains=nchains)





