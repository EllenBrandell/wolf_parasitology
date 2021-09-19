## Ellen Brandell
## September 2021

## SUCCESS MODEL

# load libraries for analysis
library(brms)      # Bayesian model
library(tidybayes) # for analysis of posterior draws of a Bayesian model
# load libraries for plotting
library(tidyverse)
library(viridis)
library(cowplot)
library(bayesplot) # for plotting posterior draws of a Bayesian model

# load data:
data <- read.csv("success_model_data.csv")
data <- data %>% mutate_at(c("geno","site","cover","month_created","period"), factor)

str(data)
summary(data)
summary(is.na(data))

## find the variables for which I can set priors:
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




