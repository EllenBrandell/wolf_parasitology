## Ellen Brandell
## September 2021

## SUCCESS MODEL

# load libraries for analysis
library(brms)      # Bayesian model
library(tidybayes) # for analysis of posterior draws of a Bayesian model
# load libraries for plotting
library(tidyverse)
# library(viridis)
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
dim(posterior) # iterations retained, chains, parameters
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
          prob=0.80,
          transformations="exp") +
  geom_vline(xintercept=1, color="grey") +
  xlim(0,50)

p1 <- mcmc_intervals(posterior, transformations="exp", 
                     prob=0.5, prob_outer=0.7, inner_size=1, point_size=3,
               pars=c("b_days_elapsed","b_cover1","b_high_temp_mean","b_total_precip"))
p2 <- mcmc_intervals(posterior, transformations="exp", 
                     prob=0.5, prob_outer=0.5, inner_size=1, point_size=3,
                     pars=c("b_periodwinter","b_periodz_den"),)

p1.1 <- p1 + scale_y_discrete(
  labels=c("exp(b_days_elapsed)"="DAYS",
           "exp(b_cover1)"="COVER (1)",
           "exp(b_high_temp_mean)"="HIGH TEMP",
           "exp(b_total_precip)"="PRECIP",
           "exp(b_periodwinter)"="WINTER",
           "exp(b_periodz_den)"="DEN" )  ) + 
  xlim(0,3) + xlab("log-odds") + geom_vline(aes(xintercept=1),color="gray60")
p1.1

p2.1 <- p2 + scale_y_discrete(
  labels=c("exp(b_periodwinter)"="WINTER",
           "exp(b_periodz_den)"="DEN" )  ) + 
  xlim(0,45) + xlab("log-odds") + geom_vline(aes(xintercept=1),color="gray60")
p2.1

plot_grid(p1.1,p2.1, nrow=1, labels=c("(A)","(B)"))
ggsave("SUCCESS_MODEL_coef_plots.jpeg", height=5, width=10, units="in", dpi=300)

