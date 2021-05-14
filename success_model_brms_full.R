## Ellen Brandell
## May 2021





################################### FULL MODEL


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




sessionInfo()
