## Ellen Brandell
## September 2021

## SUCCESS MODEL

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
data <- read.csv("success_model_data.csv")