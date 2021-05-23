## Ellen Brandell
## May 2021

## assess cortisol results, make plots

library(tidyverse)
library(cowplot)
library(viridis)
theme_set(theme_minimal(base_size=11))

setwd("~/Desktop/scat_project/DATA/data_summaries")

# individual scat data
cort <- read.csv("all_cortisol.csv")
diameter <- read.csv("scat_diameter.csv")
geno <- read.csv("all_geno_sex_results.csv")

# environmental data
meta <- read.csv("scat_metadata.csv")
location <- read.csv("locations.csv")

######################################## COMBINE DATA
data1 <- merge(geno, cort, by="scatID", all.x=T)
data1 <- merge(data1, diameter, by="scatID", all.x=T)

###### make a data set of JUST KNOWN WOLVES
wolves1 <- data1[!is.na(data1$wolfID),]

######################################## EXPLORE CORTISOL RESULTS

## ALL WOLVES
hist(data1$cortisol)
summary(data1$cortisol)

## INDIVIDUAL-LEVEL
hist(wolves1$cortisol)
summary(wolves1$cortisol)

# look at repeated samples by individual



## PACK-LEVEL



## BY AGE


## BY SEX



ggplot(data=data1, aes(x=cortisol)) + geom_histogram(fill='gray', color='black')
ggplot(data=data1, aes(y=cortisol, color=sex)) + geom_boxplot(size=1) + facet_wrap(~sex)






