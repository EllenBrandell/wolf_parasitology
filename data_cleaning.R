## Ellen Brandell
## May 2021

## combine and clean scat-related files into TWO final data set for analysis:
# 1: wolf characteristics
# 2: environmental data

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

######################################## COMBINE AND CLEAN DATA SET 1
data1 <- merge(geno, cort, by="scatID", all.x=T)
data1 <- merge(data1, diameter, by="scatID", all.x=T)
data1 <- merge(data1, meta[,c(1,4)], by="scatID", all.x=T) # I want site

## look at summary information
table(is.na(data1$wolfID))
# FALSE  TRUE 
#  60    50
60/110*100

###### AGE
ggplot(data=data1, aes(x=scat_diameter)) + geom_histogram(fill='gray', color='black')
data1$age <- ifelse(data1$scat_diameter<=2.0, "Juvenile", "Adult")
table(data1$age)
table(data1$age)/nrow(data1)
table(data1$scat_diameter>1.9 & data1$scat_diameter<2.2)
# there are quite a few of these "in between" sizes

###### make a data set of JUST KNOWN WOLVES
wolves1 <- data1[!is.na(data1$wolfID),]

###### AGE - unique wolves
# make sure age isn't changing unless the wolf ages from Juv -> Adult
# change scat 145 wolf 1228F to Juvenile

table(wolves1$age)
table(wolves1$age)/nrow(wolves1)
table(data1$scat_diameter>2 & data1$scat_diameter<=2.5)
# 34 !! A LOT of in-between values

## sex information
table(wolves1$sex, useNA="ifany")
table(wolves1$sex, useNA="ifany")/nrow(wolves1)

###### CORTISOL
hist(data1$cortisol)
hist(wolves1$cortisol)
summary(data1$cortisol)
summary(wolves1$cortisol)

ggplot(data=data1, aes(x=cortisol)) + geom_histogram(fill='gray', color='black')
ggplot(data=data1, aes(y=cortisol, color=sex)) + geom_boxplot(size=1) + facet_wrap(~sex)












