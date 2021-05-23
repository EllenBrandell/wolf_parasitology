## Ellen Brandell
## May 2021

## assess the importance of the diameter used to delineate pups/adults 

library(tidyverse)
library(cowplot)
theme_set(theme_minimal(base_size=11))

setwd("~/Desktop/scat_project/DATA/data_summaries")

# individual scat data
diameter <- read.csv("scat_diameter.csv")
geno <- read.csv("all_geno_sex_results.csv")

# environmental data
meta <- read.csv("scat_metadata.csv")
location <- read.csv("locations.csv")

######################################## COMBINE  DATA
data1 <- merge(geno, diameter, by="scatID", all.x=T)
data1 <- merge(data1, meta[,c(1,4)], by="scatID", all.x=T) # I want site

###### make a data set of JUST KNOWN WOLVES
wolves1 <- data1[!is.na(data1$wolfID),]


######################################## DIAMETER PATTERNS
ggplot(data=data1, aes(x=scat_diameter)) + geom_histogram(fill='gray', color='black')

## look at different diameter measurements by unique wolves



## look at different diameter measurements by collared wolves (known age)



######################################## EXPLORING CUTOFFS
cut <- 2.0  # less than OR equal to this is a pup
data1$age <- ifelse(data1$scat_diameter<=cut, "Pup", "Adult")
# make sure age isn't changing unless the wolf ages from Juv -> Adult


table(data1$age)
table(data1$age)/nrow(data1)
# "questionable"/"in between" scat diameters based on our cutoff
table(data1$scat_diameter>cut-0.3 & data1$scat_diameter<cut+0.3)
# there are quite a few of these "in between" sizes

###### HOW MANY MISMATCHES DO WE HAVE BASED ON CUTOFF VALUE?
# can only use resampled wolves for this


## and importantly, we want to make sure the collared wolves are correctly classified:
# pup = 1211F, 1228F, 1229F
# adult = 969F, 907F, 996M, 1005F, 1155M, 1156M


######
# change scat 145 wolf 1228F to Pup


