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
data1 <- merge(data1, location, by="scatID", all.x=T)
data1 <- merge(data1, meta[,c(1,4)], by="scatID", all.x=T) # I want site

## look at summary information
table(is.na(data1$wolfID))
# FALSE  TRUE 
#  49    50
49/110*100

###### AGE
ggplot(data=data1, aes(x=scat_diameter)) + geom_histogram(fill='gray', color='black')
data1$age <- ifelse(data1$scat_diameter<=2.0, "Juvenile", "Adult")
table(data1$age)
table(data1$age)/nrow(data1)

###### make a data set of JUST KNOWN WOLVES
wolves1 <- data1[!is.na(data1$wolfID),]

## summarize over wolves
wolves2 <- wolves1 %>% group_by(wolfID) %>% summarize(sex=sex[1])



###### AGE - unique wolves
table(wolves1$age)
table(wolves1$age)/nrow(wolves1)

## sex information
table(wolves2$sex, useNA="ifany")
table(wolves2$sex, useNA="ifany")/nrow(wolves2)

###### CORTISOL
summary(data1$cortisol)
summary(wolves1$cortisol)


######################################## CONFIRM CORRECT CLASS









