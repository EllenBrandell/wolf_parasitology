## Ellen Brandell
## September 2021

## SUCCESS MODEL DATA - make the csv and export it!

library(tidyverse)

# load data:
setwd("~/Desktop/scat_project/DATA/data_summaries")
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

data <- data %>% mutate_at(vars("scatID","cover","period","site"), as.factor)
# cover 0:open, 1:covered.
data <- data %>% mutate_at(vars("geno","days_elapsed","high_temp_mean","total_precip"), as.numeric)
# make sure there is not missing data
summary(is.na(data))
str(data)
summary(data)

setwd("~/Desktop/scat_project/GitHub/wolf_parasitology")
write.csv(data, "success_model_data.csv", row.names=F)

# USING TOKEN LOL SO DUMB