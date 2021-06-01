## Ellen Brandell
## May 2021

## getting data together and doing some summary statistics and plotting for the upcoming June 11 coauthors meetng

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

# get all infection information
para <- read.csv("sent_successful.csv")
# clean this a bit
para1 <- para[,c(1,6)]
colnames(para1) <- c("scatID","result")
para1 <- para1[!is.na(para1$result),]
para1$result <- ifelse(para1$result=="negative",0,1)
table(para1$result)
table(para1$result)/nrow(para1)


######################################## COMBINE AND CLEAN DATA SET 1
data1 <- merge(cort, diameter, by="scatID", all.x=T)
data1 <- merge(data1, location, by="scatID", all.x=T)
data1 <- merge(data1, meta[,c(1,4)], by="scatID", all.x=T) # I want site
data1 <- merge(data1, para1, by="scatID", all.x=T) # I want site

# just known wolves
data2 <- merge(geno, cort, by="scatID", all.x=T)
table(data2$e.c.)
table(data2$e.m.)
table(data2$t.m.d.)

data2$e.c. <- ifelse(data2$e.c.=="P",1,0)
data2$e.m. <- ifelse(data2$e.m.=="P",1,0)
data2$t.m.d. <- ifelse(data2$t.m.d.=="P",1,0)

unique(data2$pack)
data2$pack <- gsub("1005F Group","Crevice Lake",data2$pack)
data2$pack <- gsub("Crevice Lake/Crevice Lake","Crevice Lake",data2$pack)
data2$pack <- droplevels(as.factor(data2$pack))
table(data2$pack)

data2$inf <- data2$e.c.+data2$e.m.+data2$t.m.d. 
table(data2$inf)

## summarize over wolves
wolf.sum <- data2 %>% group_by(wolfID) %>% summarize(sex=sex[1], cort=median(cortisol), pack=pack[1],
                                                     e.c.=max(e.c.),e.m.=max(e.m.),t.m.d.=max(t.m.d.))
wolf.sum$inf <- wolf.sum$e.c.+wolf.sum$e.m.+wolf.sum$t.m.d.

## parasite species
table(wolf.sum$e.c.)
table(wolf.sum$e.c.)/nrow(wolf.sum)
table(wolf.sum$e.m.)
table(wolf.sum$e.m.)/nrow(wolf.sum)
table(wolf.sum$t.m.d.)
table(wolf.sum$t.m.d.)/nrow(wolf.sum)

## all infections
table(wolf.sum$inf)
table(wolf.sum$inf)/nrow(wolf.sum)

table(wolf.sum$pack)

### CORTISOL
summary(data1$cortisol)
summary(wolf.sum$cort)

a <- ggplot(data=data1, aes(x=cortisol)) + geom_histogram(fill='gray', color='black') +
  ggtitle("All Cortisol Measurements")
b <- ggplot(data=wolf.sum, aes(x=cort)) + geom_histogram(fill='gray', color='black') +
  ggtitle("Median Cortisol Measurements") + xlab("cortisol")
a
b

c <- ggplot(data=wolf.sum, aes(x=sex, y=cort, color=sex)) + geom_boxplot() + geom_jitter(width=0.2) +
  ggtitle("Median Cortisol Measurements") + ylab("cortisol") + xlab("sex") + theme(legend.position="none")
c

d <- ggplot(data=wolf.sum, aes(x=pack, y=cort, color=pack)) + geom_boxplot() + geom_jitter(width=0.2) +
  ggtitle("Median Cortisol Measurements") + ylab("cortisol") + xlab("pack") +theme(legend.position="none")
d
  
plot_grid(a,b,c,d, nrow=2)
ggsave("cortisol_plots.png", height=10, width=10, units="in", dpi=250)

cort.inf <- data1[!is.na(data1$result),]
e <- ggplot(data=cort.inf, aes(x=as.factor(result), y=cortisol, color=as.factor(result))) + geom_boxplot() + 
  geom_jitter(width=0.2) + ggtitle("All Cortisol Measurements") + 
  ylab("cortisol") + xlab("infection status") + theme(legend.position="none")
f <- ggplot(data=data2, aes(x=as.factor(inf), y=cortisol)) + geom_boxplot(outlier.color='white') + 
  geom_jitter(width=0.2, aes(color=wolfID)) + ggtitle("Known Wolf Cortisol Measurements") + 
  ylab("cortisol") + xlab("infection status") + theme(legend.position="right")
plot_grid(e,f, nrow=1)
ggsave("infection_cortisol_plots.png", height=5, width=11, units="in", dpi=250)






