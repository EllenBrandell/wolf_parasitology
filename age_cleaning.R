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
# data1 <- merge(data1, meta[,c(1,4)], by="scatID", all.x=T) # I want site


######################################## DIAMETER PATTERNS
ggplot(data=data1, aes(x=scat_diameter)) + geom_histogram(fill='gray', color='black')


## look at different diameter measurements by unique wolves



## look at different diameter measurements by collared wolves (known age)



######################################## EXPLORING CUTOFFS
cut1 <- 2.0  # less than OR equal to this is a pup
cut2 <- 2.5  # less than OR equal to this is a pup
data1$age <- ifelse(data1$scat_diameter<=cut1, "Pup", "Adult")
data1$age2 <- ifelse(data1$scat_diameter<=cut2, "Pup", "Adult")
# make sure age isn't changing unless the wolf ages from Juv -> Adult
ggplot(data=data1, aes(x=scat_diameter)) + geom_histogram(fill='gray', color='black') +
  geom_vline(xintercept=cut1, col='red', size=1.3)


table(data1$age)
table(data1$age)/nrow(data1)
# "questionable"/"in between" scat diameters based on our cutoff
table(data1$scat_diameter>cut1-0.3 & data1$scat_diameter<cut1+0.3)
# there are quite a few of these "in between" sizes


###### make a data set of JUST KNOWN WOLVES
wolves1 <- data1[!is.na(data1$wolfID),]

wolf.sum <- wolves1 %>% group_by(wolfID) %>% summarize(sex=sex[1], ages=length(unique(age)), age=age[1], 
                                                       age2=age2[1], pack=pack[1],
                                                     e.c.=max(e.c.),e.m.=max(e.m.),t.m.d.=max(t.m.d.),
                                                     sites=length(unique(indiv.site)), nsamples=n())
table(wolf.sum$ages)

table(wolf.sum$age)
table(wolf.sum$age)/nrow(wolf.sum)

table(wolf.sum$pack, wolf.sum$age)

a <- ggplot(data=wolves1, aes(x=scat_diameter)) + geom_histogram(fill='gray', color='black') +
  geom_vline(xintercept=cut1, col='red', size=1.3) + geom_vline(xintercept=cut2, col='red', size=1.3, linetype="dashed") +
  xlab("scat diameter") + ggtitle("Scat Diameter Known Wolves")
a
b <- ggplot(data=wolves1, aes(y=scat_diameter, x=wolfID, color=wolfID)) + geom_point() + 
  ylab("scat diameter") + ggtitle("Scat Diameter Known Wolves") + 
  theme(axis.text.x=element_text(angle=90), legend.position="none") +
  geom_hline(yintercept=cut1, col='red', size=0.5) + geom_hline(yintercept=cut2, col='red', size=0.5, linetype="dashed")
b
c <- ggplot(data=wolves1, aes(y=scat_diameter, x=age, color=age)) + geom_boxplot() + geom_jitter(width=0.1) +
  ylab("scat diameter") + ggtitle("Scat Diameter Known Wolves") +
  geom_hline(yintercept=cut1, col='red', size=0.7)
c
d <- ggplot(data=wolves1, aes(y=scat_diameter, x=age2, color=age2)) + geom_boxplot() + geom_jitter(width=0.1) +
  ylab("scat diameter") + ggtitle("Scat Diameter Known Wolves") +
  geom_hline(yintercept=cut2, col='red', size=0.7, linetype="dashed")
d

ggplot(data=wolves1, aes(y=scat_diameter, x=sex, color=sex)) + geom_boxplot() + geom_jitter(width=0.2) +
  ylab("scat diameter") + ggtitle("Scat Diameter Known Wolves") + 
  geom_hline(yintercept=cut, col='red', size=0.7) + geom_hline(yintercept=cut2, col='red', size=0.7, linetype="dashed")


plot_grid(a,b,c,d, nrow=2, ncol=2)
ggsave("scat_diameter.png", height=11, width=11, units="in", dpi=250)


###### HOW MANY MISMATCHES DO WE HAVE BASED ON CUTOFF VALUE?
# can only use resampled wolves for this


## and importantly, we want to make sure the collared wolves are correctly classified:
# pup = 1211F, 1228F, 1229F
# adult = 969F, 907F, 996M, 1005F, 1155M, 1156M


######
# change scat 145 wolf 1228F to Pup


