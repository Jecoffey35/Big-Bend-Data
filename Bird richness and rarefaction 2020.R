# Load packages -----------------------------------------------------------
library(vegan)

# #Read in data -----------------------------------------------------------
rawdat <- read.csv("ALL_BIRD_DATA.csv")
#abdat <- read.csv("BUP_mns_birds_2020.csv")
head(rawdat)


# Species Accumulation Curves----
spec <- specaccum(tot.abund, method="random",permutations = 100)
plot(spec,ylab = "Species Number")
spec17 <- specaccum(abund.mat17, method="exact")
plot(spec17,ylim=c(0,50))
#try seperate species accumulation curves for different habitats
cane <- abund.cat[which(abund.cat$category2016 =="untreated_cane"),]
new <- abund.cat[which(abund.cat$category2016 =="new_treated"),]
old <- abund.cat[which(abund.cat$category2016 =="old_treated"),]
head(cane)
#plot them
cane.sp <- specaccum(cane[,-c(1,33)], method="exact")
new.sp <- specaccum(new[,-c(1,33)], method="exact")
old.sp <- specaccum(old[,-c(1,33)], method="exact")
cane.sp
plot(x=cane.sp$sites,y=cane.sp$richness,ylab = "Species Number")
plot(new.sp$sites,new.sp$richness)
plot(old.sp$sites,old.sp$richness)

#Now use specpool to get pooled estimates of richness by management habitat----
#species pool uses incidence data, so convert
head(abund.cat)
head(tot.abund)
incid.cat <- decostand(tot.abund, method = "pa")
incid.cat16 <- decostand(abund.cat[,-32], method="pa")
incid.cat17 <- decostand(abund.cat17[,-32], method="pa")
head(incid.cat)

# Play around with different vegan richness estimators----
pool17 <- specpool(incid.cat17,abund.cat17$category2017)
pool17
specpool(incid.cat,abund.cat$category2016)

specnumber(abund.mat16)
specnumber(abund.mat17)
#The following functions estimate richness using counts of species at sites,as well as the standard errors of the estimates. 
#These only concern the number of added species, and assume that there is no variance in the observed richness.
estimated.rich <- estimateR(abund.mat17)
summary(estimated.rich)
estimated.rich
#write.csv(estimated.rich, file="estimated.rich2016.csv")
incid.mat17$site <- rownames(incid.mat17)
incid.mat17 <- incid.mat17[which(incid.mat17$site %in% Bird.clusts17$site),]
abund.mat17$site <- rownames(abund.mat17)
abund.mat17 <- abund.mat17[which(abund.mat17$site %in% Bird.clusts17$site),]
dim(abund.mat17)
incid.mat17 <- incid.mat17[-58]
head(abund.mat16)
rarefy(abund.mat16, min(rowSums(abund.mat16)))
rarefy(abund.mat17, min(rowSums(abund.mat17)))
specaccum(abund.mat17,method = "exact")
#Try to plot comparisons among species richness estimators
droplevels(abund.cat17$category2017)
pool17
observed.se <-c(0,0,0,0,0,0,0)
pool17 <- cbind(pool17, observed.se)
pool17
poolse <-c(pool17$observed.se,pool17$chao.se, pool17$jack1.se, pool17$boot.se)
pooly <- c(pool17$Species,pool17$chao, pool17$jack1, pool17$boot)
birdpool17 <- ggplot(pool17, aes(x=rownames(pool17), y=pool17$Species)) + 
  geom_bar(position=position_dodge(), stat="identity",colour="black", size=.3) + # Thinner lines
  geom_errorbar(aes(ymin=pool17$Species-pool17$observed.se, ymax=pool17$Species+pool17$observed.se),size=.3,width=.2,position=position_dodge(.9)) +
  xlab(" ") +
  ylab("Number of Species") +
  ggtitle("Comparison of Species Estimators")

birdpool17
dim(incid.mat17)
divers.17 <- diversity(incid.mat16, index = "shannon")
length(divers.17)
#write.csv(divers.17,file="Shannon_diversity_by_site_2016.csv")

H <-diversity(abund.mat16) #Shannon diversity
J <- H/log(specnumber(abund.mat16)) # evenness

