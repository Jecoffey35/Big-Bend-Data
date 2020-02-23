# #load libraries ---------------------------------------------------------
require(vegan) #community ecology stuff
require(ggplot2) #plotting
require(reshape2) #data manipulation
require(unmarked) #multispecies occupancy modelling
require(tibble) #data manip
require(tidyr) #data manip
require(tidyverse)
require(dplyr) #data manip
require(MuMIn) #model dredging, install if necessary
require(AICcmodavg) #model averaging. fitting in unmarked
require(DiversityOccupancy)
require(boot)#for inverse logit
# #Load in data -----------------------------------------------------------
birds <- read.csv("ALL_BIRD_DATA.csv", header = T, sep = ",", stringsAsFactors = TRUE)[c(2,4,7,8,9,11)]
head(birds)
bird.filt <- birds[(!grepl("FL|110", birds$distance)),]
unique(bird.filt$distance)
##Create lists of species for subsetting if desired###
final <- read.csv("Species_list.csv",sep = ",",stringsAsFactors = T)
head(final)
colnames(final)[1] <- "Species"
bird.char <- bird.filt[which(bird.filt$species %in% final$Species == TRUE), ]
unique(bird.char$species)
bird.char$species <- factor(bird.char$species)
levels(bird.char$species) #No species we were trying to remove, good.

#Now we don't need distance column anymore
head(bird.char)
bird.dat <- bird.char[,c(1,2,3,4,6)]
head(bird.dat)

#Create hybridized variable---
bird.dat$visityear <- paste(bird.dat$year, bird.dat$visit, sep = "_") 
head(bird.dat)

bird.hist <- table(bird.dat[c("visityear","site","species")],useNA="no")
bird.hist
transpose <- t(bird.hist[,,1])
str(transpose)
#Okay now deal with NAs for Co and that one SC site
na.dat <- read.csv('bird_visit_history.csv')
row.has.na <- apply(na.dat, 1, function(x){any(is.na(x))}) 
na.dat2 <- na.dat[row.has.na,] # This is a subset of any sites with NAs
na.dat2
#change colnames to match visityear. 
colnames(na.dat2) <- c("site","2016_1","2016_2","2016_3","2017_1","2017_2","2017_3")
na.dat2

na.mlt <- melt(na.dat2,id.vars = "site",measure.vars = c(2,3,4,5,6,7),variable.name = "visityear")
na.mlt
true.na <-na.mlt[is.na(na.mlt$value),]
true.na
true.na$combo <- paste(true.na$site,true.na$visityear,sep = "_")

#okay now convert matrix to df and then fill in nas by merging
bird.df <- data.frame(bird.hist)
head(bird.df)
bird.df$combo <- paste(bird.df$site,bird.df$visityear,sep = "_")
bird.df$Freq[bird.df$combo %in% true.na$combo] <- "NA"
#Check that it worked
bird.df[bird.df$site =="CO1",]
bird.df <- bird.df[,-5]# remove combo column now
head(bird.df)

# Calculate Site Species Richness -----------------------------------------
#Prior to calculating this, going to need to take the NAs in the file that represent sites that were
#where we're only using 2016 or 2017 data. So we could do that by filtering out the combination of 
#sites and years. i.e. don't use CO5 if the year is 2016, or likewise for LC.

#Wait, duh need to remove NAs first. 
bird.nose <- bird.df[!bird.df$Freq=="NA",]
bird.noos <- bird.nose[!bird.nose$Freq=="0",]
nrow(bird.noos)
bird.noos
nrow(bird.df)
head(bird.nose)

# So this will be a count of unique species per site.Need site by species matrix with inside as frequency
bird.df2 <- melt(data = bird.noos,id.vars = "site",measure.vars = "species")
head(bird.df2)

bird.rich <- bird.df2 %>% group_by(site) %>% summarise(length(unique(value)))
bird.rich

rich.xp <- as.data.frame(bird.rich)
rich.xp
colnames(rich.xp) <- c("site","richness")

#write.csv(rich.xp, file = "Bird_richness_both_yrs.csv")

#now re-matrify it
bird.mat <- acast(bird.df, site~visityear~species, value.var="Freq")
bird.mat#Success!

head(bird.mat)
rowSums(x = bird.mat)

#For convenience, just use Lars' covariate df
obsCov <- read.csv("BEVI_detection history and covariates.csv")[,8:19]
obsCov
dim(obsCov) #168 by 12 columns --> observers and years

class(bird.mat) <- "numeric"
dimnames(bird.mat)#168 by 6 cols by 25 bird matrix slices
bird.mat[,,1]
bird.mat[]
#Checking out some dimensions/observations prior to modelling
check1 <- bird.filt[which(bird.filt$species == "ATFL"),]
check2 <- birds[which(birds$species == "ATFL"),]
dimnames(bird.mat)
check2
check1
length(unique(check1$site))
unique(check1$site)

df <- unmarkedFramePCount(y=bird.mat[,,1], obsCovs = list(obs=obsCov[,1:6], year=obsCov[,7:12]))
summary(df)
  
#For now, using observer and year models for abundance estimates
#det.mods <- list(null <- pcount(~1~1, df), 	# null model
#obs <- pcount(~obs~1, df), 	#observer influenced detection prob 
#yr <- pcount(~year~1, df), #year influenced detection prob
#both <- pcount(~obs + year ~ 1, df))#year and observer influenced detection prob		

#Define models outside of list:
#null <- pcount(~1~1, df)	# null model
#obs <- pcount(~obs~1, df) 	#observer influenced detection prob 
#yr <- pcount(~year~1, df) #year influenced detection prob
both <- pcount(~obs + year ~ 1, df)
both  
p_det <- backTransform(null, "det") #can only use this with null model
p_det@estimate
g.det.prob[c(1:2),] <- p_det

ad.boot <- nonparboot(both, B = 100)#can compare parametric and nonparametric results here
ad.boot

predict(both, type="det")
SE(both, type="det") #method='nonparboot') 
hist(EBUP1) # EBUP is adjusted site-specific abundance estimates.
CI1 <- confint(re1, level = 0.95)                   # also just cols in re1
# estimate of global abundance across all sites, with EB interval (CI)

s <- nrow(dat)# number of rows           
rbind(PAO = c(Estimate = sum(EBUP1), colSums(CI1))/s)

##or model selection for best fit
#modsel.tab <- aictab(cand.set = det.mods,modnames = c("null","obs","yr","both"),sort = T)
#modsel.tab
#remove(det.mods, df, modsel.tab)

#store model selection results 
mod.results <- data.frame(matrix(nrow = 168, ncol = 290))
mod.results
colnames(mod.result) <- c(occupancy.list[i], paste0("K","AICc","Delta_AICc","Cum.Wt","LL"))
mod.results[i:5,i:8]

# empirical Bayes estimate of proportion of sites occupied
  #re1 <- ranef(null)						# estimates table
  #re2 <- ranef(obs)
  #re3 <- ranef(yr)
  re4 <- ranef(both)
  re4
  
remove(null,obs,yr,both,re1,re2,re3,re4)

# For loop across all species output --------------------------------------------------
results <- matrix(, nrow = 168, ncol = 25) #create empty dataframe for storing output

for (i in seq.int(from = 1,to = 25,by=1)){
  df <- unmarkedFramePCount(y=bird.mat[,,i], obsCovs = list(obs=obsCov[,1:6], year=obsCov[,7:12]))
  both <- pcount(~obs + year ~ 1, df)
  re <- ranef(both)
  results[,i] <- bup(re,stat = 'mean')
}
results
dimnames(bird.mat[,,c(1:25)])
bird.names <- dimnames(bird.mat)[[3]]
rownames(results) <- rownames(bird.mat[,,1])
colnames(results) <- bird.names
head(results)
#write.csv(x = results,file = "BUP_mns_birds_2020.csv")

#Goodness of fit tests
#Estimate the posterior distribution of the latent variable given the data
#and the estimates of the fixed effects (the MLEs). The mean or the mode of the estimated posterior
#distibution is referred to as the empirical best unbiased predictor (EBUP), which in unmarked can
#be obtained by applying the bup function to the estimates of the posterior distributions returned by
#the ranef function. 
EBUP <- bup(re4, stat="mean")
hist(EBUP) 
CI <- confint(re4, level=0.95)# These are kinda wide for several sites, 
#would be good to use a secondary measure (probs bootstrapping and SE) to account for uncertainty
CI

# estimate of abundance across all sites, with EB interval (CI)
rbind(PAO = c(Estimate = sum(EBUP), colSums(CI))/168) # I don't get why you'd want the estimate of abundance across all sites, not particularly useful to us here
#would be good to have CI or SE for all sites.
#rbind(EBUP, CI)

# model fit function: chi squared statistic and parametric bootstrap
chisq <- function(fm) {
  umf <- getData(fm)
  y <- getY(umf)
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}

# apply model fit function to a model (can take a minute!)
pb <- parboot(both, statistic=chisq, nsim=100, parallel=F)
pb
