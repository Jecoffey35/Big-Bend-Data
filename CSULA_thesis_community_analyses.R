# ##Read in libraries -----------------------------------------------------
require(vegan)#dist matrices, NMDS stuff
require(GGally)# ggpairs extension
require(ggplot2)#plotting
require(ggdendro)#dendrogram ggplot manips
require(tidyverse)#data shaping
require(dplyr)# data manip
require(factoextra)# stuff with ordination
require(reshape)
require(cluster)# manipfor dendro and clustering
require(indicspecies)
require(dynamicTreeCut)#for help cutting trees
require(fpc)
library(MASS)
library(RVAideMemoire)#pairwise permanova

# Read in data ------------------------------------------------------------
#species abundances
ab <- read.csv("BUP_mns_birds_category_updates.csv")# This is all species, need to subset by detection probability cutoff
head(ab)
#shortlist of species with higher detection probabilities
final <- read.csv("Species_list.csv",sep = ",",stringsAsFactors = T)
colnames(final)[1] <- "Species"
final$Species
#set site to rownames before subsetting
rownames(ab) <- ab$site
#subset colnames of ab by species in species list
abf <- ab[colnames(ab) %in% final$Species]
abf$site <- rownames(abf)
head(abf)
veg <- read.csv("env_all.csv")[,-c(1,3,12,14)]
head(veg)
#lastly need group memberships for additional subsetting
hab <- read.csv("Bird_point_treatment_categories_all.csv")[,c(1,3)]#only using 2017
head(hab)
levels(hab$category2017)
str(hab)
treated <- hab %>% filter(category2017 %in% c("untreated_cane","new_treated","old_treated"))  
nrow(treated)#66

mat <- abf[which(abf$site %in% treated$site==T),]
veg <-veg[which(veg$site %in% treated$site==T),]
nrow(mat)#subset to 66 sites

#Now let's look at some summary stats for the bird community
summary <- abf %>%
  summarise_if(is.numeric, list(median=median,mn=mean,min=min, max=max, sd=sd))
summary
#write.csv(summary, "Species_abundance_summary_subset_NEW.csv")

#Now examine correlations, prior to any standardizations
#site column to rowname
rownames(mat) <- mat$site
str(mat)
mat <- mat[,-26]
comm <- decostand(mat, method = "hellinger")
comm
# check total abundance in each sample
apply(comm, 1, sum)
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,method = "spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(veg,lower.panel =panel.cor)

#Now for ordination
comm.nmds <- metaMDS(comm, dist = "bray",k = 3,trymax = 200)
stressplot(comm.nmds)
# plot site scores as text
ordiplot(comm.nmds, display = "sites", type = "text")
ordiplot(comm.nmds, display = "species", type = "text")

levels(treated$category2017)
mds.fig <- ordiplot(comm.nmds, type = "none")
# plot just the samples, colour by habitat, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = treated$category2017 == "untreated_cane")
points(mds.fig, "sites", pch = 19, col = "red", select = treated$category2017 == "new_treated")
points(mds.fig, "sites", pch = 19, col = "orange", select = treated$category2017 == "old_treated")
# add confidence ellipses around habitat types
ordiellipse(comm.nmds, treated$category2017, conf = 0.95, label = TRUE)

# plot Species abundance. cex increases the size of bubbles.
ordisurf(comm.nmds, comm[, "COYE"], bubble = TRUE, main = "COYE abundance", cex = 3)

ordiplot(comm.nmds)
# calculate and plot environmental variable correlations with the axes use
# the subset of metadata that are environmental data
str(veg)
plot(envfit(comm.nmds, veg[, c(5,7,8,9,12)]))
