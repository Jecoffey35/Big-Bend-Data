# Libraries ---------------------------------------------------------------
#required packages: use install.package() to install if not already
require(vegan)#dist matrices, NMDS stuff
require(GGally)# ggpairs extension
require(ggplot2)#plotting
require(ggsignif)
theme_set(theme_bw())
library(ggord)
require(ggdendro)#dendrogram ggplot manips
require(tidyr)#data shaping
library(plyr) #Mostly for re-mapping factor values
require(dplyr)# data manip
require(factoextra)# stuff with ordination
require(reshape)
require(cluster)# manipfor dendro and clustering
require(indicspecies)
require(dynamicTreeCut)#for help cutting trees
require(fpc)
library(MASS)
library(SPECIES)# For richness/abundance estimation
library(RVAideMemoire)#pairwise permanova
library(stringr)
##Read in data -----------------------------------------------------------
#setwd("~/BIBE Analysis") #for veg and birds
#setwd("~/BIBE Analysis/Occupancy/Butterfly") #for veg and birds

#bird.mn.occ <- read.csv("occupancy_means.csv", header=T)[-1]#don't read in first column, just nonsense
#Burd abundance matrix-both years
bird.abund.mat <- read.csv("BUP_mns_birds_category_updates.csv")#[-1]
colnames(bird.abund.mat)[1] <- "site"
head(bird.abund.mat)
dim(bird.abund.mat)# This still has all 168 sites, so will need to exclude some 

#bt.abund <- read.csv("BUP_mns_butt_sub.csv")
#dim(bt.abund)#168
#head(bt.abund)#only 10 species

#Site groupings
remove(habitats)
habitats <- read.csv("Bird_point_treatment_categories_all.csv")[c(1,2)]
#habitats <- read.csv("Butterfly_transect_habitat_cats.csv")[,c("site","final_cat")]
head(habitats)
#colnames(habitats)[2] <- "cat_final"
#Environmental variables
veg <- read.csv("env_all.csv")[c(1:67),c(2,3,5,7,10,11)]
#veg <- read.csv("Butterfly_environmental_vars.csv")[,-1]
dim(veg)
head(veg)
treated.only <- habitats[which(habitats$cat_final == "new_treated"|habitats$cat_final =="old_treated"|habitats$cat_final =="untreated_cane"),"site"]
length(treated.only)# 56 for butterflies, 67 for birds

#OLD CODE FOR TREATING YEARS SEPARATELY
#okay, need to combine rows by site to get abundance matrix for both years----
#abund16 <- bird.abund.mat[which(bird.abund.mat$year=="2016"),]
#abund17 <- bird.abund.mat[which(bird.abund.mat$year=="2017"),]
head(abund17)
str(abund16)
abund16 <- abund16[which(abund16$site %in% abund17$site==TRUE),] 
abund17 <- abund17[which(abund17$site %in% abund16$site==TRUE),] 
rownames(abund16) <- abund16$site
rownames(abund17) <- abund17$site
abund16 <- abund16[-c(1,2)]
abund17 <- abund17[-c(1,2)]
head(abund16)
mcol <- match(colnames(abund16),colnames(abund17))#Match columns
m3 <- abund16#assign a matrix
m3[,mcol]<-m3[,mcol]+abund17#fpr the matching columns of each matrix, add together
tot.abund <- m3#rename
head(tot.abund)
#Additional exclusions for upland sites, as defined by creosote percentage greater than 5
#upland.sites <- c("BO3","SC7","TS1", "GS7","BO4","BO7","TS2","TS6","GP5","RP3","BD4","GP1")#12 sites
#Read in veg variables
#veg2017 <- read.csv("ALL_proofed_bird_veg_2017.csv",header=T)[,c("site","Total_willow","circle_radius", "Grass","herbaceous","live_arundo","arundo_treatment_regrowth","Mesquite_tree","Mesquite_shrub","salt_cedar_tree","salt_cedar_shrub")]
#veg2017 <- veg2017[which(veg2017$circle_radius=="100"),] # Use 100m because available for both years
tot.abund$site <- rownames(tot.abund)
tot.abund <- bird.abund.mat[which(bird.abund.mat$site %in% treated.only ==TRUE),]
head(tot.abund)
dim(tot.abund)# 55 by 32
#quickly get species totals only within subset sites 
colSums(tot.abund[,-32])

# Merge all environmental variables into site by variable matrix ----------
env.all <- veg[which(veg$site %in% treated.only == TRUE),]
#butt.all <- bt.abund[which(bt.abund$site %in% treated.only == TRUE),]
bird.all <- bird.abund.mat[which(bird.abund.mat$site %in% treated.only == TRUE),]
dim(bird.all)#67 by 26
head(env.all)

#Before we go any further, create new df where we merge veg with habitat data by sites, so 
#we can symbolize NMDS plot by that later
veg.hab <- merge(env.all,habitats,by = "site")
head(veg.hab) #hmm only 54?
#veg.hab <- veg.hab[which(veg.hab$site %in% c("CO1", "CO4","CO5")==F),]

# Bird dataframe ----------------------------------------------------------
#Merge with habitat vars
bird.hab <- merge(bird.all,habitats,by = "site")
head(bird.hab)
rownames(bird.hab) <- bird.hab$site
tot.abund <- bird.hab[,-c(1,8,22)]#Remove site column, as well as CLSW and VABU bc strong outliers
head(tot.abund)
colnames(tot.abund)[24] <- "Restoration_group"
tot.abund$Restoration_group <- recode_factor(tot.abund$Restoration_group,new_treated="Burned <4 yrs",old_treated="Burned >4 yrs", untreated_cane= "Unburned A.donax")
levels(tot.abund$Restoration_group)
tot.abund$Restoration_group <- droplevels(tot.abund$Restoration_group)

#butt.hab <- merge(butt.all,habitats,by = "site")
#head(butt.hab)
#nrow(butt.hab)#67 sites, subset by treatment relevant

#create function to rescale NDVI to be 0-100----
rescale <- function(x) (((100-0)/(max(x)-min(x)))*(x-max(x)))+100
#Test on small vector with similar range to NDVI
vec <- c(-1,0,.5,1)
rescale(vec)
head(veg.hab)
veg.hab$NDVI_new <- rescale(veg.hab$NDVI_new)
veg.hab$naip_text_cv <- rescale(veg.hab$naip_text_cv)
veg.hab$cane_rs <-rescale(veg.hab$cane_rs)
veg.hab$herbs <-rescale(veg.hab$herbs)
veg.hab$tree_cover <-rescale(veg.hab$tree_cover)

#for buttefly:
veg.hab$ndvi_mn <- rescale(veg.hab$ndvi_mn)
veg.hab$text_cv <- rescale(veg.hab$text_cv)
range(veg.hab$cane_rs)
#Drop out non-numerical columns, like site
head(veg.hab)
colnames(veg.hab)[7] <- "Restoration_group"
#Re-factor levels of the habitats so they look nice on plots
veg.hab$Restoration_group <- droplevels(veg.hab$Restoration_group)
veg.hab$Restoration_group <- recode_factor(veg.hab$Restoration_group,new_treated="Burned <4 yrs",old_treated="Burned >4 yrs", untreated_cane= "Unburned A.donax")
levels(veg.hab$Restoration_group)
rownames(veg.hab) <- veg.hab$site
env.final <- veg.hab[,-c(1,2)]#remove site and NDVI columns
head(env.final)
#Pairs plot----
cor(env.final[,-5])
head(env.final)
str(env.final)
#aggregate(env.final$herbaceous~as.factor(env.final$cat_final),FUN = function(x) mean(x))

hist(env.final$tree_cover)
#none of the vegetation variables are normally distributed, but can probably wait for the
#wisconsin double standardization to correct for this during nmds

# #Pairs plot to look at distributions and correlations -------------------
pairs.plot <- ggpairs(env.final[,-5],upper =list(continuous=wrap('cor',size=6,col='black')),lower=list('dot'),columnLabels = colnames(env.final[,-6]))+
             ggtitle("Bird Vegetation Variables")+
             theme_bw(base_family = "Times New Roman",base_size = 12)+
             theme(plot.title = element_text(hjust = 0.5)) 
pairs.plot
#ggsave(pairs.plot,filename = "bird_ggpairs.png",device = "png",dpi = 1200)
#Change some column names around so they'll look nicer on plots----
head(env.final)
#colnames(env.final)[1]<-"Grass"
#colnames(env.final)[2]<-"Herbaceous"
#colnames(env.final)[3]<-"Vegetation_Greenness"
colnames(env.final)<-c("Vegetation_Structure","A.donax_cover","Herbaceous",
"Tree_cover","Restoration_group")
head(env.final)

# okay now standardize and calculate distance matrix -------------------------------------------------------------

#Create and explore distance matrix---
#Can't have negative values, so need to standardize with something other than mean=0
#bcdist <- vegdist(env.all,method = "bray")#Bray Curtis distance matrix=default. This is very sensitive to weird values
head(env.final)
dim(tot.abund)
envbray <- vegdist(env.final[,-5], method="bray",binary = F)
birdbray <- vegdist(tot.abund[-24], method="bray",binary = F)
#btbray <- vegdist(tot.abund[-11], method="bray")
birdbray
hist(envbray, breaks=40, col='gray') #more peaked and multimodal than unstandardized
#Calculate Shannon diversity
bird.shan <- diversity(x = tot.abund[-24],index = 'shannon')
head(bird.shan,5)

# #Plot shannon diversity as it relates to group membership ---------------
H.dat <- read.csv("Simpson_diversity_birds.csv",sep = ",")
head(H.dat)
H.dat <- merge(H.dat,y=habitats,by = "site")
H.dat$cat_final <- droplevels(H.dat$cat_final)
H.dat$cat_final <- recode_factor(H.dat$cat_final,new_treated="Burned <4 yrs",old_treated="Burned >4 yrs", untreated_cane= "Unburned A.donax")
levels(H.dat$cat_final)

means <- aggregate(D ~  cat_final, H.dat, mean)
means$D <- round(means$D,3)
means

Q <- ggplot(H.dat, aes(x = reorder(cat_final, D, FUN = median), y = D))+ 
  geom_boxplot()+
  geom_signif(comparisons = list(c("Burned <4 yrs", "Burned >4 yrs","Unburned A.donax")), 
              map_signif_level=TRUE,)+
  stat_summary(fun.y = mean, geom="point")+
  geom_text(data = means, aes(label = D, y = D))+
  ggtitle("Shannon H by Habitat Group")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        plot.title = element_text(hjust=0.5),axis.text = element_text(size = 12))
Q

pairwise.perm.t.test(resp = H.dat$H,fact = H.dat$cat_final,p.method = "bonferroni")


# #Non-metric dimensional scaling ordination ------------------------------
#wisconsin double standardization and possibly also sqrt transformation if 
#autotransform not set to false, note if you use distance matrix as input you won't be able to 
#get species scores
#stress values less than 0.2 are considered good.The smaller the stress value,the better the match 
#between data and ordination distances and thus the better the new two-dimensional 
#configuration represents the patterns based on all the original variables.

#For butterflies-need to ensure same sites in each
#tot.abund <- tot.abund[which(rownames(tot.abund) %in% c("CO1","CO4","CO5")==F),]#missing CO1, CO4
#dim(tot.abund)
head(env.final)
remove(vegord)
vegord <- metaMDS(veg.hab[,-c(1,2,7)], try=c(50,500),k=2,wascores = T,autotransform = F,distance = "bray")
#btord <-metaMDS(tot.abund[,-11],k=3, try=c(50,500))
birdord <- metaMDS(tot.abund[,-24], k=3, try=c(50,500),wascores = T, distance = "bray")
stressplot(birdord)
stressplot(vegord)
vegord
birdord
#text(x=0.6,y = 0.015,labels = "Stress=0.141")#For bird veg
text(x=0.22,y = 0.1,labels = "Stress=0.193")#for bird ord
#text(x=0.35,y = 0.1,labels = "Stress=0.192")
#title("Butterfly Community NMDS Stressplot")

##NMDS points
veg.NMDS.data <- env.final
veg.NMDS.data$NMDS1<-vegord$points[ ,1]
veg.NMDS.data$NMDS2<-vegord$points[ ,2] 
veg.NMDS.data
bird.NMDS.data <- tot.abund
bird.NMDS.data$NMDS1<-birdord$points[ ,1]
bird.NMDS.data$NMDS2<-birdord$points[ ,2] 

#bt.NMDS.data <- tot.abund
#bt.NMDS.data$NMDS1<-btord$points[ ,1]
#bt.NMDS.data$NMDS2<-btord$points[ ,2] 

#Permutations of ordination with variables to get p-values
remove(vegfit)
vegfit <- envfit(vegord,env = env.final[,-5], display = 'sites', perm=9999, 
                 choices = c(1,2))
vegfit$vectors

birdfit <- envfit(birdord,env = env.final[,-5], display = 'sites', perm = 9999,
                  choices = c(1,2))
birdfit$vectors

#remove(btfit)
#btfit <- envfit(btord,env = env.final[,-7], display = 'sites', perm = 9999, choices = c(1,2))
#btfit$vectors

#Examine species and site scores
variableScores <- birdord$species
sampleScores <- birdord$points

# Summarizing Species Data ------------------------------------------------
#yr2sites <- c("CP10","CP14","LC4","LC5","SV5","SV6")# These are the sites that
#should have the greatest discrepancies between vegetation and categories if we just use 2017 data

# data for the envfit arrows
env.scores <- as.data.frame(scores(vegfit, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names
env.scores
# function for ellipsess - just run this, is used later
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#Pairwise adonis-----
bird.test <- anosim(birdbray,grouping = tot.abund$Restoration_group,permutations = 9999)
summary(bird.test) #anosim r is 0.23, p of 0.0004 so similar with some differences.

#and for butterflies
bt.test <- anosim(btbray,grouping = tot.abund$cat_final,permutations = 9999)
summary(bt.test) #anosim r is 0.2347, so similar with some differences. p of 0.0004
plot(bt.test)

veg.dif <- anosim(envbray,grouping = env.final$Restoration_group,permutations = 9999)
summary(veg.dif) #no transformations, bray distance--anosim r is 0.1352 so not much difference. p-value of 0.0008

#Look closer to see which variables vary the most between groups
#bt.dif <- anosim(btbray~factor(cat_final),data=tot.abund,permutations = 9999)
head(env.final)
#ord <- dbrda(birdbray ~ Vegetation_Structure + A.donax_cover +
#               Herbaceous + as.factor(Restoration_group), data = env.final)
#anova(ord, by = 'margin')

pairwise <- pairwise.perm.manova(envbray,fact=factor(tot.abund$cat_final),nperm=9999,p.method = "bonferroni",test="Wilks")
pairwise

#Plot with sites colored by burn history, add environmental vectors and species----
table(tot.abund$Restoration_group)#count of number of sites in each category

#plot ordinations----
#try plot with symbols instead of colors
#specify points as black circle, grey circle, and hollow diamond
pch_lookup <- c("Unburned Arundo" = 19, "Burned >4yrs ago" = 22, "Burned <4yrs ago" = 5)
#Vegetation plot----
plot.new()
plot(vegord, type="t",choices = c(1,2), cex=0.8)
points(vegord,display="sites",pch=pch_lookup[as.factor(env.final$Restoration_group)],bg="darkgrey",choices = c(1,2))
plot(vegfit, col='black',cex=1, choices=c(1,2),p.max = 0.05)
legend(x=0.75,y=0.75, legend=c("Unburned Arundo", "Burned >4yrs ago", "Burned <4yrs ago"), 
      pch=c(5,22,19),pt.bg="darkgrey",cex = 1,bty="o",pt.lwd = 1)
text(x=1.2,y = -0.65,labels = paste("ANOSIM R = 0.1352",paste("p-value=",veg.dif$signif),sep = "\n"))
title("Environmental Variables")

# data for labelling the ellipse
NMDS.mean=aggregate(veg.NMDS.data[ ,c("NMDS1", "NMDS2")], 
                         list(group = veg.NMDS.data$Restoration_group), mean)

# Try with ggplot. Black and white. Dots symbolized by habitat type, vectors overlayed 
remove(veg.NMDS.data)
head(env.scores)
env.scores$newx = str_wrap(env.scores$env.variables, width = 8)
env.scores
#Rename groups
NMDS.mean$group <- c("Burned <4yrs","Burned >4 yrs","Unburned A.donax")
NMDS.mean
mult <- 1
veg.NMDS.data
veg.gg <- ggplot(data = veg.NMDS.data, aes(x = NMDS1,y = NMDS2))+ #sets up the plot. brackets around the entire thing to make it draw automatically
  geom_point(aes(shape = Restoration_group,colour=Restoration_group), size = 3) +
  scale_colour_grey(start = .2, end = .75)+
  geom_segment(data = env.scores,
               aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = NMDS1, y = NMDS2, label=str_wrap(c("Vegetation Texture","A.donax","Herbaceous","Tree"),width = 10)),
            size = 5,nudge_x = -0.25,nudge_y = 0.1,
            hjust = -0.25)
veg.gg
#ggsave(plot = veg.gg, filename = "Veg_NMDS2020.jpeg",dpi = "retina",width = 9.5,height = 6,units = "in")

# Butterfly NMDS plot -----------------------------------------------------
bird.NMDS.data
bird.gg <- ggplot(data =bird.NMDS.data , aes(x = NMDS1,y = NMDS2))+ #sets up the plot
  geom_point(aes(shape = Restoration_group,colour=Restoration_group), size = 3) +
  scale_colour_grey(start = .2, end = .75)+
  geom_segment(data = env.scores,
            aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
            arrow = arrow(length = unit(0.10, "cm")), colour = "black") +  
  geom_text(data = env.scores, #labels the environmental variables
            aes(x = NMDS1, y = NMDS2, label=env.variables),hjust = -0.25)
bird.gg
#Bird plot----
head(tot.abund)#7 is CLSW 21 is VABU
plot.new()
plot(birdord, type="n",choices = c(1,2),cex.axis=1,display = "species")
points(birdord,dis="sites",pch=pch_lookup[as.factor(tot.abund$Restoration_group)],bg="darkgrey",cex=1.5,choices = c(1,2))
plot(birdfit, col='black',cex=1, choices=c(1,2),p.max = 0.05)
orditorp(x=birdord,display = "species",choices = c(1,2),cex = 1)
legend(x=0.3,y=0.3, legend=c("Unburned A.donax", "Burned >4yrs", "Burned <4yrs"), 
      pch=c(5,22,19),pt.bg="darkgrey",cex = 1,bty="n")
text(x=0.39,y = -0.15,labels = paste("ANOSIM R = 0.25",paste("p-value=",bird.test$signif),sep = "\n"))

#Or fancy ellipse/spider plots----
plot.new()
plot(vegord, display ="sites", type='p',choices=c(1,2),cex=0.7)#for quick labels, use type='t'
plot(vegfit, col='black',cex=0.7, choices=c(1,2))
orditorp(vegord,pch = 1,display="sites",col=col_palette,air=1,choices = c(1,2))#,choices = c(2,3)
ordispider(ord = vegord,groups = habitats$cat_final,label = F,spiders="median")
ordiellipse(ord = vegord,groups = habitats$cat_final,kind = "ehull",label = T)
legend('topcenter', legend=c("Never Burned", "Burned >4yrs ago", "Burned <4yrs ago"), 
       col=unique(col_palette), pch = 1,cex = .5,bty="n")
dev.off()
#par(mfrow=c(2,2))
plot(btord, display ="sites", type='p',choices=c(1,2),cex=0.7)#for quick labels, use type='t'
#orditorp(btord,pch = 1,display="sites",col=col_palette,air=1)#,choices = c(2,3)
ordiellipse(ord = btord,groups = tot.abund$cat_final,kind = "ehull",label = T)
plot(btfit, col='black',cex=0.7,choices = c(1,2))#, choices=c(2,3)
orditorp(btord,display="species",col="darkgreen",air=1,cex = 1)#,choices = c(2,3)
dev.off()

#plot variable scores with arrows, need species scores
plot(btord, display = "sites", type = "p")
#with(tot.abund, ordiellipse(btord, groups = as.factor(cat_final), col = "black",label = T))
with(tot.abund, ordispider(btord, groups = as.factor(cat_final), col = "black",label = F))
orditorp(btord,display="species",col="black",air=1,cex = 1)#,choices = c(2,3)


# Hierarchical Clustering -------------------------------------------------
#LOOK INTO RECURSIVE PARTITIONING for decision tree stuff
tree1 = rpart(factor(burned_less_4_2016) ~ ., data=env.all)
tree1_p = as.party(tree1)
tree1_p

hclust.veg <- hclust(d = envbray,method = "ward.D2")
hclust.stand <- hclust(d = birdbray,method = "ward.D2")
hclust.stand
#Explore number of groups, clustering metrics
# To determine optimal k, try
head(env.stand)
fviz_nbclust(env.all, kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2)+
  labs(subtitle = "Elbow method", title="Optimal number of clusters, bird vegetation 2016")

#Silhouette method, which
fviz_nbclust(env.all, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method",title="Optimal number of clusters, bird vegetation 2016")

#Same, but with gap statistic
set.seed(53)
#Same, but with gap statistic
birdgaps <-clusGap(env.all,FUN=kmeans,nstart=25,B = 500, K.max=8)
birdgaps
fviz_gap_stat(birdgaps,linecolor = "blue")+
  labs(subtitle = "Gap statistic method",title="Optimal number of clusters, bird vegetation 2016")

# Cut tree into 6 groups
sub_grp <- cutree(hclust.veg, k = 2)#Learned from following this through
#that non-standardized values don't do as good a job of pulling out willow, cane etc.
sub_grp

# Look at spearman's correlations of veg variables with NMDS axes

# Number of members in each cluster
table(sub_grp)
head(env.all16)
env.all[("BC6"),]
env.all16$site <- rownames(env.all16)
k6.df <- as.data.frame(sub_grp)
k6.df$site <- rownames(k6.df)
k6.df
env.group.membs <- merge(env.all16,k6.df, by="site")
env.group.membs
#write.csv(env.group.membs, file="Bird_clusters_2016.csv")
head(env.all16)
env.all16 <- env.all16[-1]
head(env.all)
#Now use these sub group memberships to visualize
fviz_cluster(list(data = env.all, cluster = sub_grp),show.clust.cent = T,
             main = "Hierarchical Clusters Bird Vegetation, k=2")

plot(hclust.stand, cex = 0.5)
rect.hclust(hclust.stand, k = 5, border = 2:5)

#Try to assign color groups based on burn year
#First, make sure burn years have same sites as env.all
env.all$site <- rownames(env.all)
burn.yrs <- years[(years$site %in% env.all$site),]
nrow(burn.yrs)

head(env.all)
env.dist = dist(env.all)
plot(silhouette(cutree(divisive.clust,4),env.dist))# This can give a measure of how strong of a structure was found
#less than 0.35 means no structure, >.5 is reasonable

#Or try with GG
hcdata <- dendro_data(hclust.stand, type="rectangle")
hcdata

#merge with burn data
hcdata$labels <- merge(x = hcdata$labels, y = burn.yrs,  by.x = "label", by.y = "site")

ggplot() +
  geom_segment(data=hcdata$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = hcdata$labels, aes(x=x, y=y, label=label, colour =as.factor(burn.yrs$year_since_burn), hjust=0), size=1) +
  #geom_point(data = hcdata$labels, aes(x=x, y=y), size=3, shape = 21) +
  coord_flip() +
  scale_y_reverse(expand=c(0.2, 0)) +
  scale_colour_brewer(palette = "Dark2") + 
  theme_dendro() 

# Indicator Species Analysis  ----------------------------------------------------------
#quickly play with anosim

ab <- anosim(birdbray,grouping = sub_grp,permutations = 999)
summary(ab)

##indicspecies package stuff
head(env.all)
env.all <- env.all[-11]
## Use group assignments from hclust to determine associations
diffs <- cbind(sub_grp_no5, sub_grp2)
diffs

## Determine sensitivity of individual species
B=strassoc(env.all, cluster=sub_grp,func="B")# sub_grp2 is k=6, fits more ecologically
A=strassoc(env.all, cluster=sub_grp,func="A")
A
B
## Select species with more than 20% of sensitivity for the first group
sel=which(B[,1]>0.1)
sel
## Run indicator analysis with species combinations for the first group
sc= indicators(X=env.stand[,sel], cluster=sub_grp_k5, group=1, verbose=TRUE, At=0.5, Bt=0.2)
## Determine the coverage of the selected set of indicators
coverage(sc)
## Plot the coverage against the threshold At
plotcoverage(sc)
plotcoverage(sc, max.order=2, add=TRUE, lty=2)

## Runs the combination analysis using IndVal.g as statistic
envpt = multipatt(env.stand, sub_grp_k5, func = "r.g",duleg=T,control = how(nperm=999))
envptInd = multipatt(env.stand, sub_grp_k5, func = "IndVal.g",duleg=T,control = how(nperm=999))
envpt
## Determines the coverage for each site group combination
envptInd
#write.csv(envpt,file="bird_env_var_group_contribs.csv")
coverage(env.all, envpt)