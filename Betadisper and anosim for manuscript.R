library(plyr)
library(ggplot2)
library(vegan)
library(ggrepel)#spacing for plot labels

bird <- read.csv("Bird_abundance_matrix_final.csv")[-1]
head(bird)
birdenv <- read.csv("Final_bird_veg_subset.csv") [,-1]
birdenv

#Butterfly data, may need different working directory "BT_for_Eric"
#setwd("~/Big-Bend-Data/BT_For_Eric")
bt.treated <- read.csv("abundance_max.csv")
colnames(bt.treated)[1] <- "site"
btenv <- read.csv("butterfly_veg_subset_2020.csv") [,-1]
head(btenv)
btenv <- btenv[which(btenv$final_cat %in% c("Burned >4 yrs","Unburned A.donax","Burned <4 yrs")),]
dim(btenv) #54 
dim(bt.treated)#56...

head(btenv)
head(bt.treated)
burns <- btenv[which(btenv$final_cat %in% c("Burned >4 yrs","Burned <4 yrs")),]#"
old.un <- btenv[which(btenv$final_cat %in% c("Unburned A.donax","Burned >4 yrs")),]
new.un <- btenv[which(btenv$final_cat %in% c("Unburned A.donax","Burned <4 yrs")),]
butts1 <- bt.treated[which(bt.treated$site %in% burns$site==T),]
butts2 <- bt.treated[which(bt.treated$site %in% old.un$site==T),]
butts3 <- bt.treated[which(bt.treated$site %in% new.un$site==T),]
butts <- bt.treated[which(bt.treated$site %in% btenv$site==T),]
btenv <- btenv[which(btenv$site %in% butts$site==T),]
dim(butts)#54

#Now create distance matrices-euclidean for veg, bray for community
ds.burns <- dist(burns[,3:8], method ="euclidean")
ds.oldun<- dist(old.un[,3:8], method ="euclidean")
ds.newun<- dist(new.un[,3:8], method ="euclidean")
ds.all<- dist(btenv[,3:8], method ="euclidean")
#Butterfly community matrices
ds.btall <- vegdist(butts[,-1],method = "bray")
ds.bt1 <- vegdist(butts1[,-1],method = "bray")
ds.bt2 <- vegdist(butts2[,-1],method = "bray")#oldun
ds.bt3 <- vegdist(butts3[,-1],method = "bray")

btb <- betadisper(d = ds.btall,group = btenv$final_cat, type = "centroid",add = T)
permutest(btb, pairwise = TRUE, permutations = 9999)
btb
remove(btb)

##########

sub1 <- birdenv[which(birdenv$cat_final %in% c("untreated_cane","new_treated")),]#
sub2 <- birdenv[which(birdenv$cat_final %in% c("untreated_cane","old_treated")),]#
sub3 <- birdenv[which(birdenv$cat_final %in% c("old_treated","new_treated")),]#
head(sub1)
bird1 <- bird[which(bird$site %in% sub1$site==T),]
bird2 <- bird[which(bird$site %in% sub2$site==T),]
bird3 <- bird[which(bird$site %in% sub3$site==T),]
head(bird3)
dim(bird3)#30
dim(sub3)
d1 <- vegdist(bird1[,-1],method = "bray")
d2 <- vegdist(bird2[,-1],method = "bray")
vd <- vegdist(sub1[,c(2:6)],method = "euclidean")
d3 <- vegdist(bird3[,-1],method = "bray")
vd
vdall <-vegdist(birdenv[,c(2:6)],method = "euclidean")
#betadisper
remove(vd,b)
b <- betadisper(d = vdall,group = birdenv$cat_final, type = "centroid",add = T)
permutest(b, pairwise = TRUE, permutations = 9999)
b

#anosim
anos <- anosim(x = d3, grouping = sub3$cat_final,distance = 'bray',permutations = 999)
anos
remove(anos,b)

####
birdord <- metaMDS(d1, try=c(50,500),k=3,wascores = T,autotransform = F, distance="bray")
vegord <- metaMDS(vd, try=c(50,500),k=3,wascores = T,autotransform = F, distance="euclidean")

#stressplot(birdord)#0.18 stress
birdord

birdfit <- envfit(birdord,env = birdenv[,2:6], display = 'sites', perm=9999, choices = c(1,2))
birdfit$vectors


#Create a dataframe with the results for displaying
sc <- as.data.frame(scores(vegord, display = "sites")) #save NMDS results into dataframe
sc <- cbind(sc, Management = birdenv$cat_final) #add grouping variable "Management" to dataframe
sc <- cbind(sc, Site = rownames(sc)) #add site names as variable if you want to display on plot
head(sc)

arrows <- as.data.frame(birdfit$vectors$arrows*sqrt(birdfit$vectors$r)/2)
#arrows
spp.scrs <- cbind(arrows, Species = rownames(arrows))
spp.scrs <- cbind(spp.scrs, pval=birdfit$vectors$pvals)
sig.scrs <- spp.scrs[spp.scrs$pval < 0.05,]
sig.scrs$Species <- c("NDVI","Image texture","A.donax","Tree")
sig.scrs
sc
#write.csv(sig.scrs,file = "Bird_ord_arrows.csv")
#revalue to match Lars' figures
sc$Management <- revalue(sc$Management, c("old_treated"="Old Burn", "new_treated" = "New Burn",
                                          "untreated_cane"="Unburned"))

#write.csv(sc,"Bird_veg_ord_axes_values.csv")

#plot

bird.mds <- ggplot(sc, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, colour = factor(Management)), size = 4) +   
  geom_segment(data = sig.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),arrow=arrow(length = unit(0.25, "cm"))) +
  geom_text(x=0.2,y=0.17,label = "(A)", size=14) +
  geom_text_repel(data = sig.scrs, aes(x = NMDS1, y = NMDS2, 
                  label = str_wrap(Species, width = 10)),nudge_x = -0.025,nudge_y = 0.015,segment.color = 'transparent') + #Gives space between labels
  scale_colour_manual(values=c("Unburned" = "#009E73", "Old Burn" = "#E69F00","New Burn" = "#D55E00" )) +
  xlab("Axis 1") + 
  ylab("Axis 2") +
  theme_bw()+
  theme(text = element_text(size=18),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Restoration Group")
bird.mds
#ggsave(plot = bird.mds,filename = "birdord_test.tiff",device = "tiff",dpi = 500)
#Lars' suggestion:
#tiff(filename, width = 1.9, height = 4, units='in', res=500, compression="lzw")
#dev.off()


#vegetation plot
vord <- metaMDS(comm = vd,k = 3,trymax = 200,autotransform = F, wascores = T)
stressplot(vord)#stress=0.03974
vc <- as.data.frame(scores(vord, display = "sites")) #save NMDS results into dataframe
vc <- cbind(vc, Management = birdenv$cat_final) #add grouping variable "Management" to dataframe
head(vc)
scaleFUN <- function(x) sprintf("%.1f", x)

vbird.mds<- ggplot(vc, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, colour = factor(Management)), size = 4) +   
  scale_colour_manual(values=c("Unburned" = "#009E73", "Old Burn" = "#E69F00","New Burn" = "#D55E00")) +
  geom_text(x=0.54,y=0.51,label = "(B)", size=14) +
  scale_y_continuous(labels=scaleFUN)+
  xlab("Axis 1") + 
  ylab("Axis 2") +
  theme_bw()+
  theme(text = element_text(size=18),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        line = element_line(linetype = "solid"), 
        panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Restoration Group")
vbird.mds

#arrange!
TADA <- gridExtra::grid.arrange(bird.mds, vbird.mds, bt.mds, btveg.mds, ncol = 2, nrow=2,
                                widths=c(3, 3))       
TADA
table(butt$Management)
#ggsave(file="NMDSplots.pdf", TADA,device = "pdf",width = 10,height = 10,units = "in",dpi = 350)
