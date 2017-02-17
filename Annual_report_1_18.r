
#Cuckoo Presence vs vegetation cover
bird.veg<-read.csv("bird_veg_data.csv", header = T, sep = ",", stringsAsFactors=FALSE)
#nrow(bird_comp)
YBCU<-read.csv("YBCU_Presence.csv", header=T, sep=",", stringsAsFactors=FALSE)[,c('point','presence')]
#head(YBCU, 5)
#merge data based on points, finds common attribute, excludes non matches
bird.comp<-merge(bird.veg,YBCU)
#nrow(bird_comp)

treatment<-read.csv("sites_by_treatment.csv", header = T, sep = ",", stringsAsFactors=FALSE)
treatment<-treatment[, colSums(is.na(treatment)) != nrow(treatment)] #remove any na values
nrow(treatment)
#head(treatment)

recent.treatment<-read.csv("2016_treated.csv", header = T, sep = ",", stringsAsFactors=TRUE)
recent.treatment<-recent.treatment[, colSums(is.na(recent.treatment)) != nrow(recent.treatment)]
#head(recent.treatment)

bosque.area<-read.csv("Bosque_area_edited.csv",header=T, sep=",", stringsAsFactors=FALSE)
bosque.rank<-bosque.area[c("bosque","area_ranking")]
#bosque_area

bird.data<-read.csv("bird_data.csv", header = T, sep = ",", stringsAsFactors=FALSE)
bird.data<-bird.data[, colSums(is.na(bird.data)) != nrow(bird.data)] #remove NA values
#head(bird_data)

bird.filtered <- bird.data[(!grepl("FL", bird.data$distance) & bird.data$distance != 110) & (!grepl("XX", bird_data$species)),]
#head(bird.filtered)

#bosque.richness<-tapply(bosque.rich$species,INDEX=bosque.rich$area_ranking, FUN=mean)
#write.table(bosque_richness, file = "bird_richness_by_bosque.csv", sep = " ",row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
#ranked<-read.csv("bird_richness_by_bosque.csv", sep=",")
#head(ranked)

bird.abundance<-aggregate(species ~ point, bird.filtered, function(x) length(x))
nrow(bird_abundance)

bird.richness<-aggregate(species ~ point, bird.filtered, function(x) length(unique(x)))
bird.richness <- bird.richness[(!grepl("CW3", bird.richness$point)) & (!grepl("CW4", bird.richness$point)) &
                               (!grepl("TS5", bird.richness$point)),]
names(bird.richness)[2] <- "richness"
#tail(bird.richness)

abund.merge<-merge(treatment, bird.abundance)
names(abund.merge)[3] <- "total_abundance"
abund.merge$mean_abundance <- (abund.merge$total_abundance / 3)
nrow(abund.merge)

bird.merge <- cbind(abund.merge,bird.richness$richness)
names(bird.merge)[5] <- "richness"
#head(bird.merge)

bird.all<-merge(bird.merge,bird.veg)
bird.all[1:5,15:25]

YBCU.merge <- merge(treatment, YBCU)
nrow(YBCU.merge)

mesquite <- aggregate(mature_mesquite ~ bosque, bird.all, FUN = mean)
#mesquite
boxplot(mature_mesquite~bosque,data=bird.all, las=2, ylab="mean percent cover, mature mesquite")

YBCU.data <- read.csv("YBCU_by_bosque.csv", header = T, sep = ",", stringsAsFactors=FALSE)
#bird.data<-bird.data[, colSums(is.na(bird.data)) != nrow(bird.data)] #remove NA values
#YBCU.data

#group and average by bosque
cacwrich<-bird.merge[grepl("CA",bird.merge$bosque)| grepl("CW", bird.merge$bosque),]
#cacwrich
#tapply(bird_merge$species,INDEX=bird_merge$bosque, FUN=function(x) mean(unique(x)))
#mean(cacwrich$species)
    
#group mean species richness by
#tapply(bird.merge$species,INDEX=bird.merge$burn, FUN=mean)
#by herbacious
#mesquite.table<-tapply(bird.merge$species,INDEX=bird.merge$mature_mesquite, FUN=mean)
#mesquite.table
#barplot(mesquite.table, main="Mesquite cover")

#plot salt cedar by YBCU
#png('scedar_YBCU.png')
library(repr)
options(repr.plot.width=5, repr.plot.height=5)
#boxplot(salt_cedar_total~presence,data=bird_comp,main="Cuckoo detections and tamarix",xlab="Cuckoo presence",ylab="Tamarix percent cover")
t.test(salt_cedar_total~presence, data=bird.comp)
#dev.off()

#boxplot arundo vs YBCU presence
#png('arundo_YBCU.png')
#boxplot(arundo_donax~presence,data=bird_comp,main="Cuckoo detections and Arundo donax",xlab="Cuckoo presence",ylab="Arundo donax percent cover")
#dev.off()
#t.test(arundo_donax~presence, data=bird.comp)

#png("tree cover v YBCU.png")
#boxplot(tree~presence,data=bird.comp,main="Cuckoo detections and shrub cover",xlab="Cuckoo presence",ylab="shrub percent cover")
#t.test(tree ~presence, data=bird.comp)
#dev.off()

#mesquite vs YBCU
#png('mesquite_YBCU.png') #create png of plot
#boxplot(mature_mesquite~presence,data=bird.comp,col=(c("white","lightblue")),main="Cuckoo detections and mesquite",xlab="Cuckoo presence",ylab="Mesquite canopy cover")
#dev.off()
t.test(mature_mesquite~presence, data=bird.comp)

#boxplot(species~treatment,data=recent.treated.data,col="white",main="Bird richness based on burn treatment",xlab="Burn Treatment",ylab="Bird Richness")

#burn vs richness
#png('Bird_richness_vs_burn.png')
boxplot(species~treatment,data=bird.merge,col=c("lightblue","white"),main="Bird richness based on burn treatment",xlab="Burn Treatment",ylab="Bird Richness")
#dev.off()
t.test(bird.merge$species~bird.merge$treatment)

#richness vs invasives
#png('Bird_richness v invasives.png')
#plot(species~invasives,data=bird.merge,main="Bird richness and invasive species",xlab="Invasive Species Percent Cover",ylab="Bird Richness")
#scatter.smooth(bird.merge$invasives,bird.merge$species)
#abline(lm(bird.merge$species~bird.merge$invasives),col="red")
#dev.off()
#note: "cannot compute exact p-value message bc some values are identical, resulting in ties when ranked
cor.test(bird.merge$species,bird.merge$invasives, method = "spearman")

#richness vs arundo
#png('Bird_richness v arundo.png')
#plot(species~arundo_donax,data=bird.merge,main="Bird richness and Arundo donax",xlab="Arundo donax Percent Cover",ylab="Bird Richness")
#abline(lm(bird_merge$species~bird_merge$arundo_donax),col="red")
#dev.off()
cor.test(bird.merge$species,bird.merge$arundo_donax, method = "spearman")

#histogram to check for normality

#hist(bird.merge$species, main="Histogram of Bird Richness", xlab="Bird Richness")
#plot residuals
#qqnorm(bird.merge$species)

#tamarix vs bird_richness
#png("bird_richness_tamarix.png")
plot(species~salt_cedar_total,data=bird.merge,main="Bird richness and Tamarix",xlab="Tamarix Percent Cover",ylab="Bird Richness")
#abline(lm(species~salt_cedar_total, data=bird.merge),col="red")
#dev.off()
cor.test(bird_merge$species,bird.merge$salt_cedar_total, method = "spearman")

#png("bird_mean_abundance_v_treatment.png")
boxplot(new~treatment,data=abund.mean,col= c("lightblue","white"),main="Bird abundance based on burn treatment",xlab="Burn Treatment",ylab="Bird Abundance")
#dev.off()
t.test(abund.mean$new~abund.mean$treatment)

#bird abundance v treament
#png("bird_abundance_v_treatment.png")
boxplot(species~treatment,data=abund.merge,col=c("lightblue", "white"),main="Bird abundance based on burn treatment",xlab="Burn Treatment",ylab="Bird Abundance")
#dev.off()
t.test(abund.merge$species~abund.merge$treatment)
