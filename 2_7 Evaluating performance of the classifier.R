library(ggplot2)
library(raster)
library(PresenceAbsence)

dataVal<-read.csv("Data/locust_data_validation.csv", header = TRUE)
fln<-"Data/Breeding_suitability_map_1000m.tif"
sp <- SpatialPoints(dataVal)
vg <- raster(x=fln)
rasValue <- raster::extract(vg, sp)

zz<-dataVal
zz$pred<-rasValue
zz<-zz[!is.na(zz$pred),]

ggplot(zz, aes(x=presence, y=pred))+
  geom_boxplot(aes(group=presence))

DATA <- data.frame(ID=1:length(zz$presence),
                   OBSERVED=zz$presence,
                   RF=zz$pred/100)
auc.roc.plot(DATA, color = TRUE, na.rm =TRUE, 
             add.legend = FALSE, main = "")

# presence.absence.hist(DATA, na.rm =TRUE, legend.cex = 1,  N.bars = 20) 
pred.prev <- predicted.prevalence(DATA = DATA, na.rm =TRUE, 
                                  threshold = 101)
pred.prev[, 2] <- round(pred.prev[, 2], digits = 2)

accu <- presence.absence.accuracy(DATA, na.rm =TRUE)
accu[, -c(1, 2)] <- signif(accu[, -c(1, 2)], digits = 2)
accu
