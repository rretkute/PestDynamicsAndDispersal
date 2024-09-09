library(ggplot2)
library(ggpubr)
require(gridExtra)
library (deSolve) 
library(zoo)
library(scales)
library(geosphere)
sf_use_s2(FALSE)
source('utilities.R')


################################
#  Hopper/bands presence data
################################

date1<-as.Date("2010-01-01") 

lcst.prs<-rbind(read.csv("Data/Hoppers.csv", header=T, stringsAsFactors=F),
          read.csv("Data/Bands.csv", header=T, stringsAsFactors=F))
lcst.prs$date<-as.Date(lcst.prs$STARTDATE)

lcst.prs<-lcst.prs[lcst.prs$COUNTRYID %in% c("ER", "DJ", "ET", "SO", "KE") &
                     lcst.prs$date>date1, ]
lcst.prs<-lcst.prs[, c(1,2,153)]
colnames(lcst.prs)<-c("lng", "lat", "date")
nrow(lcst.prs)

fig+
  geom_point(data=lcst.prs,
             aes(x=lng, y=lat), size=0.01, col="red")

################################
#  Hopper/bands absence data
################################

lcst<-rbind(read.csv("Data/Hoppers.csv", header=T, stringsAsFactors=F),
            read.csv("Data/Adults.csv", header=T, stringsAsFactors=F),
            read.csv("Data/Bands.csv", header=T, stringsAsFactors=F),
            read.csv("Data/Swarms.csv", header=T, stringsAsFactors=F))
lcst<-lcst[lcst$COUNTRYID %in% c("ER", "DJ", "ET", "SO", "KE"),]

yy<-read.csv("Data/Ecology.csv", header=T, stringsAsFactors=F)
yy$date<-as.Date(yy$STARTDATE)
yy<-yy[yy$COUNTRYID %in% c("ER", "DJ", "ET", "SO", "KE") &
         yy$date>date3,]
yy<-yy[, c(1,2, 19)]
colnames(yy)<-c("lng", "lat", "date")

ids<-c()
for(i in 1:nrow(yy)){
  dd<-distm (yy[i, 1:2], lcst[,1:2], fun = distHaversine)
  tmp<- min(dd)/1000
  if(tmp>1) ids<-c(ids, i)
}
lcst.abs<-yy[ids,]
nrow(lcst.abs)

fig+
  geom_point(data=lcst.abs,
             aes(x=lng, y=lat), size=0.01, col="blue")


################################
#  Combine
################################

lcst.prs$presence<-1
lcst.abs$presence<-0
ans<-rbind(lcst.prs[, c(1,2,4)], lcst.abs[, c(1,2,4)])

fig+
  geom_point(data=ans,
             aes(x=lng, y=lat, col=as.factor(presence)), size=0.01) +
  labs(col="Presence")

################################
#  Divide into training/validation sets
################################

set.seed(100)
indTrain<-sort(sample(1:nrow(ans), 0.8*round(nrow(ans))))
indVal<-setdiff(1:nrow(ans), indTrain)

dataTrain<-ans[indTrain, ]
dataVal<-ans[indVal, ]

fig+
  geom_point(data=dataTrain,
             aes(x=lng, y=lat, col=as.factor(presence)), size=0.01) +
  labs(col="Presence") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"))

fig+
  geom_point(data=dataVal,
             aes(x=lng, y=lat, col=as.factor(presence)), size=0.01) +
  labs(col="Presence") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"))

# write.csv(dataTrain, "Data/locust_data_training.csv", row.names = F)
# write.csv(dataVal, "Data/locust_data_validation.csv", row.names = F)


