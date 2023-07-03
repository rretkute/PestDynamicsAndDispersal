# PestDynamicsAndDispersal
A framework for modelling migratory pest population dynamics and large-scale dispersal. 
The code provided here uses synthetic data, therefore results should be used for educational or illustrative purpose, and not for guiding survellance or management.

#### Installation

```r
library(devtools)
install.packages("ggplot2")
install.packages("raster")
install.packages("gridExtra")
install.packages("ggnewscale")
install.packages("ggspatial")
install.packages("viridis")
install_github("rretkute/PestDynamicsAndDispersal")
``` 

#### Example use

```r
library(ggplot2)
library(raster)
library(gridExtra)
library(ggnewscale)
library(ggspatial)
library(viridis)
library(PestDynamicsAndDispersal)

###################
#      Functions 
###################
?suit_breeding
?egg_development
?sample_dev_period
?hopper_dev
?get_NDVI_values
?get_land_cover
?smooth_Whittaker
?fit_line
?swrm_stay_type
?swrm_stay
?get_NDVI_trend
?get_NDVI_dnst
?get_wind

###################
#       Data
###################

data(fig)  # Empty map boundaries
data(src) #  Grid for Met Office data
data(sbr)  # Breeding suitability map
data(lng.stay) # Feeding capacity
data(st.dur) # Length of stay on 2020/08/28
data(pth1)  # Simulated breeding, development and migration of DL.
data(st.dur1) # Simulated breeding, development and migration of DL.
data(snth.pth) #  Illustrate forecasting
data(lnd.pnts) #  Illustrate forecasting
data(Obs) #  Illustrate forecasting
data(LandCovercols) #  Illustrate forecasting
data(TRJ1) #  Illustrate forecasting

###################
#      Figures 
###################

# Breeding suitability
pal<-c('#2891C9', '#A0C29B', '#FAFA64', '#FB8C32', '#E80F14')
fig+
  geom_tile(data=sbr, aes(x=lng, y=lat, fill=rng)) +
  xlab("")+ ylab("") +
  scale_fill_manual(name="Breeding suitability",
                    values = pal)+
  theme(legend.position="bottom")

# Feeding capacity 
ggplot(lng.stay[!is.na(lng.stay$stay),], aes(x=tr, y=vl)) +
  geom_tile(aes(fill=stay)) +
  facet_grid(.~lc) +
  scale_fill_manual(name = "Duration",
                    values=c("#7570b3", "#1b9e77", "#d95f02"))+
  xlab("NDVI trend")+ ylab("NDVI value")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Length of stay on 2020/08/28
fig+
  geom_tile(data=st.dur, aes(x=lng, y=lat, fill=Duration)) +
  xlab("")+ ylab("") +
  scale_fill_manual(name = "Duration",
                    values=c("#7570b3", "#1b9e77", "#d95f02", "gray"))+
  theme(legend.position="none") 
  
# Simulated breeding, development and migration of DL.
fig+
  geom_point(data=pth1, 
             aes(x=lng, y=lat, group=ind, col=LandType),size=0.05) +
  scale_color_manual(values = LandCovercols) +
  new_scale_colour()+
  geom_point(data=st.dur1, 
             aes(x=lng, y=lat, col=Duration),
             pch=16, size=1.75, stroke=0.05)+
  scale_colour_manual(name = "Duration",
                      values=c("#7570b3", "#1b9e77", "#d95f02", "gray"))+
  theme(legend.position = "none") + xlab("") +ylab("")

#  Illustrate forecasting
fig1<-fig+
  geom_point(data=snth.pth, size=0.5,
             aes(x=lng, y=lat, col=LandCover)) +
  scale_color_manual(values = LandCovercols) +
  geom_point(data=lnd.pnts, aes(x=lng, y=lat, fill=LandType),
             pch=21, size=4) +
  scale_fill_manual(values = LandCovercols) +
  geom_text(data=lnd.pnts[1:12,],
            aes(x=lng, y=lat, label = ind), size=2.5)+
  geom_point(data=Obs, aes(x=lng, y=lat), pch=21,
             fill="red", size=1.5)+
  xlab("")+ ylab("") +
  guides(col="none", fill="none")

fig2<-fig+
  annotation_scale(location = "br", width_hint = 0.5)+
  geom_path(data=TRJ1[TRJ1$dd==1,], size=0.25,
            aes(x=lng, y=lat, group=id, col=flt.tm)) +
  scale_color_viridis(option = "D", direction=-1)+
  geom_path(data=snth.pth[34:51,][1:7,], size=1,
            aes(x=lng, y=lat), col='black') +
  geom_point(data=Obs[1,], aes(x=lng, y=lat), pch=21,
             fill="red", size=4)+
  geom_point(data=lnd.pnts[3,], aes(x=lng, y=lat, fill=LandCover),
             pch=21, size=8) +
  scale_fill_manual(values = LandCovercols) +
  geom_text(data=lnd.pnts[3,],
            aes(x=lng, y=lat, label = ind), size=6)+
  xlab("")+ ylab("") + ggtitle("Day of reporting")+
  guides(col="none", fill="none")  +
  coord_sf(xlim = c(45, 50),  ylim = c(6, 11), expand = FALSE)

fig3<-fig+
  annotation_scale(location = "br", width_hint = 0.5)+
  geom_path(data=TRJ1[TRJ1$dd<=2,], size=0.25,
            aes(x=lng, y=lat, group=id), col="gray") +
  geom_path(data=TRJ1[TRJ1$dd==3,], size=0.25,
            aes(x=lng, y=lat, group=id, col=flt.tm)) +
  scale_color_viridis(option = "D", direction=-1)+
  geom_path(data=snth.pth[34:51,][1:7,], size=1,
            aes(x=lng, y=lat), col='black') +
  geom_point(data=Obs[1,], aes(x=lng, y=lat), pch=21,
             fill="red", size=4)+
  geom_point(data=lnd.pnts[3,], aes(x=lng, y=lat, fill=LandCover),
             pch=21, size=8) +
  scale_fill_manual(values = LandCovercols) +
  geom_text(data=lnd.pnts[3,],
            aes(x=lng, y=lat, label = ind), size=6)+
  xlab("")+ ylab("") + ggtitle("After 3 days")+
  guides(col="none", fill="none")  +
  coord_sf(xlim = c(42, 50),  ylim = c(3, 11), expand = FALSE)

fig4<-fig+
  annotation_scale(location = "br", width_hint = 0.5)+
  geom_path(data=TRJ1, size=0.25,
            aes(x=lng, y=lat, group=id), col="gray")+
  geom_path(data=TRJ1[TRJ1$dd==6,], size=0.25,
            aes(x=lng, y=lat, group=id, col=flt.tm)) +
  scale_color_viridis(option = "D", direction=-1)+
  geom_path(data=snth.pth[34:51,][1:7,], size=1,
            aes(x=lng, y=lat), col='black') +
  geom_point(data=Obs[1,], aes(x=lng, y=lat), pch=21,
             fill="red", size=4)+
  geom_point(data=lnd.pnts[3,], aes(x=lng, y=lat, fill=LandCover),
             pch=21, size=8) +
  scale_fill_manual(values = LandCovercols) +
  geom_text(data=lnd.pnts[3,],
            aes(x=lng, y=lat, label = ind), size=6)+
  xlab("")+ ylab("") + ggtitle("After 7 days")+
  guides(col="none", fill="none")  +
  coord_sf(xlim = c(41, 50),  ylim = c(2, 11), expand = FALSE)

grid.arrange(fig1, fig2, fig3, fig4, 
             ncol = 2, nrow = 2)
``` 

