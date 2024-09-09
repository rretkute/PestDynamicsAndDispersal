library(ggplot2)
library(grid)
library(raster)
library(matlib)
library(sfsmisc)
library(Rfast)
library(geostats)
library(AMISEpi) # https://github.com/rretkute/AMISEpi
library(matrixStats)
library(rnaturalearth)
library(rnaturalearthdata)
source('utilities.R')

src<-read.csv("data/sourcelist_3860.csv", header = FALSE)
colnames(src)<-c("ID", "lng", "lat","Country")

####################
#  Plot locations
####################

loc<-data.frame(ID=c(3588, 3247, 2575), 
                lng=c(50.13281, 47.03906, 43.10156), 
                lat=c(10.546875,  6.046875,  3.046875),
                site=c("i", "ii", "iii"))
fig+
  geom_point(data=loc, aes(x=lng, y=lat), pch=21, fill='gray', size=6)+
  geom_text(data=loc, aes(x=lng, y=lat, label =site))


####################
#  Angles for all days
####################

wn.dir.d<-data.frame(date =c(), ind=c())
sim.angl<-matrix(NA, ncol=1000, nrow=0)
datesAll<-seq(as.Date("2020-01-01"), as.Date("2021-12-31"), by=1)
for(d in 1:length(datesAll)){
  cat("\n")
  date<-as.Date(datesAll[d])
  print(date)
  for(ii in 1:3){
    cat(c( ii, ""))
    aa<-get_direction_all(loc$ID[ii], date)
    wn.dir.d<-rbind(wn.dir.d, data.frame(date =date, ind=loc$ID[ii]))
    sim.angl<-rbind(sim.angl, aa$angle)
  }
}

####################
#  Fit von Mises distribution
####################


nn<-1000 #Number of trajectories
n.dat<-nrow(sim.angl) # Number of parameter sets to sample each iteration
nn.parsets<-10^3
# Max number of iterations
NN<-rep(nn.parsets,100)

#  Set up for mixture function as a multivariate Students distribution
proposal=mvtComp(df=3); mixture=mclustMix();
dprop <- proposal$d
rprop <- proposal$r

ESS.R<-100
#  log.10.beta, alpha
# Prior in Uniform[min, max]
n.param<-2
param.lim<-data.frame(min=c(0, 0), max=c(360, 700))

# Prior PDF value 
dprop0<-function(pp, param.lim){
  ans<-dunif(pp[1], min=param.lim$min[1], max=param.lim$max[1], log=FALSE)
  if(length(pp)>1){
    for(i in 2:length(pp)){
      ans<-ans*dunif(pp[i], min=param.lim$min[i], max=param.lim$max[i], log=FALSE)
    }
  }
  return(ans)
}

# Sampling from prior
rprop0<-function(param.lim){
  ans<-sapply(1:nrow(param.lim), function(a) runif(1, min=param.lim$min[a],
                                                   max=param.lim$max[a]))
  return(ans)
}

prior.dns<-log(1/param.lim$max[1])+log(1/param.lim$max[2])

Sigma <- list(NA)
Mean<-list(NA)
PP<-list(NA)
GG<-list(NA)
ess<-data.frame(it=c(), Sets=c(), ESS=c())

# Iteration 1
set.seed(1)
it<-1
cat(c("Started iteration ", it,"\n"))

param<-data.frame(m=c(), k=c())
DENS<-matrix(NA, nrow=n.dat, ncol=0)

for(ii in 1:NN[it]){
  cat(c(ii, ""))
  new.param<-rprop0(param.lim)
  m<-new.param[1]; k<-new.param[2]
  param<-rbind(param, data.frame(m=m, k=k))
  dens<-geostats::vonMises(unlist(as.list(t(sim.angl*2*pi/360))),
                           m*2*pi/360, k)
  dens<-matrix(dens, ncol=nn, byrow=TRUE)
  DENS<-cbind(DENS, rowSums(log(dens)))
}

WW<-0*DENS
for(ii in 1:n.dat){
  wl<- DENS[ii,]-prior.dns
  CC<-logSumExp(wl)
  ww<-exp(wl-CC)
  ww<-ww/sum(ww)
  WW[ii, ]<-ww
  ess<-rbind(ess, data.frame(it=it, Sets=paste("set ", ii), ESS=(sum((ww)^2))^(-1)))
}
hist(ess$ESS,10)
lst.ess<-ess$ESS

# Iterations 2+
stop<-0
while(stop==0){
  it<-it+1
  cat(c("\n Started iteration ", it,"\n"))
  wh.dt<-which(lst.ess<ESS.R)
  if(length(wh.dt)>1){
    J<-unique(sample(1:sum(NN[1:(it-1)]), 1000, prob=colMeans(WW[wh.dt,]), replace=TRUE))
  } else {
    J<-unique(sample(1:sum(NN[1:(it-1)]), 1000, prob=WW[wh.dt,], replace=TRUE))
  }
  xx<-as.matrix(param[J,1:n.param])
  clustMix <- mixture(xx)
  G <- clustMix$G
  cluster <- clustMix$cluster
  ### Components of the mixture
  ppt <- clustMix$alpha
  muHatt <- clustMix$muHat
  varHatt <- clustMix$SigmaHat
  GG[[it-1]]<-G
  G1<-0; G2<-G
  if(it>2) {
    G1<-sum(sapply(1:(it-2), function(a) GG[[a]]))
    G2<-sum(sapply(1:(it-1), function(a) GG[[a]]))
  }
  for(i in 1:G){
    Sigma[[i+G1]] <- varHatt[,,i]
    Mean[[i+G1]] <- muHatt[i,]
    PP[[i+G1]]<-ppt[i]
  }
  
  # Draw new parameters and calculate log likelihood
  while(nrow(param)<sum(NN[1:it])){
    compo <- sample(1:G,1,prob=ppt)
    x1 <- t(rprop(1,muHatt[compo,], varHatt[,,compo]))
    new.param<-as.numeric(x1)
    if(dprop0(new.param, param.lim)>0){
      m<-new.param[1]; k<-new.param[2]
      dens<-geostats::vonMises(unlist(as.list(t(sim.angl*2*pi/360))),
                               m*2*pi/360, k)
      dens<-matrix(dens, ncol=nn, byrow=TRUE)
      LL<-rowSums(log(dens))
      if(min(LL)>-Inf){
        DENS<-cbind(DENS, LL)
        param<-rbind(param, data.frame(m=m, k=k))
      }
    }
    cat(c(nrow(param), ""))
  }
  
  q <- (NN[1]/sum(NN[1:it]))*exp(prior.dns) +
    (sum(NN[2:it])/sum(NN[1:it]))* rowSums(as.matrix(sapply(1:G2, 
                                                            function(a) PP[[a]] * dprop(param[,1:n.param],mu= Mean[[a]], Sig=Sigma[[a]]))))
  
  #Calculate ESS and weighst
  WW<-0*DENS
  for(ii in 1:n.dat){
    wl<- DENS[ii,]-log(q)
    CC<-logSumExp(wl)
    ww<-exp(wl-CC)
    ww<-ww/sum(ww)
    WW[ii, ]<-ww
    ess<-rbind(ess, data.frame(it=it, Sets=paste("set ", ii), ESS=(sum((ww)^2))^(-1)))
  }
  lst.ess<-ess$ESS[ess$it==it]
  print(range(lst.ess))
  if(min(lst.ess)>=ESS.R) stop<-1
  if(it>= length(NN)) stop<-1
}
hist(lst.ess, 100)
ggplot(ess, aes(x=it, y=ESS))+
  geom_point()+
  geom_hline(yintercept=ESS.R, linetype="dashed", color = "black")+
  xlab("Iteration")+
  theme_bw()

MAP<-apply(DENS, 1, function(x) max(x, na.rm = TRUE))
dir.mp<-wn.dir.d
colnames(dir.mp)<-c("coord", "date")
dir.mp$mu<-0
dir.mp$kappa<-0
for(i in 1:nrow(dir.mp)){
  wh2<-which(DENS[i,]==MAP[i])
  dir.mp$mu[i]<-param$m[wh2]
  dir.mp$kappa[i]<-param$k[wh2]
}
dir.mp$mu1<-dir.mp$mu + 270
dir.mp$mu1<-dir.mp$mu1 %% 360

####################
#  Plot MAP
####################

ii<-1
xx<-dir.mp[which(dir.mp$coord==ii & dir.mp$date<366),]
xx$kappa[xx$kappa>=300]<-300
figC<-ggplot() + 
  geom_line(data=data.frame(date=c(1, 1), y=c(-5, 5)), aes(x=date, y=y), col='white') +
  theme_bw() + xlim(1,365) +ylim(0, 2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
figC

pal<-colorRampPalette(c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3','#e41a1c'))
for(dd in 1:365){
  cl<-pal(365)[ceiling(xx$mu1[which(xx$date==dd)])]
  figC<-figC+ 
    geom_line(data=data.frame(date=c(dd, dd), y=c(1.5,2)), aes(x=date, y=y),
              col=cl) 
}

pal<-colorRampPalette(c('#f1eef6', '#bdc9e1', '#74a9cf', '#2b8cbe',
                        '#045a8d'))
for(dd in 1:365){
  cl<-pal(300)[ceiling(xx$kappa[which(xx$date==dd)])]
  figC<-figC+  
    geom_line(data=data.frame(date=c(dd, dd), y=c(1, 1.5)), aes(x=date, y=y),
              col=cl)
}

xx<-dir.mp[which(dir.mp$coord==ii & dir.mp$date>=366),]
xx$kappa[xx$kappa>=300]<-300
xx$date<-xx$date-365

pal<-colorRampPalette(c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3','#e41a1c'))
for(dd in 1:365){
  cl<-pal(365)[ceiling(xx$mu1[which(xx$date==dd)])]
  figC<-figC+  
    geom_line(data=data.frame(date=c(dd, dd), y=c(0.5, 1)), aes(x=date, y=y),
              col=cl) 
}

pal<-colorRampPalette(c('#f1eef6', '#bdc9e1', '#74a9cf', '#2b8cbe',
                        '#045a8d'))
for(dd in 1:365){
  cl<-pal(300)[ceiling(xx$kappa[which(xx$date==dd)])]
  figC<-figC+  
    geom_line(data=data.frame(date=c(dd, dd), y=c(0, 0.5)), aes(x=date, y=y),
              col=cl) 
}
figC


ii<-2
xx<-dir.mp[which(dir.mp$coord==ii & dir.mp$date<366),]
xx$kappa[xx$kappa>=300]<-300
figC<-ggplot() + 
  geom_line(data=data.frame(date=c(1, 1), y=c(-5, 5)), aes(x=date, y=y), col='white') +
  theme_bw() + xlim(1,365) +ylim(0, 2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
figC

pal<-colorRampPalette(c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3','#e41a1c'))
for(dd in 1:365){
  cl<-pal(365)[ceiling(xx$mu1[which(xx$date==dd)])]
  figC<-figC+ 
    geom_line(data=data.frame(date=c(dd, dd), y=c(1.5,2)), aes(x=date, y=y),
              col=cl) 
}

pal<-colorRampPalette(c('#f1eef6', '#bdc9e1', '#74a9cf', '#2b8cbe',
                        '#045a8d'))
for(dd in 1:365){
  cl<-pal(300)[ceiling(xx$kappa[which(xx$date==dd)])]
  figC<-figC+  
    geom_line(data=data.frame(date=c(dd, dd), y=c(1, 1.5)), aes(x=date, y=y),
              col=cl)
}

xx<-dir.mp[which(dir.mp$coord==ii & dir.mp$date>=366),]
xx$kappa[xx$kappa>=300]<-300
xx$date<-xx$date-365

pal<-colorRampPalette(c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3','#e41a1c'))
for(dd in 1:365){
  cl<-pal(365)[ceiling(xx$mu1[which(xx$date==dd)])]
  figC<-figC+  
    geom_line(data=data.frame(date=c(dd, dd), y=c(0.5, 1)), aes(x=date, y=y),
              col=cl) 
}

pal<-colorRampPalette(c('#f1eef6', '#bdc9e1', '#74a9cf', '#2b8cbe',
                        '#045a8d'))
for(dd in 1:365){
  cl<-pal(300)[ceiling(xx$kappa[which(xx$date==dd)])]
  figC<-figC+  
    geom_line(data=data.frame(date=c(dd, dd), y=c(0, 0.5)), aes(x=date, y=y),
              col=cl) 
}
figC


ii<-3
xx<-dir.mp[which(dir.mp$coord==ii & dir.mp$date<366),]
xx$kappa[xx$kappa>=300]<-300
figC<-ggplot() + 
  geom_line(data=data.frame(date=c(1, 1), y=c(-5, 5)), aes(x=date, y=y), col='white') +
  theme_bw() + xlim(1,365) +ylim(0, 2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
figC

pal<-colorRampPalette(c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3','#e41a1c'))
for(dd in 1:365){
  cl<-pal(365)[ceiling(xx$mu1[which(xx$date==dd)])]
  figC<-figC+ 
    geom_line(data=data.frame(date=c(dd, dd), y=c(1.5,2)), aes(x=date, y=y),
              col=cl) 
}

pal<-colorRampPalette(c('#f1eef6', '#bdc9e1', '#74a9cf', '#2b8cbe',
                        '#045a8d'))
for(dd in 1:365){
  cl<-pal(300)[ceiling(xx$kappa[which(xx$date==dd)])]
  figC<-figC+  
    geom_line(data=data.frame(date=c(dd, dd), y=c(1, 1.5)), aes(x=date, y=y),
              col=cl)
}

xx<-dir.mp[which(dir.mp$coord==ii & dir.mp$date>=366),]
xx$kappa[xx$kappa>=300]<-300
xx$date<-xx$date-365

pal<-colorRampPalette(c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3','#e41a1c'))
for(dd in 1:365){
  cl<-pal(365)[ceiling(xx$mu1[which(xx$date==dd)])]
  figC<-figC+  
    geom_line(data=data.frame(date=c(dd, dd), y=c(0.5, 1)), aes(x=date, y=y),
              col=cl) 
}

pal<-colorRampPalette(c('#f1eef6', '#bdc9e1', '#74a9cf', '#2b8cbe',
                        '#045a8d'))
for(dd in 1:365){
  cl<-pal(300)[ceiling(xx$kappa[which(xx$date==dd)])]
  figC<-figC+  
    geom_line(data=data.frame(date=c(dd, dd), y=c(0, 0.5)), aes(x=date, y=y),
              col=cl) 
}
figC






