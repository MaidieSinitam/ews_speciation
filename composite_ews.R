
### Code adapted from Clements and Ozgul (2016)


rm(list=ls())
require(ggplot2)
require(gridExtra)
require(reshape2)
require(earlywarnings)
require(zoo)
require(parallel)
require(data.table)
require(pBrackets)
require(grid)
##########################################################################################
##function to calculate CV
CV <- function(x, na.rm){
  ave<-mean(x, na.rm=na.rm)
  dev<-sd(x, na.rm=na.rm)
  return((dev/ave))
}
##interpolate function
interp<-function(days, obs){
  if(length(na.omit(obs))<length(obs)){
    return(approx(days, obs, n = length(obs), xout=days, method = "linear")$y)}else{
      return(obs)}
}
##function for rolling mean
rolling_mean <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- mean(x[1:i], na.rm=T);
  }    
  return(result);
}
##function for rolling sd
rolling_sd <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- sd(x[1:i], na.rm=T);
  }    
  return(result);
}


loadAll_Clements <- function(metric){
  metrics = c("time", "me", "wreswimm", "wmaxwres", "fitall", "fitrand", "fst", "dxyall", "dxysel", "dxyneut", "ld", "nbloci", "nblocisel", "nblocineut")
  k = which(metrics == metric)
  col1 = c()
  col2 = c()
  col3=c()
  i=1
  for (n in abrupt){
    K = getMetrics(n, cut=TRUE)[,k]
    timeindex = seq(1,length(K), by=1)
    bw = floor(10 * length(K)/100)
    smYY <- ksmooth(timeindex, K, kernel = "normal", bandwidth = bw, range.x = range(timeindex), 
      x.points = timeindex)
    nsmY <- K - smYY$y
    col1 = append(col1, rep(n,each=length(K)))
    col2 = append(col2, timeindex)
    col3 = append(col3,nsmY)
    i=i+1
  }
  out = data.frame(col1,col2, col3)
  colnames(out)= c("Tube", "Day", "Count")
  return(out)
}

##########################################################################################
#tipping.points<-read.csv("Tipping points.csv")[,2:4]
tipping.points = data.frame(abrupt, specTimeRowA)
colnames(tipping.points) = c("Tube", "tipping.point")
##########################################################################################
## read in the experimental data
#dd_clem<-read.csv("Experimental_data_Clements.csv", sep=";", dec=",")


loadFST_C =loadAll_Clements("fst")
loadME_C =loadAll_Clements("me")
loadNLOCI_C =loadAll_Clements("nblocisel")
loadDXY_C =loadAll_Clements("dxysel")
loadMFIT_C =loadAll_Clements("fitrand")
loadLD_C =loadAll_Clements("ld")


dd <- loadLD_C
dd1<-dd
##add in the time of tipping points, where applicable
dd1$tipping.point<-NA
dd1$tipping.point<-tipping.points$tipping.point[match(dd1$Tube, tipping.points$Tube)]
##########################################################################################
##calculate the standardised change in each metric
split.dd<-split(dd1, list(dd1$Tube))
#x2<-split.dd[[25]]
results.ld<-do.call("rbind", mclapply(split.dd, function(x2){
  
  dat <- tail(x2,93)
  
  ##blank objects to save results
  RES<-NULL
  roll.sd<-NULL
  # roll.acf<-NULL
  # roll.body<-NULL
  roll.ar<-NULL
  roll.sk<-NULL
  roll.kurt<-NULL
  # roll.return.rate<-NULL
  # roll.density.ratio<-NULL
  
  ##looped to calculate the rolling change	
  for(i in 2:length(dat$Day)){
    #i=3
    ##subset the population of interest up until dat i
    dat.t<-subset(dat, Day<=unique(sort(dat$Day))[i])
  
    ## for ar1
    roll.ar[[i]]<-ar.ols(dat.t$Count, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)$ar[1]
    ar.t<-(ar.ols(dat.t$Count, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)$ar[1]-mean(roll.ar, na.rm=T))/sd(roll.ar, TRUE)
    ##for sd
    roll.sd[[i]]<-sd(dat.t$Count, na.rm=T)
    sd.t<-(sd(dat.t$Count, na.rm=T)-mean(roll.sk, na.rm=T))/sd(roll.sk, TRUE)
    ##for skewness
    roll.sk[[i]]<-skewness(dat.t$Count, na.rm=T)
    sk.t<-(skewness(dat.t$Count, na.rm=T)-mean(roll.sk, na.rm=T))/sd(roll.sk, TRUE)
    ##for kurtosis
    roll.kurt[[i]]<-kurtosis(dat.t$Count, na.rm=T)
    kurt.t<-(kurtosis(dat.t$Count, na.rm=T)-mean(roll.kurt, na.rm=T))/sd(roll.kurt, TRUE)
    # 
    ##save results
    RES[[i]]<-data.frame(dat.t[i,], "ar1"=ar.t, "sd"=sd.t, "sk"=sk.t, "kurt"=kurt.t)
  }
  rr<-do.call("rbind", RES)
  rr$norm.day<-(rr$Day-max(rr$Day))-1
  return(rr)	
}, mc.cores=6))
##change the results data frame from wide to long format
melt.res<-melt(results.ld, id=c(1:4, 9))
##a vector of the indicators present in the results data frame
indicators<-unique(melt.res$variable)
###########################################################################
##calculate compositie EWS for single, pairwise, triplicate, etc.combinations
##split the data by tube
split.res<-split(melt.res, list(melt.res$Tube))
#kk<-split.res[[2]]
long.run.roll.nloci<-do.call("rbind", mclapply(split.res, function(kk){
  ##save out results from loops here, using counter to specify the position in the list object finished.results
  finished.results<-NULL
  counter<-0
  
  ##loop for single predictors, double, etc
  ##calculate the number of indicators (ar1, cv, etc)
  no.pred<-length(unique(kk$variable))
  ##and their names
  indicators<-sort(unique(kk$variable))
  ##loop through from 1:the number of indicators to be included in the method
  for(h in 1:no.pred){
    #h=3
    ##made a data frame of the indicator names
    inds<-data.frame(indicators)
    ##if this method should include more than 1 indicator in it (i.e. if h>1)
    ## then loop through to the number of indicators to be included, adding a column to "inds" each time
    if(h!=1){for(o in 2:h){
      inds<-cbind(inds, indicators)}
      ##all possible combinations of inds
      combs.all<-expand.grid(inds)
      ##make a unique code for each composite indicator
      combs.all$code <- apply(combs.all[ ,1:length(combs.all)], 1, function(x) paste(sort(x), collapse = " + " ))
      ##and remove any multiple occurances of the same combination of leading indicators
      combs<-as.data.frame(do.call("rbind", mclapply(split(combs.all, list(combs.all$code)), function(j){
        j$code<-NULL
        return(j[1,])
      })))
      rownames(combs)<-NULL
    }else{combs<-inds}
    
    ##then loop through object combs
    result<-FALSE
    for(m in 1:length(combs[,1])){
      #m=20
      ##if there is only one predictor:
      if(h==1){
        ##get results for each leading indicator from the experimental data
        sub.kk<-kk[(kk$variable%in%unlist(combs[m,])),]
        ##create an id for the combined predictor
        sub.kk$id<-paste(unique(sub.kk$variable), collapse=" + ")
        ##calculate the rolling values
        roll<-colSums(do.call("rbind", mclapply(split(sub.kk, list(sub.kk$variable)), function(d){
          if(length(d[,1])>0){
            return(d$value)
          }})))
        ##make a data set of neccesary info to add the rolling value to			
        fin.kk<-sub.kk[1:(length(sub.kk[,1])/length(unique(sub.kk$variable))),c(1:3,5,6,8)]
        ##add the rolling value
        fin.kk$value<-roll
        ##calculate the sd and mean of the rolling value
        fin.kk$roll.sd<-rolling_sd(fin.kk$value)
        fin.kk$roll.mean<-rolling_mean(fin.kk$value)
        ##add in information on whether size is included
        #fin.kk$inc.size<-ifelse(length(na.omit(match(c("size", "size.sd"), unique(sub.kk$variable))))>0, "Yes", "No")
        #fin.kk$inc.size.sd<-ifelse(length(na.omit(match(c("size.sd"), unique(sub.kk$variable))))>0, "Yes", "No")
        ##has a result been generated?
        result<-TRUE
        counter<-counter+1
      }	
      ##if there are multiple predictors then make sure there arent multiples of the same predictor 
      ## in the EWS (i.e. make sure you dont include CV twice)
      if(h!=1){if(length(unique(unlist(combs[m,])))==h){
        ##get results for each leading indicator from the experimental data
        sub.kk<-kk[(kk$variable%in%unlist(combs[m,])),]
        print(sub.kk)
        ##create an id for the combined predictor
        sub.kk$id<-paste(sort(unique(sub.kk$variable)), collapse=" + ")
        ##calculate the rolling values
        roll<-colSums(do.call("rbind", mclapply(split(sub.kk, list(sub.kk$variable)), function(d){
          if(length(d[,1])>0){
            return(d$value)
          }})))
        ##make a data set of neccesary info to add the rolling value to						
        fin.kk<-sub.kk[1:(length(sub.kk[,1])/length(unique(sub.kk$variable))),c(1:3,5,6,8)]
        ##add the rolling value
        fin.kk$value<-roll
        ##calculate the sd and mean of the rolling value
        fin.kk$roll.sd<-rolling_sd(fin.kk$value)
        fin.kk$roll.mean<-rolling_mean(fin.kk$value)
        ##add in information on whether size is included
        #fin.kk$inc.size<-ifelse(length(na.omit(match(c("size", "size.sd"), unique(sub.kk$variable))))>0, "Yes", "No")
        #fin.kk$inc.size.sd<-ifelse(length(na.omit(match(c("size.sd"), unique(sub.kk$variable))))>0, "Yes", "No")
        ##has a result been generated?
        result<-TRUE
        counter<-counter+1
      }}
      ##if there is a result add it to the finished data frame, along with info on the number of predictors
      if(result==TRUE){
        fin.kk$no.of.pred<-h
        finished.results[[counter]]<-fin.kk}
    }
  }
  return(do.call("rbind",finished.results))
}))











###########################################################################

## plot part 
ar1.all <- long.run.roll.nloci[long.run.roll.nloci$id == "ar1",]
ar1.mean.roll = rep(0,each=92)
ar1.sd.roll = rep(0,each=92)
ar1.values = c()
for (n in abrupt){
  ar1.sd.roll = ar1.sd.roll + ar1.all[ar1.all$Tube == n,8]
  ar1.mean.roll = ar1.mean.roll + ar1.all[ar1.all$Tube == n,9]
  ar1.values = append(ar1.values, ar1.all[ar1.all$Tube == n,7])
}

ar1.sd.roll = ar1.sd.roll/50
ar1.mean.roll = ar1.mean.roll/50
type = rep("val", each = 92*50)
type = append(type,rep(c("mean","+1sd","+2sd","-1sd", "-2sd"), each=92))
ar1.gg = data.frame( type=type, values = c(ar1.values, ar1.mean.roll, ar1.mean.roll+ar1.sd.roll, ar1.mean.roll+2*ar1.sd.roll, ar1.mean.roll-ar1.sd.roll, ar1.mean.roll-2*ar1.sd.roll), time = rep(1:92,times=55))

AR1 <- ggplot(ar1.gg, aes(x=time, y=values, color=type, shape=type, size=type))+
  geom_point(data=ar1.gg[ar1.gg[,1]=="val", ])+
  geom_line(data=ar1.gg[ar1.gg[,1]!="val", ])+
  theme(legend.position="none")+
  scale_shape_manual(values=c(19,19,19,19,19,3))+
  scale_color_manual(values=c( "mean"="#696969","-2sd"='#A9A9A9',"+2sd"="#A9A9A9","-1sd"='#808080',"+1sd"="#808080",  "val"="#B47328"))+
  scale_size_manual(values=c(1,1,1,1,0.25,1))+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(0,25,50,75),
                     label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="ar1")



sd.all <- long.run.roll.nloci[long.run.roll.nloci$id == "sd",]
sd.mean.roll = rep(0,each=92)
sd.sd.roll = rep(0,each=92)
sd.values = c()
for (n in abrupt){
  sd.sd.roll = sd.sd.roll + sd.all[sd.all$Tube == n,8]
  sd.mean.roll = sd.mean.roll + sd.all[sd.all$Tube == n,9]
  sd.values = append(sd.values, sd.all[sd.all$Tube == n,7])
}

sd.sd.roll = sd.sd.roll/50
sd.mean.roll = sd.mean.roll/50
type = rep("val", each = 92*50)
type = append(type,rep(c("mean","+1sd","+2sd","-1sd", "-2sd"), each=92))
sd.gg = data.frame( type=type, values = c(sd.values, sd.mean.roll, sd.mean.roll+sd.sd.roll, sd.mean.roll+2*sd.sd.roll, sd.mean.roll-sd.sd.roll, sd.mean.roll-2*sd.sd.roll), time = rep(1:92,times=55))

SD <- ggplot(sd.gg, aes(x=time, y=values, color=type, shape=type, size=type))+
  geom_point(data=sd.gg[sd.gg[,1]=="val", ])+
  geom_line(data=sd.gg[sd.gg[,1]!="val", ])+
  theme(legend.position="none")+
  scale_shape_manual(values=c(19,19,19,19,19,3))+
  scale_color_manual(values=c( "mean"="#696969","-2sd"='#A9A9A9',"+2sd"="#A9A9A9","-1sd"='#808080',"+1sd"="#808080",  "val"="#B47328"))+
  scale_size_manual(values=c(1,1,1,1,0.25,1))+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(0,25,50,75),
                     label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="standard dev.")+
  ylim(-5,50)




sk.all <- long.run.roll.nloci[long.run.roll.nloci$id == "sk",]
sk.mean.roll = rep(0,each=92)
sk.sd.roll = rep(0,each=92)
sk.values = c()
for (n in abrupt){
  sk.sd.roll = sk.sd.roll + sk.all[sk.all$Tube == n,8]
  sk.mean.roll = sk.mean.roll + sk.all[sk.all$Tube == n,9]
  sk.values = append(sk.values, sk.all[sk.all$Tube == n,7])
}

sk.sd.roll = sk.sd.roll/50
sk.mean.roll = sk.mean.roll/50
type = rep("val", each = 92*50)
type = append(type,rep(c("mean","+1sd","+2sd","-1sd", "-2sd"), each=92))
sk.gg = data.frame( type=type, values = c(sk.values, sk.mean.roll, sk.mean.roll+sk.sd.roll, sk.mean.roll+2*sk.sd.roll, sk.mean.roll-sk.sd.roll, sk.mean.roll-2*sk.sd.roll), time = rep(1:92,times=55))

SK <- ggplot(sk.gg, aes(x=time, y=values, color=type, shape=type, size=type))+
  geom_point(data=sk.gg[sk.gg[,1]=="val", ])+
  geom_line(data=sk.gg[sk.gg[,1]!="val", ])+
  theme(legend.position="none")+
  scale_shape_manual(values=c(19,19,19,19,19,3))+
  scale_color_manual(values=c( "mean"="#696969","-2sd"='#A9A9A9',"+2sd"="#A9A9A9","-1sd"='#808080',"+1sd"="#808080",  "val"="#B47328"))+
  scale_size_manual(values=c(1,1,1,1,0.25,1))+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(0,25,50,75),
                     label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="skewness")




kurt.all <- long.run.roll.nloci[long.run.roll.nloci$id == "kurt",]
kurt.mean.roll = rep(0,each=92)
kurt.sd.roll = rep(0,each=92)
kurt.values = c()
for (n in abrupt){
  kurt.sd.roll = kurt.sd.roll + kurt.all[kurt.all$Tube == n,8]
  kurt.mean.roll = kurt.mean.roll + kurt.all[kurt.all$Tube == n,9]
  kurt.values = append(kurt.values, kurt.all[kurt.all$Tube == n,7])
}

kurt.sd.roll = kurt.sd.roll/50
kurt.mean.roll = kurt.mean.roll/50
type = rep("val", each = 92*50)
type = append(type,rep(c("mean","+1sd","+2sd","-1sd", "-2sd"), each=92))
kurt.gg = data.frame( type=type, values = c(kurt.values, kurt.mean.roll, kurt.mean.roll+kurt.sd.roll, kurt.mean.roll+2*kurt.sd.roll, kurt.mean.roll-kurt.sd.roll, kurt.mean.roll-2*kurt.sd.roll), time = rep(1:92,times=55))

KURT <- ggplot(kurt.gg, aes(x=time, y=values, color=type, shape=type, size=type))+
  geom_point(data=kurt.gg[kurt.gg[,1]=="val", ])+
  geom_line(data=kurt.gg[kurt.gg[,1]!="val", ])+
  theme(legend.position="none")+
  scale_shape_manual(values=c(19,19,19,19,19,3))+
  scale_color_manual(values=c( "mean"="#696969","-2sd"='#A9A9A9',"+2sd"="#A9A9A9","-1sd"='#808080',"+1sd"="#808080",  "val"="#B47328"))+
  scale_size_manual(values=c(1,1,1,1,0.25,1))+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(0,25,50,75),
                     label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="kurtosis")



ggarrange(AR1,SD,KURT,SK,
          labels = c("AR1", "SD", "KURT", "SK"),
          ncol = 2, nrow = 2)














###########################################################################
# percentage above 1 or 2 sd

ar1.sd1 = rep(0, each=92)
ar1.sd2 = rep(0, each=92)
ar1.nb.val = rep(50, each=92)
for (n in abrupt){
  dd <- ar1.all[ar1.all$Tube == n,]
  for (i in 1:92){
    if ( is.na(dd[i,7]) || is.na(dd[i,8]) || is.na(dd[i,9]) ){ ar1.nb.val[i] = ar1.nb.val[i]-1}
    else if (dd[i,7] > (dd[i,9]+2*dd[i,8])){ ar1.sd2[i] = ar1.sd2[i]+1
    ar1.sd1[i] = ar1.sd1[i]+1}
    else if (dd[i,7] > (dd[i,9]+dd[i,8])){ ar1.sd1[i] = ar1.sd1[i]+1 }
  }
}
ar1.sd1 = 100*ar1.sd1/ar1.nb.val
ar1.sd2 = 100*ar1.sd2/ar1.nb.val


sd.sd1 = rep(0, each=92)
sd.sd2 = rep(0, each=92)
sd.nb.val = rep(50, each=92)
no.run = matrix(ncol=4, nrow=92)
k = 0
for (n in abrupt){
  k = k+1
  dd <- sd.all[sd.all$Tube == n,]
  for (i in 1:92){
#    no.run[i] = c(dd[i,7],dd[i,8], dd[i,9], (dd[i,8]+2*dd[i,9]))
    if ( is.na(dd[i,7]) || is.na(dd[i,8]) || is.na(dd[i,9]) ){ sd.nb.val[i] = sd.nb.val[i]-1}
    else if (dd[i,7] > (dd[i,9]+2*dd[i,8])){ sd.sd2[i] = sd.sd2[i]+1
    sd.sd1[i] = sd.sd1[i]+1}
    else if (dd[i,7] > (dd[i,9]+dd[i,8])){ sd.sd1[i] = sd.sd1[i]+1}
  }
}
sd.sd1 = 100*sd.sd1/sd.nb.val
sd.sd2 = 100*sd.sd2/sd.nb.val



sk.sd1 = rep(0, each=92)
sk.sd2 = rep(0, each=92)
sk.nb.val = rep(50, each=92)
for (n in abrupt){
  dd <- sk.all[sk.all$Tube == n,]
  for (i in 1:92){
    if ( is.na(dd[i,7]) || is.na(dd[i,8]) || is.na(dd[i,9]) ){ sk.nb.val[i] = sk.nb.val[i]-1}
    else if (dd[i,7] > (dd[i,9]+2*dd[i,8])){ sk.sd2[i] = sk.sd2[i]+1
    sk.sd1[i] = sk.sd1[i]+1}
    else if (dd[i,7] > (dd[i,9]+dd[i,8])){ sk.sd1[i] = sk.sd1[i]+1 }
  }
}
sk.sd1 = 100*sk.sd1/sk.nb.val
sk.sd2 = 100*sk.sd2/sk.nb.val



kurt.sd1 = rep(0, each=92)
kurt.sd2 = rep(0, each=92)
kurt.nb.val = rep(50, each=92)
for (n in abrupt){
  dd <- kurt.all[kurt.all$Tube == n,]
  for (i in 1:92){
    if ( is.na(dd[i,7]) || is.na(dd[i,8]) || is.na(dd[i,9]) ){ kurt.nb.val[i] = kurt.nb.val[i]-1}
    else if (dd[i,7] > (dd[i,9]+2*dd[i,8])){ kurt.sd2[i] = kurt.sd2[i]+1
    kurt.sd1[i] = kurt.sd1[i]+1}
    else if (dd[i,7] > (dd[i,9]+dd[i,8])){ kurt.sd1[i] = kurt.sd1[i]+1 }
  }
}
kurt.sd1 = 100*kurt.sd1/kurt.nb.val
kurt.sd2 = 100*kurt.sd2/kurt.nb.val


ar1.gg2 = data.frame("time"= rep(1:92, times=2), "val"= append(ar1.sd1,ar1.sd2), "type"= rep(c("sd1", "sd2"), each=92) )
AR2 <- ggplot(ar1.gg2, aes(x=time, y=val, color=type))+
  geom_point()+
  geom_line()+
  theme(legend.position="none")+
  scale_color_manual(values=c("sd1"="#B47328", "sd2"="#624E13"))+
  scale_x_continuous(name ="Generations to speciation", 
                      breaks = c(0,25,50,75),
                      label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="% outside sd")


sd.gg2 = data.frame("time"= rep(1:92, times=2), "val"= append(sd.sd1,sd.sd2), "type"= rep(c("sd1", "sd2"), each=92) )
SD2 <- ggplot(sd.gg2, aes(x=time, y=val, color=type))+
  geom_point()+
  geom_line()+
  theme(legend.position="none")+
  scale_color_manual(values=c("sd1"="#B47328", "sd2"="#624E13"))+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(0,25,50,75),
                     label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="% outside sd")


sk.gg2 = data.frame("time"= rep(1:92, times=2), "val"= append(sk.sd1,sk.sd2), "type"= rep(c("sd1", "sd2"), each=92) )
SK2 <- ggplot(sk.gg2, aes(x=time, y=val, color=type))+
  geom_point()+
  geom_line()+
  theme(legend.position="none")+
  scale_color_manual(values=c("sd1"="#B47328", "sd2"="#624E13"))+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(0,25,50,75),
                     label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="% outside sd")


kurt.gg2 = data.frame("time"= rep(1:92, times=2), "val"= append(kurt.sd1,kurt.sd2), "type"= rep(c("sd1", "sd2"), each=92) )
KURT2 <- ggplot(kurt.gg2, aes(x=time, y=val, color=type))+
  geom_point()+
  geom_line()+
  theme(legend.position="none")+
  scale_color_manual(values=c("sd1"="#B47328", "sd2"="#624E13"))+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(0,25,50,75),
                     label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="% outside sd")




ggarrange(AR2,SD2,KURT2,SK2,
          labels = c("AR1", "SD", "KURT", "SK"),
          ncol = 2, nrow = 2)










###########################################################################
# composite index



comp.all <- long.run.roll.nloci[long.run.roll.nloci$id == "ar1 + sd + sk + kurt",]
comp.mean.roll = rep(0,each=92)
comp.sd.roll = rep(0,each=92)
comp.values = c()
for (n in abrupt){
  comp.sd.roll = comp.sd.roll + comp.all[comp.all$Tube == n,8]
  comp.mean.roll = comp.mean.roll + comp.all[comp.all$Tube == n,9]
  comp.values = append(comp.values, comp.all[comp.all$Tube == n,7])
}

comp.sd.roll = comp.sd.roll/50
comp.mean.roll = comp.mean.roll/50
type = rep("val", each = 92*50)
type = append(type,rep(c("mean","+1sd","+2sd","-1sd", "-2sd"), each=92))
comp.gg = data.frame( type=type, values = c(comp.values, comp.mean.roll, comp.mean.roll+comp.sd.roll, comp.mean.roll+2*comp.sd.roll, comp.mean.roll-comp.sd.roll, comp.mean.roll-2*comp.sd.roll), time = rep(1:92,times=55))

COMP <- ggplot(comp.gg, aes(x=time, y=values, color=type, shape=type, size=type))+
  geom_point(data=comp.gg[comp.gg[,1]=="val", ])+
  geom_line(data=comp.gg[comp.gg[,1]!="val", ])+
  theme(legend.position="none")+
  scale_shape_manual(values=c(19,19,19,19,19,3))+
  scale_color_manual(values=c( "mean"="#696969","-2sd"='#A9A9A9',"+2sd"="#A9A9A9","-1sd"='#808080',"+1sd"="#808080",  "val"="#B47328"))+
  scale_size_manual(values=c(1,1,1,1,0.25,1))+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(0,25,50,75),
                     label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="composite")+
  ylim(-5,50)



comp.sd1 = rep(0, each=92)
comp.sd2 = rep(0, each=92)
comp.nb.val = rep(50, each=92)
for (n in abrupt){
  dd <- comp.all[comp.all$Tube == n,]
  for (i in 1:92){
    if ( is.na(dd[i,7]) || is.na(dd[i,8]) || is.na(dd[i,9]) ){ comp.nb.val[i] = comp.nb.val[i]-1}
    else if (dd[i,7] > (dd[i,9]+2*dd[i,8])){ comp.sd2[i] = comp.sd2[i]+1
    comp.sd1[i] = comp.sd1[i]+1}
    else if (dd[i,7] > (dd[i,9]+dd[i,8])){ comp.sd1[i] = comp.sd1[i]+1 }
  }
}
comp.sd1 = 100*comp.sd1/comp.nb.val
comp.sd2 = 100*comp.sd2/comp.nb.val



comp.gg2 = data.frame("time"= rep(1:92, times=2), "val"= append(comp.sd1,comp.sd2), "type"= rep(c("sd1", "sd2"), each=92) )
COMP2 <- ggplot(comp.gg2, aes(x=time, y=val, color=type))+
  geom_point()+
  geom_line()+
  theme(legend.position="none")+
  scale_color_manual(values=c("sd1"="#B47328", "sd2"="#624E13"))+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(0,25,50,75),
                     label=c( "920000","670000","42000","17000" ))+
  labs(x="Generations to speciation", y="% outside sd")+
  ylim(0,100)




ggarrange(COMP,COMP2,
          labels = c(),
          ncol = 2, nrow = 1)



