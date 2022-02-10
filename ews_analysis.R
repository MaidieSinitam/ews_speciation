library(ggplot2)
library(ggpubr)
library(earlywarnings)

# Load necessary variables
# Load code names for the 6 differentiation metrics
metrics = c("me","meanfit", "fst", "dxysel", "ld", "nblocisel")
abruptN = seq(300001, 300100, by=1)
gradualN = c(300160, 300115, 300193, 300185, 300146, 300104, 300187, 300188, 300174, 300113, 300137, 300144, 
             300108, 300114, 300142, 300159, 300109, 300127, 300165, 300140, 300199, 300198, 300124, 300176, 300105, 
             300121, 300141, 300122, 300134, 300173, 300168, 300190, 300102, 300118, 300191, 300116, 300196, 300166, 
             300182, 300120, 300136, 300128, 300110, 300186, 300129, 300164, 300194, 300157, 300138, 300123)
neutralA = seq(300301, 300400, by=1) #m=0.1
neutralG = seq(300401, 300500, by=1) #m=0.01

# Load time of speciation of each simulation output from the spectime.txt file 
T = read.table("spectime.txt", sep=" ", header=TRUE)
specTimeRowAN=T$specTimeRowAN-30
timeToSpecAN = specTimeRowAN*10
specTimeRowGN=T$specTimeRowGN-30
timeToSpecGN = specTimeRowGN*10
lenAN=T$lenAN
lenGN=T$lenGN
lenNeutA = T$lenNeutA
lenNeutG=T$lenNeutG

#####################################################################################
# Detrending
# Code to make appendix S5

ts=getMetricsNew(abruptN[1],cut=TRUE, sfreq=30)
end = length(ts$time)
timeindex = seq(1,end,by=1)


#me
bwidth = c(floor(5 * end/100), floor(25 * end/100))
time = rep(ts$time, times=length(bwidth)+1)
Y = ts$me
y =c()
i=0
for (b in bwidth){
  smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = b, range.x = range(timeindex), 
                  x.points = timeindex)
  smY <- smYY$y
  y = append(y, smY)
  i=i+1
}
bwidth = c(5,25)
bwidth = apply(X=as.matrix(bwidth),1,toString)
bwidth = append(bwidth, "0 - raw")
y = append(y,Y)
bwidth = rep(bwidth, each=end)
out = data.frame(bwidth,time,y)

me = ggplot(data=out, aes(x=time, y=y, group=bwidth)) +
  geom_line(aes(color=bwidth), size=0.5)+
  scale_color_manual(values=c('#999999', "#A12823", "#17657D"))+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  scale_size_manual(values=c(0.5,0.25,0.25,0.25,0.25))+
  labs(x="Generations", y = "effective migration rate")+
  theme(legend.position = "none")




#fst
bwidth = c(floor(5 * end/100), floor(25 * end/100))
time = rep(ts$time, times=length(bwidth)+1)
Y = ts$fst
y =c()
i=0
for (b in bwidth){
  smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = b, range.x = range(timeindex), 
                  x.points = timeindex)
  smY <- smYY$y
  y = append(y, smY)
  i=i+1
}
bwidth = c(5,25)
bwidth = apply(X=as.matrix(bwidth),1,toString)
bwidth = append(bwidth, "0 - raw")
y = append(y,Y)
bwidth = rep(bwidth, each=end)
out = data.frame(bwidth,time,y)


fst = ggplot(data=out, aes(x=time, y=y, group=bwidth)) +
  geom_line(aes(color=bwidth), size=0.5)+
  scale_color_manual(values=c('#999999', "#A12823", "#17657D"))+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  scale_size_manual(values=c(0.25,0.25,0.25,0.25, 0.25))+
  labs(x="Generations", y = "FST")+
  theme(legend.position = "none")




#nblocisel
bwidth = c(floor(5 * end/100), floor(25 * end/100))
time = rep(ts$time, times=length(bwidth)+1)
Y = ts$nblocisel
y =c()
i=0
for (b in bwidth){
  smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = b, range.x = range(timeindex), 
                  x.points = timeindex)
  smY <- smYY$y
  y = append(y, smY)
  i=i+1
}
bwidth = c(5,25)
bwidth = apply(X=as.matrix(bwidth),1,toString)
bwidth = append(bwidth, "0 - raw")
y = append(y,Y)
bwidth = rep(bwidth, each=end)
out = data.frame(bwidth,time,y)


nblocisel = ggplot(data=out, aes(x=time, y=y, group=bwidth)) +
  geom_line(aes(color=bwidth), size=0.5)+
  scale_color_manual(values=c('#999999', "#A12823", "#17657D"))+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  scale_size_manual(values=c(0.25,0.25,0.25,0.25, 0.25))+
  labs(x="Generations", y = "Nb of selected variable loci")+
  theme(legend.position = "none")





#fitrand
bwidth = c(floor(5 * end/100), floor(25 * end/100))
time = rep(ts$time, times=length(bwidth)+1)
Y = ts$fitall
y =c()
i=0
for (b in bwidth){
  smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = b, range.x = range(timeindex), 
                  x.points = timeindex)
  smY <- smYY$y
  y = append(y, smY)
  i=i+1
}
bwidth = c(5,25)
bwidth = apply(X=as.matrix(bwidth),1,toString)
bwidth = append(bwidth, "0 - raw")
y = append(y,Y)
bwidth = rep(bwidth, each=end)
out = data.frame(bwidth,time,y)


fitrand = ggplot(data=out, aes(x=time, y=y, group=bwidth)) +
  geom_line(aes(color=bwidth), size=0.5)+
  scale_color_manual(values=c('#999999',"#A12823", "#17657D"))+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  scale_size_manual(values=c(0.25,0.25,0.25,0.25))+
  labs(x="Generations", y = "Mean fitness")+
  theme(legend.position = "none")





#dxysel
bwidth = c(floor(5 * end/100), floor(25 * end/100))
time = rep(ts$time, times=length(bwidth)+1)
Y = ts$dxysel
y =c()
i=0
for (b in bwidth){
  smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = b, range.x = range(timeindex), 
                  x.points = timeindex)
  smY <- smYY$y
  y = append(y, smY)
  i=i+1
}
bwidth = c(5,25)
bwidth = apply(X=as.matrix(bwidth),1,toString)
bwidth = append(bwidth, "0 - raw")
y = append(y,Y)
bwidth = rep(bwidth, each=end)
out = data.frame(bwidth,time,y)


dxy = ggplot(data=out, aes(x=time, y=y, group=bwidth)) +
  geom_line(aes(color=bwidth), size=0.5)+
  scale_color_manual(values=c('#999999', "#A12823", "#17657D"))+
  # scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  scale_size_manual(values=c(0.25,0.25,0.25,0.25, 0.25))+
  labs(x="Generations", y = "DXY")+
  theme(legend.position = "none")



#ld
bwidth = c(floor(5 * end/100), floor(25 * end/100))
time = rep(ts$time, times=length(bwidth)+1)
Y = ts$ld
y =c()
i=0
for (b in bwidth){
  smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = b, range.x = range(timeindex), 
                  x.points = timeindex)
  smY <- smYY$y
  y = append(y, smY)
  i=i+1
}
bwidth = c(5,25)
bwidth = apply(X=as.matrix(bwidth),1,toString)
bwidth = append(bwidth, "0 - raw")
y = append(y,Y)
bwidth = rep(bwidth, each=end)
out = data.frame(bwidth,time,y)


ld = ggplot(data=out, aes(x=time, y=y, group=bwidth)) +
  geom_line(aes(color=bwidth), size=0.5)+
  scale_color_manual(values=c('#999999', "#A12823", "#17657D"))+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  scale_size_manual(values=c(0.25,0.25,0.25,0.25,0.25))+
  labs(x="Generations", y = "LD")+
  theme(legend.position = "none")


ggarrange(me, fst, nblocisel, fitrand,
          labels = c("", "", "", ""),
          ncol = 2, nrow = 2)
ggarrange(dxy, ld,
          labels = c( "", ""),
          ncol = 2, nrow = 2)


#####################################################################################
# SENSITIVITY ANALYSIS 

# Modified version of the sensitivity_ews function from the earlywarnings package
# The getMetrics() function is required to run sensitivity_ews2
sensitivity_ews2 <- function(timeseries, ind, indicator = c("ar1", "sd", "acf1", "sk", "kurt", "cv", "returnrate", "densratio"), winsizerange = c(25, 75), incrwinsize = 25, 
                             detrending = c("no", "gaussian", "loess", "linear", "first-diff"), bandwidthrange = c(5, 100), spanrange = c(5, 100), degree = NULL, incrbandwidth = 20, 
                             incrspanrange = 10, logtransform = FALSE, interpolate = FALSE) {
  
  # timeseries<-ts(timeseries) #strict data-types the input data as tseries object
  # for use in later steps
  
  metrics = c("time", "me", "wreswimm", "wmaxwres", "fitall", "fitrand", "fst", "dxyall", "dxysel", "dxyneut", "ld", "nbloci", "nblocisel", "nblocineut")
  timeseries <- getMetrics(timeseries, cut=TRUE)
  ind = which(metrics == ind) #strict data-types the input data as tseries object for use in later steps
  Y <- timeseries[, ind]
  end = length(Y)
  timeindex <- seq(1,end,by=1)
  
  # Interpolation
  if (interpolate) {
    YY <- approx(timeindex, Y, n = length(Y), method = "linear")
    Y <- YY$y
  } else {
    Y <- Y
  }
  
  # Log-transformation
  if (logtransform) {
    Y <- log(Y + 1)
  }
  
  # Determine the step increases in rolling windowsize
  incrtw <- incrwinsize
  tw <- seq(floor(winsizerange[1] * length(Y)/100), floor(winsizerange[2] * length(Y)/100), 
            by = incrtw)
  twcol <- length(tw)
  low <- 6
  
  # Detrending
  detrending <- match.arg(detrending)
  if (detrending == "gaussian") {
    incrbw <- incrbandwidth
    width <- c(floor(5 * length(Y)/100), floor(10 * length(Y)/100), floor(15 * length(Y)/100), floor(20 * length(Y)/100), floor(25 * length(Y)/100))
    #width <- seq(floor(bandwidthrange[1] * length(Y)/100), floor(bandwidthrange[2] * length(Y)/100), 
    #            by = incrbandwidth)
    bwrow <- length(width)
    # Create matrix to store Kendall trend statistics
    Ktauestind <- matrix(, bwrow, twcol)
    Ktaupind <- matrix(, bwrow, twcol)
    # Estimation
    for (wi in 1:(length(width))) {
      width1 <- width[wi]
      smYY <- ksmooth(timeindex, Y, kernel = c("normal"), bandwidth = width1, 
                      range.x = range(timeindex), n.points = length(timeindex))
      nsmY <- Y - smYY$y
      for (ti in 1:length(tw)) {
        tw1 <- tw[ti]
        # Rearrange data for indicator calculation
        omw1 <- length(nsmY) - tw1 + 1  ##number of overlapping moving windows
        high <- omw1
        nMR1 <- matrix(data = NA, nrow = tw1, ncol = omw1)
        for (i in 1:omw1) {
          Ytw <- nsmY[i:(i + tw1 - 1)]
          nMR1[, i] <- Ytw
        }
        # Estimate indicator
        indicator = match.arg(indicator)
        if (indicator == "ar1") {
          indic <- apply(nMR1, 2, function(x) {
            nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, demean = TRUE, 
                           intercept = FALSE)
            nAR1$ar
          })
        } else if (indicator == "sd") {
          indic <- apply(nMR1, 2, sd)
        } else if (indicator == "sk") {
          indic <- apply(nMR1, 2, skewness)
        } else if (indicator == "kurt") {
          indic <- apply(nMR1, 2, kurtosis)
        } else if (indicator == "acf1") {
          indic <- apply(nMR1, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            nACF$acf[2]
          })
        } else if (indicator == "returnrate") {
          indic <- apply(nMR1, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            1 - nACF$acf[2]
          })
        } else if (indicator == "cv") {
          indic <- apply(nMR1, 2, function(x) {
            sd(x)/mean(x)
          })
        } else if (indicator == "densratio") {
          indic <- apply(nMR1, 2, function(x) {
            spectfft <- spec.ar(x, n.freq = omw1, plot = FALSE, order = 1)
            spectfft$spec
            spectfft$spec[low]/spectfft$spec[high]
          })
        }
        # Calculate trend statistics
        timevec <- seq(floor(0.5*length(indic)), length(indic), by=1)
        indic2 <- indic[floor(0.5*length(indic)):length(indic)]
        Kt <- cor.test(timevec, indic2, alternative = c("two.sided"), method = c("kendall"), 
                       conf.level = 0.95)
        Ktauestind[wi, ti] <- Kt$estimate
        Ktaupind[wi, ti] <- Kt$p.value
      }
    }
    
    # Output
    #out <- data.frame(Ktauestind)
    #colnames(out) <- tw
    #rownames(out) <- width
    output = list("gaussian",Ktauestind, Ktaupind, tw, indicator,width)
    return(output)
    
  } else if (detrending == "loess") {
    incrbw <- incrspanrange
    width <- seq(floor(spanrange[1] * length(Y)/100), floor(spanrange[2] * length(Y)/100), 
                 by = incrbw)
    bwrow <- length(width)
    # Create matrix to store Kendall trend statistics
    Ktauestind <- matrix(, bwrow, twcol)
    Ktaupind <- matrix(, bwrow, twcol)
    # Estimation
    if (is.null(degree)) {
      degree <- 2
    } else {
      degree <- degree
    }
    for (wi in 1:(length(width))) {
      width1 <- width[wi]
      smYY <- loess(Y ~ timeindex, span = width1, degree = degree, normalize = FALSE, 
                    family = "gaussian")
      smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
      nsmY <- Y - smY
      for (ti in 1:length(tw)) {
        tw1 <- tw[ti]
        # Rearrange data for indicator calculation
        omw1 <- length(nsmY) - tw1 + 1  ##number of overlapping moving windows
        high <- omw1
        nMR1 <- matrix(data = NA, nrow = tw1, ncol = omw1)
        for (i in 1:omw1) {
          Ytw <- nsmY[i:(i + tw1 - 1)]
          nMR1[, i] <- Ytw
        }
        # Estimate indicator
        indicator = match.arg(indicator)
        if (indicator == "ar1") {
          indic <- apply(nMR1, 2, function(x) {
            nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, demean = TRUE, 
                           intercept = FALSE)
            nAR1$ar
          })
        } else if (indicator == "sd") {
          indic <- apply(nMR1, 2, sd)
        } else if (indicator == "sk") {
          indic <- apply(nMR1, 2, skewness)
        } else if (indicator == "kurt") {
          indic <- apply(nMR1, 2, kurtosis)
        } else if (indicator == "acf1") {
          indic <- apply(nMR1, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            nACF$acf[2]
          })
        } else if (indicator == "returnrate") {
          indic <- apply(nMR1, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            1 - nACF$acf[2]
          })
        } else if (indicator == "cv") {
          indic <- apply(nMR1, 2, function(x) {
            sd(x)/mean(x)
          })
        } else if (indicator == "densratio") {
          indic <- apply(nMR1, 2, function(x) {
            spectfft <- spec.ar(x, n.freq = omw1, plot = FALSE, order = 1)
            spectfft$spec
            spectfft$spec[low]/spectfft$spec[high]
          })
        }
        # Calculate trend statistics
        timevec <- seq(1, length(indic))
        Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
                       conf.level = 0.95)
        Ktauestind[wi, ti] <- Kt$estimate
        Ktaupind[wi, ti] <- Kt$p.value
      }
    }
    #loess plot part removed here
    # Output
    #out <- data.frame(Ktauestind)
    #colnames(out) <- tw
    #rownames(out) <- width
    
    output = list("loess",Ktauestind, Ktaupind, tw, indicator,width)
    return(output)
    
  } else if (detrending == "linear") {
    nsmY <- resid(lm(Y ~ timeindex))
  } else if (detrending == "first-diff") {
    nsmY <- diff(Y)
  } else if (detrending == "no") {
    nsmY <- Y
  }
  
  # Create matrix to store Kendall trend statistics
  Ktauestind <- matrix(, twcol, 1)
  Ktaupind <- matrix(, twcol, 1)
  
  for (ti in 1:length(tw)) {
    tw1 <- tw[ti]
    # Rearrange data for indicator calculation
    omw1 <- length(nsmY) - tw1 + 1  ##number of overlapping moving windows
    high = omw1
    nMR1 <- matrix(data = NA, nrow = tw1, ncol = omw1)
    for (i in 1:omw1) {
      Ytw <- nsmY[i:(i + tw1 - 1)]
      nMR1[, i] <- Ytw
    }
    # Estimate indicator
    indicator = match.arg(indicator)
    if (indicator == "ar1") {
      indic <- apply(nMR1, 2, function(x) {
        nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, demean = TRUE, intercept = FALSE)
        nAR1$ar
      })
    } else if (indicator == "sd") {
      indic <- apply(nMR1, 2, sd)
    } else if (indicator == "sk") {
      indic <- apply(nMR1, 2, skewness)
    } else if (indicator == "kurt") {
      indic <- apply(nMR1, 2, kurtosis)
    } else if (indicator == "acf1") {
      indic <- apply(nMR1, 2, function(x) {
        nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
        nACF$acf[2]
      })
    } else if (indicator == "returnrate") {
      indic <- apply(nMR1, 2, function(x) {
        nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
        1 - nACF$acf[2]
      })
    } else if (indicator == "cv") {
      indic <- apply(nMR1, 2, function(x) {
        sd(x)/mean(x)
      })
    } else if (indicator == "densratio") {
      indic <- apply(nMR1, 2, function(x) {
        spectfft <- spec.ar(x, n.freq = omw1, plot = FALSE, order = 1)
        spectfft$spec
        spectfft$spec[low]/spectfft$spec[high]
      })
    }
    # Calculate trend statistics
    timevec <- seq(1, length(indic))
    Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
                   conf.level = 0.95)
    Ktauestind[ti] <- Kt$estimate
    Ktaupind[ti] <- Kt$p.value
  }
  
  #plot part removed here
  
  # Output out<-data.frame(tw,Ktauestind,Ktaupind) colnames(out)<-c('rolling
  # window','Kendall tau estimate','Kendall tau p-value')
  #out <- data.frame(Ktauestind)
  #Ktaupind <- data.frame(Ktaupind)
  #rownames(out) <- tw
  
  output = list(detrend="no",ktest=Ktauestind, ktp=Ktaupind, win=tw, ind=indicator)
  return(output)
} 

# Matrix to store Kendall tau trends
KAR = matrix(nrow=25000, ncol=6)
KSD = matrix(nrow=25000, ncol=6)
KKURT = matrix(nrow=25000, ncol=6)
KSK = matrix(nrow=25000, ncol=6)

# Loop to calculate Kendall trends on all abrupt replicates
k=1
for (i in metrics){
  j=1
  for (n in abrupt){
    Y = sensitivity_ews2(n, ind = i, detrending="gaussian", incrwinsize = 5, incrbandwidth = 5, bandwidthrange = c(5,25), indicator="ar1")
    tw = Y[[4]]
    bw = Y[[6]]
    
    for (p in seq(1,length(bw),by=1)){
      KAR[j:(j+length(tw)-1),k] = Y[[2]][p,]
      KSD[j:(j+length(tw)-1),k] = sensitivity_ews2(n, ind = i, incrbandwidth = 5, bandwidthrange = c(5,25), detrending="gaussian", incrwinsize = 5, indicator="sd")[[2]][p,]
      KKURT[j:(j+length(tw)-1),k] = sensitivity_ews2(n, ind = i, detrending="gaussian",bandwidthrange = c(5,25), incrwinsize = 5, indicator="kurt")[[2]][p,]
      KSK[j:(j+length(tw)-1),k] = sensitivity_ews2(n, ind = i, detrending="gaussian",bandwidthrange = c(5,25), incrwinsize = 5,indicator="sk")[[2]][p,]
      #KCV[j:(j+length(tw)-1),k] = sensitivity_ews2(n, ind = i, detrending="gaussian", incrwinsize = 5,indicator="cv")[[2]][p,]
      j=j+length(tw)
    }
  }
  k=k+1
}



#####################################################################################
### MAKING SURROGATES


generate_iAAFT = function(X, specflag) {
  
  nargin = length(as.list(match.call())) -1 
  
  if (nargin<2) {
    specflag = 1
  }
  if (length(X) == 0){
    print('input vector is empty', quote=FALSE)
    #  return
  }
  
  max_it = 500
  pp = length(X)
  
  # Initial Conditions
  rn = X[sample(pp)]
  Xsorted = sort(X)	# Desired signal distribution
  Yamp = abs(fft(X, inverse = FALSE))	# Desired amplitude spectrum
  
  prev_err = 0
  E = rep(0,max_it)
  c = 1
  prev_err = 1e10
  err = prev_err - 1
  while ((prev_err>err) & (c<max_it)) {
    
    # Match Amplitude Spec
    Yrn = fft(rn, inverse = FALSE)
    Yang = Arg(Yrn)
    sn = Re(fft(Yamp*exp(sqrt(-1+ 0i)*Yang), inverse = TRUE))/length(Yamp*exp(sqrt(-1+ 0i)*Yang))
    
    # Scale to Original Signal Distribution
    sns = sort(sn,index.return = TRUE)
    rn[sns$ix] = Xsorted
    
    # Eval Convergence
    prev_err = err
    A2 = abs(Yrn)
    err = mean(abs(A2-Yamp))
    E[c] = err
    
    c = c+1
  }
  E = E[1:c-1]
  if (specflag==1)
    Xs = sn	# Exact Amplitude Spectrum
  else
    Xs = rn	# Exact Signal Distribution
  
  # output
  out=Xs
  return(out)
}
makeSurrogates <- function(run, indic, bw){
  #100 surrogates for one run for one metric 
  
  metrics = c("time", "me", "wreswimm", "wmaxwres", "fitall", "fitrand", "fst", "dxyall", "dxysel", "dxyneut", "ld", "nbloci", "nblocisel", "nblocineut")
  indic = which(metrics == indic)
  ts <- getMetrics(run, cut=TRUE)
  end = length(ts$time)
  timeindex = seq(1,end,by=1)
  Y = ts[,indic]
  
  # detrending
  bw <- bw*end/100
  smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = bw, range.x = range(timeindex), 
    x.points = timeindex)
  nsmY <- Y - smYY$y
  smY <- smYY$y
  
  out = matrix(ncol=100, nrow=length(Y))
  
  for (i in 1:100){
    out[,i] = generate_iAAFT(Y, specflag=1) #nsmY normalement
  }
  return(out) }


# function for sensitivity analysis for surrogates
sensitivity_ews_surr <- function(timeseries, indicator = c("ar1", "sd", "acf1", "sk", "kurt", "cv", "returnrate", "densratio", "mean"), winsizerange = c(25, 75), incrwinsize = 25, 
                             detrending = c("no", "gaussian", "loess", "linear", "first-diff"), bandwidthrange = c(5, 100), spanrange = c(5, 100), degree = NULL, incrbandwidth = 20, 
                             incrspanrange = 10, logtransform = FALSE, interpolate = FALSE) {
  
  # timeseries<-ts(timeseries) #strict data-types the input data as tseries object
  # for use in later steps
  
  timeseries <- data.matrix(timeseries)
  if (dim(timeseries)[2] == 1) {
    Y = timeseries
    timeindex = 1:dim(timeseries)[1]
  } else if (dim(timeseries)[2] == 2) {
    Y <- timeseries[, 2]
    timeindex <- timeseries[, 1]
  } else {
    warning("not right format of timeseries input")
  }
  
  # Interpolation
  if (interpolate) {
    YY <- approx(timeindex, Y, n = length(Y), method = "linear")
    Y <- YY$y
  } else {
    Y <- Y
  }
  
  # Log-transformation
  if (logtransform) {
    Y <- log(Y + 1)
  }
  
  # Determine the step increases in rolling windowsize
  incrtw <- incrwinsize
  tw <- seq(floor(winsizerange[1] * length(Y)/100), floor(winsizerange[2] * length(Y)/100), 
            by = incrtw)
  twcol <- length(tw)
  low <- 6
  
  # Detrending
  detrending <- match.arg(detrending)
  if (detrending == "gaussian") {
    incrbw <- incrbandwidth
    width <- seq(floor(bandwidthrange[1] * length(Y)/100), floor(bandwidthrange[2] * 
                                                                   length(Y)/100), by = incrbw)
    bwrow <- length(width)
    # Create matrix to store Kendall trend statistics
    Ktauestind <- matrix(, bwrow, twcol)
    Ktaupind <- matrix(, bwrow, twcol)
    # Estimation
    for (wi in 1:(length(width))) {
      width1 <- width[wi]
      smYY <- ksmooth(timeindex, Y, kernel = c("normal"), bandwidth = width1, 
                      range.x = range(timeindex), n.points = length(timeindex))
      nsmY <- Y - smYY$y
      for (ti in 1:length(tw)) {
        tw1 <- tw[ti]
        # Rearrange data for indicator calculation
        omw1 <- length(nsmY) - tw1 + 1  ##number of overlapping moving windows
        high <- omw1
        nMR1 <- matrix(data = NA, nrow = tw1, ncol = omw1)
        for (i in 1:omw1) {
          Ytw <- nsmY[i:(i + tw1 - 1)]
          nMR1[, i] <- Ytw
        }
        # Estimate indicator
        indicator = match.arg(indicator)
        if (indicator == "ar1") {
          indic <- apply(nMR1, 2, function(x) {
            nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, demean = TRUE, 
                           intercept = FALSE)
            nAR1$ar
          })
        } else if (indicator == "sd") {
          indic <- apply(nMR1, 2, sd)
        } else if (indicator == "sk") {
          indic <- apply(nMR1, 2, skewness)
        } else if (indicator == "kurt") {
          indic <- apply(nMR1, 2, kurtosis)
        } else if (indicator == "mean") {
          indic <- apply(nMR1, 2, mean)
        } else if (indicator == "acf1") {
          indic <- apply(nMR1, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            nACF$acf[2]
          })
        } else if (indicator == "returnrate") {
          indic <- apply(nMR1, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            1 - nACF$acf[2]
          })
        } else if (indicator == "cv") {
          indic <- apply(nMR1, 2, function(x) {
            sd(x)/mean(x)
          })
        } else if (indicator == "densratio") {
          indic <- apply(nMR1, 2, function(x) {
            spectfft <- spec.ar(x, n.freq = omw1, plot = FALSE, order = 1)
            spectfft$spec
            spectfft$spec[low]/spectfft$spec[high]
          })
        }
        # Calculate trend statistics
        timevec <- seq(1, length(indic))
        Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
                       conf.level = 0.95)
        Ktauestind[wi, ti] <- Kt$estimate
        Ktaupind[wi, ti] <- Kt$p.value
      }
    }
    
    # Output
    #out <- data.frame(Ktauestind)
    #colnames(out) <- tw
    #rownames(out) <- width
    output = list("gaussian",Ktauestind, Ktaupind, tw, indicator,width)
    return(output)
    
  } else if (detrending == "loess") {
    incrbw <- incrspanrange
    width <- seq(floor(spanrange[1] * length(Y)/100), floor(spanrange[2] * length(Y)/100), 
                 by = incrbw)
    bwrow <- length(width)
    # Create matrix to store Kendall trend statistics
    Ktauestind <- matrix(, bwrow, twcol)
    Ktaupind <- matrix(, bwrow, twcol)
    # Estimation
    if (is.null(degree)) {
      degree <- 2
    } else {
      degree <- degree
    }
    for (wi in 1:(length(width))) {
      width1 <- width[wi]
      smYY <- loess(Y ~ timeindex, span = width1, degree = degree, normalize = FALSE, 
                    family = "gaussian")
      smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
      nsmY <- Y - smY
      for (ti in 1:length(tw)) {
        tw1 <- tw[ti]
        # Rearrange data for indicator calculation
        omw1 <- length(nsmY) - tw1 + 1  ##number of overlapping moving windows
        high <- omw1
        nMR1 <- matrix(data = NA, nrow = tw1, ncol = omw1)
        for (i in 1:omw1) {
          Ytw <- nsmY[i:(i + tw1 - 1)]
          nMR1[, i] <- Ytw
        }
        # Estimate indicator
        indicator = match.arg(indicator)
        if (indicator == "ar1") {
          indic <- apply(nMR1, 2, function(x) {
            nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, demean = TRUE, 
                           intercept = FALSE)
            nAR1$ar
          })
        } else if (indicator == "sd") {
          indic <- apply(nMR1, 2, sd)
        } else if (indicator == "sk") {
          indic <- apply(nMR1, 2, skewness)
        } else if (indicator == "kurt") {
          indic <- apply(nMR1, 2, kurtosis)
        } else if (indicator == "mean") {
          indic <- apply(nMR1, 2, mean)
        } else if (indicator == "acf1") {
          indic <- apply(nMR1, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            nACF$acf[2]
          })
        } else if (indicator == "returnrate") {
          indic <- apply(nMR1, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            1 - nACF$acf[2]
          })
        } else if (indicator == "cv") {
          indic <- apply(nMR1, 2, function(x) {
            sd(x)/mean(x)
          })
        } else if (indicator == "densratio") {
          indic <- apply(nMR1, 2, function(x) {
            spectfft <- spec.ar(x, n.freq = omw1, plot = FALSE, order = 1)
            spectfft$spec
            spectfft$spec[low]/spectfft$spec[high]
          })
        }
        # Calculate trend statistics
        timevec <- seq(1, length(indic))
        Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
                       conf.level = 0.95)
        Ktauestind[wi, ti] <- Kt$estimate
        Ktaupind[wi, ti] <- Kt$p.value
      }
    }
    #loess plot part removed here
    # Output
    #out <- data.frame(Ktauestind)
    #colnames(out) <- tw
    #rownames(out) <- width
    
    output = list("loess",Ktauestind, Ktaupind, tw, indicator,width)
    return(output)
    
  } else if (detrending == "linear") {
    nsmY <- resid(lm(Y ~ timeindex))
  } else if (detrending == "first-diff") {
    nsmY <- diff(Y)
  } else if (detrending == "no") {
    nsmY <- Y
  }
  
  # Create matrix to store Kendall trend statistics
  Ktauestind <- matrix(, twcol, 1)
  Ktaupind <- matrix(, twcol, 1)
  
  for (ti in 1:length(tw)) {
    tw1 <- tw[ti]
    # Rearrange data for indicator calculation
    omw1 <- length(nsmY) - tw1 + 1  ##number of overlapping moving windows
    high = omw1
    nMR1 <- matrix(data = NA, nrow = tw1, ncol = omw1)
    for (i in 1:omw1) {
      Ytw <- nsmY[i:(i + tw1 - 1)]
      nMR1[, i] <- Ytw
    }
    # Estimate indicator
    indicator = match.arg(indicator)
    if (indicator == "ar1") {
      indic <- apply(nMR1, 2, function(x) {
        nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, demean = TRUE, intercept = FALSE)
        nAR1$ar
      })
    } else if (indicator == "sd") {
      indic <- apply(nMR1, 2, sd)
    } else if (indicator == "sk") {
      indic <- apply(nMR1, 2, skewness)
    } else if (indicator == "kurt") {
      indic <- apply(nMR1, 2, kurtosis)
    } else if (indicator == "mean") {
      indic <- apply(nMR1, 2, mean)
    } else if (indicator == "acf1") {
      indic <- apply(nMR1, 2, function(x) {
        nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
        nACF$acf[2]
      })
    } else if (indicator == "returnrate") {
      indic <- apply(nMR1, 2, function(x) {
        nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
        1 - nACF$acf[2]
      })
    } else if (indicator == "cv") {
      indic <- apply(nMR1, 2, function(x) {
        sd(x)/mean(x)
      })
    } else if (indicator == "densratio") {
      indic <- apply(nMR1, 2, function(x) {
        spectfft <- spec.ar(x, n.freq = omw1, plot = FALSE, order = 1)
        spectfft$spec
        spectfft$spec[low]/spectfft$spec[high]
      })
    }
    # Calculate trend statistics
    timevec <- seq(1, length(indic))
    Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
                   conf.level = 0.95)
    Ktauestind[ti] <- Kt$estimate
    Ktaupind[ti] <- Kt$p.value
  }
  
  #plot part removed here
  
  # Output out<-data.frame(tw,Ktauestind,Ktaupind) colnames(out)<-c('rolling
  # window','Kendall tau estimate','Kendall tau p-value')
  #out <- data.frame(Ktauestind)
  #Ktaupind <- data.frame(Ktaupind)
  #rownames(out) <- tw
  
  output = list(detrend="no",ktest=Ktauestind, ktp=Ktaupind, win=tw, ind=indicator)
  return(output)
} 





# calculate AR1 of all surrogates
KAR_sur = matrix(nrow=81200, ncol=6)
k=0
# load surrogates
allSurr5 = read.table("allSurr_bw2_25.txt", sep=" ")
for (indi in 1:6){
  j=1
  
  for (i in 1:50){
    end =  specTimeRowA[i]
    time = seq(1,end,by=1)
    
    for (n in 1:100){
      ts = data.frame(time, allSurr5[(1:end),(k+100*(i-1)+n)])
      Y = sensitivity_ews_surr(ts, incrwinsize = 5, indicator="ar1")
      tw = Y[[4]]
      KAR_sur[j:(j+length(tw)-1),indi] = Y[[2]]
      
      j=j+length(tw)
    }
  }
  k = k + 5000
}






#####################################################################################

# AUC 

library(sfsmisc)
roc_curve <- function(real, null){
  
  Dr = density(real, from=-1, to=1) 
  Dn = density(null, from=-1, to=1)
  
  thresholds = seq(-1,1, length.out=50)
  ROC = matrix(ncol=2,nrow=length(thresholds))
  
  for (t in 1:length(thresholds)){
    Ir = integrate.xy(Dr$x,Dr$y,thresholds[t],1)
    In = integrate.xy(Dn$x,Dn$y,thresholds[t],1)
    
    ROC[t,] = c(In,Ir)
  }
  return(ROC)
}

# example : calculate the AUC and plot the ROC for the skewness of the effective migration rate
real = read.table("./data_concatenated/Abrupt_KAR_sf30.txt", sep=" ")
null = read.table("./data_concatenated/NeutralA_KAR_sf30.txt", sep=" ", header=TRUE)
roc = roc_curve(real[,1],null[,1])
auc = integrate.xy(roc[,1], roc[,2])


# plot roc curve
roc = as.data.frame(roc)
colnames(roc) = c("x","y")
pp <- ggplot(data=roc, aes(x=x, y=y)) +
  geom_line()+
  theme(legend.position="none")+
  labs(title="AR1", x="false positive", y = "true positive")


