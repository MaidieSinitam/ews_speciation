# 
#
#
#
#
#
#
abrupt = seq(206321, 206370, by=1) #abrupt outputs number


# retreive time of speciation for all 50 simulation output
getTimeToSpec <- function(files){
  
  df=c()
  
  for (n in files){
    data = read.table(paste("./Run",n,"/parameters.m",sep=""), sep=";", comment.char="%")[1]
    ltrnr = as.integer(strsplit(strsplit(as.character(data$V1[62]), "=")[[1]][2], " ")[[1]][2])
    df=append(df,ltrnr)
  }
  
  return(df)
}
# row in the timeseries corresponding to speciation for all abrupt simulation outputs
# retreived with getTimeToSpecA()
specTimeRowA = c( 168, 227, 228, 112, 143, 181, 137, 189, 142, 155, 194, 134, 149, 193, 207, 212, 183, 96, 
                  167, 104, 186, 140, 155, 218, 167, 119, 101, 152, 182, 117, 148, 144, 125, 120, 169, 136, 
                  174, 136, 127, 168, 184, 181, 119, 131, 187, 146, 167, 205, 159, 183 )
specTimeRowA = specTimeRowA-3 #




# calculate the Fst per time step for one run
cleanFst <- function(run){
  data = read.table(paste("./Run",run,"/FSTtimeSeries.txt", sep=""), sep = " ", header = FALSE, dec =".")
  colnames(data)=c("time","locus","Fst","D","E", "F")
  mean_fst = c()
  var_fst = c()
  cv = c()
  sk = c()
  
  doublons = which(duplicated(data$time))   
  timestamp = data$time[-doublons]    #Creates a list of time stamps w/o doublons
  print(timestamp)
  for (T in timestamp){
    fst = which(data$time == T)  #Take all values of Fst for one time
    mean_fst = append(mean_fst, mean(data$Fst[fst])) #Save their mean
    var_fst = append(var_fst, var(data$Fst[fst]))
    cv = append(cv, sqrt(var(data$Fst[fst]))/mean(data$Fst[fst]))
    sk = append(sk,skewness(data$Fst[fst]))
  }
  
  out = data.frame(timestamp,mean_fst, var_fst, cv,sk)
  colnames(out) = c("timeindex", "MeanFST", "VarFST", "CvFST", "SkFST")
  return(out)
}

# calculate Wmax/Wres and Wres/Wimm for one run
fitnessRatios <- function(run,cut){
  data = read.table(paste("./Run",run,"/EffectiveMigrationRates.txt", sep=""), sep = " ", header = FALSE, dec = ".")
  colnames(data) = c("Time","m1","m2","Nloci1","Nres1","Nimm1","Wres1",
                     "Wimm1","Wmaxres1","Wminres1","Wrandom1","Wtot1","Nres2","Nimm2",
                     "Wres2","Wimm2","Wmaxres2","Wminres2","Wrandom2","Wtot2")
  if (missing(cut)){
    
  } else {
    cut_row = which(data$Time == cut)
    data = data[1:cut_row,]
  }
  
  WresWimm = c()
  WmaxWres = c()
  
  for (i in 1:nrow(data)){
    WresWimm = append(WresWimm, log(data[i,7]/data[i,8]))
    WmaxWres = append(WmaxWres, log(data[i,9]/data[i,7]))
  }
  
  out = data.frame(data$Time, WresWimm, WmaxWres)
  colnames(out) = c("timeindex","WresWimm", "WmaxWres")
  return(out)
}

# calculate mean fitness of residents for one run
randomFitness <- function(file){
  data = read.table( paste("./Run",file,"/FitnessTSdeme0.txt", sep=""), sep=" ")
  time = data$V1
  n = ncol(data)-6
  data = data[,2:n]
  meanfit = c()
  varfit = c()
  cvfit = c()
  skfit = c()
  
  for ( t in (1:nrow(data)) ) {
    T = as.numeric(data[t,])
    meanfit = append(meanfit, mean(T) )
    varfit = append(varfit, var(T) )
    cvfit = append(cvfit, sqrt(var(T))/mean(T)) 
    skfit = append(skfit, skewness(T) ) 
  }
  
  out = data.frame(time, meanfit, varfit, cvfit, skfit) 
  colnames(out) = c("time", "mean", "var", "cv", "sk")
  return(out)
}

# store values of fitness ratios in a file to gain time
df=matrix(ncol=201, nrow=404)
df[,1]=seq(1000,404000,by=1000)
i = 2
for (n in abrupt){
  T = fitnessRatios(n)
  resim = T$WresWimm
  maxres = T$WmaxWres
  L = 404-nrow(T)
  diff = rep(NA, each=L)
  resim = c(T$WresWimm,diff)
  maxres = c(T$WmaxWres,diff)
  df[,i] = resim
  df[,i+1] = maxres
  i = i+2
}
#write.table(df, "fit_ratios.txt", sep=" ", dec=".", row.names=FALSE, col.names=FALSE)


# get a data frame with a column for each metric for one run
getMetrics <-function(n, cut=FALSE){
  
  if (n %in% abrupt){
    ext = "A"
    N = which(abrupt == n)
    if (cut == FALSE){
      end = lenA[N]
      cut = lenA[N]
    } else {
      end = specTimeRowA[N]
      cut = specTimeRowA[N]
    }
    Nu = 2+(N-1)*4
    Nu3 = 2+(N-1)*2
    
  } else if (n %in% gradual) {
    ext = "G"
    N = which(gradual == n)
    if (cut == FALSE){
      end = lenG[N]
      cut = lenG[N]
    } else {
      end = timeToSpecRowG[N]
      cut = timeToSpecRowG[N]
    }
    Nu = 2+(N-1)*4
    Nu3 = 2+(N-1)*2
  }
  
  time = read.table(paste("./Run",n,"/EffectiveMigrationRates.txt", sep=""), sep=" ")[1:cut,1]
  # me
  f2 = read.table(paste("./Run",n,"/EffectiveMigrationRates.txt", sep=""), sep=" ")[1:cut,2]
  # wres/wimm
  f3 = read.table( paste("fitnessRatios",ext,".txt", sep=""), sep=" ")[1:end,Nu3]
  # wmax/wres
  f4 = read.table( paste("fitnessRatios",ext,".txt", sep=""), sep=" ")[1:end,Nu3+1]
  # mean fit all res
  f5 = read.table(paste("./Run",n,"/EffectiveMigrationRates.txt", sep=""), sep=" ")[1:cut,7]
  # mean fit 200 random indiv
  f6 = read.table( paste("randomFitness",ext,".txt", sep=""), sep=" ")[1:end,Nu]
  # Fst
  f7 = read.table( paste("cleanedFst",ext,".txt", sep=""), sep=" ")[1:end,Nu]
  #f7 = read.table( paste("Fst_sel",ext,".txt", sep=""), sep=" ")[1:end,N]
  # DXY
  f8 = read.csv(paste("./Run",n,"/DXYtimeSeries.csv", sep=""), sep=",", comment.char = "T", header=FALSE)[1:cut,2]
  f9 = read.csv(paste("./Run",n,"/DXYtimeSeries.csv", sep=""), sep=",", comment.char = "T", header=FALSE)[1:cut,3]
  f10 = read.csv(paste("./Run",n,"/DXYtimeSeries.csv", sep=""), sep=",", comment.char = "T", header=FALSE)[1:cut,4]
  # Correlation allelic states for sel sites
  f11 = read.table(paste("./Run",n,"/EffectiveMigrationRates.txt", sep=""), sep=" ")[1:cut,25]
  # Nb variable selected loci
  f12 = read.table(paste("./Run",n,"/NumberVariableLoci.txt", sep=""), sep=" ")[1:cut,2]
  f13 = read.table(paste("./Run",n,"/NumberVariableLoci.txt", sep=""), sep=" ")[1:cut,3]
  f14 = read.table(paste("./Run",n,"/NumberVariableLoci.txt", sep=""), sep=" ")[1:cut,4]
  
  #f8 = f8*f12
  
  out = data.frame(time, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14)
  colnames(out) = c("time", "me", "wreswimm", "wmaxwres", "fitall", "fitrand", "fst", "dxyall", "dxysel", "dxyneut", "ld", "nbloci", "nblocisel", "nblocineut")
  return(out)
}


# function to calculate the Spearman correlation between two metrics of the same simulation output
correlationtestA <-function(runs, cut = TRUE, data1 = c("me", "fst", "ld", "dxy", "nbloci", "maxres", "meanfit", "lociratio", "resimm"), 
                            data2 = c("me", "fst", "ld", "dxy", "nbloci", "maxres", "meanfit", "lociratio", "resimm")){
  
  df=c()
  
  metrics = c("time", "me", "wreswimm", "wmaxwres", "fitall", "fitrand", "fst", "dxyall", "dxysel", "dxyneut", "ld", "nbloci", "nblocisel", "nblocineut")
  for (n in runs){
    f1 <- getMetrics(n, cut=TRUE)[,which(metrics == data1)]
    f2 <- getMetrics(n, cut=TRUE)[,which(metrics == data2)]
    if (data1 == "lociratio"){
      f = getMetrics(n, cut=TRUE)[,which(metrics == "nbloci")]
      g = getMetrics(n, cut=TRUE)[,which(metrics == "nblocineut")]
      f1 =f/g
    }
    if (data2 == "lociratio"){
      f = getMetrics(n, cut=TRUE)[,which(metrics == "nbloci")]
      g = getMetrics(n, cut=TRUE)[,which(metrics == "nblocineut")]
      f2 =f/g
    }
    df = append(df, cor(f1,f2, method = "spearman"))
  }
  return(df)
}



A1 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "me", data2 = "fst")
A2 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "me", data2 = "nblocisel")
A3 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "me", data2 = "dxysel")
A4 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "me", data2 = "ld")
A5 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "me", data2 = "wreswimm")
A6 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "fst", data2 = "nblocisel")
A7 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "fst", data2 = "dxysel")
A8 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "fst", data2 = "ld")
A9 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "fst", data2 = "wreswimm")
A10 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "nblocisel", data2 = "dxysel")
A11 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "nblocisel", data2 = "ld")
A12 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "nblocisel", data2 = "wreswimm")
A13 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "dxysel", data2 = "ld")
A14 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "dxysel", data2 = "wreswimm")
A15 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "wreswimm", data2 = "ld")
A16 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "me", data2 = "fitrand")
A17 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "fst", data2 = "fitrand")
A18 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "nblocisel", data2 = "fitrand")
A19 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "dxysel", data2 = "fitrand")
A20 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "ld", data2 = "fitrand")
A21 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "wreswimm", data2 = "fitrand")
A22 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "lociratio", data2 = "me")
A23 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "lociratio", data2 = "nblocisel")
A24 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "lociratio", data2 = "dxysel")
A25 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "lociratio", data2 = "ld")
A26 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "lociratio", data2 = "wreswimm")
A27 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "lociratio", data2 = "fst")
A28 = correlationtestA(runs = abrupt, cut=TRUE, data1 = "lociratio", data2 = "fitrand")



df1 = data.frame(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12,
                 A13, A14, A15, A16, A17, A18, A19, A20, A21, A22, 
                 A23, A24, A25, A26, A27, A28)
colnames(df1) = c("me/fst" ,"me/nbloci" ,"me/dxy" ,"me/ld" ,"me/resimm","fst/nbloci", "fst/dxy" ,
                  "fst/ld", "fst/resimm" ,"nbloci/dxy" ,"nbloci/ld" ,"nbloci/resimm", "dxy/ld" ,"dxy/resimm",
                  "resimm/ld","me/meanfit" ,"fst/meanfit" ,"nbloci/meanfit", "dxy/meanfit" ,"ld/meanfit" ,"resimm/meanfit",
                  "lociratio/me" ,"lociratio/nbloci", "lociratio/dxy" ,"lociratio/ld" ,"lociratio/resimm" ,"lociratio/fst","lociratio/meanfit")

# mean Spearman correlation
mean.corr = apply(df1,2,mean)
sd.corr = apply(df1,2,sd)




# Spearman correlation with rolling window
correlMW <- function(n, cut = FALSE, data1 = c("me", "fst", "ld", "dxy", "nbloci", "maxres", "meanfit", "lociratio", "resimm"), data2 = c("me", "fst", "ld", "dxy", "nbloci", "maxres","meanfit", "lociratio", "resimm"), winsize=50){
  
  
  # the part after sets the following parameters
  # useful to get data from the stored files
  # end = length of timeseries
  # Nu and Nu3 are the nb of the column we extract from the files
  
  N = which(abrupt == n)
  Nu = 2+(N-1)*4
  Nu3=2*N
  if (cut == TRUE){
    end = specTimeRowA[N]
  } else { end = lenA[N] }
  
  timeindex = read.table(paste("./Run",n,"/NumberVariableLoci.txt", sep=""), sep=" ")$V1[1:end]
  #Preparing the data
  if (data1 == "nbloci"){
    Y = read.table(paste("./Run",n,"/NumberVariableLoci.txt", sep=""), sep=" ")[1:end,3]
  } else if (data1 == "me"){
    Y = read.table(paste("./Run",n,"/EffectiveMigrationRates.txt", sep=""), sep=" ")[1:end,2]
  } else if (data1 == "maxres"){
    Y = read.table("fitnessRatiosA.txt", sep=" ")[1:end,Nu3+1]
  } else if (data1 == "resimm"){
    Y = read.table("fitnessRatiosA.txt", sep=" ")[1:end,Nu3]
  } else if (data1 == "fst"){
    Y = read.table("cleanedFstA.txt", sep=" ", dec=".")[1:end,Nu]
  } else if (data1 == "dxy"){ 
    Y = read.csv(paste("./Run",n,"/DXYtimeSeries.csv", sep=""), sep=",", comment.char = "T", header=FALSE)[1:end,2]
  } else if (data1 == "ld") {
    Y = read.table(paste("./Run",n,"/EffectiveMigrationRates.txt", sep=""), sep=" ")[1:end,25]
  } else if (data1 == "meanfit"){
    Y = read.table("randomFitnessA.txt", sep=" ")[1:end,Nu]
  } else if (data1 == "lociratio"){
    A = read.table(paste("./Run",n,"/NumberVariableLoci.txt", sep=""), sep=" ")[1:end,]
    totA = A$V2
    neutA = A$V4
    Y = totA/neutA
  } else {
    return("Wrong data1")
  }
  
  if (data2 == "me"){
    Yd = read.table(paste("./Run",n,"/EffectiveMigrationRates.txt", sep=""), sep=" ")[1:end,2]
  } else if (data2 == "maxres"){
    Yd = read.table( "fitnessRatiosA.txt", sep=" ")[1:end,Nu3+1]
  } else if (data2 == "resimm"){
    Yd = read.table("fitnessRatiosA.txt", sep=" ")[1:end,Nu3]
  } else if (data2 == "fst"){
    Yd = read.table("cleanedFstA.txt", sep=" ", dec=".")[1:end,Nu]
  } else if (data2 == "dxy"){ 
    Yd = read.csv(paste("./Run",n,"/DXYtimeSeries.csv", sep=""), sep=",", comment.char = "T", header=FALSE)[1:end,2]
  } else if (data2 == "ld") {
    Yd = read.table(paste("./Run",n,"/EffectiveMigrationRates.txt", sep=""), sep=" ")[1:end,25]
  } else if (data2 == "nbloci"){
    Yd = read.table(paste("./Run",n,"/NumberVariableLoci.txt", sep=""), sep=" ")[1:end,3]
  } else if (data2 == "meanfit"){
    Yd = read.table("randomFitnessA.txt", sep=" ")[1:end,Nu]
  } else if (data2 == "lociratio"){
    A = read.table(paste("./Run",n,"/NumberVariableLoci.txt", sep=""), sep=" ")[1:end,]
    totA = A$V2
    neutA = A$V4
    Yd = totA/neutA
  } 
  
  
  
  # Rearrange data for indicator calculation
  mw <- round(length(Y) * winsize/100)
  omw <- length(Y) - mw + 1  ##number of moving windows
  low <- 2
  high <- mw
  nMR <- matrix(data = NA, nrow = mw, ncol = omw)
  nMRd <- matrix(data = NA, nrow = mw, ncol = omw)
  x1 <- 1:mw
  for (i in 1:omw) {
    Ytw <- Y[i:(i + mw - 1)]
    Ytwd <- Yd[i:(i + mw - 1)]
    nMR[, i] <- Ytw
    nMRd[, i] <- Ytwd
  }
  
  # Calculate indicator
  correl <- numeric()
  for (i in 1:ncol(nMR)) {
    correl[i] <- cor(nMR[, i], nMRd[, i], method="spearman")
  }
  
  out = data.frame(time = timeindex[mw:length(Y)], corr = correl)
  return(out)
}

# calculate moving correlations for each pair and each run
# one loop calculates the correlation for all pairs of one single run and stores it in a file
for (i in abrupt){
  #A0 = correlMW(i, cut=TRUE, data1 = "me", data2 = "nbloci", winsize=45)$time
  A1 = correlMW(i, cut=TRUE, data1 = "me", data2 = "fst", winsize=25)$corr
  A2 = correlMW(i, cut=TRUE, data1 = "me", data2 = "nbloci", winsize=25)$corr
  A3 = correlMW(i, cut=TRUE, data1 = "me", data2 = "dxy", winsize=25)$corr
  A4 = correlMW(i, cut=TRUE, data1 = "me", data2 = "ld", winsize=25)$corr
  A5 = correlMW(i, cut=TRUE, data1 = "me", data2 = "resimm", winsize=25)$corr
  
  A6 = correlMW(i, cut=TRUE, data1 = "fst", data2 = "nbloci", winsize=25)$corr
  A7 = correlMW(i, cut=TRUE, data1 = "fst", data2 = "dxy", winsize=25)$corr
  A8 = correlMW(i, cut=TRUE, data1 = "fst", data2 = "ld", winsize=25)$corr
  A9 = correlMW(i, cut=TRUE, data1 = "fst", data2 = "resimm", winsize=25)$corr
  
  A10 = correlMW(i, cut=TRUE, data1 = "nbloci", data2 = "dxy", winsize=25)$corr
  A11 = correlMW(i, cut=TRUE, data1 = "nbloci", data2 = "ld", winsize=25)$corr
  A12 = correlMW(i, cut=TRUE, data1 = "nbloci", data2 = "resimm", winsize=25)$corr
  
  A13 = correlMW(i, cut=TRUE, data1 = "dxy", data2 = "ld", winsize=25)$corr
  A14 = correlMW(i, cut=TRUE, data1 = "dxy", data2 = "resimm", winsize=25)$corr
  A15 = correlMW(i, cut=TRUE, data1 = "resimm", data2 = "ld", winsize=25)$corr
  A16 = correlMW(i, cut=TRUE, data1 = "me", data2 = "meanfit", winsize=25)$corr
  A17 = correlMW(i, cut=TRUE, data1 = "fst", data2 = "meanfit", winsize=25)$corr
  A18 = correlMW(i, cut=TRUE, data1 = "nbloci", data2 = "meanfit", winsize=25)$corr
  A19 = correlMW(i, cut=TRUE, data1 = "dxy", data2 = "meanfit", winsize=25)$corr
  A20 = correlMW(i, cut=TRUE, data1 = "ld", data2 = "meanfit", winsize=25)$corr
  A21 = correlMW(i, cut=TRUE, data1 = "resimm", data2 = "meanfit", winsize=25)$corr
  
  A22 = correlMW(i, cut=TRUE, data1 = "lociratio", data2 = "me", winsize=25)$corr
  A23 = correlMW(i, cut=TRUE, data1 = "lociratio", data2 = "fst", winsize=25)$corr
  A24 = correlMW(i, cut=TRUE, data1 = "lociratio", data2 = "nbloci", winsize=25)$corr
  A25 = correlMW(i, cut=TRUE, data1 = "lociratio", data2 = "dxy", winsize=25)$corr
  A26 = correlMW(i, cut=TRUE, data1 = "lociratio", data2 = "ld", winsize=25)$corr
  A27 = correlMW(i, cut=TRUE, data1 = "lociratio", data2 = "resimm", winsize=25)$corr
  A28 = correlMW(i, cut=TRUE, data1 = "lociratio", data2 = "meanfit", winsize=25)$corr
  
  A29 = correlMW(i, cut=TRUE, data1 = "maxres", data2 = "lociratio", winsize=25)$corr
  A30 = correlMW(i, cut=TRUE, data1 = "maxres", data2 = "me", winsize=25)$corr
  A31 = correlMW(i, cut=TRUE, data1 = "maxres", data2 = "fst", winsize=25)$corr
  A32 = correlMW(i, cut=TRUE, data1 = "maxres", data2 = "nbloci", winsize=25)$corr
  A33 = correlMW(i, cut=TRUE, data1 = "maxres", data2 = "dxy", winsize=25)$corr
  A34 = correlMW(i, cut=TRUE, data1 = "maxres", data2 = "ld", winsize=25)$corr
  A35 = correlMW(i, cut=TRUE, data1 = "maxres", data2 = "resimm", winsize=25)$corr
  A36 = correlMW(i, cut=TRUE, data1 = "maxres", data2 = "meanfit", winsize=25)$corr
  
  df1 = data.frame( A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, 
                    A12, A13, A14, A15, A16, A17, A18, A19, A20, A21, A22, A23, A24, A25, A26, A27, A28, A29, A30, A31, A32, A33, A34, A35, A36)
  write.table(df1, paste("Acormw_",i,"_cut25.txt", sep=""), sep=" ", dec=".", row.names=FALSE, col.names=FALSE)
  
}





# calculate Kendall tau of each moving correlation
allK = matrix(nrow=50,ncol=36)
j=1
for (n in abrupt){
  k=1
  for (i in 1:36){
    correl = read.table(paste("./Acormw_",n,"_cut25.txt", sep=""), sep=" " )[,i]
    timevec = seq(1, length(correl))
    Kt =  cor.test(timevec, correl, alternative = c("two.sided"), method = c("kendall"), 
                   conf.level = 0.95)
    allK[j,i] = round(Kt$estimate, digits = 3)
  }
  j=j+1
}






#########################################################################################

# plot Figure 3 : Correlation plots
library(ggplot2)
library(ggpubr)


#### upper panels : plot one metric against the other

T = getMetrics(206321, cut=TRUE)
col1 = c()
col2 = c()
col3 = c()
col4 = c()
col5 = c()

i=1
for (n in abrupt){
  T = getMetrics(n, cut=TRUE)
  col1 = append(col1, rep(i, each=specTimeRowA[i]))
  col2 = append(col2, T$dxysel)
  col3 = append(col3, T$ld)
  
  col4 = append(col4, T$nblocisel)
  col5 = append(col5, T$nbloci/T$nblocineut)
  i = i+1
}
out1 = data.frame(col1,col2,col3)
colnames(out1) = c("run", "x", "y")
out2 = data.frame(col1,col4,col5)
colnames(out2) = c("run", "x", "y")

A <- ggplot(out1, aes(x=x, y=y, fill=run))+
  geom_point(aes(color=run), size=0.5, shape=3)+
  labs(x="DXY", y="LD")+
  theme(legend.position="none")+
  scale_color_gradient2(midpoint=25, mid="white", low="blue",
                        high="red",  space = "Lab")+
  geom_smooth(method=loess, se=FALSE)

B <- ggplot(out2, aes(x=x, y=y, fill=run))+
  geom_point(aes(color=run), size=0.5, shape=3)+
  labs(x="nbselloci ", y="loci ratio")+
  theme(legend.position="none")+
  scale_color_gradient2(midpoint=25, mid="white", low="blue",
                        high="red",  space = "Lab")+
  geom_smooth(method=loess, se=FALSE)






##### Lower panel : moving pairwise Spearman correlations


# names of all the files to load them later
filesAc = c()
for (i in abrupt){
  filesAc = append(filesAc, paste("./Acormw_cut50/Acormw_",i,"_cut50.txt", sep=""))
}

# length of timeseries of the moving correlations
lenRW = c( 84, 113, 114, 56, 71, 90, 68, 94, 70, 77, 96, 66, 74, 96, 103, 106, 91, 48, 83,
           52, 92, 70, 77, 108, 83, 59, 50, 76, 90, 58, 74, 72, 62, 60, 84, 68, 86, 68, 63, 
           84, 92, 90, 59, 65, 93, 72, 83, 102, 79, 91)

col1 = c()
col2 = c()
col3 = c()
col4 = c()
col5 = c()
for (n in 1:50){
  T = read.table(filesAc[n], sep=" " )
  col1 = append(col1, rep(n, each=lenRW[n]) )
  x = seq(114000-(lenRW[n]-1)*1000,114000,by=1000) # adjust time vector to have the time to speciation
  col2 = append(col2, T$V25) #nloci vs loci ratio
  col3 = append(col3, x)
}
out3 = data.frame(col1,col2,col3)
colnames(out3) = c("run", "x", "y")

col1 = c()
col2 = c()
col3 = c()
col4 = c()
col5 = c()
for (n in 1:50){
  T = read.table(filesAc[n], sep=" " )
  col1 = append(col1, rep(n, each=lenRW[n]) )
  x = seq(114000-(lenRW[n]-1)*1000,114000,by=1000)
  col2 = append(col2, T$V14) #dxy vs ld
  col3 = append(col3, x)
}
out4 = data.frame(col1,col2,col3)
colnames(out4) = c("run", "x", "y")



C <- ggplot(out3, aes(x=y, y=x))+
  geom_point(aes(color=run), size=0.5, shape=3)+
  labs(x="Generations to speciation", y="Correlation")+
  theme(legend.position="none")+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(4000, 34000,74000,114000),
                     label=c("114000", "74000", "34000", "4000"))+
  ylim(0,1)+
  scale_color_gradient2(midpoint=25,  low="blue", mid="white",
                        high="red",  space = "Lab")

D <- ggplot(out4, aes(x=y, y=x))+
  geom_point(aes(color=run), size=0.5, shape=3)+
  labs(y="Correlation")+
  theme(legend.position="none")+
  scale_x_continuous(name ="Generations to speciation", 
                     breaks = c(4000, 34000,74000,114000),
                     label=c("114000", "74000", "34000", "4000"))+
  ylim(0,1)+
  scale_color_gradient2(midpoint=25, mid="white", low="blue",
                        high="red",  space = "Lab")


ggarrange(A,B,D,C,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)




###################################################################################################

# Fig 5

F = getMetrics(206321)
T = timeToSpecA[1]

F2 = getMetrics(206321, cut=TRUE)

A <- ggplot(data.frame( "x" = F$time, "y" = F$me), aes(x=x, y=y))+
  geom_line(color = "#CD5C5C")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="Effective migration rate", x="Generations", y="eff mig rate")
A2 <- ggplot(data.frame( "x" = F2$time, "y" = F2$me), aes(x=x, y=y))+
  geom_line(color = "#CD5C5C")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="  ", x="Generations", y="eff mig rate")

B <- ggplot(data.frame( "x" = F$time, "y" = F$fst), aes(x=x, y=y))+
  geom_line(color = "#FF8C00")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="Fst", x="Generations", y="Fst")
B2 <- ggplot(data.frame( "x" = F2$time, "y" = F2$fst), aes(x=x, y=y))+
  geom_line(color = "#FF8C00")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="  ", x="Generations", y="Fst")

C = ggplot(data.frame( "x" = F$time, "y" = F$fitrand), aes(x=x, y=y))+
  geom_line(color = "#BDB76B")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="Mean fitness", x="Generations", y="mean fitness")
C2 = ggplot(data.frame( "x" = F2$time, "y" = F2$fitrand), aes(x=x, y=y))+
  geom_line(color = "#BDB76B")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="  ", x="Generations", y="Mean fitness")+
  ylim(0,4)

D <- ggplot(data.frame( "x" = F$time, "y" = F$dxysel), aes(x=x, y=y))+
  geom_line(color = "#4B0082")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="DXY", x="Generations", y="DXY")
D2 <- ggplot(data.frame( "x" = F2$time, "y" = F2$dxysel), aes(x=x, y=y))+
  geom_line(color = "#4B0082")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="  ", x="Generations", y="DXY")

E <- ggplot(data.frame( "x" = F$time, "y" = F$ld), aes(x=x, y=y))+
  geom_line(color = "#2E8B57")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="LD", x="Generations", y="LD")
E2 <- ggplot(data.frame( "x" = F2$time, "y" = F2$ld), aes(x=x, y=y))+
  geom_line(color = "#2E8B57")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  labs(title="  ", x="Generations", y="LD")

G <- ggplot(data.frame( "x" = F$time, "y" = F$nblocisel), aes(x=x, y=y))+
  geom_line(color = "#1E90FF")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  labs(title="Nb of selected variable loci", x="Generations", y="Nb loci")+
  ylim(-50, 1000)
G2 <- ggplot(data.frame( "x" = F2$time, "y" = F2$nblocisel), aes(x=x, y=y))+
  geom_line(color = "#1E90FF")+
  geom_vline(xintercept = T, linetype="dotted", color = "#585858")+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  labs(title="  ", x="Generations", y="Nb loci")


ggarrange(A,A2,B,B2,
          labels = c("A", "", "B", ""),
          ncol = 2, nrow = 2)
ggarrange(C,C2,D,D2,
          labels = c("C", "", "C", ""),
          ncol = 2, nrow = 2)
ggarrange(E,E2,G,G2,
          labels = c("E", "", "F", ""),
          ncol = 2, nrow = 2)








