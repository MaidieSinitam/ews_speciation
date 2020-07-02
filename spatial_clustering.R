#
#
#
# Spatial clustering of variable loci
# This code calculates Ripley's K for variable loci
#
# The code relies on FSTtimeSeries.txt files to get the spatial position of loci
#


library(ggplot2)
library(ggpubr)





# calculate null expectation
ripleyK_nullexp <- function(Ln, incr,n_null = 1000){
  # calculates null expectation and 2.5% and 97.5% quantiles # this function is called by ripleyKt() below
  # Ln = nb of variable sites
  # incr = increment of step for the abscissa vector
  # n_null = number of null datasets to produce
  
  # output = data frame with the following columns (null expectation, 2.5% quantile, 97.5% quantile, standard deviation)
  
  lambda = Ln/25
  den = (Ln - 1)*lambda
  
  data = matrix(runif(Ln*n_null, min=0, max=25),Ln)
  
  ix=seq(1,n_null,by=1)
  
  test=lapply(ix,function(x,data){
    K=c()
    D = dist(data[,x])
    for (i in seq(0,25,by=incr)){
      Dx = length(which(D < i))
      Kx = Dx/den
      K = append(K, Kx)
    }
    return(K)
  },data=data)
  
  dfres2=as.data.frame(do.call(rbind,test))
  
  nullexp = colMeans(dfres2)
  sd = lapply(dfres2, var)
  sd = as.data.frame(sd)
  sd = lapply(sd,sqrt)
  sd = as.data.frame(sd)
  sd = t(sd)
  
  Q = as.data.frame(do.call(rbind,lapply((dfres2), quantile, probs=c(2.5,97.5)/100)))
  lowQ =  Q[1]
  highQ = Q[2]
  
  out = data.frame(NUL = (nullexp), LOW = lowQ, HIGH = highQ, SD = sd)
  return(out) }

# calculate Ripley's K for one chromosome at one time step on one run
ripleyKt <- function(run.nb, time, chr, incr){
  # run.nb = nb of the simulation output to analyse (ex. 206321)
  # time = timestep
  # chr = chromosome nb
  # incr = increment of step for the abscissa vector
  
  # output = data frame with the following columns (null expectation, 2.5% quantile,
  # 97.5% quantile, Ripley's K)
  K = c()
  
  data = read.table( paste("./Run",run.nb,"/FSTtimeSeries.txt", sep=""), sep=" ") # change directory if needed
  
  t = which(data$V1 == time)
  data = data[t,]
  
  chrow = which(data$V7 == chr)
  if (length(chrow) == 0){return(paste("no var loci for chromosome", chr, sep=" "))}
  data = data[chrow,]
  Ln = nrow(data)
  lambda = Ln/25
  den = (Ln - 1)*lambda
  
  Xvect = seq(0,25,by=incr)
  D = dist(data$V8)
  
  for (x in Xvect){
    Dl = length(which(D < x))
    Kx = Dl/den
    K = append(K, Kx)
  }
  
  NL = ripleyK_nullexp(Ln, incr, n_null=1000)
  out = data.frame(X = Xvect, nul = NL$NUL, low = NL$X2.5., high = NL$X97.5., k = K)

  return(out)
}






### plot ripley's K for one run, chromosome 0, at 4 different time steps

T = ripleyKt(206370,180000,0,1)
col3 = c()
col2 = c()
col1 = rep( c("null", "K", "Q 97.5%", "Q 2.5% "), each= length(T$X))
col3 = T$nul
col3 = append(col3, T$k)
col3 = append(col3, T$high)
col3 = append(col3, T$low)
col2 = rep(T$X, times=4)
out7 = data.frame( col1,col2,col3 )
colnames(out7) = c("type", "x", "y")
A <- ggplot(out7, aes(x=x, y=y, fill=type))+
  theme(legend.position="none")+
  geom_line(aes(color=type), size=0.5)+
  scale_color_manual(values=c( "#FF0000",'#000000', "#696969", "#696969"))+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  labs(x="x", y="K")


T = ripleyKt(206370,181000,0,1)
col3 = c()
col2 = c()
col1 = rep( c("null", "K", "Q 97.5%", "Q 2.5% "), each= length(T$X))
col3 = T$nul
col3 = append(col3, T$k)
col3 = append(col3, T$high)
col3 = append(col3, T$low)
col2 = rep(T$X, times=4)
out8 = data.frame( col1,col2,col3 )
colnames(out8) = c("type", "x", "y")
B <- ggplot(out8, aes(x=x, y=y, fill=type))+
  geom_line(aes(color=type), size=0.5)+
  theme(legend.position="none")+
  scale_color_manual(values=c( "#FF0000",'#000000', "#696969", "#696969"))+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  labs(x="x", y="K")


T = ripleyKt(206370,182000,0,1)
col3 = c()
col2 = c()
col1 = rep( c("null", "K", "Q 97.5%", "Q 2.5% "), each= length(T$X))
col3 = T$nul
col3 = append(col3, T$k)
col3 = append(col3, T$high)
col3 = append(col3, T$low)
col2 = rep(T$X, times=4)
out9 = data.frame( col1,col2,col3 )
colnames(out9) = c("type", "x", "y")
C <- ggplot(out9, aes(x=x, y=y, fill=type))+
  geom_line(aes(color=type), size=0.5)+
  theme(legend.position="none")+
  scale_color_manual(values=c( "#FF0000",'#000000', "#696969", "#696969"))+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  labs(x="x", y="K")

T = ripleyKt(206370,183000,0,1)
col3 = c()
col2 = c()
col1 = rep( c("null", "K", "Q 97.5%", "Q 2.5% "), each= length(T$X))
col3 = T$nul
col3 = append(col3, T$k)
col3 = append(col3, T$high)
col3 = append(col3, T$low)
col2 = rep(T$X, times=4)
out10 = data.frame( col1,col2,col3 )
colnames(out10) = c("type", "x", "y")
D <- ggplot(out10, aes(x=x, y=y, fill=type))+
  geom_line(aes(color=type), size=0.5)+
  theme(legend.position="none")+
  scale_color_manual(values=c( "#FF0000",'#000000', "#696969", "#696969"))+
  #  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotdash", "dashed"))+
  labs(x="x", y="K")

ggarrange(A,B,D,C,
          labels = c("t=18000", "t=18100", "t=18200", "t=18300"),
          ncol = 2, nrow = 2)


