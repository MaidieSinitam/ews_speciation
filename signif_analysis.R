#####################################################################################
# Calculate the ROC curves and AUC 

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


#####################################################################################
# Plot Figures 6, 7 and 8

# load stored ROC curves
ar = read.table("ROC_AR.txt", sep=" ", header=TRUE)
sd = read.table("ROC_SD.txt", sep=" ", header=TRUE)
sk = read.table("ROC_SK.txt", sep=" ", header=TRUE)
kurt = read.table("ROC_KURT.txt", sep=" ", header=TRUE)

x = seq(0,1, by=0.5)
one = data.frame( "metric" = rep("zone", times=3), "x" = x, "y" = x)
ar = rbind(ar,one)
sd = rbind(one,sd)
kurt = rbind(kurt,one)
sk = rbind(sk,one)

# Plot ROC curves

AR = ggplot(data=ar, aes(x=x, y=y, group=metric, fill=metric)) +
  geom_line(aes(color=metric, linetype=metric))+
  scale_colour_manual(name = "AUC", labels = c("DXY", "Fst", "LD","Eff mig rate","Mean fitness", "Nb variable \n selected loci", "zone"), values = c("dxysel"="#9CCFE6","fst"="#0079B9", "ld"="#A5E27F", "me"="#399D19", "meanfit"="#F5BDCB", "nblocisel"="#C8181B", "zone" = "#000000"))+
  scale_linetype_manual(values=c("solid", "solid", "solid","solid","solid","solid","dotted"))+
  labs(title="AUTOCORRELATION", x="false positive rate", y = "true positive rate")

SD = ggplot(data=sd, aes(x=x, y=y, group=metric)) +
  geom_line(aes(color=metric, linetype=metric))+
  scale_colour_manual(name = "AUC", labels = c("DXY", "Fst", "LD","Eff mig rate","Mean fitness", "Nb variable \n selected loci", "zone"), values = c("dxysel"="#9CCFE6","fst"="#0079B9", "ld"="#A5E27F", "me"="#399D19", "meanfit"="#F5BDCB", "nblocisel"="#C8181B", "zone" = "#000000"))+
  scale_linetype_manual(values=c("dotted", "solid", "solid","solid","solid","solid","solid"))+
  labs(title="STANDARD DEVIATION", x="false positive rate", y = "true positive rate")


KURT = ggplot(data=kurt, aes(x=x, y=y, group=metric)) +
  geom_line(aes(color=metric, linetype=metric))+
  scale_colour_manual(name = "AUC", labels = c("DXY", "Fst", "LD","Eff mig rate","Mean fitness", "Nb variable \n selected loci", "zone"), values = c("dxysel"="#9CCFE6","fst"="#0079B9", "ld"="#A5E27F", "me"="#399D19", "meanfit"="#F5BDCB", "nblocisel"="#C8181B", "zone" = "#000000"))+
  scale_linetype_manual(values=c("solid", "solid", "solid","solid","solid","solid","dotted"))+
  labs(title="KURTOSIS", x="false positive rate", y = "true positive rate")
SK = ggplot(data=sk, aes(x=x, y=y, group=metric)) +
  geom_line(aes(color=metric, linetype=metric))+
  scale_colour_manual(name = "AUC", labels = c("DXY", "Fst", "LD","Eff mig rate","Mean fitness", "Nb variable \n selected loci", "zone"), values = c("dxysel"="#9CCFE6","fst"="#0079B9", "ld"="#A5E27F", "me"="#399D19", "meanfit"="#F5BDCB", "nblocisel"="#C8181B", "zone" = "#000000"))+
  scale_linetype_manual(values=c("solid", "solid", "solid","solid","solid","solid","dotted"))+
  labs(title="SKEWNESS", x="false positive rate", y = "true positive rate")

ggarrange(AR,SD,KURT,SK,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 2)

# Plot violins

#Load stored Kendall tau
# Each file should contain the Kendall tau trends for all replicates of one dataset type
# for all differentiation metrics and all bandwidths and rolling window sizes
# The loaded files should have one column for each differentiation metric and one line for each Kendall tau

#KAR = read.table("__.txt", sep=" ")
#KSD = read.table("__.txt",, sep=" ")
#KSK = read.table("__.txt",, sep=" ")
#KKURT = read.table("__.txt",, sep=" ")

#KAR2 = read.table("__.txt",, sep=" ")
#KSD2 = read.table("__.txt",, sep=" ")
#KSK2 = read.table("__.txt",, sep=" ")
#KKURT2 = read.table("__.txt",, sep=" ")


old.par <- par(mfrow=c(2,2), mar=c(3,3,3,3))
vioplot(KAR[,1], KAR[,2], KAR[,3], KAR[,4], KAR[,5], KAR[,6], names=metrics, col="#F0AA56", plotCentre = "line", side="left", ylim=c(-1,1))
vioplot(KAR2[,1], KAR2[,2], KAR2[,3], KAR2[,4], KAR2[,5], KAR2[,6], names=metrics, col="#C4335C", plotCentre = "line", add=TRUE, side="right", ylim=c(-1,1))
title("Autocorrelation", cex=0.5)

vioplot(KSD[,1], KSD[,2], KSD[,3], KSD[,4], KSD[,5], KSD[,6], names=metrics, col="#F0AA56", plotCentre = "line", side="left", ylim=c(-1,1))
vioplot(KSD2[,1], KSD2[,2], KSD2[,3], KSD2[,4], KSD2[,5], KSD2[,6], names=metrics, col="#C4335C", plotCentre = "line", side="right",add=TRUE, ylim=c(-1,1))
title("Standard deviation", cex=0.5)

vioplot(KKURT[,1], KKURT[,2], KKURT[,3], KKURT[,4], KKURT[,5], KKURT[,6], names=metrics, col="#F0AA56", plotCentre = "line", side="left", ylim=c(-1,1))
vioplot(KKURT2[,1], KKURT2[,2], KKURT2[,3], KKURT2[,4], KKURT2[,5], KKURT2[,6], names=metrics, col="#C4335C", plotCentre = "line", side="right",add=TRUE, ylim=c(-1,1))
title("Kurtosis", cex=0.5)

vioplot(KSK[,1], KSK[,2], KSK[,3], KSK[,4], KSK[,5], KSK[,6], names=metrics, col="#F0AA56", plotCentre = "line", side="left", ylim=c(-1,1))
vioplot(KSK2[,1], KSK2[,2], KSK2[,3], KSK2[,4], KSK2[,5], KSK2[,6], names=metrics, col="#C4335C", plotCentre = "line", side="right", add=TRUE, ylim=c(-1,1))
title("Skewness", cex=0.5)
par(old.par)


old.par <- par(mfrow=c(2,2), mar=c(3,3,3,3))
vioplot(KAR[,1], KAR[,2], KAR[,3], KAR[,4], KAR[,5], KAR[,6], names=metrics, col="#F0AA56", plotCentre = "line", side="left", ylim=c(-1,1))
vioplot(KAR2[,1], c(0,0), KAR2[,3], KAR2[,4], KAR2[,5], KAR2[,6], names=metrics, col="#A8CC6E", plotCentre = "line", add=TRUE, side="right", ylim=c(-1,1))
title("Autocorrelation", cex=0.5)

vioplot(KSD[,1], KSD[,2], KSD[,3], KSD[,4], KSD[,5], KSD[,6], names=metrics, col="#F0AA56", plotCentre = "line", side="left", ylim=c(-1,1))
vioplot(KSD2[,1], c(0,0), KSD2[,3], KSD2[,4], KSD2[,5], KSD2[,6], names=metrics, col="#A8CC6E", plotCentre = "line", side="right",add=TRUE, ylim=c(-1,1))
title("Standard deviation", cex=0.5)

vioplot(KKURT[,1], KKURT[,2], KKURT[,3], KKURT[,4], KKURT[,5], KKURT[,6], names=metrics, col="#F0AA56", plotCentre = "line", side="left", ylim=c(-1,1))
vioplot(KKURT2[,1], c(0,0), KKURT2[,3], KKURT2[,4], KKURT2[,5], KKURT2[,6], names=metrics, col="#A8CC6E", plotCentre = "line", side="right",add=TRUE, ylim=c(-1,1))
title("Kurtosis", cex=0.5)

vioplot(KSK[,1], KSK[,2], KSK[,3], KSK[,4], KSK[,5], KSK[,6], names=metrics, col="#F0AA56", plotCentre = "line", side="left", ylim=c(-1,1))
vioplot(KSK2[,1], c(0,0), KSK2[,3], KSK2[,4], KSK2[,5], KSK2[,6], names=metrics, col="#A8CC6E", plotCentre = "line", side="right", add=TRUE, ylim=c(-1,1))
title("Skewness", cex=0.5)
par(old.par)


