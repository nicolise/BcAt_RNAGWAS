#Nicole E Soltis
#021018
#test Methods for BcAt_RNAGWAS: compare thresholds & lsmeans effect sizes from Tomato T4 bigRR to z-scaled effect sizes

#---------------------------------------------------------------------------
rm(list=ls())
library(tidyr); library(ggplot2)
#setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
setwd("~/Projects/BcSolGWAS/data/GWAS_files/04_bigRRoutput/trueMAF20_20NA/")
#Import data
I.plotdata <- read.csv("SlBc_12plants_trueMAF20_20NA.HEM.PlotFormat.csv") 
I.thresh <- read.csv("SlBc_12plants_trueMAF20_20NA.HEM.Thresh.csv")
I.plotdata <- I.plotdata[,-c(1,17,18,19)]
#could add domestication traits
#D.plotdata <- read.csv()
#D.thresh
#mydat <- cbind(I.plotdata,D.plotdata)
I.thresh <- I.thresh[,-c(1)]
#now z-scale the phenotypes
mydat <- I.plotdata[,c(1:3,16,4:15)]
mydat_z <- mydat[,c(1:4)]
#4:228
for (i in c(5:16)){
  mydat[,i+12] <- scale(mydat[,i], center = TRUE, scale = TRUE)
  names(mydat)[i+12] <- paste(names(mydat)[i], "_z", sep="")
}

setwd("~/Projects/BcAt_RNAGWAS/")
#5:16
for (j in c(5:16)){
  my.95.thr.p <- I.thresh[1,j-3]
  my.95.thr.n <- I.thresh[5,j-3]
  my.99.thr.p <- I.thresh[3,j-3]
  my.99.thr.n <- I.thresh[7,j-3]
  my.999.thr.p <- I.thresh[4,j-3]
  my.999.thr.n <- I.thresh[8,j-3]
jpeg(paste("plots/testMethods/tomatoZscale/scatter_", names(mydat)[j], ".jpg", sep=""), width=7.5, height=5, units='in', res=600)
plot(ggplot(mydat, aes(x=mydat[,j], y=mydat[,j+12]))+
       theme_bw()+
       scale_x_continuous(name=paste("Effect size",names(mydat[j])))+
       scale_y_continuous(name=paste("Z-scaled Effect size",names(mydat[j])))+
       geom_point(alpha=1/2)+
       geom_vline(xintercept=my.95.thr.p, lty=2)+
       geom_vline(xintercept=my.95.thr.n, lty=2)+
       geom_vline(xintercept=my.99.thr.p, lty=1)+
       geom_vline(xintercept=my.99.thr.n, lty=1)+
       geom_vline(xintercept=my.999.thr.p, lty=3)+
       geom_vline(xintercept=my.999.thr.n, lty=3))
dev.off()
}

#then: get list of z-scaled values at thresholds
#linear regression of each trait
#empty data frame
my.linreg <- mydat[1,5:16]
row.names(my.linreg)[1] <- "Slope"
for (y in c(5:16)){
  #format y~x
  my.linreg.reg <- lm(formula = mydat[,y+12] ~ mydat[,y], data=mydat)
  #slope
  my.linreg[1,y-4] <- coef(my.linreg.reg)[2]
  #then extract y value at x = threshold
  #y = mx + b
  my.y.95.p <- coef(my.linreg.reg)[2]*I.thresh[1,y-3]+coef(my.linreg.reg)[1]
  my.linreg[2,y-4] <- my.y.95.p[[1]]
  my.y.95.n <- coef(my.linreg.reg)[2]*I.thresh[5,y-3]+coef(my.linreg.reg)[1]
  my.linreg[3,y-4] <- my.y.95.n[[1]]
  my.y.99.p <- coef(my.linreg.reg)[2]*I.thresh[3,y-3]+coef(my.linreg.reg)[1]
  my.linreg[4,y-4] <- my.y.99.p[[1]]
  my.y.99.n <- coef(my.linreg.reg)[2]*I.thresh[7,y-3]+coef(my.linreg.reg)[1]
  my.linreg[5,y-4] <- my.y.99.n[[1]]
  my.y.999.p <- coef(my.linreg.reg)[2]*I.thresh[4,y-3]+coef(my.linreg.reg)[1]
  my.linreg[6,y-4] <- my.y.999.p[[1]]
  my.y.999.n <- coef(my.linreg.reg)[2]*I.thresh[8,y-3]+coef(my.linreg.reg)[1]
  my.linreg[7,y-4] <- my.y.999.n[[1]]
}
row.names(my.linreg)[2] <- "z_95.thr.p"
row.names(my.linreg)[3] <- "z_95.thr.n"
row.names(my.linreg)[4] <- "z_99.thr.p"
row.names(my.linreg)[5] <- "z_99.thr.n"
row.names(my.linreg)[6] <- "z_999.thr.p"
row.names(my.linreg)[7] <- "z_999.thr.n"

hist(as.numeric(my.linreg[7,]),breaks=10)
write.csv(my.linreg, "data/allreadsGWAS/testMethods/BcTomato_zscale_permutThrs.csv")