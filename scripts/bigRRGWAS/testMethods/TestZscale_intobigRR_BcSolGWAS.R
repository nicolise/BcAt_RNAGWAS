#Nicole E Soltis
#021018
#test Methods for BcAt_RNAGWAS: compare thresholds & lsmeans effect sizes from Tomato T4 bigRR to z-scaled effect sizes

#---------------------------------------------------------------------------
rm(list=ls())
library(tidyr); library(ggplot2)
#setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
setwd("~/Projects/BcAt_RNAGWAS/data/allreadsGWAS/testMethods/BO5tomato/03_bigRRout")
#Import data
mydat_perm <- read.csv("tomato_BO510_bigRR.csv")

my.files <- list.files(pattern = "B05.10")

#rename all files
for(i in 1:length(my.files)) {
  #read only top row
  my.file <- read.csv(my.files[i], nrows=1)
  my.name <- names(my.file)[3]
  file.rename(from=file.path(my.files[i]), to=file.path(paste(my.name,".csv",sep="")))
}

#now merge them
mydat_zscale <- read.csv(my.files[1])
my.files <- list.files(pattern = "HEM")
#rename all files
for(i in 2:length(my.files)) {
  #read only top row
  my.file <- read.csv(my.files[i])
  mydat_zscale <- cbind(mydat_zscale, my.file[,3])
  names(mydat_zscale)[i+2] <- names(my.file)[3]
}

#can also z-scale the effect estimates
# mydat <- I.plotdata[,c(1:3,16,4:15)]
# mydat_z <- mydat[,c(1:4)]
# #4:228
# for (i in c(5:16)){
#   mydat[,i+12] <- scale(mydat[,i], center = TRUE, scale = TRUE)
#   names(mydat)[i+12] <- paste(names(mydat)[i], "_z", sep="")
# }

setwd("~/Projects/BcAt_RNAGWAS/")
#extract thresholds mydat_perm
#first remove first 4 rows (threshold data)
mydat_perm <- mydat_perm[,-c(1)]
my_thresh <- mydat_perm[1:8,]
mydat_perm <- mydat_perm[-c(1:8),]
names(mydat_perm)[1] <- "Location"

mydat_zscale <- mydat_zscale[,-c(1)]
#order mydat_zscale to match mydat_perm
mydat_zscale <- mydat_zscale[,c("Domesticated.HEM","Wild.HEM","DmWoD.HEM","LA1547.HEM","LA1589.HEM","LA1684.HEM","LA2093.HEM","LA2176.HEM","LA2706.HEM","LA3008.HEM","LA3475.HEM","LA410.HEM","LA4345.HEM","LA4355.HEM","LA480.HEM")]
plotdat <- cbind(mydat_perm, mydat_zscale[,2:16])

setwd("~/Projects/BcAt_RNAGWAS/plots/testMethods/tomatoZscale/zscale_phenos")
#2:16
for (j in c(2:16)){
  my.95.thr.p <- my_thresh[1,j]
  my.95.thr.n <- my_thresh[5,j]
  my.99.thr.p <- my_thresh[3,j]
  my.99.thr.n <- my_thresh[7,j]
  my.999.thr.p <- my_thresh[4,j]
  my.999.thr.n <- my_thresh[8,j]
jpeg(paste("scatter_", names(mydat_perm)[j], ".jpg", sep=""), width=7.5, height=5, units='in', res=600)
plot(ggplot(plotdat, aes(x=plotdat[,j], y=plotdat[,j+16]))+
       theme_bw()+
       scale_x_continuous(name=paste("Effect size",names(plotdat[j])))+
       scale_y_continuous(name=paste("Z-scaled Effect size",names(plotdat[j])))+
       geom_point(alpha=1/8)+
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