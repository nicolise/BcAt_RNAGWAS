#Michelle Stitzer Script

library(data.table)


## read in the positions of each SNP
snppos=read.table('/home/berdeja/germination_gwas/fastphase/clump1_d11_4C_todo.fastPHASE.positions', skip=1)

## set up a data frame to fill in with haplotype assignments
clusters=data.frame(snppos=snppos)

for (i in 1:1135){
  a=read.table(paste0('/home/berdeja/germination_gwas/fastphase/fastphase_E-Cluster-memberships.txt_', i))
  map=sapply(1:nrow(a), function(x) which.max(a[x,]))  ## find max haplotype assignment
  clusters[,paste0('E', i)]=map  ## add to output data frame

}

dim(clusters)

# get ecotype IDs
eco=fread('/home/berdeja/germination_gwas/fastphase/clump1_d11_4C_todo.vcf.recode.vcf', skip=7, nrows=1)
eco=as.character((eco[,10:ncol(eco)]))
eco=paste0('E', eco)
colnames(clusters)=c('snppos', eco)

write.table(clusters, 'fastphase_clusters.ecotypeCol.snpRow.txt', row.names=F, col.names=T, sep='\t', quote=F)


## sort by focal snp (i'm guessing on position here!)

## can pick the specific row that is actually the focal snp!
clusters.sorted=clusters[,order(clusters[2500,],clusters[2000,],clusters[3000,])]  ## note that this puts the snp position at the end, but it's still called V1

c15=palette(rainbow(15)) 

png('fastphase_haplos.alejandra.png', 1000,1000)  ## png here because there are a lot of points
plot(clusters.sorted[,'snppos'], seq(0,1140,length.out=length(clusters.sorted[,1])), col='white')
for(i in 1:1135){
#  points(clusters.sorted[,'snppos'], rep(i, length(clusters.sorted[,1])), col=c15[clusters.sorted[,i]], cex=0.1)
  segments(x0=c(0,clusters.sorted[-nrow(clusters.sorted),'snppos']), y0=rep(i, length(clusters.sorted[,1])),x1=clusters.sorted[,'snppos'], y1=rep(i, length(clusters.sorted[,1])), col=c15[clusters.sorted[,i]], cex=0.1)
  
  }
dev.off()

