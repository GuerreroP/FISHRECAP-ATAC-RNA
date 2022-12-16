##### BiocManager package installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Mfuzz")


###### Load libraries
library(Mfuzz)
library("Biobase")

###### Import data. It is advisable to remove genes that don't vary enough.
final<-as.matrix(read.table("normcount.csv" , header=TRUE, sep=";", row.names=1, as.is=TRUE))
final2<-log(final+1)
minimalSet <- ExpressionSet(assayData=final2)
filteredSet <- filter.NA(minimalSet, thres=0.5)
filledSet <- fill.NA(filteredSet)
pdf("tmp.pdf")
tmp <- filter.std(filledSet,min.std=0)
dev.off()
standardisedSet <- standardise(filledSet)

##### Calculate fuzziness parameter:

m1 <- mestimate(standardisedSet)

##### Decide cluster numbers.
pdf("cselection.pdf")
tmp3<-cselection(standardisedSet, m=m1, crange=seq(2,40,2), repeats=5, visu=TRUE)
dev.off()
tmp3<-cselection(standardisedSet, m=m1, crange=seq(2,50,2), repeats=1, visu=TRUE)
pdf("Dmin.pdf")
tmp_4  <- Dmin(standardisedSet,m=m1,crange=seq(10,40,10),repeats=30,visu=TRUE)
dev.off()
tmp_4
pdf("Dmin.2.pdf")
tmp_4  <- Dmin(standardisedSet,m=m1,crange=seq(5,40,5),repeats=30,visu=TRUE)
dev.off()

ncluster <- 30

##### Mfuzz cluster plot based on decided cluster number
cl <- mfuzz(standardisedSet,c=ncluster,m=m1)
mfuzz.plot(standardisedSet,cl=cl,mfrow=c(4,5))
pdf("file.pdf", height=10, width=10)
dev.off()

####Save output files

write.table( cl$centers, file="Center.txt", sep="\t")
write.table( cl$membership, file="Membership.txt", sep="\t")
write.table( cl$size, file="Size.txt", sep="\t")
write.table( cl$cluster, file="Cluster.txt", sep="\t")
savehistory("History.txt")
