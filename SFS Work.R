setwd("~/Documents/GitHub/Ecological-Genomics/myresults/ANGSD")
list.files()

SFS <- scan("XDS_outFold.sfs")
sumSFS <- sum(SFS)
sumSFS

pctPoly = 100*(1-(SFS[1]/sumSFS))
plotSFS <- SFS[-c(1,length(SFS))]
barplot(plotSFS,main="Folded SFS XDS",xlab="minor allele frequency", ylab="number of sites",col = "orange")

div <- read.table("XDS_folded_allsites.thetas.idx.pestPG")
colnames(div) = c("window","chrname","wincenter","tW","tP","tF","tH","tL","tajD","fulif","fuliD","fayH","zengsE","numSites")
print(div)

div$tWpersite= div$tW/div$numSites
div$tPpersite= div$tP/div$numSites

pdf("XDS_diversity_stats.pdf")
par(mfrow=c(2,2))

hist(div$tWpersite,col="gray",xlab="Theta-W",main="")
hist(div$tPpersite,col="gray",xlab="Theta-Pi",main="")
hist(div$tajD,col="gray",xlab="Tajima's D",main="")
barplot(plotSFS)
dev.off()

summary(div)


