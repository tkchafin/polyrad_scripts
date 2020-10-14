library(VariantAnnotation)
library(ggplot2)
library(polyRAD)

myVCF <- "polyrad.vcf"
mybg <- bgzip(myVCF)
indexTabix(mybg, format = "vcf")

#possible ploidies: 1(2) diploid; 2(4) autotetraploid, 3(2,2) allotetraploid
mydata <- VCF2RADdata(myVCF, possiblePloidies = list(2, 4, c(2,2)),
                      expectedLoci = 50000, expectedAlleles = 12000,
                      min.ind.with.reads=48)

myhindhe <- HindHe(mydata)
myhindheByLoc <- colMeans(myhindhe, na.rm = TRUE)

pdf("hindheByLoc.pdf")
hist(myhindheByLoc, col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by locus")
abline(v = 0.5, col = "blue", lwd = 2)
dev.off()

#NOTE: Manually identify peak below 0.5 in the above histogram
peak <- 0.35
coeff<- InbreedingFromHindHe(peak, ploidy = 2)

#simulate what 'good' data might look like if diploid & with mild overdispersion
exp<- ExpectedHindHe(mydata, inbreeding = coeff, ploidy = 2)

#get 95% lower bound, representing loci w/ problems 
lb <- mean(exp, na.rm=T)-2*(sd(exp, na.rm=T))

#subset loci 
keeploci <- names(myhindheByLoc)[myhindheByLoc >= lb]
mydata <- SubsetByLocus(mydata, keeploci)

mydataPopStruct <- IteratePopStruct(mydata, nPcsInit = 8, tol = 5e-03,
                                    overdispersion = 9)

pdf("pca.pdf")
myallele <- 1
freqcol <- heat.colors(101)[round(mydataPopStruct$alleleFreqByTaxa[,myallele] * 100) + 1]
plot(mydataPopStruct, pch = 21, bg = freqcol)
dev.off()

#looking at inheritance mode fits
pdf("dip_auto.pdf")
plot(mydataPopStruct$ploidyChiSq[1,], mydataPopStruct$ploidyChiSq[2,], 
     xlab = "Chi-squared for diploid model",
     ylab = "Chi-squared for autotetraploid model", log = "xy")
abline(a = 0, b = 1, col = "blue", lwd = 2)
dev.off()
pdf("dip_allo.pdf")
plot(mydataPopStruct$ploidyChiSq[1,], mydataPopStruct$ploidyChiSq[3,], 
     xlab = "Chi-squared for diploid model",
     ylab = "Chi-squared for allotetraploid model", log = "xy")
abline(a = 0, b = 1, col = "blue", lwd = 2)
dev.off()
pdf("auto_allo.pdf")
plot(mydataPopStruct$ploidyChiSq[2,], mydataPopStruct$ploidyChiSq[3,], 
     xlab = "Chi-squared for autotetraploid model",
     ylab = "Chi-squared for allotetraploid model", log = "xy")
abline(a = 0, b = 1, col = "blue", lwd = 2)
dev.off()

#
pdf("dip_auto_chi.pdf")
myChiSqRat <- mydataPopStruct$ploidyChiSq[1,] / mydataPopStruct$ploidyChiSq[2,]
myChiSqRat <- tapply(myChiSqRat, mydataPopStruct$alleles2loc, mean)
allelesPerLoc <- as.vector(table(mydataPopStruct$alleles2loc))
ggplot(mapping = aes(x = myhindheByLoc[GetLoci(mydata)], y = myChiSqRat, fill = as.factor(allelesPerLoc))) +
  geom_point(shape = 21, size = 3) +
  labs(x = "Hind/He", y = "Ratio of Chi-squared values, diploid to autotetraploid",
       fill = "Alleles per locus") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0.5) +
  scale_fill_brewer(palette = "YlOrRd")
dev.off()
pdf("dip_allo_chi.pdf")
myChiSqRat <- mydataPopStruct$ploidyChiSq[1,] / mydataPopStruct$ploidyChiSq[3,]
myChiSqRat <- tapply(myChiSqRat, mydataPopStruct$alleles2loc, mean)
allelesPerLoc <- as.vector(table(mydataPopStruct$alleles2loc))
ggplot(mapping = aes(x = myhindheByLoc[GetLoci(mydata)], y = myChiSqRat, fill = as.factor(allelesPerLoc))) +
  geom_point(shape = 21, size = 3) +
  labs(x = "Hind/He", y = "Ratio of Chi-squared values, diploid to allootetraploid",
       fill = "Alleles per locus") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0.5) +
  scale_fill_brewer(palette = "YlOrRd")
dev.off()
pdf("auto_allo_chi.pdf")
myChiSqRat <- mydataPopStruct$ploidyChiSq[2,] / mydataPopStruct$ploidyChiSq[3,]
myChiSqRat <- tapply(myChiSqRat, mydataPopStruct$alleles2loc, mean)
allelesPerLoc <- as.vector(table(mydataPopStruct$alleles2loc))
ggplot(mapping = aes(x = myhindheByLoc[GetLoci(mydata)], y = myChiSqRat, fill = as.factor(allelesPerLoc))) +
  geom_point(shape = 21, size = 3) +
  labs(x = "Hind/He", y = "Ratio of Chi-squared values, autotetraploid to autotetraploid",
       fill = "Alleles per locus") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0.5) +
  scale_fill_brewer(palette = "YlOrRd")
dev.off()

myHindHe <- HindHe(mydata)
TotDepthT <- rowSums(mydata$locDepth)
myHindHeByInd <- rowMeans(myHindHe, na.rm = TRUE)

pdf("ind_hindhe.pdf")
ggplot(data.frame(Depth = TotDepthT, HindHe = myHindHeByInd), 
       mapping = aes(x = Depth, y = HindHe)) +
  geom_point() +
  scale_x_log10() +
  labs(x = "Read Depth", y = "Hind/He")
dev.off()

#the above can be used to filter out individuals

pdf("loc_hindhe.pdf")
myHindHeByLoc <- colMeans(myHindHe, na.rm = TRUE)
hist(myHindHeByLoc, breaks = 50, xlab = "Hind/He",
     main = "Distribution of Hind/He among loci",
     col = "lightgrey")
abline(v = 0.5, col = "blue", lwd = 2)
dev.off()

#filter out bad individuals
ind_thresh<-mean(TotDepthT)-2*sd(TotDepthT)
print(paste0("ind threshold:", ind_thresh))
goodinds <- TotDepthT[TotDepthT >= mean(TotDepthT)-2*sd(TotDepthT)]

mydata <- SubsetByTaxa(mydata, names(goodinds))

#filter out loci above 0.5 HindHe
goodloci50 <- colnames(myHindHe)[myHindHeByLoc < 0.5]
goodloci75 <- colnames(myHindHe)[myHindHeByLoc < 0.75]

mydata_75 <- SubsetByLocus(mydata, goodloci50)
mydata_50 <- SubsetByLocus(mydata, goodloci75)

RADdata2VCF(mydata_75, file = "polyrad_out_75.vcf")
RADdata2VCF(mydata_50, file = "polyrad_out_50.vcf")




