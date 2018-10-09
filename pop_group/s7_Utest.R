library(Matrix)
readSC <- function(cellLine){
  geCounts <- readMM(paste0(cellLine,"/GRCh38/matrix.mtx"))
  gNames <- read.delim(paste0(cellLine,"/GRCh38/genes.tsv"), header = FALSE, stringsAsFactors = FALSE)
  cellSum <- apply(geCounts,2,sum)
  geCountsN <- ceiling(t(t(geCounts)*(mean(cellSum)/cellSum)))
  #geCountsN <- geCountsN[,!colSums(geCountsN) %in% boxplot.stats(colSums(geCountsN))$out]
  selectedGenes <- apply(geCountsN,1,sum) > 0
  geCountsN <- geCountsN[selectedGenes,]
  gNames <- gNames[selectedGenes,]
  out <- list(sc = geCountsN, names = gNames)
  return(out)
}
GM18502 <- readSC("GM18502")
mix <- readSC("mix/")
GM12878 <- readSC("GM12878")


shared <-  GM12878$names[,1] %in% GM18502$names[,1]
GM12878$sc <- GM12878$sc[shared,]
GM12878$names <- GM12878$names[shared,]

shared <- GM18502$names[,1] %in% GM12878$names[,1]
GM18502$sc <- GM18502$sc[shared,]
GM18502$names <- GM18502$names[shared,]

shared <- mix$names[,1] %in% GM12878$names[,1]
mix$sc <- mix$sc[shared,]
mix$names <- mix$names[shared,]

GM12878$sc <- as.matrix(GM12878$sc)
GM18502$sc <- as.matrix(GM18502$sc)
mix$sc <- as.matrix(mix$sc)

library(NMF)
scM <- apply(mix$sc,1,var)
scM <- scM > quantile(scM,0.60)
sum(scM)
out <- nmf(as.matrix(mix$sc[scM,]), rank = 2, .options='v3')
H <- out@fit@H
H <- round(t(t(H)/apply(H,2,sum)),2)
H[H == 0.5] <- 0
H <- round(H)
EUR <- which(H[1,] == 1)
AFR <- which(H[2,] == 1)

x <- c(length(EUR), length(AFR))
x / sum(x)
# 
# pGM12878 <- read.csv("scGM12878NB_Parameters.tsv", sep = "\t", stringsAsFactors = FALSE)
# pGM18502 <- read.csv("scGM18502NB_Parameters.tsv", sep = "\t", stringsAsFactors = FALSE)
# pMix<- read.csv("scmixNB_Parameters.tsv", sep = "\t", stringsAsFactors = FALSE)
# pMix$var <- (pMix$mu + pMix$mu^2/pMix$size)
# pMix$ff <- pMix$var / pMix$mu
# 
# #s <- order(pMix$ff, decreasing = TRUE)[1:1000]
# #dMix <- mix$sc[mix$names[,2] %in% pMix$GENE[s],]
# dMix <- sapply(1:1000,function(x){apply(mix$sc[,sample(seq_len(ncol(mix$sc)),100)],1,mean)})
# dMix <- dMix[apply(dMix==0,1,mean) < 1,]
# out <- nmf(as.matrix(dMix),2, .options='v3')
# H <- round(t(t(out@fit@H)/apply(out@fit@H,2,sum)))
# apply(H,1,mean)
# 
# 
pGM12878 <- pGM12878[pGM12878$ENSEMBL %in% pGM18502$ENSEMBL,]
pGM18502 <- pGM18502[pGM18502$ENSEMBL %in% pGM12878$ENSEMBL,]
rownames(pGM12878) <- pGM12878$ENSEMBL
rownames(pGM18502) <- pGM18502$ENSEMBL
markers <- apply(cbind(pGM12878$mu,pGM18502$mu),1,function(x){log2(max(x)/min(x)) > 5})
# sum(markers)
markerData <- mix$sc[mix$names[,1] %in% pGM18502$ENSEMBL[markers],]
rownames(markerData) <- mix$names[mix$names[,1] %in% pGM18502$ENSEMBL[markers],2]
cellSum <- colSums(markerData)
markerData <- t(t(markerData)/cellSum) * mean(cellSum)

library(NMF)
out <- nmf(as.matrix(markerData), rank = 2, .options='v3')
H <- out@fit@H
H <- round(t(t(H)/apply(H,2,sum)),2)
H[H == 0.5] <- 0
H <- round(H)
which(grepl("IGLC2", rownames(markerData)))
EUR <- which(H[1,] == 1)
AFR <- which(H[2,] == 1)
mean(markerData[108,EUR])
mean(markerData[108,AFR])
cellID <- numeric(ncol(markerData))
cellID[EUR] <- "EUR"
cellID[AFR] <- "AFR"
cellID <- cbind(seq_len(ncol(markerData)),cellID)
write.table(cellID, quote = FALSE, row.names = FALSE, col.names = FALSE, file = "cellLineage.tsv", sep = "\t")
# write.table(x = markerData, file = "mix_U.tsv", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
# 
# # install.packages("NMF")
# acc <- sapply(1:1000,function(x){
#   initialP <- (round(runif(n = 1, min = 0, max = 50)))/100
#   knowProportion <- cbind(GM12878$sc[,sample(seq_len(ncol(GM12878$sc)), 1000*initialP)],GM18502$sc[,sample(seq_len(ncol(GM18502$sc)), 1000*(1-initialP))])
#   cellSums <- apply(knowProportion,2,sum)
#   knowProportion <- (t(t(knowProportion)/cellSums)) * mean(cellSums)
#   sMarkers <- GM12878$names[,1] %in% pGM18502$ENSEMBL[markers]
#   knowProportion <- knowProportion[sMarkers,]
#   nonEmpty <- apply(knowProportion,1,sum) != 0
#   sMarkers <- sMarkers[nonEmpty]
#   knowProportion <- knowProportion[nonEmpty,]
#   
#   
#   library(NMF)
#   out <- nmf(knowProportion, rank = 2, .options='v3')
#   H <- out@fit@H
#   H <- round(t(t(H)/apply(H,2,sum)),2)
#   H[H > 0.5] <- 1
#   H[H < 0.5] <- 0
#   c(min(apply(H,1,mean)),min(c(initialP,(1-initialP))))
# })
# acc <- round(t(acc),2)[,c(2,1)]
# plot(acc, xlab="Expected", ylab="Observed")
# abline(lm(acc[,2]~acc[,1]), col = "red", lwd = 2)
