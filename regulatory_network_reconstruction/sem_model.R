library(Matrix)
library(lavaan)
library(bnlearn)
library(semPlot)

# Reading data
GM12878 <- readMM("GM12878/GRCh38/matrix.mtx")
rownames(GM12878) <- read.csv("GM12878/GRCh38/genes.tsv", sep = "\t", 
                              header = FALSE, stringsAsFactors = FALSE)[,2]
colnames(GM12878) <- readLines("GM12878/GRCh38/barcodes.tsv")

# Filtering genes
sGenes <- c("REL","RELA","PRDM1","AICDA","PAX5","IRF4","BACH2","BCL6")
gValues <- as.matrix(t(GM12878[sGenes,]))
gValues <- as.data.frame(scale(gValues))
rownames(gValues) <- NULL

# Correlation matrix
cValues <- cor(gValues)

# Definition of the model
model <- c(
  "IRF4 ~ PRDM1 + RELA + REL",
  "PAX5 ~ PRDM1",
  "BACH2 ~ PAX5 + REL + PRDM1",
  "BCL6 ~ PAX5 + REL + PRDM1",
  "REL ~ PRDM1",
  "PRDM1 ~ BACH2 + BCL6 + IRF4 + RELA",
  "AICDA ~ REL + RELA + PAX5 + IRF4"
)

# Fitting the model
out <- lavaan(model = model, data = data.frame(gValues), std.lv=TRUE, auto.var=TRUE, auto.cov.lv.x=TRUE, se="bootstrap")

# Result
summary(out, fit.measures = TRUE)

# Figure
png("oModel.png", width = 1800, height = 1800, res = 300)
semPaths(out,"std","est", layout = "circle",whatLabels = "std",style = "lisrel", residuals = FALSE,nCharNodes = 5)
dev.off()
