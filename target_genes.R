## Script para determinar los genes dianas de un factor de transcripción
## a partir del fichero narrowPeak generado por MaCS2.

## Instalar chipseeker y paquete de anotación de Arabidopsis thaliana

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("clusterProfiler")
BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
BiocManager::install("org.At.tair.db")

library(ChIPseeker)
library(ChIPpeakAnno)
library(clusterProfiler)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(org.At.tair.db)

txdb <- TxDb.Athaliana.BioMart.plantsmart28
annotation.atha<-org.At.tair.db


###########################################
## En ChIPpeakAnno:

## Leer fichero de picos

? toGRanges

?system.file

bed <- system.file("extdata", "MACS_output.bed", package="ChIPpeakAnno")
gr1 <- toGRanges(bed, format="narrowPeak", header=FALSE) 


###########################################
## En ChIPseeker:

## Leer fichero de picos
prr5.peaks <- readPeakFile(peakfile = "prr5_summits.bed",header=FALSE)

covplot(prr5.peaks,weightCol = "V5")

## Definir la región que se considera promotor entorno al TSS

promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, 
                         downstream=1000)

## Anotación de los picos
prr5.peakAnno <- annotatePeak(peak = prr5.peaks, 
                             tssRegion=c(-1000, 1000),
                             TxDb=txdb)

plotAnnoPie(prr5.peakAnno,ndigit=0.001)
plotDistToTSS(prr5.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

## Convertir la anotación a data frame
prr5.annotation <- as.data.frame(prr5.peakAnno)

promoter.prr5.annotation <- subset(prr5.annotation,annotation=="Promoter")
regulome.prr5 <- promoter.prr5.annotation$geneId

target.genes <- prr5.annotation$geneId[prr5.annotation$annotation == "Promoter"]
write(x = target.genes,file = "prr5_target_genes.txt")


genes.arabidopsis <- as.data.frame(genes(txdb))
genes.arabidopsis.chr1 <- subset(genes.arabidopsis, seqnames==1)
genes.arabidopsis.chr1.names <- rownames(genes.arabidopsis.chr1)

## Enriquecimiento de terminos de GO

library(clusterProfiler)
ego <- enrichGO(gene          = regulome.prr5,
                universe      = genes.arabidopsis.chr1.names,
                OrgDb         = org.At.tair.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType       = "TAIR")

head(ego)

## Visualizacion de enriquecimiento de terminos de GO

library(DOSE)
library(enrichplot)

barplot(ego, showCategory=15,main="GO Enrichment Analysis")
dotplot(ego, showCategory=20)

## Cluster-profiling de términos de GO

edox <- setReadable(ego, 'org.At.tair.db', 'TAIR')
p1 <- cnetplot(edox, foldChange=geneList)

## Análisis de las palabras reconocidas por el factor
## de transcripción mediante HOMER
## (http://homer.ucsd.edu/homer/ngs/peakMotifs.html)


