################################################################################
### WGCNA analysis
### Erliang Zeng
### University of South Dakota
### 07/25/2017
###############################################################################

## Read and pre-process DS data

library(WGCNA)
options(stringsAsFactors = FALSE)
sdata<-read.table("DS.counts.txt", header=T, sep="\t", row.names=1)*1
datExpr <- t(sdata)

## Determine appropriate power for DS network construction

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## Build DS network and identify modules

net = blockwiseModules(datExpr, power = 9, minModuleSize = 10, reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE,
                       saveTOMFileBase = "SigGene_padj_0001", verbose = 5)
TOM <- TOMsimilarityFromExpr(datExpr, power = 9)
dissTOM = 1-TOM
colors1 = labels2colors(net$colors)
moduleLabels = net$colors

## Save DS network for Cytoscape visualization. 

cyt = exportNetworkToCytoscape(TOM, edgeFile = paste("CytoscapeInput-edges- SnRNA_Drought_stress_threshold.05", ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-SnRNA_Drought_stress_threshold.05", ".txt", sep=""),
                               weighted = TRUE, threshold = 0.05,nodeNames = rownames(dat1), nodeAttr = colors1)

## Read and pre-process WW data

wdata<-read.table("WW.counts.txt", header=T, sep="\t", row.names=1)*1
datExpr <- t(wdata)

## Determine appropriate power for DS network construction

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## Build WW network and identify modules

net = blockwiseModules(datExpr0, power = 9, minModuleSize = 10, reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "SigGene_padj_0001", verbose = 5)
TOM <- TOMsimilarityFromExpr(datExpr, power = 9)
dissTOM = 1-TOM
colors2 = labels2colors(net$colors)
moduleLabels = net$colors

## Save WW network for Cytoscape visualization. 

cyt = exportNetworkToCytoscape(TOM, edgeFile = paste("CytoscapeInput-edges-Interaction_SnRNA_Water_File_threshold.05", ".txt", sep=""), nodeFile = paste("CytoscapeInput-nodes-Interaction_SnRNA_Water_File_threshold.05", ".txt", sep=""), weighted = TRUE, threshold = 0.05,nodeNames = rownames(dat2), nodeAttr = colors2)

## Preparation for network preservation analysis. 

setLabels = c("color1", "color2")
multiExpr = list(Drought = list(data = t(sdata)), Water = list(data = t(wdata)))
multiColor = list(Drought = colors1,Water=colors2)
system.time( { mp = modulePreservation(multiExpr, multiColor,referenceNetworks = c(1:2), loadPermutedStatistics = FALSE, nPermutations = 200, verbose = 3)} )
save(mp, file = "SnRNAModulePreservation.RData")
load(file = "SnRNAModulePreservation.RData")

## Network preservation analysis (WW modules against DS network). 

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
plotMods = !(modColors %in% c("grey", "gold"));
text = modColors[plotMods];
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], 
mp$preservation$Z[[ref]][[test]][, 2])
mains = c("Preservation Median rank", "Preservation Zsummary");
sizeGrWindow(10, 5);
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))

for (p in 1:2) {
 
  min = min(plotData[, p], na.rm = TRUE); max = max(plotData[, p], na.rm = TRUE);
  if (p==2){
     if (min -max/10) min = -max/10;
     ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)) 
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min));
    plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,main = mains[p], cex = 2.4, ylab = mains[p], xlab = "Module size", log = "x", ylim = ylim, xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4);
    labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
 
 if (p==2) {
    abline(h=0);
    abline(h=2, col = "blue", lty = 2);
    abline(h=10, col = "darkgreen", lty = 2);
 }
}

## Network preservation analysis (DS modules against WW network). 

ref = 2
test = 1
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])


modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
plotMods = !(modColors %in% c("grey", "gold"));
text = modColors[plotMods];
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
mains = c("Preservation Median rank", "Preservation Zsummary");
sizeGrWindow(10, 5);
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))

for (p in 1:2) {
 
  min = min(plotData[, p], na.rm = TRUE); max = max(plotData[, p], na.rm = TRUE);
  if (p==2){
     if (min -max/10) min = -max/10;
     ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)) 
  } else
     ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min));
     plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,main = mains[p], cex = 2.4, ylab = mains[p], xlab = "Module size", log = "x", ylim = ylim, xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4);labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
 
 if (p==2) {
    abline(h=0);
    abline(h=2, col = "blue", lty = 2);
    abline(h=10, col = "darkgreen", lty = 2);
 }
}
