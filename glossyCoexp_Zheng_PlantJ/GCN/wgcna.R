### INITIAL WGCNA PACKAGE INSTALL CODE ###
WD = "/data1/home/liu3zhen/B73coexpNetwork/genes/0-WGCNA.sanzhen"
setwd(WD)
options(stringsAsFactors = FALSE)

#source("http://bioconductor.org/biocLite.R")
#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
#install.packages("WGCNA")
library(WGCNA)

### FORMATTING DATA ###
genedata.raw = read.csv("../fromDan/Network_Initial_Data/genedata_M4_CV1.GCN.subset.csv")
genedata = genedata.raw
genedata = genedata[, -c(1:4)]
rownames(genedata) = genedata.raw[, 1]
genedata = t(genedata)
dim(genedata)
### NEW GENE NETWORKS WITH VARIANT PARAMETERS ###
net = blockwiseModules(genedata, power=5,
                               TOMType = "unsigned", minModuleSize = 75,
                               reassignThreshold = 0, mergeCutHeight = 0.1,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs = FALSE,verbose = 3, deepSplit = 2)

### CREATING COLOR ASSIGNMENT DENDROGRAM ###
sizeGrWindow(12, 9)
moduleColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

### TOM CALCULATION AND PLOTTING
TOM = TOMsimilarityFromExpr(genedata, power = 5)
dissTOM = 1 - TOM
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
plotTOM = dissTOM^6
diag(plotTOM) = NA
TOMplot(plotTOM, geneTree, moduleColors, main = "Network Heatmap Plot")

### EXPORTING CYTOSCAPE FILES ###
moduleColors = labels2colors(net$colors)
cytdata = exportNetworkToCytoscape(TOM, edgeFile = paste("CytoscapeInput-edges-P5M75", ".txt", sep = ""),
                                   nodeFile = paste("CytoscapeInput-nodes-P5M75", ".txt", sep = ""),
                                   weighted = TRUE, threshold = 0.1681, nodeNames = colnames(genedata),
                                   nodeAttr = moduleColors)

