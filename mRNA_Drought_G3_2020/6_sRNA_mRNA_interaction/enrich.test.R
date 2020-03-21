enrich.test <- function (gene2feature, gene2cat, gene2bias=NULL, test.feature=NULL, method="Sampling",
                         sampling=2000, pval.cutoff=0.01, enrich.cutoff=2, min.freq = 1, fdr.cutoff=NULL,
                         outsave=T, outpath=".", outfile=NULL, help=F) {
### this is the modified version of goseq
### only the sampling method is implemented in this version
### Sanzhen Liu
### 10/2/2014
### current version:
ver <- "v0.01"

  ###################################################################################################
  ### help information
  help <- function (x) {
    cat("Current version is", ver)
    cat("This script is designed to perform the enrichment test using weighted sampling approach\n")
    cat("that is similar to the method implemented in the package of GOSeq\n")
    cat("  gene2feature MUST be a two-column dataframe: 1st gene, 2nd features\n")
    cat("  gene2cat MUST a two-column dataframe: 1st gene, 2nd category\n")
    cat("  gene2bias is a two-column dataframe: 1st gene, 2nd bias data, if it is not NULL\n")
    cat("      if gene2bias is NULL, completely random method will be used\n")
    cat("  Currently only the sampling method is implemented\n")
  }
  ###################################################################################################
  message("current version is: ", ver)

###################################################################################################
goseqmakespline <- function (x, y, newX = NULL, nKnots = 6, lower_bound = 10^-3) {
  ### this function was copied from the package of goseq
  library(mgcv)  ### for the function of gam
  
  ww = order(x)
  size = ceiling(length(y)/10)
  low = sum(y[ww][1:size])
  hi = sum(y[ww][(length(y) - size):length(y)])
  if (hi <= low) {
    reflectionFactor = 10^10
    x = reflectionFactor - x
  }
  nKnots <- round(nKnots)
  if (is.null(newX)) 
    newX <- x
  if (hi > low) {
    x = c(0, x)
    y = c(0, y)
  }
  f.ug <- gam(y ~ s(x, k = nKnots, bs = "cr"), family = binomial())
  dat <- data.frame(x = x, y = y)
  sm <- smoothCon(s(x, k = nKnots, bs = "cr"), dat, knots = NULL)[[1]]
  if (length(sm$xp) < 6) 
    warning("Few than 6 nKnots were specified!\n Possible inaccurate fitting!\n")
  F <- mono.con(sm$xp, TRUE)
  F$A = rbind(F$A, c(1, rep(0, ncol(F$A) - 1)))
  F$b = c(F$b, 0)
  G <- list(X = sm$X, C = matrix(0, 0, 0), sp = f.ug$sp, p = sm$xp, 
            y = y, w = y * 0 + 1, Ain = F$A, bin = F$b, S = sm$S, 
            off = 0)
  if (any((G$Ain %*% G$p - G$bin) < 0)) {
    show(G$p)
    stop("Constraints for spline fit not satisfied by initial parameters")
  }
  p <- pcls(G)
  fv <- Predict.matrix(sm, data.frame(x = newX)) %*% p
  fv <- as.vector(fv)
  if (min(fv) < lower_bound) 
    fv = fv - min(fv) + lower_bound
  return(fv)
}
###################################################################################################

### check gene2feature
  if (!is.data.frame(gene2feature)) {
    stop("gene2feature was expected to be a dataframe of mapping categories to genes.\nCheck gene2feature input and try again.")
  } else {
    if (ncol(gene2feature) != 2) {
      stop("gene2feature was expected to be two columns of which gene on the 1st and cat on the 2nd")
    }
  }

### check gene2cat
  if (!is.data.frame(gene2cat)) {
    stop("gene2cat was expected to be a dataframe of mapping categories to genes.\nCheck gene2cat input and try again.")
  } else {
    if (ncol(gene2cat) != 2) {
      stop("gene2cat was expected to be two columns of which gene on the 1st and cat on the 2nd")
    }
  }

### check gene2bias
  nosplinefit <- 0
  if (is.null(gene2bias)) {
    nosplinefit <- 1
	gene2bias <- data.frame(Gene=gene2feature[, 1], Bias=rep(1, nrow(gene2feature)))
  } else {
    if (!is.data.frame(gene2bias)) {
      stop("gene2bias was expected to be a dataframe of mapping categories to genes.\nCheck gene2feature input and try again.")
    } else {
      if (ncol(gene2feature) != 2) {
        stop("gene2feature was expected to be two columns of which gene on the 1st and cat on the 2nd")
      }
    }
  }

### check method
  if (method != "Sampling") {
    stop("Invalid calculation method selected. Valid options are Sampling for this version.")
  }

  #####################################  
  ### remove NA data from a table
  rm.na.from.table <- function (x) {
    for (i in 1:ncol(x)) {
      x <- x[!is.na(x[,i]), ]
    }
    x
  }
  #####################################
  
  ### organizing input data  
  colnames(gene2feature) <- c("Gene", "Feature")
  gene2cat <- rm.na.from.table(gene2cat)
  colnames(gene2cat) <- c("Gene", "Cat")
  gene2cat <- gene2cat[gene2cat$Gene %in% gene2feature$Gene, ]
  gene2bias <- rm.na.from.table(gene2bias)
  colnames(gene2bias) <- c("Gene", "Bias")

  if (is.null(test.feature)) {
    allfeatures <- unique(gene2feature$Feature)
  	allfeatures <- allfeatures[!is.na(allfeatures)]
  } else {
    allfeatures <- test.feature
  }

  ### test for all the features
  out.list <- vector("list", length(allfeatures))
  for (eachfeature in allfeatures) {
    message("Analyzing ", eachfeature, " ...\n")
    gene2feature$Feature01 <- as.numeric(gene2feature$Feature == eachfeature)  ### 0 or 1
    featured.genes <- gene2feature$Gene[gene2feature$Feature01==1]
    nfeatures <- length(featured.genes)  ### number of featured genes
    featured.categories <- gene2cat$Cat[gene2cat$Gene %in% featured.genes]
    cat.stat <- as.data.frame(table(featured.categories))
    colnames(cat.stat) <- c("Cat", "Obs")
    
    gfb <- merge(gene2feature, gene2bias, by="Gene")
    
	if (nosplinefit) {
		gfb$Weight <- 1
	} else {
		gfb$Weight <- goseqmakespline(gfb$Bias, gfb$Feature01)
    }

    selected.cats <- as.character(cat.stat$Cat)
    
    ### simulations
    for (i in 1:sampling) {
      weighted.sampling.order <- order(runif(nrow(gfb))^(1/gfb$Weight), decreasing=T)
      sfeatured.genes <- gfb$Gene[weighted.sampling.order][1:nfeatures]  ### simulated featured genes
      ### category
      sfeatured.categories <- gene2cat$Cat[gene2cat$Gene %in% sfeatured.genes]
      
      sim.cat.stat <- table(sfeatured.categories)
      insect.sim.cat.stat <- as.character(sim.cat.stat[selected.cats])
      
	    old.colnames <- colnames(cat.stat)
      cat.stat <- cbind(cat.stat, insect.sim.cat.stat)
      colnames(cat.stat) <- c(old.colnames, i)
      
      ### report running status
      if ((i %% 100) == 0) {
        message(paste(i, "out of", sampling, "was finished\n"))
      }
    }
    
    ### pvalues for enrichment for each category
    cal.pval <- function(x) {
      all <- x[-1]
      all <- as.numeric(as.character(all))
      all[is.na(all)] <- 0
      obs <- all[1]
      p <- sum(all>=obs)/length(all)
      return(p)
    }
    
    ### calculate row means with NA=0
    cal.rowMeans.na <- function(x) {
      x <- as.numeric(as.character(x))
      x[is.na(x)] <- 0
      xmean <- mean(x)
      return(xmean)
    }
    
	######################
	#print(head(cat.stat[, 1:2]))
	######################

    ### determine pvalues
    all.pvals <- apply(cat.stat, 1, cal.pval)
    
	all.cats <- table(gene2cat$Cat)
    cat_ALL_gene_num <- all.cats[as.character(cat.stat$Cat)]
	cat_ALL_gene_num <- as.numeric(as.character(cat_ALL_gene_num)) 
 	######################
	#print(head(cat_ALL_gene_num))
	#print(length(cat_ALL_gene_num))
	#print(nrow(cat.stat))
	######################


    in_feature <- paste0("Cat_in_", eachfeature, "_gene_num")
    expected_in_feature <- paste0("Cat_expected_in_", eachfeature, "_gene_num")
    infeature.values <- cat.stat[, 2]
    expected_in_feature.values <- apply(cat.stat[, -1], 1, cal.rowMeans.na)
    expected_in_feature.values <- round(expected_in_feature.values, 1)
    out <- data.frame(cat.stat[, 1], cat_ALL_gene_num, infeature.values, expected_in_feature.values)
    colnames(out) <- c("Cat", "Cat_all_gene_num", in_feature, expected_in_feature)
    out$enrich <- round(out[, in_feature]/out[, expected_in_feature],2)
    out$pvals <- all.pvals
    out$qvals <- p.adjust(all.pvals, method="BH")
    
    if (is.null(fdr.cutoff)) {
      out.filter <- (out$pvals<=pval.cutoff & out$enrich>=enrich.cutoff & out$in_feature >= min.freq )
    } else {
      out.filter <- (out$pvals<=pval.cutoff & out$enrich>=enrich.cutoff & out[, in_feature] >= min.freq & out$qvals<=fdr.cutoff)
    }

    #print(head(out[out.filter, ]))
    out2 <- out[out.filter, ]
    #print(head(out2[order(out2$pvals, decreasing = 1), ]))
    
    out.list[[eachfeature]] <- out
    
    cat.db <- merge(gene2cat, gene2feature[gene2feature$Feature==eachfeature, ], by="Gene")
    colnames(cat.db) <- c("Gene", "Cat", "Feature", eachfeature)
    cat.db <- cat.db[, c("Gene", "Cat", "Feature")]
    cat.db2 <- merge(cat.db, out[out.filter, c("Cat", "enrich", "pvals")], by="Cat")

    ### output files
    if (outsave & nrow(cat.db2)>0) {
      if (is.null(outfile)) {
        outfile <-  paste0(eachfeature, ".sig.enriched.catogories.genes.txt")
      }
      write.table(cat.db2, paste0(outpath, "/", outfile), quote=F, row.names=F, sep="\t")
    }
  }
  
  ### output
  invisible(out.list)
}

