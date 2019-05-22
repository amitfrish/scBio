#dependencies
#library(foreach)
#library(doSNOW)
#library(parallel)
#library(raster)
#library(fields)
#library(limma)
#library(LiblineaR)
#library(sp)

########## create sub-reference for the deconvolution
#' @keywords internal
createSpecificRef <- function(currRefference, modelSize, neighborhoodSize, genePercents, chosenCells, chosenNeigCells){
  specificRef = do.call(cbind,lapply(chosenCells,function(chosenCell){
    currRefference[,chosenNeigCells[[chosenCell]]]
  }))
  colnames(specificRef) = unlist(lapply(chosenCells,function(cell){
    rep(paste(colnames(currRefference)[cell],cell,sep="_"),neighborhoodSize)
  }))
  row.names(specificRef) = row.names(currRefference)
  specificRef = specificRef[sample(1:dim(currRefference)[1],round(genePercents*dim(currRefference)[1])),]
  list(ref = specificRef, chosenCells = chosenCells)
}

########## Remove cell duplications from reference
#' @keywords internal
createNoCellDupeReference <- function(refference){
  currCellNames = colnames(refference)
  refferenceNoDups = as.data.frame(matrix(0,nrow = dim(refference)[1], ncol=length(unique(currCellNames))))
  for (i in 1:length(unique(currCellNames))){
    if (length(which(currCellNames == unique(currCellNames)[i]))!=1){
      refferenceNoDups[,i] = rowMeans(refference[,which(currCellNames == unique(currCellNames)[i])])
    }else{
      refferenceNoDups[,i] = refference[,which(currCellNames == unique(currCellNames)[i])]
    }
  }
  row.names(refferenceNoDups) = row.names(refference)
  colnames(refferenceNoDups) = unique(currCellNames)
  return(refferenceNoDups)
}

########## Score genes using ANOVA
#' @keywords internal
GeneBasedAnova <- function(specificRefference, nonZeroRatio = NULL){
  cellNamesForAnova = colnames(specificRefference)
  #genes_to_take = row.names(specificRefference)[!(grepl("Rpl", row.names(specificRefference)) | grepl("Rps", row.names(specificRefference)) | grepl("mt-", row.names(specificRefference)) | grepl("Gm[0-9]", row.names(specificRefference)))]
  genes_to_take = row.names(specificRefference)
  genes_to_take = names(which(rowMeans(specificRefference[genes_to_take,])!=0))
  if(!is.null(nonZeroRatio)){
    genes_to_take = unlist(apply(as.data.frame(genes_to_take),1,function(gene){
      if((length(which(as.numeric(specificRefference[gene,which(cellNamesForAnova==1)])!=0))/length(as.numeric(specificRefference[gene,which(cellNamesForAnova==1)])))>nonZeroRatio){
        gene
      }else{
        NULL
      }
    }))
  }

  dat = cbind(rep(0,length(genes_to_take)),specificRefference[genes_to_take,])
  group = c("",cellNamesForAnova)
  #dat = specificRefference[genes_to_take,]
  #group = c(cellNamesForAnova)
  #print(cellNamesForAnova)
  dmat <- stats::model.matrix(~ group)
  # fit the ANOVA model
  fit <- limma::lmFit(dat, dmat)
  fit = fit[,-1]
  fit <- limma::eBayes(fit)
  fitF = fit$F
  res = as.data.frame(cbind(gsub("group","",colnames(fit$coefficients)[apply(fit$coefficients,1,function(x){order(x,decreasing = T)[1]})]),fitF))
  colnames(res) = c("group", "score")

  listOfGenes = apply(as.data.frame(unique(cellNamesForAnova)),1,function(cellGroup){
    selectedIndexes = which(as.character(res$group)==as.character(cellGroup))
    (genes_to_take[selectedIndexes])[order(res$score[selectedIndexes],decreasing = T)]
  })
  return(listOfGenes)
}

########## Select final gene sets for reference using the kappa function
#' @keywords internal
selectGenesUsingKappa <- function(refferenceNoDups, allGenes){
  bestKappa = Inf
  bestG = 0
  mul = 1
  maxNumberOfGenesPerCell = 50
  bestGenes = c()
  indexRange = 2:maxNumberOfGenesPerCell
  for (i in indexRange){
    selectedGenes = unique(as.character(unlist(lapply(1:length(allGenes), function(listIndex){
      unlist(allGenes[listIndex])[which(!is.na(unlist(allGenes[listIndex])[1:as.numeric(i*mul)]))]
    }))))
    currRefferenceNoDups = refferenceNoDups[match(selectedGenes, row.names(refferenceNoDups)),]
    newKappa = kappa(currRefferenceNoDups)
    if (newKappa<bestKappa){
      bestKappa = newKappa
      bestG = i
      bestGenes = unique(selectedGenes)
    }
  }
  finalRefference = refferenceNoDups[which(row.names(refferenceNoDups) %in% bestGenes),]

  return(list(reference = finalRefference, G = bestG, kappa = bestKappa))
}

########## Check the variance of none zero genes
#' @keywords internal
checkVariableGenes = function(a, ratio) {
  count_nonZeros = length(which(a > min(a)))
  if (count_nonZeros/length(a) > ratio) {
    var(a)/ mean(a)
  } else {
    0
  }
}

########## Run svr based deconvolution
#' @keywords internal
runLibLinear = function(ref_matrix, sample_matrix){
  X <- data.matrix(ref_matrix)
  X <- X[apply(X,1,sd)!=0,]
  Y <- data.matrix(sample_matrix)
  Y = Y[match(row.names(X),row.names(Y)),]

  #intersect genes
  Y <- Y[row.names(Y) %in% row.names(X),]
  X <- X[row.names(X) %in% row.names(Y),]

  #standardize sig matrix
  #X <- (X - mean(X)) / sd(as.vector(X))
  X <- t(apply(X,1,function(rowX){
     (rowX - mean(rowX)) / sd(as.vector(rowX))
  }))

  C = LiblineaR::heuristicC(X)

  #iterate through mixtures
  predictionMatrix = do.call(rbind,lapply(1:dim(Y)[2],function(index){
    y <- Y[,index]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    (LiblineaR::LiblineaR(data = X, target = y, type = 11, cost = C)$W)[1:dim(X)[2]]
  }))
  colnames(predictionMatrix) = colnames(X)
  row.names(predictionMatrix) = colnames(Y)
  predictionMatrix
}

########## Select cells for each run and caculate the desired number of runs
#' @keywords internal
choseCellsForRuns = function(XY, refNames, modelSize, minSelection, neighborhoodSize){
  #### cell selection
  k = floor(modelSize/length(unique(refNames)))
  if(k==0){
    k=1
  }

  minValueToReduceTo = 10^-10

  initialGrids = lapply(unique(refNames), function(currCluster){
    clusterIndexes = which(refNames==currCluster)
    nbins = max(k,length(clusterIndexes)/neighborhoodSize)
    if(is.null(dim(XY))){
      currXY = XY[clusterIndexes]
      breaks = seq(min(currXY)-10^-7,max(currXY)+10^-7, (max(currXY)-min(currXY)+2*10^-7)/ceiling(nbins))
      grid <- rep(NA,ceiling(nbins))
      cellLocationOnGrid = rep(NA,length(currXY))
      for(currBreakIndex in 2:length(breaks)){
        cellLocationOnGrid[which(currXY>breaks[currBreakIndex-1] & currXY<breaks[currBreakIndex])] = currBreakIndex-1
      }
      tab <- table(cellLocationOnGrid)
      grid[as.numeric(names(tab))] <- tab
    }else{
      currXY = XY[clusterIndexes,]
      ch <- grDevices::chull(currXY)
      coords <- currXY[c(ch, ch[1]), ]
      # ch = geometry::convhulln(currXY[,1:3])
      # coords = t(apply(ch,1,function(currCh){
      #   unlist(lapply(1:length(currCh),function(currIndex){
      #     currXY[currCh[currIndex],currIndex]
      #   }))
      # }))
      poly = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), "x")))
      grid <- raster::raster(raster::extent(poly), nrows = ceiling(sqrt(nbins)), ncols= ceiling(sqrt(nbins)))
      sp::proj4string(grid)<-sp::proj4string(poly)

      cellLocationOnGrid = raster::cellFromXY(grid, currXY)
      tab <- table(cellLocationOnGrid)
      grid[as.numeric(names(tab))] <- tab
    }
    list(grid = grid, clusterIndexes = clusterIndexes, cellLocationOnGrid = cellLocationOnGrid, nFullbins = length(tab), maxBinSize = max(tab))
  })

  numOfRuns = ceiling(minSelection*max(unlist(lapply(initialGrids,function(clusterData){ clusterData$nFullbins * clusterData$maxBinSize  }))) / k)

  meanDistMatrix = rep(1,length(refNames))
  chosenCellList = lapply(1:numOfRuns, function(runNum){
    chosenCells = as.numeric(unlist(lapply(unique(refNames),function(currCluster){
      initialGrid = initialGrids[[which(unique(refNames)==currCluster)]]
      clusterIndexes = initialGrid$clusterIndexes
      grid = initialGrid$grid
      cellLocationOnGrid = initialGrid$cellLocationOnGrid
      kToUse = k
      if(k>length(which(!is.na(grid[])))){
        kToUse = length(which(!is.na(grid[])))
      }
      gridCellsToUse = sample(which(!is.na(grid[])),kToUse,replace = F)
      chosenCellsForCluster = clusterIndexes[unlist(lapply(gridCellsToUse, function(currCell){
        chosenCell = which(cellLocationOnGrid==currCell)
        if(length(chosenCell)>1){
          chosenCell = sample(chosenCell,1,prob = meanDistMatrix[clusterIndexes[chosenCell]])
        }
        chosenCell
      }))]
      chosenCellsForCluster
    })))
    cellsToReduce = chosenCells[which(meanDistMatrix[chosenCells]>minValueToReduceTo)]
    meanDistMatrix[cellsToReduce] <<- meanDistMatrix[cellsToReduce]/10
    chosenCells
  })

  chosenNeigList = lapply(1:length(refNames),function(cellIndex){
      selectedCellType = refNames[cellIndex]
      selectedCellIndexes = which(refNames == selectedCellType)
      cellXY = XY[cellIndex,]
      cellDist = fields::rdist(t(as.matrix(cellXY)),XY[selectedCellIndexes,])
      chosenRepeats = order(as.numeric(cellDist),decreasing = F)[1:neighborhoodSize]
      chosenRepeats = chosenRepeats[!is.na(chosenRepeats)]
      selectedCellIndexes[chosenRepeats]
  })
  list(chosenCellList = chosenCellList, chosenNeigList = chosenNeigList, numOfRuns = numOfRuns)
}

########## CPM main part
#' @keywords internal
CPMMain = function(refference,refferenceNames, Y, chosenCellList, chosenCellNeigList ,numOfRuns, modelSize, neighborhoodSize,
                     no_cores, genePercents, quantifyTypes, typeTransformation, calculateCI){
  YReduced = Y[row.names(Y) %in% row.names(refference),]

  ##### Revome genes low in reference data  #####
  geneVarianceRef = apply(refference,1,function(gene){checkVariableGenes(as.numeric(as.matrix(gene)),0.1)})
  geneVarianceFinalRef = sort(geneVarianceRef[geneVarianceRef>0],decreasing = T)
  mutualGenes = names(geneVarianceFinalRef)[names(geneVarianceFinalRef) %in% row.names(YReduced)]

  YReduced = YReduced[mutualGenes,]
  refferenceSmaller = refference[mutualGenes,]

  ##### Main algorithm runs #####
  if(is.null(no_cores)){
    no_cores = max(1, parallel::detectCores() - 1)
  }
  cl<-parallel::makeCluster(no_cores)
  parallel::clusterExport(cl=cl, varlist=c("refferenceNames", "refferenceSmaller", "YReduced","neighborhoodSize","modelSize",
                                 "createSpecificRef","GeneBasedAnova", "chosenCellList", "chosenCellNeigList" ,"createNoCellDupeReference", "selectGenesUsingKappa", "runLibLinear", "genePercents"),
               envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(min = 1, max = numOfRuns, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #t = proc.time()
  `%dopar2%` <- foreach::`%dopar%`
  runNumber = NULL
  resultSmallMatrixes <- foreach::foreach(runNumber = 1:numOfRuns, .options.snow = opts) %dopar2% {
  #for(runNumber in 1:numOfRuns) {
    print(runNumber)
    completeSpecificRefBefore = createSpecificRef(refferenceSmaller, modelSize, neighborhoodSize, genePercents, chosenCellList[[runNumber]], chosenCellNeigList)
    completeSpecificRef = completeSpecificRefBefore$ref

    #list of genes between clusters
    clusterNamesVector = rep("", dim(completeSpecificRef)[2])
    for(cluster in unique(refferenceNames)){
      selectedCellsForCluster = completeSpecificRefBefore$chosenCells[which(refferenceNames[completeSpecificRefBefore$chosenCells] == cluster)]
      selectedNamesForCluster = paste(colnames(refferenceSmaller)[selectedCellsForCluster],selectedCellsForCluster,sep="_")
      clusterNamesVector[!is.na(match(colnames(completeSpecificRef),selectedNamesForCluster))] = cluster
    }

    allGenes = c()

    #list of genes inside clusters
      allGenesList = lapply(unique(refferenceNames), function(cluster){
        specificClusterRef = completeSpecificRef[,which(clusterNamesVector == cluster)]
        colnames(specificClusterRef) = colnames(completeSpecificRef)[clusterNamesVector == cluster]
        GeneBasedAnova(specificClusterRef)
      })

    for (list in allGenesList){
      allGenes = c(allGenes, list)
    }

      specificRefNoDupes = createNoCellDupeReference(completeSpecificRef)

      results = selectGenesUsingKappa(specificRefNoDupes, allGenes)
      X = results$reference
      X = refferenceSmaller[row.names(X),completeSpecificRefBefore$chosenCells]
      X = X[rowSums(X)!=0,]

      YRefinedReduced = YReduced[row.names(YReduced) %in% row.names(X),]
      PBSReductionData = YRefinedReduced

    setTxtProgressBar(pb, runNumber)

    resMatrix = t(runLibLinear(X, PBSReductionData))
    row.names(resMatrix) = chosenCellList[[runNumber]]
    resMatrix
  }
  #print(proc.time() - t)
  parallel::stopCluster(cl)
  close(pb)

  ##### Combining cell predictions #####
  print("Combining CPM iterations")
  predictedCells = matrix(0, nrow = dim(YReduced)[2], ncol = dim(refferenceSmaller)[2])
  predictedCellsCounts = matrix(0, nrow = dim(YReduced)[2], ncol = dim(refferenceSmaller)[2])

  for(resultMatrix in resultSmallMatrixes){
    completeResultMatrix = matrix(0, nrow = dim(resultMatrix)[2], ncol = dim(refferenceSmaller)[2])
    completeResultMatrix[,as.numeric(as.matrix(row.names(resultMatrix)))] = t(resultMatrix)
    predictedCells = predictedCells + completeResultMatrix
    predictedCellsCounts = predictedCellsCounts + abs(sign(completeResultMatrix))
  }
  predictedCellsFinal = predictedCells/predictedCellsCounts

  ##### Smoothing #####
  print("Smoothing")
  allClusterMeansMatrix = t(do.call(rbind,lapply(1:length(refferenceNames),function(cell){
    rowMeans(predictedCellsFinal[,chosenCellNeigList[[cell]]])
  })))
  colnames(allClusterMeansMatrix) = colnames(refference)
  row.names(allClusterMeansMatrix) = colnames(Y)

  cellTypeRes = NULL
  seRes = NULL
  confMatrix = NULL

  #### Cell type prediction ####
  if(quantifyTypes){
    print("Calculating cell type quantities")
    allClusterMeansMatrixForCellTypes = allClusterMeansMatrix
    if(typeTransformation){
      allClusterMeansMatrixForCellTypes = t(apply(t(allClusterMeansMatrixForCellTypes),2,function(x){
        x-min(x)
      }))
    }

    # cellTypeRes = do.call(cbind,lapply(unique(refferenceNames),function(currCluster){
    #   apply(allClusterMeansMatrixForCellTypes,1,function(x){
    #     #ks.test(x,x[currCluster==refferenceNames],alternative = "less")$p.value
    #     df = as.data.frame(cbind(sample(x,length(which(currCluster==refferenceNames))),x[currCluster==refferenceNames]))
    #     lm(V2~V1+0, data = df)$coefficients[1]
    #   })
    # }))

    cellTypeRes = do.call(cbind,lapply(unique(refferenceNames),function(currCluster){
      rowMeans(allClusterMeansMatrixForCellTypes[,currCluster==refferenceNames])
    }))
    colnames(cellTypeRes) = unique(refferenceNames)

    if(typeTransformation){
      cellTypeRes = t(apply(t(cellTypeRes),2,function(x){
        #x = (x-min(x))
        x/sum(x)
      })
      )
    }
  }

  #### Standard error prediction ####
  if(calculateCI){
    print("Calculating the confidence interval matrix")

    resultOriginalSizeMatrixes = lapply(resultSmallMatrixes, function(resultSmallMatrix){
      completeResultMatrix = matrix(NA, nrow = dim(resultSmallMatrix)[2], ncol = dim(refferenceSmaller)[2])
      completeResultMatrix[,match(colnames(allClusterMeansMatrix)[as.numeric(as.matrix(row.names(resultSmallMatrix)))],colnames(refferenceSmaller))] = t(resultSmallMatrix)
      completeResultMatrix
    })

    seRes <- do.call(rbind,lapply(colnames(YReduced), function(sample){
      sampleMatrix = do.call(rbind, lapply(resultOriginalSizeMatrixes,function(currRes){
        currRes[which(colnames(YReduced)==sample),]
      }))
      apply(sampleMatrix,2,function(x){
        sd(x[!is.na(x)])/sqrt(length(which(!is.na(x))))
      })
    }))

    seResNorm = t(do.call(rbind,lapply(1:length(refferenceNames),function(cell){
      rowMeans(seRes[,chosenCellNeigList[[cell]]])
    })))

    confMatrix = matrix(paste(allClusterMeansMatrix-1.96*seResNorm,allClusterMeansMatrix+1.96*seResNorm,sep = " <-> "),ncol = dim(allClusterMeansMatrix)[2])

    colnames(seRes) = colnames(confMatrix) = colnames(refference)
    row.names(seRes) = row.names(confMatrix) = colnames(Y)
  }

  print("Done")
  list(predictions = allClusterMeansMatrix, cellTypePredictions = cellTypeRes, sePredictions = seRes, confMatrix = confMatrix)
}

########## CPM
#' The Cellular Population Mapping (CPM) algorithm.
#'
#' This function initiate the Cellular Population Mapping (CPM) algorithm - a deconvolution algorithm in which single-cell genomics is required in only one or a few samples, where in other samples of the same tissue, only bulk genomics is measured and the underlying fine resolution cellular heterogeneity is inferred.
#' CPM predicts the abundance of cells (and cell types) ranging monotonically from negative to positive levels. Using a relative framework these values correspond to decrease and increase in cell abundance levels, respectively. On the other hand, in an absolute framework lower values (including negatives) correspond to lower abundances and vise versa. These values are comparable between samples.
#'
#' @param SCData A matrix containing the single-cell RNA-seq data. Each row corresponds to a certain gene and each column to a certain cell. Importantly, CPM relies on many iterative processes and therefore might take a long running time. For extremely large single cell datasets, we suggest to use only a portion of the data, using random sampling of the cells.
#' @param SCLabels A vector containing the labels of each of the cells.
#' @param BulkData A matrix containing heterogenous RNA-seq data for one or more samples. Each row corresponds to a certain gene and each column to a certain sample.
#' @param cellSpace The cell state space corresponding to the single-cell RNA-seq data. It can be a vector for a 1-dim space or a 2D matrix for a two space where each column represents a different dimension. The cell space should incorporate the similarities of cells within cell types. Similarities between cells from different cell types, based on the cell space, are not taken into account in CPM.
#' @param no_cores A number for the amount of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.
#' @param neighborhoodSize Cell neighborhood size which will be used for the analysis. This should be lower than the number of cells in the smallest cell type. The defalt is 10.
#' @param modelSize The reference subset size in each iteration of CPM. This should be lower than the total number of cells. The defalt is 50.
#' @param minSelection The minimum number of times in which each reference cell is selected. Increasing this value might have a large effect on the algorithm's running time. The defalt is 5.
#' @param quantifyTypes A boolean parameter indicating whether the prediction of cell type quantities is needed. This is recommended only in the case of homogenicity within cell types. Cell types with high inner cellular variability will recieve less reliabe values. The default is FALSE.
#' @param typeTransformation This parameter will have an effect only if quantifyTypes = TRUE. A boolean parameter indicating whether cell type deconvolution should be provided in fractions. This is done by substracting all cell types by values of the minimal cell type and dividing in their sum. This is not recommended, since it reduces the comparability between sample. The default is FALSE.
#' @param calculateCI A boolean parameter indicating whether the calculation of confidence itervals is needed. The default is FALSE.
#' @return A list including:
#' \item{predicted}{CPM predicted cell abundance matrix. Each row represents a sample and each column a single cell.}
#' \item{cellTypePredictions}{CPM predicted cell-type abundance matrix. Each row represnts a sample and each column a single cell-type. This is calculated if quantifyTypes = TRUE. }
#' \item{confIntervals}{A matrix containing the confidence iterval for each cell and sample. Each row represnts a sample and each column a single cell. This is calculated if calculateCI = TRUE.}
#' \item{numOfRuns}{The number of deconvolution repeats preformed by CPM. }
#' @references
#' Frishberg, A., Peshes-Yaloz, N., Cohn, O., Rosentul, D., Steuerman, Y., Valadarsky, L., Yankovitz, G., Mandelboim, M., Iraqi, F.A., Amit, I. et al. (2019) Cell composition analysis of bulk genomics using single-cell data. Nature Methods, 16, 327-332.
#' @examples
#' data(SCLabels)
#' data(SCFlu)
#' data(BulkFlu)
#' data(SCCellSpace)
#'
#' # Creating relative bulk data (Infleunza infection compared to PBS):
#' BulkFluReduced = BulkFlu - rowMeans(BulkFlu[,grep("pbs",colnames(BulkFlu))])
#' BulkFluReduced = BulkFluReduced[,grep("flu",colnames(BulkFluReduced))]
#'
#' # Running CPM using only a single cell-type:
#' oneCellTypeIndexes = which(SCLabels == "MPS")
#' res = CPM(SCData = SCFlu[,oneCellTypeIndexes], SCLabels = SCLabels[oneCellTypeIndexes],
#'           BulkData = BulkFluReduced, cellSpace = SCCellSpace[oneCellTypeIndexes,], no_cores = 2)
#'
#' \dontrun{
#'
#' # Running CPM using a variety of cell-types:
#' res = CPM(SCFlu, SCLabels, BulkFluReduced, SCCellSpace, no_cores = 2)
#' ### Full multi-threading : CPM(SCFlu, SCLabels, BulkFluReduced, SCCellSpace)
#' }
#' @export
#' @importFrom "utils" "setTxtProgressBar"
#' @importFrom "stats" "sd" "var"
#' @importFrom "grDevices" "chull"
CPM = function(SCData, SCLabels, BulkData, cellSpace, no_cores = NULL, neighborhoodSize = 10, modelSize = 50, minSelection = 5, quantifyTypes = F, typeTransformation = F, calculateCI = F){
  genePercents = 0.4
  if(min(table(SCLabels))<neighborhoodSize){
    neighborhoodSize = min(table(SCLabels))
    print(paste("Neighborhood size was switched to:",neighborhoodSize,sep=" "))
  }
  if(length(SCLabels)<modelSize){
    modelSize = length(SCLabels)
    print(paste("Model size was switched to:",modelSize,sep=" "))
  }
  if(!is.null(SCData) & !is.null(SCLabels) & !is.null(BulkData) & !is.null(cellSpace)){
    print("Selecting cells for each iteration")
  }
  cellSelection = choseCellsForRuns(cellSpace, SCLabels, modelSize, minSelection,neighborhoodSize)
  numOfRunsToUse = cellSelection$numOfRuns
  print(paste("Number of iteration:",numOfRunsToUse,sep=" "))
  cellSelectionList = cellSelection$chosenCellList
  cellNeigSelectionList = cellSelection$chosenNeigList
  print("Running CPM, this may take a few minutes")
  deconvolutionRes = CPMMain(SCData, SCLabels,BulkData, cellSelectionList, cellNeigSelectionList, numOfRunsToUse,modelSize, neighborhoodSize, no_cores, genePercents, quantifyTypes, typeTransformation, calculateCI)
  list(predicted = deconvolutionRes$predictions, cellTypePredictions = deconvolutionRes$cellTypePredictions, confIntervals = deconvolutionRes$confMatrix, numOfRuns = numOfRunsToUse)
}

#' Gene expression profiles of flu and pbs sample.
#'
#' A dataset containing the RNA-seq profiles of colaborative-cross (CC) mice 2 days after the infection with either the flu virus or PBS.
#'
#' @format A matrix with 1858 rows (genes) and 74 columns (samples).
"BulkFlu"

#' Gene expression profiles of lung cells after influenza infection.
#'
#' A dataset containing the RNA-seq profiles of lung cells from multiple cell types, taken from two mice 2 days after the infection with either the flu virus or PBS.
#'
#' @format A matrix with 1858 rows (genes) and 349 columns (cells).
"SCFlu"

#' Single-cell classification into cell types.
#'
#' A dataset containing the classification of each cell in the SCFlu dataset to a specific cell type.
#'
#' @format A vector with 349 values.
"SCLabels"

#' Single-cell cell space.
#'
#' A dataset containing a 2-dim cell space of all single-cells in the SCFlu dataset.
#'
#' @format A matrix with 349 rows (cells) and 2 columns (dimensions).
"SCCellSpace"
