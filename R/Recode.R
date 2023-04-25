#' Perform ROMA on a datasets
#'
#' @param ExpressionMatrix matrix, a numeric matrix containing the gene expression information. Columns indicate samples and rows indicate genes.
#' @param ModuleList list, gene module list
#' @param UseWeights logical, should the weights be used for PCA calculation?
#' @param ExpFilter logical, should the samples be filtered?
#' @param MinGenes integer, the minimum number of genes reported by a module available in the expression matrix to process the module
#' @param MaxGenes integer, the maximum number of genes reported by a module available in the expression matrix to process the module
#' @param nSamples integer, the number of randomized genes sampled (per module)
#' @param ApproxSamples integer (between 0 and 100), the approximation parameter to reuse samples. This is the minimal percentage variation to reuse samples.
#' For example, 5 means that samples are re recalculated only if the number of genes in the geneset has increased by at least 5\%.
#' @param OutGeneNumber scalar, number of median-absolute-deviations away from median required for the total number of genes expressed in a sample to be called an outlier
#' @param Ncomp integer, number of principal components used to filter samples in the gene expression space
#' @param OutGeneSpace scalar, number of median-absolute-deviations away from median required for in a sample to be called
#' an outlier in the gene expression space. If set to NULL, the gene space filtering will not be performed.
#' @param GeneOutDetection character scalar, the algorithm used to filter genes in a module. Possible values are
#' \itemize{
#' \item 'L1OutVarPerc': percentage variation relative to the median variance explained supported by a leave one out approach
#' \item 'L1OutVarDC': dendrogram clustering statistics on variance explained supported by a leave one out approach
#' \item 'L1OutExpOut': number of median-absolute-deviations away from median explained variance
#' \item 'L1OutSdMean': Number of standard deviations away from the mean
#' }
#' The option "L1OutExpOut" requires the scater package.
#' @param GeneOutThr scalar, threshold used by gene filtering algorithm in the modules. It can represent maximum size of filtered cluster ("L1OutVarDC"),
#' minimal percentage variation (L1OutVarPerc) or the number of median-absolute-deviations away from median ("L1OutExpOut")
#' @param GeneSelMode character scalar, mode used to sample genes: all available genes ("All") or genes not present in the module ("Others")
#' @param centerData logical, should the gene expression values be centered over the samples?
#' @param MoreInfo logical, should detailed information on the computation be printed?
#' @param SampleFilter logical, should outlier detection be applied to sampled data as well?
#' @param PCADims integer, the number of PCA dimensions to compute. Should be >= 2. Note that the value 1 is allowed,
#' but is not advisable under normal circumstances.
#' Larger values decrease the error in the estimation of the explained variance but increase the computation time.
#' @param DefaultWeight integer scalar, the default weight to us if no weight is specified by the modole file and an algorithm requiring weights is used
#' @param PCSignMode character scalar, the modality to use to determine the direction of the principal components. The following options are currently available:
#' \itemize{
#' \item 'none' (The direction is chosen at random)
#' \item 'PreferActivation': the direction is chosen in such a way that the sum of the projections is positive
#' \item 'UseAllWeights': as 'PreferActivation', but the projections are multiplied by the weights, missing weights are set to DefaultWeight
#' \item 'UseKnownWeights': as 'UseAllWeights', but missing weights are set to 0
#' \item 'CorrelateAllWeightsByGene': the direction is chosen in such a way as to maximize the positive correlation between the expression of genes with positive (negative) weights
#' and the (reversed) PC projections, missing weights are set to DefaultWeight
#' \item 'CorrelateKnownWeightsByGene': as 'CorrelateAllWeights', but missing weights are set to 0
#' \item 'CorrelateAllWeightsBySample': the direction is chosen in such a way as to maximize the positive correlation between the expression of genes and the PC corrected weight
#' (i.e., PC weights are multiplied by gene weights), missing weights are set to DefaultWeight
#' \item 'CorrelateKnownWeightsBySample': as 'CorrelateAllWeightsBySample', but missing weights are set to 0
#' \item 'UseMeanExpressionAllWeights': Only looking at the x% of most extreme gene weights. Then multiplying by corresponding mean expression across samples (scaled) and sum results for each gene. If negative, change sign.
#' \item 'UseMeanExpressionKnownWeights': as 'UseMeanExpressionAllWeights', but missing weights are set to 0
#' }
#' If 'CorrelateAllWeights', 'CorrelateKnownWeights', 'CorrelateAllWeightsBySample' or 'CorrelateKnownWeightsBySample' are used
#' and GroupPCSign is TRUE, the correlations will be computed on the groups defined by Grouping.
#' @param PCSignThr numeric scalar, a quantile threshold to limit the projections (or weights) to use, e.g., if equal to .9
#' only the 10\% of genes with the largest projections (or weights) in terms of absolute value will be considered.
#' @param UseParallel boolean, should a parallel environment be used ? Note that using a parallel environment will increase the memory usage as a
#' copy of the gene expression matrix is needed for each core
#' @param nCores integer, the number of cores to use if UseParallel is TRUE. Set to NULL for auto-detection
#' @param ClusType string, the cluster type to use. The default value ("PSOCK") should be available on most systems, unix-like environments also support "FORK",
#' which should be faster.
#' @param SamplingGeneWeights named vector, numeric. Weight to use when correcting the sign of the PC for sampled data.
#' @param FillNAMethod name list, additional parameters to pass to the mice function
#' @param Grouping name vector, the groups associated with the sample.
#' @param FullSampleInfo boolean, should full PC information be computed and saved for all randomized genesets?
#' @param SampleSign boolean, should we try to reorient sampled genesets ? Usefull only if FullSampleInfo is TRUE
#' @param GroupPCSign boolean, should grouping information be used to orient PCs?
#' @param CorMethod character string indicating which correlation coefficient is to be used
#' to orient the principal components. Can be "pearson", "kendall", or "spearman".
#' @param SuppressWarning boolean, should warnings be displayed? This option can be ignored in non-interactive sessions.
#' @param ShowParallelPB boolean, should the progress bar be displayed when using parallel processing ? Note that the
#' progress bar is displayed via the pbapply package. This may slow down the computation, especially when using FORK clusters.
#' @param OutlierRarelyFoundThr numeric scalar, the minimum number of pathways in which a gene has to be found to potentially be considered as an outlier
#' @param OutlierFisherThr numeric scalar, the fisher p value to use as threshold to detect outliers.
#' @param OutliersPerc numeric, how much of the data should be removed based on outliers ? (eg 0.05 for a maximum of 5% of outliers)
#' @return
#' @export
#'
#' @examples
rRoma.R <- function(ExpressionMatrix,
                           ModuleList,
                           centerData = TRUE,
                           ExpFilter=FALSE,
                           UseWeights = FALSE,
                           DefaultWeight = 1,
                           MinGenes = 10,
                           MaxGenes = 1000,
                           ApproxSamples = 5,
                           nSamples = 100,
                           OutGeneNumber = 5,
                           Ncomp = 10,
                           OutGeneSpace = NULL,
                           GeneOutDetection = "L1OutExpOut",
                           GeneOutThr = 5,
                           GeneSelMode = "All",
                           SampleFilter = TRUE,
                           MoreInfo = FALSE,
                           PCADims = 2,
                           PCSignMode ='none',
                           PCSignThr = 0.90,
                           UseParallel = FALSE,
                           nCores = NULL,
                           ClusType = "PSOCK",
                           SamplingGeneWeights = NULL,
                           FillNAMethod = list(),
                           Grouping = NULL,
                           FullSampleInfo = TRUE,
                           SampleSign = FALSE,
                           GroupPCSign = FALSE,
                           CorMethod = "pearson",
                           SuppressWarning = FALSE,
                           ShowParallelPB = FALSE,
                           OutlierRarelyFoundThr = 5,
                           OutlierFisherThr = 0.05,
                           OutliersPerc = 0.05) {
  
  if(PCADims < 1){
    stop("PCADims should be >= 1")
  }
  
  if(PCADims == 1){
    print("PCADims = 1. Be aware that this is not advisable under normal circumstances")
    if(!SuppressWarning){
      readline("Press any key")
    }
  }
  
  if(is.null(SamplingGeneWeights)){
    SamplingGeneWeights = rep(DefaultWeight, nrow(ExpressionMatrix))
    names(SamplingGeneWeights) <- rownames(ExpressionMatrix)
  }
  
  if(any(is.na(ExpressionMatrix))){
    
    if(!requireNamespace("mice", quietly = TRUE)){
      stop("Unable to load mice. Impossible to proceed")
    } else {
      print("Filling NA with mice")
      library(mice)
    }
    
    imp <- do.call(what = mice, args = append(list(data = ExpressionMatrix), FillNAMethod))
    ExpressionMatrix <- complete(imp)
    
  }
  
  ExpressionMatrix <- data.matrix(ExpressionMatrix)
  
  SAMPLE_WARNING <- 3
  
  AllGenesModule <- unique(unlist(lapply(ModuleList, "[[", "Genes")))
  AllGenesMatrix <- rownames(ExpressionMatrix)
  
  if(sum(AllGenesMatrix %in% AllGenesModule) < 3){
    stop("Not enough module genes found in the matrix. Impossible to proceed")
  }
  
  if(any(is.na(AllGenesMatrix) | is.null(AllGenesMatrix))){
    print("Missing gene name in following line(s):")
    print(which(is.na(AllGenesMatrix) | is.null(AllGenesMatrix)))
    print("Please fix this before proceding.")
    stop("Impossible to proceed")
  }
  
  if(ncol(ExpressionMatrix) <= SAMPLE_WARNING & interactive()){
    print(paste("Only", ncol(ExpressionMatrix), "sample found"))
    print("The number of samples may be too small to guarantee a reliable analysis")
    if(!SuppressWarning){
      Ans <- readline("Do you want to continue anyway? (y/n)")
      if(Ans != "y" & Ans != "Y"){
        return(NULL)
      }
    }
  }
  
  if(any(AllGenesMatrix[duplicated(AllGenesMatrix)] %in% AllGenesModule)){
    stop("Module genes are not unique in the matrix. Impossible to proceed.")
  }
  
  if(any(duplicated(AllGenesMatrix)) & interactive()){
    print("Duplicated gene names detected. This may create inconsistencies in the analysis. Consider fixing this problem.")
    if(!SuppressWarning){
      readline("Press any key")
    }
  }
  
  if(SampleSign & interactive()){
    print("PC projections and weights will be computed and reoriented for sampled genesets. This is potentially very time consuming.")
    if(!SuppressWarning){
      Ans <- readline("Are you sure you want to do that? (y/n)")
      if(Ans != "y" & Ans != "Y"){
        SampleSign <- FALSE
        print("PC projections and weights will NOT be reoriented for sampled genesets.")
      } else {
        print("PC projections and weights will be computed and reoriented for sampled genesets.")
      }
    }
  }
  
  # Cleanup data
  
  if(ExpFilter){
    
    print("Filtering samples with abnormal numbers of genes detected")
    
    GeneDetectesOUT <- scater::isOutlier(colSums(ExpressionMatrix==0), nmads = OutGeneNumber)
    
    print(paste(sum(GeneDetectesOUT), "sample(s) will be filtered:"))
    if(sum(GeneDetectesOUT)>0){
      print(names(which(GeneDetectesOUT)))
    }
    
    print("Filtering samples in the gene space")
    
    if(is.null(Ncomp) | Ncomp > ncol(ExpressionMatrix)){
      Ncomp <- ncol(ExpressionMatrix)
    }
    
    if(!is.null(OutGeneSpace)){
      GeneSpaceOUT <- scater::isOutlier(rowSums(prcomp(t(ExpressionMatrix), center = TRUE, retx = TRUE)$x[,1:Ncomp]^2),
                                        nmads = OutGeneSpace)
      
      print(paste(sum(GeneSpaceOUT), "sample(s) will be filtered:"))
      if(sum(GeneSpaceOUT)>0){
        print(names(which(GeneSpaceOUT)))
      }
      
      ExpressionMatrix <- ExpressionMatrix[, !GeneDetectesOUT & !GeneSpaceOUT]
    }
    
  }
  
  # Look at groups
  
  if(!is.null(Grouping)){
    GrpNames <- unique(Grouping)
    GrpCol <- rainbow(length(GrpNames))
    names(GrpCol) <- GrpNames
    
    NotFoundSampNames <- which(!(colnames(ExpressionMatrix) %in% names(Grouping)))
    FoundSampNames <- which((colnames(ExpressionMatrix) %in% names(Grouping)))
    
    if(length(NotFoundSampNames)>1){
      print(paste("The following samples do not have an associated group:"))
      print(colnames(ExpressionMatrix)[NotFoundSampNames])
      print("They will NOT be considered for group associated analysis")
    }
    
  } else {
    FoundSampNames <- NULL
    Grouping <- rep("N/A", ncol(ExpressionMatrix))
    names(Grouping) <- colnames(ExpressionMatrix)
  }
  
  OrgExpMatrix = ExpressionMatrix
  
  if(centerData){
    print("Centering gene expression over samples (genes will have 0 mean)")
    ExpressionMatrix <- t(scale(t(ExpressionMatrix), center = TRUE, scale = FALSE))
    GeneCenters <- attr(ExpressionMatrix, "scaled:center")
    attr(ExpressionMatrix, "scaled:center") <- NULL
  } else {
    print("Gene expression will not be centered over samples (genes will not have 0 mean)")
    GeneCenters = rep(0, nrow(ExpressionMatrix))
  }
  
  ModuleSummary <- list()
  ModuleMatrix <- NULL
  
  SampleMatrix <- NULL
  WeightList <- list()
  
  PVVectMat <- NULL
  
  OutLiersList <- list()
  UsedModules <- NULL
  
  OutliersRank <- list()
  nGenes <- rep(NA, length(ModuleList))
  
  AllPCA1 <- list()
  
  # Filter genes for compatibility with the expression matrix
  
  for(i in 1:length(ModuleList)){
    Preserve <- ModuleList[[i]]$Genes %in% rownames(ExpressionMatrix)
    nGenes[i] <- sum(Preserve)
    
    ModuleList[[i]]$Genes <- ModuleList[[i]]$Genes[Preserve]
    ModuleList[[i]]$Weights <- ModuleList[[i]]$Weights[Preserve]
    
  }
  
  # Filter genesets depending on the number of genes
  ToFilter <- (nGenes > MaxGenes | nGenes < MinGenes)
  ToUse <- !ToFilter
  
  if(sum(ToFilter)>0){
    print("The following geneset(s) will be ignored due to the number of usable genes being outside the specified range")
    print(
      paste(unlist(lapply(ModuleList[ToFilter], "[[", "Name")), "/", nGenes[ToFilter], "gene(s)")
    )
    
    if(sum(ToUse) == 0){
      print("No geneset available for the analisis. The analysis cannot proceed.")
      return(NULL)
    }
    
    ModuleList <- ModuleList[ToUse]
    nGenes <- nGenes[ToUse]
  } else {
    print("All the genesets will be used")
  }
  
  ExpressionMatrix <- scale(ExpressionMatrix, center = TRUE, scale = FALSE)
  
  if(UseParallel){
    
    no_cores <- parallel::detectCores()
    
    if(is.null(nCores)){
      # Calculate the number of cores
      nCores = no_cores - 1
    }
    
    if(nCores > no_cores){
      nCores = no_cores - 1
      print(paste("Too many cores selected!", nCores, "will be used"))
    }
    
    if(nCores == no_cores){
      print("Using all the available cores. This will likely render the computer unresponsive during the analysis.")
    }
    
    if(ClusType == "PSOCK"){
      # Initiate sock cluster
      cl <- parallel::makePSOCKcluster(nCores)
      
      parallel::clusterExport(cl=cl, varlist=c("SampleFilter", "GeneOutDetection", "GeneOutThr",
                                               "ModulePCACenter", "ExpressionMatrix", "DetectOutliers",
                                               "PCADims", "OrgExpMatrix", "FullSampleInfo", "FoundSampNames", 
                                               "GroupPCSign"),
                              envir = environment())
    }
    
    if(ClusType == "FORK"){
      # Initiate sock cluster
      cl <- parallel::makeForkCluster(nnodes = nCores)
    }
    
  }
  
  ModuleOrder <- order(unlist(lapply(lapply(ModuleList, "[[", "Genes"), length)))
  ModuleList <- ModuleList[ModuleOrder]
  
  OldSamplesLen <- 0
  
  # Detecting outliers for each Module
  Genes_pathways_counts <- list()
  Genes_outliers_counts <- list()
  ModuleLengths <- c()
  
  for(i in 1:length(ModuleList)){
    CompatibleGenes <- ModuleList[[i]]$Genes
    ModuleLengths <- c(ModuleLengths, length(CompatibleGenes))
    for(j in CompatibleGenes){
      if(!j %in% names(Genes_pathways_counts)){
        Genes_pathways_counts[[j]] <- 1
        Genes_outliers_counts[[j]] <- 0
        
      }
      else{
        Genes_pathways_counts[[j]] <- Genes_pathways_counts[[j]]+  1
      }
    }

    print(ModuleList[[i]]$Name)
    
    if(UseParallel){
      SelGenes <- DetectOutliers(GeneOutDetection = GeneOutDetection, GeneOutThr = GeneOutThr, ModulePCACenter = FALSE,
                                        CompatibleGenes = CompatibleGenes, ExpressionData = ExpressionMatrix[CompatibleGenes, ],
                                        ModuleName = ModuleList[[i]]$Name, Mode = 3,
                                        ClusType = ClusType, cl = cl)
    } else {
      SelGenes <- DetectOutliers(GeneOutDetection = GeneOutDetection, GeneOutThr = GeneOutThr, ModulePCACenter = FALSE,
                                        CompatibleGenes = CompatibleGenes, ExpressionData = ExpressionMatrix[CompatibleGenes, ],
                                        ModuleName = ModuleList[[i]]$Name, Mode = 1)
    }
    
    OutLiersList[[i]] <- SelGenes[["FiltGenes"]]
    OutliersRank[[i]] <- SelGenes[["FiltRank"]]
    AllPCA1[[i]] <- SelGenes[["AllPCA1"]]
    
    Genes_outliers_counts[OutLiersList[[i]]] <- lapply(Genes_outliers_counts[OutLiersList[[i]]], function(x){x+1})
    
  }
  
  #Adjusting for rarely found genes
  genes_to_keep_counts <- c()
  genes_to_keep_fisher <- c()
  genes_to_remove <- c()
  tot_pathways_counts <- sum(sapply(Genes_pathways_counts, "[[", 1))
  tot_outliers_counts <- sum(sapply(Genes_outliers_counts, "[[", 1))
  
  for(i in names(Genes_pathways_counts)){
    if (Genes_pathways_counts[[i]] < OutlierRarelyFoundThr){
      genes_to_keep_counts <- c(genes_to_keep_counts, i)
    }
    if (fisher.test(matrix(c(tot_outliers_counts - Genes_outliers_counts[[i]], Genes_outliers_counts[[i]], 
                             tot_pathways_counts - Genes_pathways_counts[[i]], 
                             Genes_pathways_counts[[i]]), nrow = 2))$p.value >= OutlierFisherThr){
      genes_to_keep_fisher <- c(genes_to_keep_fisher, i)
    }
  }
  
  genes_to_remove <- names(Genes_pathways_counts)[! names(Genes_pathways_counts) %in% union(genes_to_keep_counts, genes_to_keep_fisher)]
  OriginalOutliers <- OutLiersList
  for(i in c(1:length(OutLiersList))){
    ranks <- rank(OutliersRank[[i]][OutLiersList[[i]] %in% genes_to_remove]) <= max(c(round(ModuleLengths[i]*OutliersPerc), 1))
    OutLiersList[[i]] <- OutLiersList[[i]][OutLiersList[[i]] %in% genes_to_remove][ranks]
  }
  
  KeptGenes <- list()
  for(i in 1:length(ModuleList)){
    KeptGenes[[i]] <- setdiff(ModuleList[[i]]$Genes, OutLiersList[[i]])
  }
  
  SampleCenters <- attr(ExpressionMatrix, "scaled:center")
  attr(ExpressionMatrix, "scaled:center") <- NULL
  
  # Filter genesets depending on the number of genes
  nGenes <- sapply(KeptGenes, length)
  ToFilter <- (nGenes > MaxGenes | nGenes < MinGenes)
  ToUse <- !ToFilter
  
  if(sum(ToFilter)>0){
    print("The following geneset(s) will be ignored due to the number of usable genes avec outlier detection being outside the specified range")
    print(
      paste(unlist(lapply(ModuleList[ToFilter], "[[", "Name")), "/", nGenes[ToFilter], "gene(s)")
    )
    
    if(sum(ToUse) == 0){
      print("No geneset remaining for the analisis. The analysis cannot proceed.")
      return(NULL)
    }
    
  } else {
    print("All the remainibg genesets will be used")
  }
  
  
  if(MoreInfo){
    print("The following genesets will be used:")
    paste(unlist(lapply(ModuleList, "[[", "Name")), "/", nGenes, "gene(s)")
  }
  
  sampled_sets_pc1_mean <- c()
  modules_pc1_mean <- c()
  
  for(i in 1:length(ModuleList)){
    
    print(Sys.time())
    gc()
    
    print(paste("[", i, "/", length(ModuleList), "] Working on ", ModuleList[[i]]$Name, " - ", ModuleList[[i]]$Desc, sep = ""))
    
    SelGenes <- KeptGenes[[i]]
    
    if(!ToUse[i]){
      print("Number of selected genes outside the specified range")
      print("Skipping module")
      next()
    } else {
      UsedModules <- c(UsedModules, i)
    }
    
    print(paste(length(SelGenes), "genes used for analysis"))
    
    if(MoreInfo){
      print("The following genes will be used:")
      print(SelGenes)
    }
    
    # Working first on the unfiltered data (only for reference)
    
    CompatibleGenes <- ModuleList[[i]]$Genes
    
    if(UseWeights){
      print("Using weights")
      Correction <- ModuleList[[i]]$Weights
      names(Correction) <- CompatibleGenes
      Correction[!is.finite(Correction)] <- DefaultWeight
    } else {
      print("Not using weights for PCA computation")
      Correction <- rep(1, length(CompatibleGenes))
      names(Correction) <- CompatibleGenes
    }
    
    BaseMatrix <- Correction[CompatibleGenes]*ExpressionMatrix[CompatibleGenes, ]
    
    if(PCADims > min(dim(ExpressionMatrix[CompatibleGenes, ])/3)){
      PCBase <- prcomp(x = BaseMatrix, center = FALSE, scale. = FALSE, retx = TRUE)
    } else {
      PCBase <- irlba::prcomp_irlba(x = BaseMatrix, n = PCADims,
                                    work = min(PCADims+7, min(dim(ExpressionMatrix[CompatibleGenes, ]))),
                                    center = FALSE, scale. = FALSE, maxit = 10000, retx = TRUE)
    }
    
    
    ExpVar <- apply(PCBase$x[,1:2], 2, var)/sum(apply(scale(BaseMatrix, center = FALSE, scale = FALSE), 2, var))
    
    PC1MeanUnf <- median(PCBase$x[, 1])
    PCBaseUnf <- PCBase
    ExpVarUnf <- ExpVar
    
    MedianExp <- median(OrgExpMatrix[CompatibleGenes, ])
    
    print("Pre-filter data")
    print(paste("L1 =", ExpVar[1], "L1/L2 =", ExpVar[1]/ExpVar[2], "PC1mean =", PC1MeanUnf))
    print(paste("Median expression (uncentered):", MedianExp))
    print(paste("Median expression (global centered/weighted):", median(BaseMatrix)))
    
    # Computing PC on the filtered data
    
    
    BaseMatrix <- Correction[SelGenes]*ExpressionMatrix[SelGenes, ]
    
    if(PCADims > min(dim(ExpressionMatrix[SelGenes, ])/3)){
      PCBase <- prcomp(x = BaseMatrix, center = FALSE, scale. = FALSE, retx = TRUE)
    } else {
      PCBase <- irlba::prcomp_irlba(x = BaseMatrix, n = PCADims,
                                    work = min(PCADims+7, min(dim(ExpressionMatrix[SelGenes, ]))),
                                    center = FALSE, scale. = FALSE, maxit = 10000, retx = TRUE)
    }
    
    ExpVar <- apply(PCBase$x[,1:2], 2, var)/sum(apply(scale(BaseMatrix, center = FALSE, scale = FALSE), 2, var))
    
    PC1Mean <- median(PCBase$x[, 1])
    
    modules_pc1_mean <- c(modules_pc1_mean, PC1Mean)
    
    MedianExp <- median(OrgExpMatrix[SelGenes, ])
    
    print("Post-filter data")
    print(paste("L1 =", ExpVar[1], "L1/L2 =", ExpVar[1]/ExpVar[2], "PC1mean =", PC1Mean))
    print(paste("Median expression (uncentered):", MedianExp))
    print(paste("Median expression (global centered/weighted):", median(BaseMatrix)))
    
    print(paste("Previous sample size:", OldSamplesLen))
    print(paste("Current sample size:", length(CompatibleGenes)))
    
    # Comparison with sample genesets
    if(nSamples > 0){
      
      if((length(CompatibleGenes)*(100-ApproxSamples)/100 <= OldSamplesLen) & (OldSamplesLen > 0)){
        
        if(length(CompatibleGenes) == OldSamplesLen){
          print("Reusing previous sampling (Same metagene size)")
        } else {
          print("Reusing previous sampling (Comparable metagene size)")
        }
        
      } else {
        
        # Define base analysis function
        
        TestGenes <- function(Gl, UpdatePB = FALSE){
          # library(irlba)
          
          if(SampleFilter){
            SampleSelGenes <- setdiff(Gl, DetectOutliers(GeneOutDetection = GeneOutDetection, GeneOutThr = GeneOutThr, ModulePCACenter = FALSE,
                                                    CompatibleGenes = Gl, ExpressionData = ExpressionMatrix[Gl, ],
                                                    ModuleName = '', PrintInfo = FALSE)[["FiltGenes"]])
            if(length(SampleSelGenes)<PCADims){
              warning(paste("Size of filtered sample geneset too small (",  length(SampleSelGenes), "). This may cause inconsitencies. Increase MinGenes or GeneOutThr to prevent the problem from happening"))
            }
          } else {
            SampleSelGenes <- Gl
          }
          
          if(length(SampleSelGenes) <= 1){
            SampMedian <- median(ExpressionMatrix[SampleSelGenes])
            warning(paste("Size of filtered sample geneset extremely small (",  length(SampleSelGenes), "). This may cause inconsitencies. Increase MinGenes or GeneOutThr to prevent the problem from happening"))
            return(list("ExpVar"=c(1, rep(0, PCADims - 1)), "MedianExp"= SampMedian))
          }
          
          
          BaseMatrix <- ExpressionMatrix[SampleSelGenes, ]
          
          SampMedian <- median(OrgExpMatrix[SampleSelGenes, ])
          
          if(PCADims <= min(dim(BaseMatrix)/3)){
            PCSamp <- irlba::prcomp_irlba(x = BaseMatrix, n = PCADims,
                                          center = FALSE, scale. = FALSE, retx = TRUE)
            TotVar <- PCSamp$totalvar
          } else {
            PCSamp <- prcomp(x = BaseMatrix, center = FALSE, scale. = FALSE, retx = TRUE)
            TotVar <- sum(PCSamp$sdev **2)
            
          }
          
          VarVect <- apply(PCSamp$x[,1:2], 2, var)
          
          
          PC1Mean <- median(PCSamp$x[, 1])
          
          # if(VarVect[2] > VarVect[1]){
          # print(sqrt(VarVect)/PCSamp$sdev)
          # }
          
          if(length(VarVect)>PCADims){
            VarVect <- VarVect[1:PCADims]
          }
          
          if(length(VarVect)<PCADims){
            VarVect <- c(VarVect, rep(0, PCADims - length(VarVect)))
          }
          
          if(UpdatePB){
            pb$up(pb$getVal() + 1)
          }
          
          
          if(FullSampleInfo & SampleSign){
            
            if(GroupPCSign){
              GroupPCsVect <- Grouping
            } else {
              GroupPCsVect <- NULL
            }
            
            ExpMat <- NULL
            
            if(PCSignMode %in% c("CorrelateAllWeightsByGene", "CorrelateKnownWeightsByGene",
                                 "CorrelateAllWeightsBySample", "CorrelateKnownWeightsBySample", "UseMeanExpressionKnownWeights", "UseMeanExpressionAllWeights")){
              ExpMat <- OrgExpMatrix[SampleSelGenes, ]
            }
            
            GeneScore1 = PCSamp$x[,1]
            SampleScore1 = PCSamp$rotation[,1]
            
            CorrectSign1 <- FixPCSign(GeneScore = GeneScore1, SampleScore = SampleScore1,
                                      Wei = SamplingGeneWeights[SampleSelGenes],
                                      Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                                      Grouping = GroupPCsVect, ExpMat = ExpMat, CorMethod = CorMethod)
            
            if(PCADims > 1){
              
              GeneScore2 = PCSamp$x[,2]
              SampleScore2 = PCSamp$rotation[,2]
              
              CorrectSign2 <- FixPCSign(GeneScore = GeneScore2, SampleScore = SampleScore2,
                                        Wei = SamplingGeneWeights[SampleSelGenes],
                                        Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                                        Grouping = GroupPCsVect, ExpMat = ExpMat, CorMethod = CorMethod)
            } else {
              GeneScore2 = NULL
              SampleScore2 = NULL
              CorrectSign2 = NULL
            }
            
            return(list("ExpVar" = VarVect/sum(apply(scale(BaseMatrix, center = TRUE, scale = FALSE), 2, var)),
                        "MedianExp"= SampMedian,
                        "GenesWei"= cbind(CorrectSign1*GeneScore1, CorrectSign2*GeneScore2),
                        "SampleScore"= cbind(CorrectSign1*SampleScore1, CorrectSign2*SampleScore2),
                        "PC1Mean" = PC1Mean))
            
            
          } else if (FullSampleInfo){
            return(list("ExpVar" = VarVect/sum(apply(scale(BaseMatrix, center = TRUE, scale = FALSE), 2, var)),
                        "MedianExp"= SampMedian,
                        "GenesWei"= cbind(PCSamp$x[,1], PCSamp$x[,2]),
                        "SampleScore"= cbind(PCSamp$rotation[,1], PCSamp$rotation[,2]),
                        "PC1Mean" = PC1Mean))
          } else{
          
            
            return(list("ExpVar"=VarVect/sum(apply(scale(BaseMatrix, center = TRUE, scale = FALSE), 2, var)),
                        "MedianExp"=SampMedian, "GenesWei" = NULL, "SampleScore" = NULL,
                        "PC1Mean" = PC1Mean))
            
          }
          
        }
        
        
        if(SampleFilter & !UseParallel){
          # pb <- txtProgressBar(min = 0, max = nSamples, initial = 0, style = 3)
          GeneToSample <- length(SelGenes)
        } else {
          GeneToSample <- length(CompatibleGenes)
        }
        
        Sampling<- function(i){
          set.seed(i)
          return(sample(x = rownames(ExpressionMatrix), size = GeneToSample, replace = FALSE))
        }
        
        if(GeneSelMode == "All"){
          SampledsGeneList <- lapply(as.list(1:nSamples), Sampling)
        }
        
        if(GeneSelMode == "Others"){
          SampledsGeneList <- lapply(as.list(1:nSamples), Sampling)
        }
        
        if(!exists("SampledsGeneList")){
          stop("Incorrect sampling mode")
        }
        
        
        if(!UseParallel){
          print("Computing samples")
          pb <- txtProgressBar(min = 0, max = nSamples, initial = 0, style = 3)
          tictoc::tic()
          SampledExp <- lapply(SampledsGeneList, TestGenes, UpdatePB = TRUE)
          cat("\n")
          tictoc::toc()
          pb$kill()
        } else {
          if(ShowParallelPB){
            print("Computing samples via parallel execution with pbapply progress bar (If this is too slow try setting ShowParallelPB = FALSE)")
            pbo <- pbapply::pboptions(nout = 20)
            tictoc::tic()
            SampledExp <- pbapply::pblapply(SampledsGeneList, TestGenes, UpdatePB = FALSE, cl = cl)
            tictoc::toc()
            pbo <- pbapply::pboptions(pbo)
          } else {
            print("Computing samples (no progress bar will be shown)")
            tictoc::tic()
            SampledExp <- parallel::parLapply(cl, SampledsGeneList, TestGenes, UpdatePB = FALSE)
            tictoc::toc()
          }
        }
        
        SampleExpVar <- sapply(SampledExp, "[[", "ExpVar")
        SampleMedianExp <- sapply(SampledExp, "[[", "MedianExp")
        SamplePC1Mean <- sapply(SampledExp, "[[", "PC1Mean")
        
        sampled_sets_pc1_mean <- c(sampled_sets_pc1_mean, SamplePC1Mean)
        
        if(PCADims >= 2){
          SampleExpVar <- rbind(SampleExpVar[1,], SampleExpVar[1,]/SampleExpVar[2,], SampleExpVar[2,])
        } else {
          SampleExpVar <- rbind(SampleExpVar, rep(NA, length(SampleExpVar)), rep(NA, length(SampleExpVar)))
        }
        
        rownames(SampleExpVar) <- c("Sampled L1", "Sampled L1/L2", "Sampled L2")
        
        OldSamplesLen <- length(CompatibleGenes)
        
      }
      
    }
    
    L1Vect <- SampleExpVar[1,] - ExpVar[1]
    L1Vect <- L1Vect[is.finite(L1Vect)]
    
    L1L2Vect <- SampleExpVar[2,] - ExpVar[1]/ExpVar[2]
    L1L2Vect <- L1L2Vect[is.finite(L1L2Vect)]
    
    MedianVect <- SampleMedianExp - MedianExp
    MedianVect <- MedianVect[is.finite(MedianVect)]
    
    PC1MeanVect <- SamplePC1Mean - PC1Mean
    PC1MeanVect <- PC1MeanVect[is.finite(PC1MeanVect)]
    
    PVVect <- rep(NA, 6)
    
    if(length(L1Vect) > 5){
      PVVect[1] <- wilcox.test(L1Vect, alternative = "less")$p.value
      PVVect[2] <- wilcox.test(L1Vect, alternative = "greater")$p.value
    }
    
    if(length(L1L2Vect) > 5){
      PVVect[3] <- wilcox.test(L1L2Vect, alternative = "less")$p.value
      PVVect[4] <- wilcox.test(L1L2Vect, alternative = "greater")$p.value
    }
    
    if(length(PC1MeanVect) > 5){
      PVVect[5] <- wilcox.test(PC1MeanVect, alternative = "less")$p.value
      PVVect[6] <- wilcox.test(PC1MeanVect, alternative = "greater")$p.value
    }
    
    PVVectMat <- rbind(PVVectMat, PVVect)
    
    ModuleMatrix <- rbind(ModuleMatrix,
                          c(ExpVar[1],
                            median(SampleExpVar[1,], na.rm = TRUE),
                            sum(sign(SampleExpVar[1,] - ExpVar[1])==1, na.rm = TRUE)/(nSamples-sum(is.na(sign(SampleExpVar[1,] - ExpVar[1])==1))),
                            ExpVar[1]/ExpVar[2],
                            median(SampleExpVar[2,], na.rm = TRUE),
                            sum(sign(SampleExpVar[2,] - ExpVar[1]/ExpVar[2])==1, na.rm = TRUE)/(nSamples-sum(is.na(sign(SampleExpVar[2,] - ExpVar[1]/ExpVar[2])==1))),
                            PC1Mean,
                            sum(sign(abs(SamplePC1Mean) - abs(PC1Mean))==1)/nSamples))
    
    # Compute the sign correction
    
    if(GroupPCSign){
      GroupPCsVect <- Grouping
    } else {
      GroupPCsVect <- NULL
    }
    
    
    
    ExpMat <- NULL
    if(PCSignMode %in% c("CorrelateAllWeightsByGene", "CorrelateKnownWeightsByGene",
                         "CorrelateAllWeightsBySample", "CorrelateKnownWeightsBySample", "UseMeanExpressionKnownWeights", "UseMeanExpressionAllWeights")){
      ExpMat <- OrgExpMatrix[CompatibleGenes, ]
    }
    
    GeneScore1 = PCBaseUnf$x[,1]
    SampleScore1 = PCBaseUnf$rotation[,1]
    
    GeneScore1Unf <- GeneScore1
    SampleScore1Unf <- SampleScore1
    
    CorrectSignUnf <- FixPCSign(GeneScore = GeneScore1, SampleScore = SampleScore1,
                                Wei = ModuleList[[i]]$Weights[ModuleList[[i]]$Genes %in% CompatibleGenes],
                                Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                                Grouping = GroupPCsVect, ExpMat = ExpMat, CorMethod = CorMethod)
    
    
    
    ExpMat <- NULL
    if(PCSignMode %in% c("CorrelateAllWeightsByGene", "CorrelateKnownWeightsByGene",
                         "CorrelateAllWeightsBySample", "CorrelateKnownWeightsBySample", "UseMeanExpressionKnownWeights", "UseMeanExpressionAllWeights")){
      ExpMat <- OrgExpMatrix[SelGenes, ]
    }
    
    GeneScore1 = PCBase$x[,1]
    SampleScore1 = PCBase$rotation[,1]
    
    CorrectSign1 <- FixPCSign(GeneScore = GeneScore1, SampleScore = SampleScore1,
                              Wei = ModuleList[[i]]$Weights[ModuleList[[i]]$Genes %in% SelGenes],
                              Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                              Grouping = GroupPCsVect, ExpMat = ExpMat, CorMethod = CorMethod)
    
    
    if(PCADims >= 2){
      
      GeneScore2 = PCBase$x[,2]
      SampleScore2 = PCBase$rotation[,2]
      
      CorrectSign2 <- FixPCSign(GeneScore = GeneScore2, SampleScore = SampleScore2,
                                Wei = ModuleList[[i]]$Weights[ModuleList[[i]]$Genes %in% SelGenes],
                                Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                                Grouping = GroupPCsVect, ExpMat = ExpMat, CorMethod = CorMethod)
      
    } else {
      GeneScore2 = NULL
      SampleScore2 = NULL
      CorrectSign2 = NULL
    }
    
    # Correct the sign of the first PC projections
    names(GeneScore1) <- SelGenes
    names(SampleScore1) <- colnames(ExpressionMatrix)
    
    SampleMatrix <- rbind(SampleMatrix, SampleScore1 * CorrectSign1)
    
    WeightList[[length(WeightList)+1]] <- GeneScore1 * CorrectSign1
    
    ModuleSummary[[length(ModuleSummary)+1]] <-
      list(ModuleName = ModuleList[[i]]$Name,
           ModuleDesc = ModuleList[[i]]$Desc,
           OriginalGenes = CompatibleGenes,
           UsedGenes = SelGenes,
           SampledGenes = SampledsGeneList,
           PCABase = PCBase,
           PCBaseUnf = PCBaseUnf,
           CorrectSignUnf = CorrectSignUnf,
           CorrectSign1 = CorrectSign1,
           CorrectSign2 = CorrectSign2,
           ExpVarBase = ExpVar,
           ExpVarBaseUnf = ExpVarUnf,
           SampledExp = SampledExp,
           SampleScore = SampleScore1,
           SampleScoreUnf = SampleScore1Unf,
           PC1Mean = PC1Mean,
           GeneWeight = GeneScore1,
           GeneWeightUnf = GeneScore1Unf,
           GMTWei = ModuleList[[i]]$Weights[ModuleList[[i]]$Genes %in% SelGenes])
    
    
  }
  
  modules_pc1_mean <- sapply(modules_pc1_mean, abs)
  sampled_sets_pc1_mean <- sapply(sampled_sets_pc1_mean, abs)
  
  ModuleMatrix <- cbind(ModuleMatrix, rep(NA, nrow(ModuleMatrix)), rep(NA, nrow(ModuleMatrix)), rep(NA, nrow(ModuleMatrix)))
  
  #Ici ajouter le calcul de q value BH
  ModuleMatrix[, 9] <- p.adjust(ModuleMatrix[, 3], "BH")
  ModuleMatrix[, 10] <- p.adjust(ModuleMatrix[, 6], "BH")
  
  for(i in c(1:length(ModuleSummary))){
    pc1_mean <- abs(ModuleMatrix[i, 7])
    ModuleMatrix[i, 11] <- min(c(1, (sum(sampled_sets_pc1_mean > pc1_mean) / length(sampled_sets_pc1_mean)) / (sum(modules_pc1_mean > pc1_mean) / (length(modules_pc1_mean) -1))))
  }
  
  ModuleMatrix[is.na(ModuleMatrix[, 11]), 11] <- min(ModuleMatrix[!is.na(ModuleMatrix[, 11]), 11])
  
  #Make sure q values are ordered in the same way as p values
  
  j<- 1
  for(i in order(ModuleMatrix[, 8])){
    equal <- which(ModuleMatrix[, 8] == ModuleMatrix[, 8][i])
    ModuleMatrix[, 11][i] <- min(ModuleMatrix[unique(c(equal, order(ModuleMatrix[, 8])[j:nrow(ModuleMatrix)])), 11])
    j <- j+1
  }
  
  if(UseParallel){
    
    # Stop cluster
    parallel::stopCluster(cl)
    
  }
  
  # Makes sure ModuleMatrix is treated as a matrix
  if(length(ModuleMatrix) == 11){
    dim(ModuleMatrix) <- c(1, 11)
  }
  
  colnames(ModuleMatrix) <- c("L1", "Median L1", "ppv L1", "L1/L2", "Median L1/L2", "ppv L1/L2", "Median Exp", "ppv Median Exp", "q L1", "q L1/L2", "q Median Exp")
  rownames(ModuleMatrix) <- unlist(lapply(ModuleList, "[[", "Name"))[UsedModules]
  
  # Makes sure PVVectMat is treated as a matrix
  if(length(PVVectMat) == 6){
    dim(PVVectMat) <- c(1, 6)
  }
  
  colnames(PVVectMat) <- c("L1 WT less pv", "L1 WT greater pv", "L1/L2 WT less pv", "L1/2 WT greater pv",
                           "PC1 Mean less pv", "PC1 Mean greater pv")
  rownames(PVVectMat) <- unlist(lapply(ModuleList, "[[", "Name"))[UsedModules]
  
  # Makes sure SampleMatrix is treated as a matrix
  if(length(SampleMatrix) == ncol(ExpressionMatrix)){
    dim(SampleMatrix) <- c(1, ncol(ExpressionMatrix))
  }
  
  colnames(SampleMatrix) <- colnames(ExpressionMatrix)
  rownames(SampleMatrix) <- unlist(lapply(ModuleList, "[[", "Name"))[UsedModules]
  
  
  InputParList <- list(centerData = centerData, ExpFilter = ExpFilter, ModuleList = ModuleList,
                       UseWeights = UseWeights, DefaultWeight = DefaultWeight, MinGenes = MinGenes,
                       MaxGenes = MaxGenes, ApproxSamples = ApproxSamples, nSamples = nSamples,
                       OutGeneNumber = OutGeneNumber, Ncomp = Ncomp, OutGeneSpace = OutGeneSpace,
                       GeneOutDetection = GeneOutDetection, GeneOutThr = GeneOutThr,
                       GeneSelMode = GeneSelMode, SampleFilter = SampleFilter, MoreInfo = MoreInfo,
                       PCSignMode = PCSignMode, PCSignThr = PCSignThr,
                       UseParallel = UseParallel, nCores = nCores, ClusType = ClusType,
                       SamplingGeneWeights = SamplingGeneWeights, FillNAMethod = FillNAMethod,
                       Grouping = Grouping, FullSampleInfo = FullSampleInfo, GroupPCSign = GroupPCSign,
                       CorMethod = CorMethod)
  
  if(length(UsedModules)>1){
    ReorderIdxs <- order(ModuleOrder[UsedModules])
    
    return(list(ModuleMatrix = ModuleMatrix[ReorderIdxs,], SampleMatrix = SampleMatrix[ReorderIdxs,], ModuleSummary = ModuleSummary[ReorderIdxs],
                WeightList = WeightList[ReorderIdxs], PVVectMat = PVVectMat[ReorderIdxs,], OutLiersList = OutLiersList[ReorderIdxs], 
                genes_to_keep_counts = genes_to_keep_counts, genes_to_keep_fisher = genes_to_keep_fisher, AllPCA1 = AllPCA1[ReorderIdxs],
                OriginalOutliers = OriginalOutliers[ReorderIdxs], GeneCenters = GeneCenters, SampleCenters = SampleCenters, InputPars = InputParList))
  } else {
    return(list(ModuleMatrix = ModuleMatrix, SampleMatrix = SampleMatrix, ModuleSummary = ModuleSummary,
                WeightList = WeightList, PVVectMat = PVVectMat, OutLiersList = OutLiersList,
                genes_to_keep_counts = genes_to_keep_counts, genes_to_keep_fisher = genes_to_keep_fisher, OriginalOutliers = OriginalOutliers,
                AllPCA1 = AllPCA1, GeneCenters = GeneCenters, SampleCenters = SampleCenters, InputPars = InputParList))
  }
  
}

