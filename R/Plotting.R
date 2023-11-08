#' View genesets L1 and median PC1 compared to null distribution. 
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param Plot string, which parameter should be visualized ? Can be "L1", "PC1Median" or "both"
#'
#' @return
#' @export
#'
#' @examples

Plot.Genesets.vs.Sampled <- function(RomaData, Selected = NULL, Plot = "both"){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Geneset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "genesets selected"))
  }
  
  if(Plot %in% c("both", "L1")){
    B <- boxplot(lapply(RomaData$ModuleSummary[Selected], function(x){
      tSampleExp <- sapply(x$SampledExp, "[[", "ExpVar")
      if(!is.null(ncol(tSampleExp))){
        tSampleExp[1,]
      } else {
        tSampleExp
      }
    }), at = 1:length(Selected), las = 2, ylab = "Explained variance", main = "Selected genesets",
    names = unlist(lapply(RomaData$ModuleSummary[Selected], "[[", "ModuleName")),
    ylim = c(0,1))
    points(x = 1:length(Selected), y = RomaData$ModuleMatrix[Selected,1], pch = 20, col="red", cex = 2)
  }
  
  if(Plot %in% c("both", "PC1Median")){
    PlotData <- lapply(RomaData$ModuleSummary[Selected], function(x){
      sapply(x$SampledExp, "[[", "PC1Mean")
    })
  
    B <- boxplot(PlotData, at = 1:length(Selected), las = 2, ylab = "Median expression PC1", main = "Selected genesets",
                names = unlist(lapply(RomaData$ModuleSummary[Selected], "[[", "ModuleName")),
                ylim = range(c(unlist(PlotData), RomaData$ModuleMatrix[Selected,7]), na.rm = TRUE))
    points(x = 1:length(Selected), y = RomaData$ModuleMatrix[Selected,7], pch = 20, col="red", cex = 2)
  }
}



#' Plot genesets information. Creates a heatmap of sample scores for all selected gene sets. If groups and functions
#' are specified, creates other heatmaps with data aggregated by groups using the functions provided.
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param GenesetMargin scalat, integer. The number of rows used to draw the geneset names.
#' @param SampleMargin scalat, integer. The number of rows used to draw the sample names.
#' @param ColorGradient vector, string. The colors used for the heatmap.
#' @param cluster_cols boolean, should the samp^le be reordered according to the dendrogram?
#' @param GroupInfo vector, character. A vector describing the group association of each sample.
#' @param HMTite scalar, string. The title of the heatmap
#' @param AggByGroupsFL list, string. A list of function names (as strings) that will be used to aggregate the
#' geneset weights and produce additional heatmaps.
#' @param Normalize boolean, shuold weights be normalized to c(-1, 1) for each geneset
#' @param Transpose boolean, should the samples by plotted on the rows instead of the columns?
#' @param ZeroColor string, the color to use to mark the points closed to 0 (e.g., "#FFFFFF")
#'
#' @return A list of matrices containing the aggregated data used to produce the heatmaps.
#' @export
#'
#' @examples
#' 

RomaData <- rRoma.output
Selected <- shifted.modules
GenesetMargin <-  4
SampleMargin <- 4
Fixed <- 60
# ColorGradient <- colorRamps::blue2red(50)
cluster_cols <- TRUE
GroupInfo <- Type
HMTite <- "Selected Genesets"
Normalize <- FALSE
Transpose <- FALSE
ZeroColor <- NULL

Plot.Genesets.Samples <- function(RomaData, 
                                  Selected = NULL,
                                  GenesetMargin = 4, 
                                  SampleMargin = 4,
                                  Fixed = NULL,
                                  # ColorGradient = colorRamps::blue2red(50),
                                  cluster_cols = FALSE, 
                                  GroupInfo = NULL,
                                  HMTite = "Selected Genesets", 
                                  AggByGroupsFL = list(),
                                  Normalize = FALSE, 
                                  Transpose = FALSE, 
                                  ZeroColor = NULL){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(Selected) < 2){
    print("Heatmap is not clustering module activities when less than 2 modules are selected")
  }
  # if((length(ColorGradient) %% 2) != 0){
  #   stop("the length of ColorGradient MUST be a multiple of 2")
  # }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Geneset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "genesets selected"))
  }
  
  op <- par(mar=c(GenesetMargin, 5, 4, 2))
  
  par(op)
  
  if(!is.null(GroupInfo)){
    
    FoundSamp <- intersect(colnames(RomaData$SampleMatrix), names(GroupInfo))
    print(paste(length(FoundSamp), "samples have an associated group"))
    
    GroupInfo <- GroupInfo[FoundSamp]
    
    if(length(FoundSamp) > 2){
      if(!is.factor(GroupInfo)){
        GroupInfo <- as.factor(GroupInfo)
      }
      AddInfo <- data.frame(Groups = GroupInfo)
    } else {
      AddInfo = NULL
    }
    
  } else {
    AddInfo = NULL
  }
  
  PlotMat <- RomaData$SampleMatrix[Selected,]
  
  if(Normalize){
    
    PlotMat <- t(apply(PlotMat, 1, function(x){
      x[x>0] <- x[x>0]/max(x[x>0])
      x[x<0] <- -x[x<0]/min(x[x<0])
      x
    }))
    
  }
  
  if(is.null(dim(PlotMat))){
    dim(PlotMat) <- c(1, length(PlotMat))
    colnames(PlotMat) <- names(RomaData$SampleMatrix[Selected,])
    rownames(PlotMat) <- rownames(RomaData$SampleMatrix)[Selected]
    tClusCol <- FALSE
    tClusRow <- FALSE
  } else {
    tClusCol <- cluster_cols
    tClusRow <- TRUE
  }
  
  # Compute breaks
  
  if (is.null(Fixed)){
    
    # Classic palette BuPu, with 4 colors
    coul <- brewer.pal(9, "RdBu")
    ColorGradient <- colorRampPalette(coul)(60)
    
    MatRange <- range(PlotMat)
  } else {
    
    coul <- brewer.pal(9, "RdBu")
    ColorGradient <- colorRampPalette(coul)(Fixed)
    
    MatRange <- c(-1, 1)
  }
  
  tColorGradient <- ColorGradient
  BrkPoints <- rep(NA, length(ColorGradient) + 1)
  
  if(prod(MatRange) > 0){
    # Only positive or only negative values are observed
    
    if(MatRange[1]>0){
      # Only positive values are available
      if(is.null(ZeroColor)){
        nBrkPoints <- length(ColorGradient)/2 + 1
        tColorGradient <- ColorGradient[(length(ColorGradient)/2 +1):length(ColorGradient)]
        BrkPoints <- seq(to = MatRange[2], from = 0, length.out =  nBrkPoints)
      } else {
        nBrkPoints <- length(ColorGradient)/2 + 2
        tColorGradient <- c(ZeroColor, ColorGradient[(length(ColorGradient)/2 +1):length(ColorGradient)])
        BrkPoints <- seq(to = MatRange[2], from = 0, length.out = nBrkPoints)
      }
    } else {
      # Only negative values are available
      if(is.null(ZeroColor)){
        nBrkPoints <- length(ColorGradient)/2 + 1
        tColorGradient <- ColorGradient[1:(length(ColorGradient)/2)]
        BrkPoints <- seq(from = MatRange[1], to = 0, length.out = nBrkPoints)
      } else {
        nBrkPoints <- length(ColorGradient)/2 + 2
        tColorGradient <- c(ColorGradient[1:(length(ColorGradient)/2)], ZeroColor)
        BrkPoints <- seq(from = MatRange[1], to = 0, length.out = nBrkPoints)
      }
    }
    
  } else {
    # Both positive and negative values are observed
    if(is.null(ZeroColor)){
      nBrkPoints <- length(ColorGradient)/2 + 1
      LowBrks <- seq(from = MatRange[1], to = 0, length.out = nBrkPoints)
      UpBrks <- seq(to = MatRange[2], from = 0, length.out = nBrkPoints)
      BrkPoints <- c(LowBrks, UpBrks[-1])
    } else {
      nBrkPoints <- length(ColorGradient)/2 + 1
      LowBrks <- seq(from = MatRange[1], to = .5*MatRange[1]/nBrkPoints, length.out = nBrkPoints)
      UpBrks <- seq(to = MatRange[2], from = .5*MatRange[2]/nBrkPoints, length.out = nBrkPoints)
      BrkPoints <- c(LowBrks, UpBrks)
      tColorGradient <- c(
        ColorGradient[1:(length(ColorGradient)/2)],
        ZeroColor,
        ColorGradient[(length(ColorGradient)/2+1):length(ColorGradient)])
    }
    
    
  }
  
  Groups <- AddInfo$Groups
  GroupsColor <- brewer.pal(8, "Set2")[1:length(unique(AddInfo$Groups))]
  # Groups <- as.character(Groups)
  names(GroupsColor) <- unique(AddInfo$Groups)
  anno_colors <- list(Groups = GroupsColor)
  
  if(Transpose){
    pheatmap::pheatmap(t(PlotMat),
                       color = tColorGradient, 
                       breaks = BrkPoints,
                       main = HMTite,
                       cluster_rows = tClusCol, 
                       cluster_cols = tClusRow,
                       annotation_row = AddInfo)
  } else {
    pheatmap::pheatmap(PlotMat,
                       color = tColorGradient, 
                       breaks = BrkPoints,
                       main = HMTite,
                       cluster_cols = tClusCol, 
                       cluster_rows = tClusRow,
                       border_color = "NA",
                       annotation_col = AddInfo,
                       annotation_colors =anno_colors)
  }
  
  
  newCols <- colorRampPalette(grDevices::rainbow(length(unique(AddInfo$Groups))))
  mycolors <- newCols(length(unique(AddInfo$Groups)))
  names(mycolors) <- unique(AddInfo$Groups)
  mycolors <- list(category = mycolors)
  
  
  if(length(AggByGroupsFL)>0 & !is.null(AddInfo)){

    if(length(Selected) == 1){
      tDF <- data.frame(PlotMat[,FoundSamp])
      colnames(tDF) <- rownames(RomaData$SampleMatrix)[Selected]
      SplitData <- split(tDF, f=AddInfo$Groups)
    } else {
      SplitData <- split(data.frame(t(PlotMat[,FoundSamp])), f=AddInfo$Groups)
    }
    
    RetData <- list()
    
    for(i in 1:length(AggByGroupsFL)){
      
      Aggmat <- sapply(SplitData, function(x) {
        apply(x, 2, get(AggByGroupsFL[[i]]))
        })
      
      if(is.null(dim(Aggmat))){
        dim(Aggmat) <- c(1, length(Aggmat))
        colnames(Aggmat) <- names(SplitData)
        rownames(Aggmat) <- rownames(RomaData$SampleMatrix)[Selected]
        tClusCol <- FALSE
        tClusRow <- FALSE
      } else {
        tClusCol <- cluster_cols
        tClusRow <- TRUE
      }
      
      
      # Compute breaks
      
      MatRange <- range(Aggmat)
      BrkPoints <- rep(NA, length(ColorGradient) + 1)
      
      tColorGradient <- ColorGradient
      
      if(prod(MatRange) > 0){
        # Only positive or negative values are observed
        
        # Only positive or negative values are observed
        if(MatRange[1]>0){
          # Only positive values are available
          if(is.null(ZeroColor)){
            nBrkPoints <- length(ColorGradient)/2 + 1
            tColorGradient <- ColorGradient[(length(ColorGradient)/2 +1):length(ColorGradient)]
            BrkPoints <- seq(to = MatRange[2], from = 0, length.out =  nBrkPoints)
          } else {
            nBrkPoints <- length(ColorGradient)/2 + 2
            tColorGradient <- c(ZeroColor, ColorGradient[(length(ColorGradient)/2 +1):length(ColorGradient)])
            BrkPoints <- seq(to = MatRange[2], from = 0, length.out = nBrkPoints)
          }
        } else {
          # Only negative values are available
          if(is.null(ZeroColor)){
            nBrkPoints <- length(ColorGradient)/2 + 1
            tColorGradient <- ColorGradient[1:(length(ColorGradient)/2)]
            BrkPoints <- seq(from = MatRange[1], to = 0, length.out = nBrkPoints)
          } else {
            nBrkPoints <- length(ColorGradient)/2 + 2
            tColorGradient <- c(ColorGradient[1:(length(ColorGradient)/2)], ZeroColor)
            BrkPoints <- seq(from = MatRange[1], to = 0, length.out = nBrkPoints)
          }
        }
        
      } else {
        # Both positive and negative values are observed
        if(is.null(ZeroColor)){
          nBrkPoints <- length(ColorGradient)/2 + 1
          LowBrks <- seq(from = MatRange[1], to = 0, length.out = nBrkPoints)
          UpBrks <- seq(to = MatRange[2], from = 0, length.out = nBrkPoints)
          BrkPoints <- c(LowBrks, UpBrks[-1])
        } else {
          nBrkPoints <- length(ColorGradient)/2 + 1
          LowBrks <- seq(from = MatRange[1], to = .5*MatRange[1]/nBrkPoints, length.out = nBrkPoints)
          UpBrks <- seq(to = MatRange[2], from = .5*MatRange[2]/nBrkPoints, length.out = nBrkPoints)
          BrkPoints <- c(LowBrks, UpBrks)
          tColorGradient <- c(
            ColorGradient[1:(length(ColorGradient)/2)],
            ZeroColor,
            ColorGradient[(length(ColorGradient)/2+1):length(ColorGradient)])
        }
        
        
      }
      
      
      if(Transpose){
        pheatmap::pheatmap(t(Aggmat),
                           color = tColorGradient, breaks = BrkPoints,
                           main = paste(HMTite, "/", AggByGroupsFL[[i]]),
                           cluster_rows = tClusCol, cluster_cols = tClusRow)
      } else {
        pheatmap::pheatmap(Aggmat,
                           color = tColorGradient, breaks = BrkPoints,
                           main = paste(HMTite, "/", AggByGroupsFL[[i]]),
                           cluster_cols = tClusCol, cluster_rows = tClusRow)
      }
      
      RetData[[length(RetData)+1]] <- Aggmat

    }
    
    names(RetData) <- AggByGroupsFL

    return(RetData)
  }
  
}








RomaData <- Okuda.ROMA_output
ExpressionMatrix <- as.data.frame(Saint_Criq@RNA_matrix)
LogExpression <- FALSE
Selected <- NULL
PlotGenes <- 40
PlotWeightSign <- FALSE

#' Plot gene weight across selected samples
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param PlotGenes scalar, numeric. The number of genes to plot
#' @param ExpressionMatrix matrix, numeric. The expression matrix used to produce gene expression boxplot. If NULL (default), no gene expression information is reported
#' @param LogExpression boolean, should gene expression be logtransformed?
#' @param PlotWeightSign boolean, should the sign of the genes weight be used to color the plots?
#'
#' @return
#' @export
#'
#' @examples
PlotGeneWeight <- function(RomaData, PlotGenes = 40,
                           ExpressionMatrix = NULL, LogExpression = TRUE,
                           Selected = NULL, PlotWeightSign = FALSE){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Geneset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset(s) selected"))
  }
  
  
  for(i in Selected){
    
    FiltGenes <- RomaData$ModuleSummary[[i]]$GeneWeight * RomaData$ModuleSummary[[i]]$CorrectSign1
    names(FiltGenes) <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    GeneNames <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    nGenes <- min(PlotGenes, length(GeneNames))
    
    GMTWei = sign(RomaData$ModuleSummary[[i]]$GMTWei)
    names(GMTWei) = names(FiltGenes)
    GMTWei[is.na(GMTWei)] <- 0
    
    DF <- data.frame(cbind(names(FiltGenes), FiltGenes,
                           GMTWei),
                     row.names = NULL)
    colnames(DF) <- c("Gene", "Value", "Wei")
    
    DF$Value <- as.numeric(as.character(DF$Value))
    DF <- DF[order(abs(DF$Value), decreasing = TRUE)[1:nGenes],]
    DF <- DF[order(DF$Value), ]
    
    DF$Gene <- factor(as.character(DF$Gene), levels = DF$Gene)
    DF$Wei <- factor(as.character(DF$Wei), levels = c("-1", "0", "1"))
    
    if(PlotWeightSign){
      if(any(as.character(DF$Wei)=="0")){
        p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Gene, x = Value, color = Wei)) +
          ggplot2::scale_color_manual(values = c("red", "black", "blue"), guide=FALSE)
      } else {
        p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Gene, x = Value, color = Wei)) +
          ggplot2::scale_color_manual(values = c("red", "blue"), guide=FALSE)
      }
      
      
    } else {
      p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Gene, x = Value))
    }
    
     p <- p + ggplot2::geom_point() +
      # ggplot2::scale_y_continuous(breaks = DF$Position[1:nGenes], labels = DF$Gene[1:nGenes]) +
      # ggplot2::scale_y_discrete(breaks=levels(DF$Gene)) +
      ggplot2::labs(x = "Gene weight", y = "Gene") +
      ggplot2::theme(panel.grid.minor = ggplot2::element_line(colour = NA))
    
    if(is.null(ExpressionMatrix)){
      
      print(p + ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                                       "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                                       "L1/L2=", signif(RomaData$ModuleMatrix[i,4], 3))))
      
    } else {
      Data_to_reshape <- ExpressionMatrix
      # Data_to_reshape$Gene <- rownames(ExpressionMatrix)
      ReshapedData <- reshape::melt(Data_to_reshape[as.character(DF$Gene),])
      colnames(ReshapedData) <- c("X1", "X2", "value")
      
      # ReshapedData <- as.data.frame(tidyr::pivot_longer(ExpressionMatrix[as.character(DF$Gene), ], names_to = "X1", values_to = "X2", cols = colnames(ExpressionMatrix)))
      # rownames(ReshapedData) <- 
      ReshapedData$X1 <- factor(as.character(ReshapedData$X1), levels = levels(DF$Gene))
      ReshapedData <- cbind(ReshapedData, GMTWei[as.character(ReshapedData$X1)])
      colnames(ReshapedData)[4] <- "Wei"
      ReshapedData$Wei <- factor(as.character(ReshapedData$Wei), levels = c("-1", "0", "1"))
      
      if(PlotWeightSign){
        if(any(as.character(ReshapedData$Wei)=="0")){
          p1 <- ggplot2::ggplot(ReshapedData, ggplot2::aes(x=X1, y=value, fill=Wei)) +
            ggplot2::scale_fill_manual(values = c("red", "white", "blue"), guide=FALSE)
        } else {
          p1 <- ggplot2::ggplot(ReshapedData, ggplot2::aes(x=X1, y=value, fill=Wei)) +
            ggplot2::scale_fill_manual(values = c("red", "blue"), guide=FALSE)
        }
      } else {
        p1 <- ggplot2::ggplot(ReshapedData, ggplot2::aes(x=X1, y=value))
      }
      
      p1 <- p1 + ggplot2::geom_boxplot() + ggplot2::coord_flip() + 
        ggplot2::labs(x = "Gene", y = "Expression")
      
      if(LogExpression){
        p1 <- p1 + ggplot2::scale_y_log10()
      }
      
      gridExtra::grid.arrange(p, p1, ncol=2, top=paste(RomaData$ModuleSummary[[i]]$ModuleName,
                                                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                                                               "L1/L2=", signif(RomaData$ModuleMatrix[i,4], 3)))
      
      # boxplot(t(ExpressionMatrix[as.character(DF$Gene[1:nGenes])[order(DF$Position[1:nGenes])], ]))
      
    }
    
  }
  
}

















#' Plot sample score across selected samples
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param PlotSamples scalar, numeric. The number of samples to plot
#' @param ExpressionMatrix matrix, numeric. The expression matrix used to produce gene expression boxplot. If NULL (default), no gene expression information is reported
#' @param LogExpression boolean, should gene expression be logtransformed?
#' @param Selected vector, integer. The position of the genesets to plot
#' @param FullExpDist boolean, should the all the genes be used when showing gene expression
#' distribution (TRUE) or only the genes of the geneset under consideration (FALSE) ?
#'
#' @return
#' @export
#'
#' @examples
PlotSampleProjections <- function(RomaData, PlotSamples = 40,
                                ExpressionMatrix = NULL, LogExpression = TRUE,
                                Selected = NULL, FullExpDist = FALSE){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Geneset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset(s) selected"))
  } 
  
  for(i in Selected){
    
    FiltSamp <- RomaData$ModuleSummary[[i]]$SampleScore * RomaData$ModuleSummary[[i]]$CorrectSign1
    # names(FiltSamp) <- colnames(RomaData$SampleMatrix)
    
    GeneNames <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    nSamples <- min(PlotSamples, length(FiltSamp))
    
    DF <- data.frame(cbind(names(FiltSamp),
                           FiltSamp),
                     row.names = NULL)
    colnames(DF) <- c("Samples", "Value")
    
    DF$Value <- as.numeric(as.character(DF$Value))
    DF <- DF[order(abs(DF$Value), decreasing = TRUE)[1:nSamples],]
    DF <- DF[order(DF$Value), ]
    
    DF$Samples <- factor(as.character(DF$Samples), levels = DF$Samples)
    
    p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Samples, x = Value)) + ggplot2::geom_point() +
      # ggplot2::scale_y_continuous(breaks = DF$Position[1:nGenes], labels = DF$Gene[1:nGenes]) +
      # ggplot2::scale_y_discrete(breaks=levels(DF$Gene)) +
      ggplot2::labs(x = "Sample Score", y = "Sample") +
      ggplot2::theme(panel.grid.minor = ggplot2::element_line(colour = NA))
    
    if(is.null(ExpressionMatrix)){
      
      print(p + ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                                       "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                                       "L1/L2=", signif(RomaData$ModuleMatrix[i,4], 3))))
      
    } else {
      
      if(FullExpDist){
        ReshapedData <- reshape::melt(ExpressionMatrix[, as.character(DF$Samples)])
        colnames(ReshapedData) <- c("X2", "value")
      } else {
        
        ReshapedData <- reshape::melt(ExpressionMatrix[GeneNames, as.character(DF$Samples)])
        colnames(ReshapedData) <- c("X2", "value")
        
      }
      
      
      ReshapedData$X2 <- factor(as.character(ReshapedData$X2), levels = levels(DF$Samples))
      
      p1 <- ggplot2::ggplot(ReshapedData, ggplot2::aes(x=X2, y=value)) + 
        ggplot2::geom_boxplot() + ggplot2::coord_flip() +
        ggplot2::labs(x = "Sample", y = "Expression")
      
      if(LogExpression){
        p1 <- p1 + ggplot2::scale_y_log10()
      }
      
      gridExtra::grid.arrange(p, p1, ncol=2, top=paste(RomaData$ModuleSummary[[i]]$ModuleName,
                                                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                                                               "L1/L2=", signif(RomaData$ModuleMatrix[i,4], 3)))
      
      # boxplot(t(ExpressionMatrix[as.character(DF$Gene[1:nGenes])[order(DF$Position[1:nGenes])], ]))
      
    }
    
  }
    
}



  
#' Plot PC projections score across selected samples
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param PlotPCProj vector, string. The plotting modality of projections. It can containing any combination of the following strings: 'Points', 'Density', or 'Bins'.
#' Any other value will result in projections not being plotted
#'
#' @return
#' @export
#'
#' @examples
PlotPCProjections <- function(RomaData, Selected = NULL, PlotPCProj = 'Points'){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Geneset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset(s) selected"))
  } 
  
  for(i in Selected){
    
    PrjList <- lapply(RomaData$ModuleSummary[[i]]$SampledExp, "[[", "GenesWei")
    
    
    if(any(sapply(PrjList, is.null)) | !all(PlotPCProj != "none")){
      print(paste("Information on reference distribution not available. Skipping ", RomaData$ModuleSummary[[i]]$ModuleName, ". Consider changing FullSampleInfo to TRUE in rRoma.R", sep = "" ))
      next()
    }
    
    PC1Sample <- unlist(lapply(PrjList, function(x){if(all(!is.na(x))) x[,1]}), use.names = FALSE)
    PC2Sample <- unlist(lapply(PrjList, function(x){if(all(!is.na(x))) x[,2]}), use.names = FALSE)
    
    PC1Data <- RomaData$ModuleSummary[[i]]$PCABase$x[,1]*RomaData$ModuleSummary[[i]]$CorrectSign1
    PC2Data <- RomaData$ModuleSummary[[i]]$PCABase$x[,2]*RomaData$ModuleSummary[[i]]$CorrectSign2
    
    DF <- data.frame(PC1 = c(PC1Sample, PC1Data), PC2 = c(PC2Sample, PC2Data),
                     Source = c(rep("Sampling", length(PC1Sample)),
                                rep("Data", length(PC1Data))
                     )
    )
    
    XLims <- quantile(PC1Sample, c(.01, .99))
    XLims[1] <- min(XLims[1], min(PC1Data))
    XLims[2] <- max(XLims[2], max(PC1Data))
    
    YLims <- quantile(PC2Sample, c(.01, .99))
    YLims[1] <- min(YLims[1], min(PC2Data))
    YLims[2] <- max(YLims[2], max(PC2Data))
    
    if(any(PlotPCProj == 'Points')){
      p <- ggplot2::ggplot(DF, ggplot2::aes(x = PC1, y = PC2, colour = Source, alpha=Source)) + ggplot2::geom_point() +
        ggplot2::scale_alpha_manual(values=c(Data=1, Sampling=.2), guide=FALSE) +
        ggplot2::scale_x_continuous(limits = XLims) +
        ggplot2::scale_y_continuous(limits = YLims) +
        ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                               "L1/L2=", signif(RomaData$ModuleMatrix[i,4], 3))) +
        ggplot2::geom_hline(yintercept=0) + ggplot2::geom_vline(xintercept=0)
      
      print(p)
    }
    
    
    if(any(PlotPCProj == 'Density')){
      p <- ggplot2::ggplot(DF[DF$Source=="Sampling",], ggplot2::aes(x = PC1, y = PC2)) +
        ggplot2::stat_density_2d(ggplot2::aes(fill = ..density..), geom="raster", contour = FALSE, n = 250) +
        ggplot2::scale_x_continuous(limits = XLims) +
        ggplot2::scale_y_continuous(limits = YLims) +
        ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                               "L1/L2=", signif(RomaData$ModuleMatrix[i,4], 3))) +
        ggplot2::geom_hline(yintercept=0) + ggplot2::geom_vline(xintercept=0) +
        ggplot2::geom_point(data = DF[DF$Source=="Data",], mapping = ggplot2::aes(x = PC1, y = PC2), color="red")
      
      print(p)
    }
    
    
    
    if(any(PlotPCProj == 'Bins')){
      p <- ggplot2::ggplot(DF[DF$Source=="Sampling",], ggplot2::aes(x = PC1, y = PC2)) +
        ggplot2::geom_bin2d(bins=75) +
        ggplot2::scale_x_continuous(limits = XLims) +
        ggplot2::scale_y_continuous(limits = YLims) +
        ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                               "L1/L2=", signif(RomaData$ModuleMatrix[i,4], 3))) +
        ggplot2::geom_hline(yintercept=0) + ggplot2::geom_vline(xintercept=0) +
        ggplot2::geom_point(data = DF[DF$Source=="Data",], mapping = ggplot2::aes(x = PC1, y = PC2), color="red")
      
      print(p)
    }
    
  }
  
}


#' Plot several plots to investigate a specific gene
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param GeneName string, the name of the gene to study
#' @param ExpressionMatrix matrix, the input data used to compute rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param GroupInfo vector, character. A vector describing the group association of each sample.
#' @param PlotExpVar boolean. Do you want to plot the expression variance of the gene compared to all other genes in the dataset ?
#' @param PlotExpMean boolean. Do you want to plot the expression mean of the gene compared to all other genes in the dataset ?
#' @param PlotExpDist boolean. Do you want to plot the distribution of the expression of the gene in all samples ? 
#' If Groups are supplied, also shows potential significant differences between groups for this gene 
#' @param PlotWeight boolean. Do you want to plot the gene weight for the considered gene sets, compared to all other gene weights ?
#' @param PlotExpSampleScore boolean. Do you want to plot gene expression against sample score for the considered gene sets ?
#'
#' @return
#' @export
#'
#' @examples
ExploreGeneProperties <- function(
  RomaData,
  GeneName,
  ExpressionMatrix,
  Selected = NULL,
  GroupInfo = NULL,
  PlotExpVar = TRUE,
  PlotExpMean= TRUE,
  PlotExpDist = TRUE,
  PlotWeight = FALSE,
  PlotExpSampleScore = FALSE
){
  
  if(PlotExpVar){
    VarVect <- apply(ExpressionMatrix, 1, var)
    
    if(!is.null(GroupInfo)){
      ExpMatByGroup <- split(data.frame(t(ExpressionMatrix)),
                             f=GroupInfo[colnames(ExpressionMatrix)])
      
      VarByGroup <- sapply(ExpMatByGroup, function(x) {
        apply(x, 2, var)
      })
      
      rownames(VarByGroup) <- names(VarVect)
      
      VarByGroup <- cbind(VarByGroup, Total = VarVect)
      
    } else {
      
      VarByGroup <- data.frame(X1 = names(VarVect), Total = VarVect)
      
    }
    
   
    VarByGroup.Melt <- reshape::melt(VarByGroup)
    colnames(VarByGroup.Melt) <- c("Gene", "Condition", "value")
    
    if(is.factor(VarByGroup.Melt$Condition)){
      VarByGroup.Melt$Condition <- relevel(VarByGroup.Melt$Condition, "Total")
    } else {
      VarByGroup.Melt$Condition <- relevel(factor(VarByGroup.Melt$Condition), "Total")
    }
    
    VarByGroup.Melt.Selected <- VarByGroup.Melt[VarByGroup.Melt$Gene == GeneName,]
    
    p <- ggplot2::ggplot(data = VarByGroup.Melt, mapping = ggplot2::aes(x=value, fill = Condition)) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(data = VarByGroup.Melt.Selected,
                          mapping = ggplot2::aes(xintercept = value), color="black") +
      ggplot2::scale_x_log10() +
      ggplot2::facet_wrap(~Condition) +
      ggplot2::labs(x = "Expression variance (line indicates the target gene)", title = GeneName) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::guides(fill = "none")
    
    print(p)
  }
  
  
  
  
  
  
  
  
  
  
  if(PlotExpMean){
    MeanVect <- apply(ExpressionMatrix, 1, mean)
    
    if(!is.null(GroupInfo)){
      ExpMatByGroup <- split(data.frame(t(ExpressionMatrix)),
                             f=GroupInfo[colnames(ExpressionMatrix)])
      
      MeanByGroup <- sapply(ExpMatByGroup, function(x) {
        apply(x, 2, mean)
      })
      
      rownames(MeanByGroup) <- names(MeanVect)
      
      MeanByGroup <- cbind(MeanByGroup, Total = MeanVect)
      
    } else {
      
      MeanByGroup <- data.frame(X1 = names(MeanVect), Total = MeanVect)
      
    }
    
    
    MeanByGroup.Melt <- reshape::melt(MeanByGroup)
    colnames(MeanByGroup.Melt) <- c("Gene", "Condition", "value")
    
    if(is.factor(MeanByGroup.Melt$Condition)){
      MeanByGroup.Melt$Condition <- relevel(MeanByGroup.Melt$Condition, "Total")
    } else {
      MeanByGroup.Melt$Condition <- relevel(factor(MeanByGroup.Melt$Condition), "Total")
    }
    
    
    MeanByGroup.Melt.Selected <- MeanByGroup.Melt[MeanByGroup.Melt$Gene == GeneName,]
    
    p <- ggplot2::ggplot(data = MeanByGroup.Melt, mapping = ggplot2::aes(x=value, fill = Condition)) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(data = MeanByGroup.Melt.Selected,
                          mapping = ggplot2::aes(xintercept = value), color="black") +
      ggplot2::facet_wrap(~Condition) +
      ggplot2::labs(x = "Mean expression (line indicates the target gene)", title = GeneName) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::guides(fill = "none")
    
    print(p)
  }
  
  
  
  
  
  
  
  
  if(PlotExpDist){
    ExpData <- data.frame(
      Sample = colnames(ExpressionMatrix),
      Exp = as.numeric(ExpressionMatrix[GeneName,]),
      Condition = "Total",
      stringsAsFactors = FALSE
      )
    
    ToDisplay <- list()
    
    if(!is.null(GroupInfo)){
      ExpData <- rbind(
        ExpData,
        data.frame(
          Sample = colnames(ExpressionMatrix),
          Exp = as.numeric(ExpressionMatrix[GeneName,]),
          Condition = GroupInfo[colnames(ExpressionMatrix)],
          stringsAsFactors = FALSE
        )
      )
      ExpData <- ExpData[!is.na(ExpData$Condition), ]
      if(is.factor(ExpData$Condition)){
        ExpData$Condition <- relevel(ExpData$Condition, "Total")
      } else {
        ExpData$Condition <- relevel(factor(ExpData$Condition), "Total")
      }
      PWComp <- pairwise.wilcox.test(ExpData$Exp, ExpData$Condition, p.adjust.method = "BH")$p.value
      ExtPWComp <- matrix(rep(NA, length(unique(ExpData$Condition))^2), nrow = length(unique(ExpData$Condition)))
      colnames(ExtPWComp) <- sort(unique(ExpData$Condition))
      rownames(ExtPWComp) <- colnames(ExtPWComp)
      
      ExtPWComp[rownames(PWComp), colnames(PWComp)] <- PWComp
      
      ToDisplay <- apply(which(ExtPWComp <= .05, arr.ind = TRUE), 1, list)
      ToDisplay <- lapply(ToDisplay, function(x){as.integer(unlist(x))})
      ToDisplay <- ToDisplay[sapply(ToDisplay, function(x){all(x != 1)})] # To prevent comparison with total
    }
    
  
  
    
    
    p <- ggplot2::ggplot(data = ExpData, mapping = ggplot2::aes(y = Exp, x = Condition, fill=Condition)) + 
      ggplot2::geom_boxplot() +
      ggplot2::labs(y = "Expression", title = GeneName) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        axis.text.x = ggplot2::element_blank()) +
      ggsignif::geom_signif(comparisons = ToDisplay,
                            map_signif_level=TRUE, test = "wilcox.test", step_increase = .1)
    
    print(p)
  }
  
  
  
  if(PlotWeight){
    if(is.null(Selected)){
      Selected <- 1:nrow(RomaData$SampleMatrix)
    }
    
    ModulesGenes <- lapply(RomaData$ModuleSummary[Selected], "[[", "UsedGenes")
    names(ModulesGenes) <- sapply(RomaData$ModuleSummary[Selected], "[[", "ModuleName")
    
    WhichModules <- sapply(ModulesGenes, function(x){
      GeneName %in% x
    })
    
    Found <- sapply(RomaData$ModuleSummary, "[[", "ModuleName") %in%
      names(which(WhichModules))
    
    print(paste("Gene found in", sum(Found), "modules"))
    
    if(sum(Found) < 1){
      return()
    }
    
    WeiList <- lapply(RomaData$ModuleSummary[Found], function(x){
      x$GeneWeight * x$CorrectSign1
    })
    
    names(WeiList) <- sapply(RomaData$ModuleSummary[Found], "[[", "ModuleName")
    
    WeiList <- lapply(1:length(WeiList), function(i) {
      x <- WeiList[[i]]
      return(
        data.frame(Wei = c(x, abs(x)),
                   Name = rep(names(x), 2),
                   Order = c(rank(x, ties.method = "average"),
                             rank(abs(x), ties.method = "average"))/length(x),
                   Type = rep(c("Signed", "Absolute"), each = length(x)),
                   Module = rep(names(WeiList)[i], 2*length(x)),
                   stringsAsFactors = FALSE)
      )
    })
    
    
    Wei.DF <- do.call(rbind, WeiList)
    Wei.DF.Selected <- Wei.DF[Wei.DF$Name == GeneName,]
    
    
    p <- ggplot2::ggplot(data = Wei.DF,
                         mapping = ggplot2::aes(y=Wei, x="")) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_grid(Module ~ Type, scales = "free") + 
      ggplot2::geom_point(data = Wei.DF.Selected,
                          mapping = ggplot2::aes(y=Wei, x=""),
                          inherit.aes = FALSE,
                          color = "red",
                          size = 3) +
      ggplot2::labs(y = "Weight", x = "", title = GeneName) + 
      ggplot2::geom_text(data = Wei.DF.Selected,
                         mapping = ggplot2::aes(x = "", y = Wei,
                                                label = paste(signif(100*Order, 2), "%")
                         ),
                         hjust = 0, nudge_x = 0.05, size = 7,
                         inherit.aes = FALSE) + 
      ggplot2::geom_text(data = Wei.DF.Selected,
                         mapping = ggplot2::aes(x = "", y = Wei,
                                                label = signif(Wei, 4)
                         ),
                         hjust = 1, nudge_x = -0.05, size = 7,
                         inherit.aes = FALSE) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    print(p)
  }
  
  
  
  if(PlotExpSampleScore){
    GeneExp <- ExpressionMatrix[GeneName,]
    
    ModScoreList <- lapply(RomaData$ModuleSummary[Found], function(x){
      x$SampleScore * x$CorrectSign1
    })
    names(ModScoreList) <- sapply(RomaData$ModuleSummary[Found], "[[", "ModuleName")
    
    ModMat <- do.call(cbind, ModScoreList)
    Mod.DF <- data.frame(ModMat,
                         Sample = rownames(ModMat),
                         stringsAsFactors = FALSE)
    
    Mod.DF <- reshape::melt(Mod.DF, id.vars = "Sample")
    
    if(!is.null(GroupInfo)){
      
      Mod.DF <- data.frame(Mod.DF,Exp =  as.numeric(GeneExp[Mod.DF$Sample]), Groups = GroupInfo[Mod.DF$Sample])
      
      SplDS <- split(Mod.DF, Mod.DF$variable)
      
      SpeCor <- sapply(SplDS, function(x){
        cor(x$value, x$Exp, method = "spe")
      })
      
      Mod.DF$Lab <- paste(as.character(Mod.DF$variable), signif(SpeCor[as.character(Mod.DF$variable)], 3), sep = " / Cor=")
      
      p <- ggplot2::ggplot(data = Mod.DF, mapping = ggplot2::aes(x = value, y = Exp)) +
        ggplot2::geom_smooth(color = "black") +
        ggplot2::geom_point(mapping = ggplot2::aes(color = Groups)) + 
        ggplot2::facet_wrap(~Lab) +
        ggplot2::labs(x = "Sample Score", y = "Gene expression", title = GeneName) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      
    } else {
      
      Mod.DF <- data.frame(Mod.DF, Exp = GeneExp[Mod.DF$Sample])
      
      SplDS <- split(Mod.DF, Mod.DF$variable)
      
      SpeCor <- sapply(SplDS, function(x){
        cor(x$value, x$Exp, method = "spe")
      })
      
      Mod.DF$Lab <- paste(as.character(Mod.DF$variable), signif(SpeCor[as.character(Mod.DF$variable)], 3), sep = " / Cor=")
      
      p <- ggplot2::ggplot(data = Mod.DF, mapping = ggplot2::aes(x = value, y = Exp)) +
        ggplot2::geom_smooth() + 
        ggplot2::geom_point() + 
        ggplot2::facet_wrap(~Lab) +
        ggplot2::labs(x = "Sample Score", y = "Gene expression", title = GeneName) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
    }
    
    print(p)
  }
  
}
                      
#' Plot several plots to investigate outliers
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param ExpressionMatrix matrix, the input data used to compute rRoma
#' @param Selected vector, integer. The position of the genesets to study
#' @param Plot_L1 boolean, to you want to plot the amount of variance explained by PC1 when a specific gene is removed ?
#' @param Plot_Weights boolean, do you want to plot the distribution of gene weights ?
#' @param Plot_Projection boolean, do you want to plot the outlier genes in the PC1/PC2 space ?
#' 
#' @return
#' @export
#'
#' @examples
PlotOutliers <- function(
  RomaData,
  ExpressionMatrix,
  Selected = NULL,
  Plot_L1 = TRUE,
  Plot_Weights = FALSE,
  Plot_Projection = FALSE){

  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(Plot_L1){
    for(i in Selected){
      to_plot <- data.frame("Genes" = RomaData$ModuleSummary[[i]]$OriginalGenes, 
                            "L1_Data" = RomaData$AllPCA1[[i]], "Status" = rep("Kept_Base", length(RomaData$ModuleSummary[[i]]$OriginalGenes)))
      to_plot[to_plot$Genes %in% RomaData$OriginalOutliers[[i]], c("Status")] <- "Kept_%_outliers"
      to_plot[to_plot$Genes %in% RomaData$OutLiersList[[i]], c("Status")]<- "Outlier"
      to_plot[to_plot$Genes %in% intersect(RomaData$genes_to_keep_fisher, RomaData$OriginalOutliers[[i]]), c("Status")] <- "Kept_Fisher"
      to_plot[to_plot$Genes %in% intersect(RomaData$genes_to_keep_counts, RomaData$OriginalOutliers[[i]]), c("Status")] <- "Kept_Counts"
      
      
      B <- boxplot(x = to_plot$L1_Data, at = 1, horizontal = FALSE, ylab = "Variance explained by PC1", main = RomaData$ModuleSummary[[i]]$ModuleName)
      if(sum(to_plot$Status == "Outlier") > 0){
        points(y= to_plot[to_plot$Status == "Outlier", c("L1_Data")], x=rep(1, length(to_plot[to_plot$Status == "Outlier", c("L1_Data")])), col='red', pch=20)
        text(y= to_plot[to_plot$Status == "Outlier", c("L1_Data")], x = rep(1,length(to_plot[to_plot$Status == "Outlier", c("L1_Data")])), 
            labels = to_plot[to_plot$Status == "Outlier", c("Genes")], pos = 4)
      }
      
      if(sum(to_plot$Status == "Kept_Counts") > 0){
        points(y= to_plot[to_plot$Status == "Kept_Counts", c("L1_Data")], x=rep(1, length(to_plot[to_plot$Status == "Kept_Counts", c("L1_Data")])), col='blue', pch=20)
        text(y= to_plot[to_plot$Status == "Kept_Counts", c("L1_Data")], x = rep(1,length(to_plot[to_plot$Status == "Kept_Counts", c("L1_Data")])), 
            labels = to_plot[to_plot$Status == "Kept_Counts", c("Genes")], pos = 4)
      }
      
      if(sum(to_plot$Status == "Kept_Fisher") > 0){
        points(y= to_plot[to_plot$Status == "Kept_Fisher", c("L1_Data")], x=rep(1, length(to_plot[to_plot$Status == "Kept_Fisher", c("L1_Data")])), col='green', pch=20)
        text(y= to_plot[to_plot$Status == "Kept_Fisher", c("L1_Data")], x = rep(1,length(to_plot[to_plot$Status == "Kept_Fisher", c("L1_Data")])), 
            labels = to_plot[to_plot$Status == "Kept_Fisher", c("Genes")], pos = 4)
      }
      
      if(sum(to_plot$Status == "Kept_%_outliers") > 0){
        points(y= to_plot[to_plot$Status == "Kept_%_outliers", c("L1_Data")], x=rep(1, length(to_plot[to_plot$Status == "Kept_%_outliers", c("L1_Data")])), col='purple', pch=20)
        text(y= to_plot[to_plot$Status == "Kept_%_outliers", c("L1_Data")], x = rep(1,length(to_plot[to_plot$Status == "Kept_%_outliers", c("L1_Data")])), 
            labels = to_plot[to_plot$Status == "Kept_%_outliers", c("Genes")], pos = 4)
      }
      
      legend("left", pch = c(20, 20, 20, 20), col=c('red', 'blue', 'green', "purple"), legend = c("Outlier(s)", "Kept Counts", "Kept Fisher", "Kept % Outliers"))
      
    }
  }
  
  if(Plot_Weights){
    for(i in Selected){
      to_plot <- data.frame("Genes" = RomaData$ModuleSummary[[i]]$OriginalGenes, 
                            "Weights" = RomaData$ModuleSummary[[i]]$GeneWeightUnf, "Status" = rep("Kept_Base", length(RomaData$ModuleSummary[[i]]$OriginalGenes)))
      to_plot[to_plot$Genes %in% RomaData$OriginalOutliers[[i]], c("Status")] <- "Kept_%_outliers"
      to_plot[to_plot$Genes %in% RomaData$OutLiersList[[i]], c("Status")]<- "Outlier"
      to_plot[to_plot$Genes %in% intersect(RomaData$genes_to_keep_fisher, RomaData$OriginalOutliers[[i]]), c("Status")] <- "Kept_Fisher"
      to_plot[to_plot$Genes %in% intersect(RomaData$genes_to_keep_counts, RomaData$OriginalOutliers[[i]]), c("Status")] <- "Kept_Counts"
      
      
      B <- boxplot(x = to_plot$Weights, at = 1, horizontal = FALSE, ylab = "Gene Weight", main = RomaData$ModuleSummary[[i]]$ModuleName)
      if(sum(to_plot$Status == "Outlier") > 0){
        points(y= to_plot[to_plot$Status == "Outlier", c("Weights")], x=rep(1, length(to_plot[to_plot$Status == "Outlier", c("Weights")])), col='red', pch=20)
        text(y= to_plot[to_plot$Status == "Outlier", c("Weights")], x = rep(1,length(to_plot[to_plot$Status == "Outlier", c("Weights")])), 
             labels = to_plot[to_plot$Status == "Outlier", c("Genes")], pos = 4)
      }
      
      if(sum(to_plot$Status == "Kept_Counts") > 0){
        points(y= to_plot[to_plot$Status == "Kept_Counts", c("Weights")], x=rep(1, length(to_plot[to_plot$Status == "Kept_Counts", c("Weights")])), col='blue', pch=20)
        text(y= to_plot[to_plot$Status == "Kept_Counts", c("Weights")], x = rep(1,length(to_plot[to_plot$Status == "Kept_Counts", c("Weights")])), 
             labels = to_plot[to_plot$Status == "Kept_Counts", c("Genes")], pos = 4)
      }
      
      if(sum(to_plot$Status == "Kept_Fisher") > 0){
        points(y= to_plot[to_plot$Status == "Kept_Fisher", c("Weights")], x=rep(1, length(to_plot[to_plot$Status == "Kept_Fisher", c("Weights")])), col='green', pch=20)
        text(y= to_plot[to_plot$Status == "Kept_Fisher", c("Weights")], x = rep(1,length(to_plot[to_plot$Status == "Kept_Fisher", c("Weights")])), 
             labels = to_plot[to_plot$Status == "Kept_Fisher", c("Genes")], pos = 4)
      }
      
      if(sum(to_plot$Status == "Kept_%_outliers") > 0){
        points(y= to_plot[to_plot$Status == "Kept_%_outliers", c("Weights")], x=rep(1, length(to_plot[to_plot$Status == "Kept_%_outliers", c("Weights")])), col='purple', pch=20)
        text(y= to_plot[to_plot$Status == "Kept_%_outliers", c("Weights")], x = rep(1,length(to_plot[to_plot$Status == "Kept_%_outliers", c("Weights")])), 
             labels = to_plot[to_plot$Status == "Kept_%_outliers", c("Genes")], pos = 4)
      }
      
      legend("left", pch = c(20, 20, 20, 20), col=c('red', 'blue', 'green', "purple"), legend = c("Outlier(s)", "Kept Counts", "Kept Fisher", "Kept % Outliers"))
      
    }
  }
  
  if(Plot_Projection){
    for(i in Selected){
      
      to_plot <- data.frame("Genes" = RomaData$ModuleSummary[[i]]$OriginalGenes, 
                            "PC1" = RomaData$ModuleSummary[[i]]$PCBaseUnf$x[,1]*RomaData$ModuleSummary[[i]]$CorrectSignUnf, 
                            "PC2" = RomaData$ModuleSummary[[i]]$PCBaseUnf$x[,2]*RomaData$ModuleSummary[[i]]$CorrectSignUnf,
                            "Status" = rep("Kept_Base", length(RomaData$ModuleSummary[[i]]$OriginalGenes)))
      to_plot$Color <- "gray"
      
      to_plot[to_plot$Genes %in% RomaData$OriginalOutliers[[i]], c("Status")] <- "Kept_%_outliers"
      to_plot[to_plot$Genes %in% RomaData$OriginalOutliers[[i]], c("Color")] <- "purple"
      
      to_plot[to_plot$Genes %in% RomaData$OutLiersList[[i]], c("Status")]<- "Outlier"
      to_plot[to_plot$Genes %in% RomaData$OutLiersList[[i]], c("Color")]<- "red"
      
      to_plot[to_plot$Genes %in% intersect(RomaData$genes_to_keep_fisher, RomaData$OriginalOutliers[[i]]), c("Status")] <- "Kept_Fisher"
      to_plot[to_plot$Genes %in% intersect(RomaData$genes_to_keep_fisher, RomaData$OriginalOutliers[[i]]), c("Color")] <- "green"
      
      to_plot[to_plot$Genes %in% intersect(RomaData$genes_to_keep_counts, RomaData$OriginalOutliers[[i]]), c("Status")] <- "Kept_Counts"
      to_plot[to_plot$Genes %in% intersect(RomaData$genes_to_keep_counts, RomaData$OriginalOutliers[[i]]), c("Color")] <- "blue"
  
      
      XLims <- quantile(to_plot$PC1, c(.01, .99))
      XLims[1] <- min(XLims[1], min(to_plot$PC1))
      XLims[2] <- max(XLims[2], max(to_plot$PC1))
      
      YLims <- quantile(to_plot$PC2, c(.01, .99))
      YLims[1] <- min(YLims[1], min(to_plot$PC2))
      YLims[2] <- max(YLims[2], max(to_plot$PC2))
      
      to_plot_unique <- unique(to_plot[, c("Status", "Color")])
      
      p <- ggplot2::ggplot(to_plot, ggplot2::aes(x = PC1, y = PC2, colour = Status, alpha=Status, label = Status, colour = Status)) + ggplot2::geom_point() +
        ggplot2::scale_alpha_manual(values=c(Data=1, Sampling=.2), guide=FALSE) +
        ggplot2::scale_x_continuous(limits = XLims) +
        ggplot2::scale_y_continuous(limits = YLims) +
        ggplot2::scale_color_manual(breaks = to_plot_unique$Status, values = as.character(to_plot_unique$Color)) +
        ggplot2::ggtitle(RomaData$ModuleSummary[[i]]$ModuleName) +
        ggplot2::geom_hline(yintercept=0) + ggplot2::geom_vline(xintercept=0) +
        ggplot2::geom_text(ggplot2::aes(label=ifelse(Status %in% c("Kept_%_outliers", "Outlier", "Kept_Fisher", "Kept_Counts"),as.character(Genes),'')),hjust=0,vjust=0)
        
        print(p)

    }
  }

}





