#' Compare ROMA scores across samples from different populations
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Groups string vector, a vector of group identifiers. Must contain the names of the samples.
#' @param Selected Genesets used to perform the analysis 
#' @param TestMode string, the type of statistical methodology to assess sample difference. Currently only ANOVA ("Aov") is available
#' @param PlotDiag boolean, should diagnostic plot be displayed?
#' @param Plot string, whether difference should be plotted by group or by sample. Possible values are "group", "sample", "both"
#'
#' @return
#' @export
#'
#' @examples
GlobalCompareAcrossSamples <- function(RomaData, Groups, Selected = NULL,
                                 TestMode = "Aov", PlotDiag = FALSE, Plot = "both") {
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Geneset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset(s) selected"))
  }
  
  tMat <- RomaData$SampleMatrix[Selected,]
#  names(Groups) <- colnames(tMat)  Groups doit déjà contenir les noms des samples
  
  MeltData <- reshape::melt(tMat)
  MeltData <- cbind(MeltData, Groups[as.character(MeltData$X2)])
  
  
  # MeltData <- data.frame(MeltData)
  
  colnames(MeltData) <- c("GeneSet", "Sample", "Value", "Group") 
  
  MeltData <- MeltData[!is.na(MeltData$Group), ]
  
  if(TestMode == "Aov"){
    
    print("Performing Type III AOV (R default)")
    
    AOVFitTypeI <- aov(formula = Value ~ Group, data = MeltData)
    
    if(PlotDiag){
      plot(AOVFitTypeI)
    }
    
    SumData <- summary(AOVFitTypeI)
    print(SumData)
    
    if(Plot %in% c("both", "group")){
      p <- ggplot2::ggplot(MeltData, ggplot2::aes(y=Value, x=Group, fill=Group)) +
        ggplot2::geom_boxplot() + ggplot2::guides(fill = "none") +
        ggsignif::geom_signif(comparisons = GetComb(unique(MeltData$Group)),
                              map_signif_level=TRUE, test = "wilcox.test", step_increase = .1) +
        ggplot2::labs(y="Sample score", x="Groups", title = "Groups")
    
      print(p)
    }
    
    if(Plot %in% c("both", "sample")){
      p <- ggplot2::ggplot(MeltData, ggplot2::aes(y=Value, x=Sample, fill=Group)) +
        ggplot2::geom_boxplot() + ggplot2::coord_flip() +
        ggplot2::labs(y="Sample score", x="Samples", title = "Groups") +
        ggplot2::theme(axis.text.y = ggplot2::element_blank())
        
      print(p)
    }
    
  }

}


#' Compare ROMA sample scores across samples from different populations, for each selected geneset
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Groups string vector, a vector of group identifiers. Must contain the names of the samples.
#' @param Selected Genesets used to perform the analysis 
#' @param TestMode string, the type of statistical methodology to assess sample difference. Currently only ANOVA ("Aov") is available
#' @param PlotDiag boolean, should diagnostic plot be displayed?
#'
#' @return
#' @export
#'
#' @examples
GlobalCompareAcrossSamplesGenesets <- function(RomaData, Groups, Selected = NULL,
                                       TestMode = "Aov", PlotDiag = FALSE) {
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Geneset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset(s) selected"))
  }
  
  tMat <- RomaData$SampleMatrix[Selected,]
  #  names(Groups) <- colnames(tMat)  Groups doit déjà contenir les noms des samples
  
  MeltData <- reshape::melt(tMat)
  MeltData <- cbind(MeltData, Groups[as.character(MeltData$X2)])
  
  
  # MeltData <- data.frame(MeltData)
  
  colnames(MeltData) <- c("GeneSet", "Sample", "Value", "Group") 
  
  MeltData <- MeltData[!is.na(MeltData$Group), ]
  
  if(TestMode == "Aov"){
    
    print("Performing Type III AOV (R default) for each geneset")
    
    AOVFitTypeI <- aov(formula = Value ~ Group/GeneSet, data = MeltData)
    
    if(PlotDiag){
      plot(AOVFitTypeI)
    }
    
    SumData <- summary(AOVFitTypeI)
    print(SumData)
    
    Sep <- seq(from = 1, to = length(unique(MeltData$GeneSet)), by = 4)
    if(max(Sep) < length(unique(MeltData$GeneSet))+1) Sep <- c(Sep, length(unique(MeltData$GeneSet))+1)
    
    for(i in 2:length(Sep)){
      
      p <- ggplot2::ggplot(MeltData[as.integer(MeltData$GeneSet) %in% Sep[i-1]:(Sep[i]-1),],
                           ggplot2::aes(y=Value, x=Group, fill=Group)) + ggplot2::geom_boxplot() +
        ggsignif::geom_signif(comparisons = GetComb(unique(MeltData$Group)),
                              map_signif_level=TRUE, test = "wilcox.test", step_increase = .1) +
        ggplot2::labs(y="Sample score", x="Groups", title = paste("Geneset VS Groups - Part", i-1)) +
        ggplot2::facet_wrap( ~ GeneSet, scales = "free_y", ncol = 2) + ggplot2::theme(strip.text.x = ggplot2::element_text(size=6, face = "bold")) +
        ggplot2::guides(fill = "none")
      
      print(p)
    }
    
  }
  
}










#' Select genesets accoding to specific conditions
#'
#' @param RomaData list, the analysis returned by rRoma 
#' @param VarThr numeric between 0 and 1, the threshold PV to select significantly over- or under-dispersed genesets
#' @param VarMode string, the test to use to select over- or under-dispersed genesets.
#' Currently it can be either 'Wil' (Wilcoxon test) or 'PPV' (permutation base p-value)
#' @param VarType string, the type of statistical difference to select. Currently it can be either 'Over' (overdispersed) or 'Under' (underdispersed)
#' @param MedThr numeric between 0 and 1, the threshold PV to select significantly over- or under-expressed genesets, using median expression
#' @param MedMode string, the test to use to select over- or under-expressed genesets, using median expression.
#' Currently it can be either 'Wil' (Wilcoxon test) or 'PPV' (permutation base p-value)
#' @param MedType string, the type of statistical difference to select. Currently it can be either 'Over' (overexpressed) or 'Under' (underexpressed)
#' @param PValAdjust string, the p-value adjustment methods. Can be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none"
#'
#' @return
#' @export
#'
#' @examples
SelectGeneSets <- function(RomaData,
                           VarThr = 1e-3, VarMode = "Wil", VarType = "Over",
                           RatThr = NULL, RatMode = "Wil", RatType = "Over",
                           MedThr = NULL, MedMode = "Wil", MedType = "Over",
                           PValAdjust = "none") {
  
  
  if(VarType %in% c("Over", "Under") & !is.null(VarThr) & VarMode %in% c("Wil", "PPV")){ 
    
    if(VarMode == 'Wil' & VarType == "Over"){
      print(paste("Using genesets overdispersed according to Wilcoxon test. VarThr =", VarThr))
      SelectedVar <- which(p.adjust(RomaData$PVVectMat[,1], method = PValAdjust)<VarThr)
    }
    
    if(VarMode == 'Wil' & VarType == "Under"){
      print(paste("Using genesets underdispersed according to Wilcoxon test. VarThr =", VarThr))
      SelectedVar <- which(p.adjust(RomaData$PVVectMat[,2], method = PValAdjust)<VarThr)
    }
    
    if(VarMode == 'PPV' & VarType == "Over"){
      print(paste("Using genesets overdispersed according to pseudo pv. VarThr =", VarThr))
      SelectedVar <- which(p.adjust(RomaData$ModuleMatrix[,3], method = PValAdjust)<VarThr)
    }
    
    if(VarMode == 'PPV' & VarType == "Under"){
      print(paste("Using genesets underdispersed according to pseudo pv. VarThr =", VarThr))
      SelectedVar <- which(p.adjust(1-RomaData$ModuleMatrix[,3], method = PValAdjust)<VarThr)
    }
    
  } else {
    print("No dispersion filter selected")
    SelectedVar <- 1:nrow(RomaData$ModuleMatrix)
  }
  
  
  
  
  if(!is.null(RatThr) & RatType %in% c("Over", "Under") & RatMode %in% c("Wil", "PPV")){
    
    if(RatMode == 'Wil' & RatType == "Over"){
      print(paste("Using genesets overcoordinated according to Wilcoxon test. RatThr =", RatThr))
      SelectedRat <- which(p.adjust(RomaData$PVVectMat[,3], method = PValAdjust)<RatThr)
    }
    
    if(RatMode == 'Wil' & RatType == "Under"){
      print(paste("Using genesets undercoordinated according to Wilcoxon test. RatThr =", RatThr))
      SelectedRat <- which(p.adjust(RomaData$PVVectMat[,4], method = PValAdjust)<RatThr)
    }
    
    if(RatMode == 'PPV' & RatType == "Over"){
      print(paste("Using genesets overcoordinated according to pseudo pv. RatThr =", RatThr))
      SelectedRat <- which(p.adjust(RomaData$ModuleMatrix[,6], method = PValAdjust)<RatThr)
    }
    
    if(RatMode == 'PPV' & RatType == "Under"){
      print(paste("Using genesets undercoordinated according to pseudo pv. RatThr =", RatThr))
      SelectedRat <- which(p.adjust(1-RomaData$ModuleMatrix[,6], method = PValAdjust)<RatThr)
    }

  } else {
    print("No coordination filter selected")
    SelectedRat <- 1:nrow(RomaData$ModuleMatrix)
  }
  
  
  
  if(!is.null(MedThr) & MedType %in% c("Over", "Under") & VarMode %in% c("Wil", "PPV")){
    
    if(MedMode == 'Wil' & MedType == "Over"){
      print(paste("Using genesets overexpressed according to Wilcoxon test. MedThr =", MedThr))
      SelectedMed <- which(p.adjust(RomaData$PVVectMat[,5], method = PValAdjust)<MedThr)
    }
    
    if(MedMode == 'Wil' & MedType == "Under"){
      print(paste("Using genesets underexpressed according to Wilcoxon test. MedThr =", MedThr))
      SelectedMed <- which(p.adjust(RomaData$PVVectMat[,6], method = PValAdjust)<MedThr)
    }
    
    if(MedMode == 'PPV' & MedType == "Over"){
      print(paste("Using genesets overexpressed according to pseudo pv. MedThr =", MedThr))
      SelectedMed <- which(p.adjust(RomaData$ModuleMatrix[,8], method = PValAdjust)<MedThr)
    }
    
    if(MedMode == 'PPV' & MedType == "Under"){
      print(paste("Using genesets underexpressed according to pseudo pv. MedThr =", MedThr))
      SelectedMed <- which(p.adjust(1-RomaData$ModuleMatrix[,8], method = PValAdjust)<MedThr)
    }
    
  } else {
    print("No expression filter selected")
    SelectedMed <- 1:nrow(RomaData$ModuleMatrix)
  }
  
  Select <- intersect(intersect(SelectedMed, SelectedRat), SelectedVar)
  
  if(length(Select) > 0){
    return(Select)
  } else {
    return(NULL)
  }

}
