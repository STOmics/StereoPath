#' Perform celltype annotation
#'
#' This function performs celltype annotation, including vascular signature scoring, cell type assignment, and spatial visualization.
#'
#' @param seurat_obj A Seurat object containing spatial transcriptome data.
#' @param cutoff A cutoff value of AUC Score to determine endothelial cells (default: 0.1).
#'
#' @return A modified Seurat object with additional cell type annotations.
#'
#' @export
#'
#' @import Seurat
#' @import SingleR
#' @import AUCell
#' @import ggplot2
#' @import tidyverse
#' @import futile.logger
STanno <- function(seurat_obj, cutoff = 0.1){
  flog.info(sprintf("\n\n\tSTEP 1: Endothelial annotation\n"))
  flog.info("Normalization")
  adata = seurat_obj
  adata = NormalizeData(adata, verbose = F)

  ####################################################################################################
  # ---------------- vascular signature -----------------
  flog.info("Vascular signature")
  genes = c('KIAA1671', 'ADGRF5', 'ADGRL4', 'AQP1', 'ARHGAP29', 'CD300LG', 'CDH13', 'CDH5', 'CLDN5', 'CLU', 'COL4A1', 'COL4A2', 'CRIM1', 'CYYR1', 'ECE1', 'ECSCR', 'EGFL7', 'EMCN', 'EPAS1', 'ESAM', 'ESM1', 'FABP4', 'FLT1', 'FLT4', 'GPIHBP1', 'IGFBP3', 'INHBB', 'KDR', 'LDB2', 'LRG1', 'LTBP4', 'MGLL', 'MMRN2', 'MYCT1', 'PCDH17', 'PECAM1', 'PLXNA2', 'PODXL', 'PREX2', 'PRSS23', 'PTPRB', 'RAMP2', 'RAPGEF5', 'RASIP1', 'S1PR1', 'SORBS2', 'STC1', 'TIE1', 'TP53I11')
  # genes <- nichenetr::convert_mouse_to_human_symbols(symbols = genes)

  genes_o <- intersect(genes, rownames(adata))

  geneSets <- list(vacular=genes_o)
  exprMatrix = GetAssayData(adata, slot = 'data')

  flog.info("Cal AUC score")
  cells_AUC <- AUCell_run(exprMatrix, geneSets)
  save(cells_AUC, file = 'cells_AUC.Rdata')


  set.seed(123)
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=4, assign=TRUE)
  save(cells_assignment, file = "cells_assignment.Rdata")

  # -------------------------- plot auc score -------------------------
  flog.info("Spatial visualization of AUC score")
  auc_matrix <- AUCell::getAUC(cells_AUC)

  adata$auc = auc_matrix[1, ]

  f <- FeaturePlot(adata, features = 'auc', reduction = 'spatial', order = T, raster = F) + ggtitle('EC AUC Score')
  print(f)
  pdf('EC_AUC_score.pdf')
  print(f)
  dev.off()

  flog.info("Binary AUC score")
  if (is.null(cutoff)) {
    on_cells <- cells_assignment$vacular$assignment
  }else{
    on_cells <- colnames(adata)[adata$auc > cutoff]
  }

  save(on_cells, file = 'on_cells.Rdata')

  flog.info("Spatial visualization of Endothelial")
  adata$EC = NA
  adata$EC[colnames(adata) %in% on_cells] = 'EC'
  # table(adata$EC, useNA='always')

  d <- DimPlot(adata, group.by = 'EC', cols = "firebrick", na.value = 'lightgray', reduction = 'spatial', order = T, raster = F)
  print(d)
  pdf("EC.pdf", height = 7, width = 7)
  print(d)
  dev.off()

  flog.info(sprintf("\n\n\tSTEP 2: Epithelial annotation\n"))
  flog.info("Start round1 anno to identify epithelial")
  adata_rmec = subset(adata, EC == 'EC', invert = T)

  hpca.comb = hpca.comb[, hpca.comb$label.main != 'Endothelial_cells']


  # first anno
  Matrix_for_SingleR = GetAssayData(adata_rmec, slot = 'data')

  pred <- SingleR(test = Matrix_for_SingleR,
                  ref = list(Blue.comb, hpca.comb),
                  labels = list(Blue.comb$label.main, hpca.comb$label.main))

  save(pred, file = 'round1_singler.Rdata')

  # load("round1_singler.Rdata")

  adata_rmec$ec_singler_epcam_celltype = NA
  adata_rmec$ec_singler_epcam_celltype[colnames(adata_rmec) %in% rownames(pred)] = pred$labels

  flog.info(sprintf("\n\n\tSTEP 3: other celltype annotation\n"))
  flog.info("Start round2 anno to identify other celltypes")
  # round2
  hpca.comb = hpca.comb[, hpca.comb$label.main != 'Epithelial_cells']

  epi_spot = WhichCells(adata_rmec, expression = EPCAM > 0 & ec_singler_epcam_celltype == 'Epithelial_cells')
  length(epi_spot)
  save(epi_spot, file = 'epi_spot.Rdata')
  # load('./epi_spot.Rdata')

  cols = c('#1f77b4','#ff7f0e','#279e68','#d62728','#aa40fc','#8c564b','#e377c2','#b5bd61','#17becf','#aec7e8','#ffbb78','#98df8a')

  adata_rmec_epi = subset(adata_rmec, cells = epi_spot, invert = T)

  Matrix_for_SingleR = GetAssayData(adata_rmec_epi, slot = 'data')
  pred <- SingleR(test = Matrix_for_SingleR,
                  ref = list(Blue.comb, hpca.comb),
                  labels = list(Blue.comb$label.main, hpca.comb$label.main))

  save(pred, file = 'round2_singler.Rdata')

  # load("./round2_singler.Rdata")

  adata$ec_singler_epcam_celltype = NA
  adata$ec_singler_epcam_celltype[colnames(adata) %in% on_cells] <- 'Endothelial_cells'
  adata$ec_singler_epcam_celltype[colnames(adata) %in% epi_spot] <- 'Epithelial_cells'
  adata$ec_singler_epcam_celltype[colnames(adata) %in% rownames(pred)] <- pred$labels
  table(adata$ec_singler_epcam_celltype, useNA='ifany')
  # table(adata$ec_singler_epcam_celltype, adata$singler_epcam_celltype, useNA = 'ifany')

  flog.info("Spatial visualization of celltype")
  d <- DimPlot(adata, group.by = 'ec_singler_epcam_celltype', split.by = 'ec_singler_epcam_celltype', reduction = 'spatial', ncol = 4, cols = cols, raster = F)
  print(d)
  png("sep_celltype_split.png", width = 1500, height = 1500)
  print(d)
  dev.off()

  d <- DimPlot(adata, reduction = 'spatial', group.by = 'ec_singler_epcam_celltype', cols = cols, pt.size = 0.8, raster = F)
  pdf('spatial_celltype.pdf')
  print(d)
  dev.off()

  adata$show = NA
  adata$show[adata$ec_singler_epcam_celltype == 'Epithelial_cells'] <- 'Epithelial_cells'

  dd <- DimPlot(adata, group.by = 'show', reduction = 'spatial', raster = F, order = T, cols = '#8c564b', na.value = 'lightgray')
  ggsave(filename = 'spatial_epi.pdf', dd, dpi = 300, height = 7, width = 7)

  write.csv(adata@meta.data, file = "meta_celltype.csv")

  return(adata)
}









