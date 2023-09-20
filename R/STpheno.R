#' Calculate Score, P-Value, and Call for a Single Gene
#'
#' This function calculates a score, p-value, and call for a single gene based on its expression value.
#'
#' @param x The expression value of the gene.
#' @param m The mean expression value of the gene.
#' @param sd The standard deviation of the gene's expression values.
#' @param p.thresh The p-value threshold for the call (default: 0.01).
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{score}{The calculated score for the gene, representing its deviation from the mean, scaled by the standard deviation implied by the interquartile range (IQR).}
#'   \item{pval}{The p-value associated with the gene's score.}
#'   \item{call}{A binary call (0 or 1) indicating whether the gene's p-value is below the specified threshold.}
#' }
#'
#' @export
mmr_loss_score = function(x,m,sd,p.thresh=0.01)
{
  ### calc linear score: diff from center, scaled by sd implied by iqr:
  score = (x-m)/sd
  ### calc pval:
  pval = pnorm(score)
  # call:
  call = 1*(pval < p.thresh)

  # save results:
  names(score) = names(call) = names(pval) = names(x)
  out = list(score = score, pval=pval,call=call)
  return(out)
}

#' Perform MMR Loss Score Analysis for Multiple Seurat Objects
#'
#' This function performs MMR (Mismatch Repair) Loss Score analysis for multiple Seurat objects, specifically related to colorectal cancer (CRC) phenotypes.
#'
#' @param seurat_obj_list A list of Seurat objects, each representing spatial transcriptome data for different samples.
#' @param sample_id A character vector of sample IDs corresponding to the Seurat objects.
#' @param phenotype_id A character vector of phenotype labels for each sample.
#'
#' @return This function generates various plots and does not have a formal return value.
#'
#' @export
#'
#' @import Seurat
#' @import mclust
#' @import ggplot2
#' @import ggpubr
#' @import ggrepel
#' @import tidyverse
#' @import patchwork
#' @import reshape2
#' @import ggsci
STpheno <- function(seurat_obj_list, sample_id, phenotype_id){
  # ----------------------------- filter genes ---------------------------------------
  gene_name = c(rownames(seurat_obj_list[[1]]))
  for (i in seurat_obj_list) {
    gene_name = intersect(gene_name, rownames(i))
  }
  # ------------------------------------------------------

  # ------------------------------- create mean expression profile in malignant cells ----------------------
  mean_df = NULL
  for (i in 1:length(seurat_obj_list)){
    tmp = GetAssayData(seurat_obj_list[[i]], slot = 'data')

    print(paste0(names(seurat_obj_list)[i], ': ', length(colnames(subset(seurat_obj_list[[i]], malignant == 'Malignant cells')))))
    tmp = tmp[gene_name, colnames(subset(seurat_obj_list[[i]], malignant == 'Malignant cells')) ]

    tmp = rowMeans(tmp)

    mean_df = cbind(mean_df, tmp)
  }
  dim(mean_df)

  colnames(mean_df) = sample_id

  meta = data.frame(sample = sample_id, phenotype = phenotype_id, rownames = sample_id)
  #######################################################

  #######-------------------------- MMR loss score -------------------------------------
  key4genes = c("MLH1","MSH2","MSH6","PMS2") #,"EPCAM")

  normalized_counts = t(mean_df)
  centers = sds = matrix(NA,1,4,dimnames = list(names(normalized_counts),key4genes))

  ######################################################

  for(gene in key4genes) {
    # for the genes where it mistakenly gets 2 clusters, force G to 1:
    G=1:2

    # model-based clustering to get classification:
    set.seed(0)
    mc = Mclust(normalized_counts[, gene],G=G)
    print(mc$classification)
    # identify upper class:
    top = which(mc$parameters$mean==max(mc$parameters$mean))
    # save mean of the upper cluster
    centers[1,gene]=mc$parameters$mean[top]

    # unbiased procedure: get SD in MSS samples:
    tempmsi = as.character(meta$phenotype)
    sds[1,gene] = sd(normalized_counts[(tempmsi=="pMMR")|(is.na(tempmsi)), gene] )
  }

  names(mc$classification) <- sample_id


  scores = list()
  scores[['CRC']] = normalized_counts[,key4genes]*NA

  for (gene in key4genes) {
    temp = mmr_loss_score(x = normalized_counts[, gene], m=centers[1,gene],sd=sds[1,gene],p.thresh=0.005)
    scores[['CRC']][,gene] = temp$score
  }

  tempscore = c(apply(scores[['CRC']],1,min))
  # sort(tempscore)

  tempdf = data.frame(score = tempscore, phenotype = phenotype_id)

  g = ggplot(tempdf, aes(x = score, y = 0.1)) + geom_point(aes(color = phenotype, size = 10)) + geom_label_repel(aes(score, 0.1, label = rownames(tempdf)), segment.colour="black", min.segment.length = 0) + guides(size = 'none', color = guide_legend(title = 'Phenotype')) +
    theme_bw() + theme(text = element_text(size = 13)) + ylab('') + scale_y_discrete(labels = '') + xlab('MMR Loss score') + ggtitle('Stereo_CRC') + scale_color_manual(values = c("#f391bd","#8bc0ee")) + theme(panel.grid = element_blank())
  pdf('scatter_dMMR_pMMR_classification.pdf')
  print(g)
  dev.off()

  # --------------------- boxplot
  tempdf$patient = rownames(tempdf)
  test.method = "wilcox" # t.test, wilcox
  g <- ggplot(tempdf, aes(x = phenotype, y = score)) + geom_boxplot(aes(fill = phenotype)) +
    geom_point(aes(x = phenotype, y = score, color = patient, size = 4)) +
    theme(legend.background = element_blank(), legend.key = element_blank(), text = element_text(size = 25), axis.text = element_text(size = 25, colour = 'black'), panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black", size = 0.5)) + ylab('MMR Loss Score') + xlab('') + stat_compare_means(method = test.method, size = 8) + scale_color_brewer(palette = 'Paired') +
    guides( color=guide_legend(title="Patient", override.aes = list(size=3)), size = F, fill = guide_legend(title = 'Group')) + scale_fill_manual(values = c("#f391bd","#8bc0ee"))
  # png("boxplot_mmr_loss_score.png")
  # print(g)
  # dev.off()

  ggsave(filename = 'boxplot_MMR_loss_score.pdf', g, dpi = 300)
}
