#' Perform spatial transcriptome Copy Number Variation (CNV) Analysis
#'
#' This function conducts a comprehensive analysis of spatial transcriptome copy number variations (CNV) using Seurat and InferCNV.
#'
#' @param seurat_obj A Seurat object containing spatial transcriptome data.
#' @param test_adata A Seurat object containing test spatial transcriptome data.
#' @param output_dir The directory where the output files and results will be saved.
#' @param gene_order_file A .Rdata file containing gene order information.
#' @param metadata A vector specifying group variables to identify test and ref.
#' @param patient_id A unique identifier for the patient or sample being analyzed.
#' @param ref_group_names A character vector of reference group names for CNV estimation.
#' @param cutoff Cut-off for the min average read counts per gene among reference cells (default: 0.01).
#' @param cluster_by_groups If observations are defined according to groups (ie. patients), each group of cells will be clustered separately. (default=FALSE, instead will use k_obs_groups setting) (default: TRUE).
#' @param denoise A logical indicating whether to perform data denoising (default: TRUE).
#' @param HMM A logical indicating whether to use Hidden Markov Models (HMM) for CNV inference (default: FALSE).
#' @param num_threads The number of CPU threads to use for parallel processing (default: 15).
#' @param tumor_subcluster_partition_method The method used for partitioning tumor subclusters (default: 'qnorm').
#' @param analysis_mode The mode of analysis ('samples' or 'subclusters' or 'cells'; default: 'samples').
#' @param k_nn The number of nearest neighbors for clustering (default: 30).
#' @param leiden_resolution The resolution parameter for Leiden clustering (default: 1).
#' @param tumor_subcluster_pval The p-value threshold for tumor subcluster identification (default: 0.05).
#' @param iter_max The maximum number of iterations for k-means clustering (default: 100).
#' @param file_name file name to save the heatmap plot.
#'
#' @return A Seurat object with additional metadata columns (cnv_score, malignant/non-malignant).
#'
#' - CNV_Score.txt: CNV scores in spatial.
#' - kmeans_df_s.txt: Results of k-means clustering.
#' - spatial_CNV_score.png: Spatial plot of CNV scores.
#' - spatial_kmeans.png: Spatial plot of k-means clustering results.
#' - vln_cnv_score_kmeans.png: Violin plot of CNV scores by k-means clusters.
#' - spatial_malignant.pdf: Spatial plot of malignant cell identification.
#' - vln_cnv_score_malig.png: Violin plot of CNV scores between malignant/non-malignant.
#' @export
#'
#' @import Seurat
#' @import tidyverse
#' @import infercnv
#' @import ggsci
#' @import patchwork
STCNV <- function(seurat_obj, test_adata, output_dir, gene_order_file, metadata, patient_id, ref_group_names,
                  cutoff=0.01, cluster_by_groups=T, denoise=T, HMM=F, num_threads=15, tumor_subcluster_partition_method='qnorm', analysis_mode='samples',
                  k_nn=30, leiden_resolution=1, tumor_subcluster_pval=0.05, iter_max=100, file_name = "D8_rm_leiden3_cnv.pdf"){

  flog.info(sprintf("\n\n\tMODULE 1: CNV ESTIMATION\n"))
  # get raw count
  raw_count = GetAssayData(seurat_obj, slot = 'counts')

  anno <- FetchData(seurat_obj, vars = metadata)

  # --------------------------------------------- Step1: start run infercnv
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(raw_count),
                                      annotations_file=anno,
                                      delim="\t",
                                      gene_order_file=geneOrder,
                                      ref_group_names=ref_group_names)

  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=cutoff,  # use 1 for smart-seq, 0.1 for 10x-genomics, 0.01 for stereo-seq
                               out_dir=output_dir,  # dir is auto-created for storing outputs

                               #  多个样本混合的话，test是否按照病人进行分组计算
                               cluster_by_groups=cluster_by_groups,   # cluster

                               # 降噪
                               denoise=denoise,

                               #  tumor subclone
                               analysis_mode  = analysis_mode,

                               tumor_subcluster_partition_method = tumor_subcluster_partition_method,
                               k_nn = k_nn,
                               leiden_resolution = leiden_resolution,


                               # if tumor_subcluster_partition_method = "random_trees" or "qnorm"
                               tumor_subcluster_pval = tumor_subcluster_pval,

                               # infer region state
                               HMM=HMM,

                               # 并行计算
                               num_threads = num_threads,

                               useRaster = FALSE
  )

  flog.info(sprintf("\n\n\tMODULE 2: Malignant identification\n"))
  # --------------------------------------------- Step2: kmeans cluster
  # ------------- 1.CNV score
  infercnv_obj = readRDS(paste0('./', output_dir, "/run.final.infercnv_obj"))

  cnv_test = infercnv_obj@expr.data[, unlist(infercnv_obj@observation_grouped_cell_indices)]
  # cnv_test = read.table(paste0(opt$workdir, '/', opt$output, '/infercnv.observations.txt'), header = T, check.names = F)

  expr1=cnv_test-1
  expr2=expr1 ^ 2
  CNV_score=as.data.frame(colSums(expr2))
  colnames(CNV_score)="CNV_score"

  write.table(CNV_score, file = "CNV_Score.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)

  seurat_obj_test = subset(test_adata, ec_singler_epcam_celltype == 'Epithelial_cells')
  # seurat_obj_test = subset(seurat_obj, cells=colnames(seurat_obj)[seurat_obj@meta.data[, metadata] %in% c("Epi_19N", "Epi_34N")], invert = T)
  seurat_obj_test = RenameCells(seurat_obj_test, add.cell.id = patient_id)
  test_adata = RenameCells(test_adata, add.cell.id = patient_id)

  fil_seurat_obj_test = subset(seurat_obj_test, cells = rownames(CNV_score))
  fil_seurat_obj_test$cnv_score = CNV_score[colnames(fil_seurat_obj_test), 1]

  test_adata$cnv_score = NA
  test_adata$cnv_score[colnames(fil_seurat_obj_test)] = CNV_score[colnames(fil_seurat_obj_test), 1]

  f <- FeaturePlot(fil_seurat_obj_test, features = 'cnv_score', reduction = 'spatial')
  png('spatial_CNV_score.png')
  print(f)
  dev.off()

  # ------------- 2.kmeans聚类
  expr <- infercnv_obj@expr.data
  # dim(expr)
  # min(expr)
  # max(expr)

  set.seed(20210418)
  kmeans.result <- kmeans(t(expr), iter.max = iter_max, 2)

  # 构建dataframe
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  table(kmeans_df)

  kmeans_df$CB=rownames(kmeans_df)
  head(kmeans_df)

  # 取出正常spot对应的名字
  normal_loc <- unlist(infercnv_obj@reference_grouped_cell_indices)

  # 取出test spot对应的名字
  test_loc <- unlist(infercnv_obj@observation_grouped_cell_indices)

  anno.df=data.frame(
    CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
    class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
  )
  head(anno.df)

  kmeans_df=kmeans_df %>% inner_join(anno.df,by="CB") #合并
  kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
  rownames(kmeans_df_s)=kmeans_df_s$CB
  kmeans_df_s$CB=NULL
  kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class)

  head(kmeans_df_s)
  table(kmeans_df_s$kmeans_class, kmeans_df_s$class)

  write.table(kmeans_df_s, file = "kmeans_df_s.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)

  need_kmean_df_s = kmeans_df_s[colnames(fil_seurat_obj_test), ]
  head(need_kmean_df_s)

  all(rownames(need_kmean_df_s) == colnames(fil_seurat_obj_test))

  fil_seurat_obj_test$kmean_cluster = as.character(need_kmean_df_s$kmeans_class)

  my_cols = pal_npg()(2)

  v <- VlnPlot(fil_seurat_obj_test, group.by = 'kmean_cluster', features = 'cnv_score', pt.size = 0,
          cols = if(tapply(fil_seurat_obj_test$cnv_score, list(fil_seurat_obj_test$kmean_cluster), mean)[1]>tapply(fil_seurat_obj_test$cnv_score, list(fil_seurat_obj_test$kmean_cluster), mean)[2]){
            my_cols
            }else{
              rev(my_cols)
              })
  png("vln_cnv_score_kmeans.png")
  print(v)
  dev.off()

  d <- DimPlot(fil_seurat_obj_test, group.by = 'kmean_cluster', reduction = 'spatial', split.by = 'kmean_cluster', pt.size = 0.4, cols = if(tapply(fil_seurat_obj_test$cnv_score, list(fil_seurat_obj_test$kmean_cluster), mean)[1]>tapply(fil_seurat_obj_test$cnv_score, list(fil_seurat_obj_test$kmean_cluster), mean)[2]){
    my_cols
  }else{
    rev(my_cols)
  })
  png('spatial_kmeans.png', width = 700)
  print(d)
  dev.off()

  table(fil_seurat_obj_test$kmean_cluster, dnn = list(col.name = c('Malignant', 'Non-malignant')))

  # --------------------------------------------- Step3: Identification of malignant cells
  # add malignant info to meta.data
  test_adata$malignant = NA

  test_adata$malignant[colnames(test_adata) %in% rownames(kmeans_df_s)[kmeans_df_s$kmeans_class == names(tapply(fil_seurat_obj_test$cnv_score, list(fil_seurat_obj_test$kmean_cluster), mean))[which.max(tapply(fil_seurat_obj_test$cnv_score, list(fil_seurat_obj_test$kmean_cluster), mean))]]] <- 'Malignant cells'

  test_adata$malignant[colnames(test_adata) %in% rownames(kmeans_df_s)[kmeans_df_s$kmeans_class == names(tapply(fil_seurat_obj_test$cnv_score, list(fil_seurat_obj_test$kmean_cluster), mean))[which.min(tapply(fil_seurat_obj_test$cnv_score, list(fil_seurat_obj_test$kmean_cluster), mean))]]] <- 'Non-Malignant cells'

  table(test_adata$malignant, useNA = 'ifany')

  tmp = subset(test_adata, malignant %in% c('Malignant cells', 'Non-Malignant cells'))

  d1 <- DimPlot(tmp, group.by = 'malignant', reduction = 'spatial', cols = if(tapply(tmp$cnv_score, list(tmp$malignant), mean)[1]>tapply(tmp$cnv_score, list(tmp$malignant), mean)[2]){ my_cols }else{ rev(my_cols) }) + NoLegend()
  d2 <- DimPlot(tmp, group.by = 'malignant', reduction = 'spatial', split.by = 'malignant', cols = if(tapply(tmp$cnv_score, list(tmp$malignant), mean)[1]>tapply(tmp$cnv_score, list(tmp$malignant), mean)[2]){ my_cols }else{ rev(my_cols) })
  pdf('spatial_malignant.pdf', width = 18)
  d1 + d2 + plot_layout(widths = c(1, 2))
  dev.off()

  v <- VlnPlot(tmp, group.by = 'malignant', features = 'cnv_score', pt.size = 0, cols = if(tapply(tmp$cnv_score, list(tmp$malignant), mean)[1]>tapply(tmp$cnv_score, list(tmp$malignant), mean)[2]){my_cols}else{rev(my_cols)}) + NoLegend() + theme(axis.title = element_blank())
  png("vln_cnv_score_malig.png", width =  400)
  print(v)
  dev.off()

  saveRDS(seurat_obj_test, file = 'seurat_malignant.rds')

  CNV_heatmap(test_adata, output_dir = output_dir, file_name = file_name, patient_id = patient_id)

  return(test_adata)
}

#' Create a CNV Heatmap
#'
#' This function generates a CNV (Copy Number Variation) heatmap using input data.
#'
#' @param adata An input data object containing information for CNV analysis.
#' @param output_dir The directory where the output heatmap file will be saved.
#' @param file_name The name of the output heatmap file (default is "D8_rm_leiden3_cnv.pdf").
#' @param patient_id A unique identifier for the patient or sample being analyzed.
#'
#' @return This function doesn't return a value but saves the CNV heatmap as a PDF file.
#'
#' @export
#'
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import ggsci
#' @import patchwork
#' @import infercnv
#' @import tidyverse
#' @import Seurat
CNV_heatmap <- function(adata, output_dir, file_name = "D8_rm_leiden3_cnv.pdf", patient_id){

  malignant_cells = colnames(adata)[which(as.character(adata$malignant) == 'Malignant cells')]
  non_malignant_cells = colnames(adata)[which(as.character(adata$malignant) == 'Non-Malignant cells')]

  tumor_cell_order = c(malignant_cells, non_malignant_cells)

  infercnv_obj=readRDS(paste0("./", output_dir, "/run.final.infercnv_obj"))
  # plot_cnv(infercnv_obj,
  #          k_obs_groups=1,
  #          cluster_by_groups=T,
  #          cluster_references=T,
  #          out_dir="../results/Fig3C/",
  #          x.center=1,
  #          x.range="auto",
  #          title="D8_inferCNV",
  #          output_filename="D8_infercnv",
  #          output_format="png",
  #          write_expr_matrix=FALSE,
  #          png_res=300,
  #          useRaster=T)

  dt_tumor = infercnv_obj@expr.data[, unlist(infercnv_obj@observation_grouped_cell_indices)] %>% t()
  dt_normal = infercnv_obj@expr.data[, unlist(infercnv_obj@reference_grouped_cell_indices)] %>% t()

  dt_tumor = dt_tumor[tumor_cell_order, ]


  my_cols = pal_npg()(2)

  df_tumor = data.frame(Epi_P54_HE_ST = c(rep('Malignant cells', length(malignant_cells)), rep('Non-Malignant cells', length(non_malignant_cells))))
  row_ha_tumor <- rowAnnotation(df = df_tumor, col = list(Epi_P54_HE_ST = c('Malignant cells' = my_cols[1], 'Non-Malignant cells' = my_cols[2])), show_annotation_name = F)


  h_t <- Heatmap(as.matrix(dt_tumor),name="Expression",col =colorRamp2(c(0.9, 1, 1.1), c("darkblue", "white", "darkred")),
                 border_gp=gpar(col="black"), show_row_names = F,show_column_names = F,cluster_rows =F,cluster_columns = F,row_dend_width = unit(1.5, "cm"),
                 column_split = infercnv_obj@gene_order$chr, column_gap = unit(0, "mm"), column_title_side = 'bottom', column_title_rot = 90, column_title_gp = gpar(fontsize = 12),
                 row_title = NULL,
                 row_title_gp = gpar(fontsize = 18, fontface = 'bold'), use_raster = T,
                 left_annotation = row_ha_tumor)
  # column_labels =rep("          ",ncol(dt)),column_names_side = c("top"),

  df = data.frame(Ref = c(rep('Epi_P19N', length(infercnv_obj@reference_grouped_cell_indices$Epi_19N)), rep('Epi_P34N', length(infercnv_obj@reference_grouped_cell_indices$Epi_34N))))

  row_ha = rowAnnotation(df = df, col = list(Ref = c('Epi_P19N' = '#3B4992FF', 'Epi_P34N' = "#008280FF")), show_annotation_name = F)


  h_n <- Heatmap(as.matrix(dt_normal),name="Expression1",
                 col =colorRamp2(c(0.9, 1, 1.1), c("darkblue", "white", "darkred")),
                 border_gp=gpar(col="black"),
                 show_row_names = F,show_column_names = F,cluster_rows =F,cluster_columns = F, row_dend_width = unit(0, "cm"),
                 column_split = infercnv_obj@gene_order$chr, column_gap = unit(0, "mm"), column_title_side = 'bottom', column_title_rot = 90, column_title_gp = gpar(fontsize = 12),
                 # row_title = 'References',
                 row_title = NULL,
                 row_title_gp = gpar(fontsize = 18, fontface = 'bold'), column_title = NULL, show_heatmap_legend = F,
                 left_annotation = row_ha, use_raster = T,
                 row_split = df, row_gap = unit(0, "mm"), height = unit(3, "cm"))

  ht_list  <- h_n %v% h_t

  pdf(file_name, width = 13)
  draw(ht_list, main_heatmap = 'Expression', padding = unit(c(5, 2, 2, 2), "mm"))
  dev.off()
}




