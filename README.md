# StereoPath
StereoPath is a toolkit that provides functions to repeat traditional pathology of CRC using stereo-seq data. It enables the identification of malignant cells, analysis of vasculature structure, estimation of proliferative proportion, and classification of dMMR/pMMR phenotype.

## Installation

StereoPath R package can be easily installed from Github using devtools:  

```
devtools::install_github("STOmics/StereoPath")
```

or

```
git clone https://github.com/STOmics/StereoPath.git

# In R
install.packages('StereoPath', repo=NULL)
```

## Usage

```
# celltype annotation
adata_anno <- STanno(adata, cutoff = 0.1)

# Identification of malignant cells
adata_cnv <- STCNV(seurat_obj, test_adata = adata, output_dir = 'output_dir', gene_order_file = geneOrder, metadata = 'orig', patient_id = 'Epi_D8', ref_group_names = c("Epi_19N", "Epi_34N"), cutoff=0.01, cluster_by_groups = F, denoise = T, HMM=F, num_threads = 20, iter_max = 100)

# classfication of dMMR/pMMR
STpheno(obj, sample_id = names(obj), phenotype_id = c(rep('dMMR', 3), rep('pMMR', 9)))
```
