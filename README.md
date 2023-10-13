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
library(StereoPath)
library(dplyr)
library(grid)

# celltype annotation
adata_anno <- STanno(adata, cutoff = 0.1)

# Identification of malignant cells
adata_cnv <- STCNV(seurat_obj, test_adata = adata, output_dir = 'output_dir', gene_order_file = geneOrder, metadata = 'orig', patient_id = 'Epi_D8', ref_group_names = c("Epi_19N", "Epi_34N"), cutoff=0.01, cluster_by_groups = F, denoise = T, HMM=F, num_threads = 20, iter_max = 100)

# classfication of dMMR/pMMR
STpheno(obj, sample_id = names(obj), phenotype_id = c(rep('dMMR', 3), rep('pMMR', 9)))
```


```
# Estimation of proliferation proportion in python
# import function
# wget https://github.com/STOmics/StereoPath/archive/refs/heads/main.zip and unzip main.zip
import sys  
sys.path.append("/StereoPath-main/")  
from stereotyping import *

# read orig IF registered image
IF = read_IF_img("/jdfssz2/ST_BIOINTEL/P22Z10200N0618/huanghuaqiang/st_crc/pathology/34_IF_ST/IF/KI67_IF_SS200000979TL_E6_regist_32bit.tif")

# pooling IF images to be consistent with bin50 stereo sequence data
# you can use different pooling func, such as sum, max, min, mean, q75 etc..
q75_wide = bin_IF_img(IF, func=percentile(75), file_fix='q75', bin_size=50)

# you can compare transform image with original image via kl divergence
kl_div_plot(I=IF, J=q75_wide, func='q75')

# read Stereo-seq data
adata_34_if_st = sc.read("/jdfssz2/ST_BIOINTEL/P22Z10200N0618/huanghuaqiang/st_crc/stereotyping/adata_34_if_st.h5ad")
adata_34_if_st.obs['spot'] = adata_34_if_st.obs.index

# To predict proliferative proportion, you should provide ground truth label. You can use ImageJ color threshold to binary IF image signal to get label.
# read ImageJ mean color threshold binary image
IF_binary = read_IF_img("/jdfssz2/ST_BIOINTEL/P22Z10200N0618/huanghuaqiang/st_crc/pathology/34_IF_ST/IF/cutout_if_bin50_percentile_mean_cutoff_binary.tif")
if_q75_mean_binary_df = array2df(IF_binary, index = list(range(adata_34_if_st.obsm['spatial'][:,1].min(), adata_34_if_st.obsm['spatial'][:,1].max()+1)), columns = list(range(adata_34_if_st.obsm['spatial'][:,0].min(), adata_34_if_st.obsm['spatial'][:,0].max()+1)), func = 'q75_cutout_mean')

adata_34_if_st.obs = adata_34_if_st.obs.merge(if_q75_mean_binary_df, how='left', on='spot').set_index(adata_34_if_st.obs.index)

adata_34_if_st.obs['mean_ki67_positive'] = np.where(adata_34_if_st.obs['value'] == 255, 'Proliferative cells', 'Non-Proliferative cells')

# Now, you have a stereo-seq data matrix and label for each bin. You should do above process workflow to data you want to use for training.


# ---------------
# read the merged data you used for training
adata = sc.read("/jdfssz2/ST_BIOINTEL/P22Z10200N0618/huanghuaqiang/st_crc/stereotyping/adata_total.h5ad")

# to estimation of proliferative proportion, you should extract epithelial cells
adata_malig = adata[adata.obs.malignant.isin(['Malignant cells'])]

# and you should provide genes used for training. I only keep genes that overlap in all data, including training dataset, validation dataset, test dataset
used_all_genes = pd.read_csv("/jdfssz2/ST_BIOINTEL/P22Z10200N0618/huanghuaqiang/st_crc/pathology/PCM/used_all_genes.csv", index_col = 0)

# Now, get the matrix for training
data = sc.get.obs_df(adata_malig, keys=used_all_genes['0'].tolist())
# get the ground truth label
label = sc.get.obs_df(adata_malig, keys=['mean_ki67_positive'])
label['mean_ki67_positive'] = np.where(label['mean_ki67_positive']=='Proliferative cells', 1, 0)

# split data to train and test data
X_train, X_test, y_train, y_test = train_test_split(data, label, test_size=.3, random_state=1234, stratify=label)

# do the feature selection to remove redundancy gene and improve the model performance
Feature_Selector = feature_select(X_train, y_train, tree_method = 'gpu_hist', gpu_id = 1, n_jobs=15, objective = 'binary:logistic', pvalue=0.05,
                   random_state=132, n_trials=50, percentile=100, plot_res = True, plot_filename = "borutashap_feature_importance.pdf")

# read features we select
borutashap_genes = pd.read_csv("/jdfssz2/ST_BIOINTEL/P22Z10200N0618/huanghuaqiang/st_crc/stereopath/borutashap_genes.csv", index_col=0)

# Now, start to train xgboost model
data = sc.get.obs_df(adata_malig, keys=borutashap_genes['0'].tolist())
label = sc.get.obs_df(adata_malig, keys=['mean_ki67_positive'])
label['mean_ki67_positive'] = np.where(label['mean_ki67_positive']=='Proliferative cells', 1, 0)
X_train, X_test, y_train, y_test = train_test_split(data, label, test_size=.3, random_state=1234, stratify=label)

xgboost_bo = model_fine_tune(X_train, y_train, X_test, y_test, gpu_id=3, n_splits=5, scoring='roc_auc', n_jobs=15, save_model='clf_fine_tune.json',
                    title='borutashap_genes')

# after training, we can use this model to predict proliferative cells and calculate proliferative proportion
# read anndata that you want to predict
adata_29 = sc.read("/jdfssz2/ST_BIOINTEL/P22Z10200N0618/huanghuaqiang/st_crc/stereotyping/adata_29.h5ad")

adata_29 = model_predict(adata_29, xgboost_bo, gene_list=xgboost_bo.get_booster().feature_names, pt_size = 3, save_col = 'res', sub_col = 'malignant', sub_group='Malignant cells')

# finally, we can explore the feature importance in model, which can help us to explain the model
shap_genes_plot(xgboost_bo, X_test, feature_perturbation = "tree_path_dependent", approximate = True, filename = "shap_genes.pdf")
```

