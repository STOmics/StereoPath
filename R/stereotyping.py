# -*- coding: utf-8 -*-
# @Author: huanghuaqiang
# @Time: 2023-08-11
# @Version = "1"

import pandas as pd
import numpy as np
import seaborn as sns
import os, sys, glob, re
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sp
import pickle
import cv2
import tifffile as tifi
from PIL import Image
import imagecodecs
import scanpy as sc
from xgboost import XGBClassifier
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn import metrics
import shap
from BorutaShap import BorutaShap
import logging
from hyperopt import tpe, STATUS_OK, Trials, hp, fmin, STATUS_OK, space_eval
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn import metrics
import anndata

format="%(asctime)s-%(levelname)s-%(message)s"
logging.basicConfig(format=format, level=logging.INFO)

def read_IF_img(img_path:str=None) -> np.ndarray:
    """
    Read an IF image from the specified path.

    Parameters
    ----------
    img_path : str
        The path of the IF image.

    Returns
    -------
    numpy.ndarray
        Numpy array representing the IF image.
    """
    IF = tifi.imread(img_path)
    return IF

def percentile(n:int=75) -> int:
    """
    Calculate the n-quantile.

    Parameters
    ----------
    n : int
        Number of quantile.

    Returns
    -------
    function
        A function that calculates the specified quantile 'n' of an input array.
    """
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_

def bin_IF_img(IF:np.ndarray, func, file_fix:str=None, bin_size:int=50) -> np.ndarray:
    """
    Pool IF signal by window, perform aggregation, and save results.

    Parameters
    ----------
    IF : numpy.ndarray
        Numpy array representing the IF image.
    func
        The function to aggregate IF signal.
    file_fix : str, optional
        Differentiator for output file name, by default None.
    bin_size : int, optional
        The size of the n*n window to aggregate the signal, by default 50.

    Returns
    -------
    np.ndarray
        Numpy array of bin-transformed IF image.

    Save
    ----
    .csv
        IF image data in long format with bin coordinates.
    .tif
        Bin-transformed IF image.
    """
    IF_df = pd.DataFrame(IF)
    IF_df.index.name = 'y'
    IF_df.reset_index(inplace=True)
    IF_df_long = IF_df.melt(id_vars='y')
    IF_df_long = IF_df_long.rename(columns={'variable': 'x'})

    IF_df_long['x_bin' + str(bin_size)] = (IF_df_long['x'] // bin_size) + 1
    IF_df_long['y_bin' + str(bin_size)] = (IF_df_long['y'] // bin_size) + 1

    tmp = IF_df_long.groupby(['x_bin' + str(bin_size), 'y_bin' + str(bin_size)]).agg({"value":func}).reset_index()
    tmp['spot'] = tmp['x_bin' + str(bin_size)].astype('str') + '_' + tmp['y_bin' + str(bin_size)].astype('str')
    tmp.to_csv('if_bin' + str(bin_size) + '_' + file_fix + '.csv')

    tmp_wide = np.array(tmp.pivot_table(index = 'y_bin' + str(bin_size), columns = 'x_bin' + str(bin_size), values = 'value', fill_value=0), dtype='float32')

    img = Image.fromarray(tmp_wide)
    img.save("if_bin" + str(bin_size) + "_" + file_fix + ".tif", compression="tiff_lzw")

    # plt.figure(figsize=(8, 8))
    # plt.imshow(img, cmap = 'gray')
    # plt.axis('off')

    return tmp_wide

def prob_dist(I:np.ndarray) -> np.ndarray:
    """
    Estimate the probability distribution of IF image.

    Parameters
    ----------
    I
        Numpy.ndarray of IF image.

    Returns
    -------
    The probability distribution of IF image.
    """
    return np.histogramdd(np.ravel(I), bins = 256)[0] / I.size

def kl_divergence(I: np.ndarray, J: np.ndarray) -> float:
    """
    Calculate the Kullback-Leibler (KL) divergence between two IF images.

    Parameters
    ----------
    I : numpy.ndarray
        Numpy array representing the first IF image.
    J : numpy.ndarray
        Numpy array representing the second IF image.

    Returns
    -------
    float
        The computed KL divergence between the two IF images.
    """
    epsilon = 1e-10
    P = prob_dist(I) + epsilon
    Q = prob_dist(J) + epsilon
    return np.where(P != 0, P * np.log2(P / Q), 0).sum()

def kl_div_plot(I:np.ndarray, J:np.ndarray, func, filename='KL_orig_bin50.pdf'):
    """
    Generate a KL divergence plot comparing probability distributions.

    Parameters
    ----------
    I : numpy.ndarray
        Numpy array representing the first IF image.
    J : numpy.ndarray
        Numpy array representing the second IF image.
    func : str
        The function used to aggregate IF signal for the second image.
    filename : str, optional
        The name of the output plot file, by default 'KL_orig_bin50.pdf'.

    Returns
    -------
    None

    Examples
    --------
    kl_div_plot(IF_df, max_wide, func='max')
    kl_div_plot(IF_df, q75_wide, func='q75')
    kl_div_plot(IF_df, mean_wide, func='mean')
    kl_div_plot(IF_df, median_wide, func='median')
    kl_div_plot(IF_df, min_wide, func='min')
    kl_div_plot(IF_df, sum_wide, func='sum')
    """
    ax = pd.DataFrame(prob_dist(I) + 1e-10).plot()
    ax.plot(prob_dist(J)+1e-10)
    ax.legend(['orig', func])
    ax.text(s='KL_div:'+"{:.2f}".format(kl_divergence(I, J)), x=50, y=0.12)
    plt.savefig(filename)

def array2df(IF:np.ndarray, index:list=None, columns:list=None, func:str=None, filename:str='if_binary.csv', save:bool=False) -> pd.DataFrame:
    """
    Convert a 2D array to a DataFrame, reshape it, and optionally save.

    Parameters
    ----------
    IF : numpy.ndarray
        2D array to be converted to a DataFrame.
    index : list or array-like, optional
        Index for the DataFrame, by default None.
    columns : list or array-like, optional
        Columns for the DataFrame, by default None.
    func : str, optional
        Function name to use in the saved filename, by default None.
    filename : str, optional
        Name of the output CSV file, by default 'if_binary.csv'.
    save : bool, optional
        Whether to save the DataFrame to a CSV file, by default False.

    Returns
    -------
    pandas.DataFrame
        The reshaped DataFrame.
    """
    df = pd.DataFrame(IF, index = index, columns = columns)
    df.index.name = 'y'
    df.reset_index(inplace = True)
    df = df.melt(id_vars = 'y')
    df = df.rename(columns = {'variable': 'x'})
    df['spot'] = df['x'].astype('str') + '_' + df['y'].astype('str')
    if save:
        df.to_csv(func + '_' + filename)
    return df

def feature_select(X_train: pd.DataFrame, y_train: pd.DataFrame, tree_method: str = 'gpu_hist', gpu_id: int = 0, n_jobs: int = 20,
                   objective: str = 'binary:logistic', pvalue: float = 0.05, random_state: int = 132,
                   n_trials: int = 50, percentile: int = 100, plot_res: bool = True,
                   plot_filename: str = "borutashap_feature_importance.pdf") -> BorutaShap:
    """
    Perform feature selection using BorutaShap method.

    Parameters
    ----------
    X_train : pandas.DataFrame
        Training data matrix.
    y_train : pandas.Series
        Training data labels.
    tree_method : str, optional
        The tree method for XGBoost, by default 'gpu_hist'.
    gpu_id : int, optional
        GPU ID to use, by default 0.
    n_jobs : int, optional
        Number of parallel jobs, by default 20.
    objective : str, optional
        The learning task objective, by default 'binary:logistic'.
    pvalue : float, optional
        P-value threshold for feature selection, by default 0.05.
    random_state : int, optional
        Random seed, by default 132.
    n_trials : int, optional
        Number of trials for BorutaShap, by default 50.
    percentile : int, optional
        Percentile threshold for feature selection, by default 100.
    plot_res : bool, optional
        Whether to plot the BorutaShap results, by default True.
    plot_filename : str, optional
        Filename for the plot, by default "borutashap_feature_importance.pdf".

    Returns
    -------
    BorutaShap.BorutaShap
        The trained BorutaShap feature selector.
    """
    params = {'tree_method': tree_method, 'gpu_id': gpu_id, 'n_jobs': n_jobs}

    model = XGBClassifier(objective='binary:logistic', seed = 12345, **params)

    Feature_Selector = BorutaShap(model=model, importance_measure='shap', classification=True, percentile=percentile, pvalue=pvalue)

    Feature_Selector.fit(X=X_train, y=y_train.values, n_trials=50, random_state=132)

    with open('feature_selector.pkl', 'wb') as file:
        pickle.dump(Feature_Selector, file)

    pd.DataFrame(Feature_Selector.accepted).to_csv('borutashap_genes.csv')

    if plot_res:
        plot_borutashap(Feature_Selector, plot_filename)

    return Feature_Selector

def plot_borutashap(Feature_Selector: BorutaShap, filename: str= "borutashap_feature_importance.pdf"):
    """
    Plot the feature importance results from BorutaShap.

    Parameters
    ----------
    Feature_Selector : BorutaShap
        The BorutaShap feature selector.
    filename : str, optional
        Filename for the plot, by default "borutashap_feature_importance.pdf".

    Returns
    -------
    None
    """
    data = Feature_Selector.history_x.iloc[1:]
    data['index'] = data.index
    data = pd.melt(data, id_vars='index', var_name='Methods')
    decision_mapper = Feature_Selector.create_mapping_of_features_to_attribute(maps=['Tentative','Rejected','Accepted', 'Shadow'])
    data['Decision'] = data['Methods'].map(decision_mapper)
    data.drop(['index'], axis=1, inplace=True)
    options = { 'accepted' : Feature_Selector.filter_data(data,'Decision', 'Accepted'),
                'tentative': Feature_Selector.filter_data(data,'Decision', 'Tentative'),
                'rejected' : Feature_Selector.filter_data(data,'Decision', 'Rejected'),
                'all' : data
                }
    Feature_Selector.check_if_which_features_is_correct('accepted')
    data = options['accepted'.lower()]
    data = data[(data['Methods'] != 'Min_Shadow') & (data['Methods'] != 'Max_Shadow')]
    minimum = data['value'].min()
    data['value'] += abs(minimum) + 0.01
    order = data.groupby(by=["Methods"])["value"].mean().sort_values(ascending=False).index
    my_palette = Feature_Selector.create_mapping_of_features_to_attribute(maps= ['yellow','red','green','blue'])
    # Use a color palette
    plt.figure(figsize=(12,8))
    ax = sns.boxplot(x=data["Methods"], y=data["value"],
                order=order, palette=my_palette)
    ax.set(yscale="log")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, size=8)
    ax.set_title('Feature Importance')
    ax.set_ylabel('Z-Score')
    ax.set_xlabel('Features')
    plt.savefig(filename)

def model_fine_tune(X_train: pd.DataFrame, y_train: pd.DataFrame, X_test: pd.DataFrame, y_test: pd.DataFrame, gpu_id: int=0, n_splits: int=5, scoring: str='roc_auc', n_jobs: int=15, save_model: str='clf_fine_tune.json',
                    title: str='borutashap_genes') -> xgb.sklearn.XGBClassifier:
    """
    Fine-tunes an XGBoost model using Bayesian optimization and evaluates it.

    Parameters
    ----------
    X_train : pandas.DataFrame
        Training data matrix.
    y_train : pandas.DataFrame
        Training data labels.
    X_test : pandas.DataFrame
        Testing data matrix.
    y_test : pandas.DataFrame
        Testing data labels.
    gpu_id : int, optional
        GPU identifier for XGBoost, by default 0.
    n_splits : int, optional
        Number of splits for cross-validation, by default 5.
    scoring : str, optional
        Scoring metric for cross-validation, by default 'roc_auc'.
    n_jobs : int, optional
        Number of parallel jobs for cross-validation, by default 15.
    save_model : str, optional
        File name to save the trained model, by default 'clf_fine_tune.json'.
    title : str, optional
        Title for the evaluation plots, by default 'borutashap_genes'.

    Returns
    -------
    xgboost_bo : XGBClassifier
        The trained XGBoost model.
    """
    # Space
    space = {
        'learning_rate': hp.choice('learning_rate', [0.01, 0.05, 0.1, 0.15, 0.2]),
        'max_depth' : hp.choice('max_depth', np.arange(5,101,5).tolist()),
        'gamma' : hp.choice('gamma', [0, 0.25, 0.5, 1]),
        'colsample_bytree' : hp.choice('colsample_bytree', [0.5, 0.8, 1]),     
        'reg_alpha' : hp.choice('reg_alpha', [1e-5, 1e-2, 0, 0.1, 1, 10, 100]), 
        'reg_lambda' : hp.choice('reg_lambda', [1e-5, 1e-2, 0.1, 1, 10, 100]),
        'n_estimators' : hp.choice('n_estimators', np.arange(50, 2001, 50).tolist()),
        "scale_pos_weight": hp.choice('scale_pos_weight', [1, 2.4, 3, 5]),
        "subsample": hp.choice('subsample', [0.5, 0.6, 0.7, 0.8, 0.9, 1]),
        'min_child_weight': hp.choice('min_child_weight', [0, 1, 2, 5, 11])
    }
    # Set up the k-fold cross-validation
    kfold = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)

    # Objective function
    def objective(params):
        
        xgboost = XGBClassifier(seed=0, **params, tree_method='gpu_hist', gpu_id=gpu_id, max_bin = 100)
        score = cross_val_score(estimator=xgboost,  
                                X=X_train, 
                                y=y_train, 
                                cv=kfold, 
                                scoring=scoring, 
                                n_jobs=n_jobs).mean()
        # Loss is negative score
        loss = - score
        # Dictionary with information for evaluation
        return {'loss': loss, 'params': params, 'status': STATUS_OK}
    # Optimize
    best = fmin(fn = objective, space = space, algo = tpe.suggest, max_evals = 100, trials = Trials())
    # Print the values of the best parameters
    print(space_eval(space, best))

    # Train model using the best parameters
    xgboost_bo = XGBClassifier(seed=0, 
                            colsample_bytree=space_eval(space, best)['colsample_bytree'], 
                            gamma=space_eval(space, best)['gamma'], 
                            learning_rate=space_eval(space, best)['learning_rate'], 
                            max_depth=space_eval(space, best)['max_depth'], 
                            reg_alpha=space_eval(space, best)['reg_alpha'],
                            reg_lambda=space_eval(space, best)['reg_lambda'],
                            n_estimators=space_eval(space, best)['n_estimators'],
                            scale_pos_weight=space_eval(space, best)['scale_pos_weight'],
                            subsample=space_eval(space, best)['subsample'],
                            min_child_weight=space_eval(space, best)['min_child_weight'],
                            n_jobs = 15, tree_method='gpu_hist', gpu_id=gpu_id).fit(X_train, y_train)
    xgboost_bo.save_model(save_model)
    
    bayesian_opt_predict_prob = xgboost_bo.predict_proba(X_test)[:,1]
    plotROC_curve(y_test, bayesian_opt_predict_prob, pos_label=1, title = title)
    plotPRC_curve(y_test, bayesian_opt_predict_prob, pos_label=1, title = title)
    confusion_matrix_plot(y_test, bayesian_opt_predict_prob, pos_label=1, title = title)
    return xgboost_bo

def model_predict(obj: anndata, model: xgb.sklearn.XGBClassifier, auc_cutoff: float=None, pt_size: float=3, gene_list: list=None, save_col: str='res', sub_col: str='malignant', sub_group: str='Malignant cells') -> anndata:
    """
    Predicts cell labels using a trained model and updates the object's annotations.

    Parameters
    ----------
    obj : AnnData
        Annotated data object.
    model : object
        Trained machine learning model.
    auc_cutoff : float, optional
        AUC cutoff threshold.
    pt_size : int, optional
        Size of points in visualization, by default 3.
    gene_list : list, optional
        List of genes for prediction.
    save_col : str, optional
        Column name for saving prediction results, by default 'res'.
    sub_col : str, optional
        Column name for subgroup annotation, by default 'malignant'.
    sub_group : str, optional
        Subgroup label to use for prediction, by default 'Malignant cells'.

    Returns
    -------
    obj : AnnData
        Annotated data object with updated predictions.
    """
    if (sub_col is not None) & (sub_group is not None):
        tmp_obj = obj[obj.obs[sub_col].isin([sub_group])]
        tmp_data = sc.get.obs_df(tmp_obj, gene_list)
        pred_tmp = model.predict_proba(tmp_data)
        binary_tmp = np.where(pred_tmp[:, 1]>auc_cutoff, 'Positive cells', 'Negative cells')
        tmp_obj.obs['predict_'+save_col] = binary_tmp
        tmp_obj.obs['predict_prob_'+save_col] = pred_tmp[:, 1]
        obj.obs = obj.obs.merge(tmp_obj.obs[['predict_'+save_col, 'predict_prob_'+save_col]], how = 'left', left_index = True, right_index=True)
        print(obj.obs['predict_'+save_col].value_counts())
        sc.pl.embedding(tmp_obj, color = 'predict_'+save_col, basis = 'spatial', size=pt_size, save='predict.pdf')
    else:
        tmp_data = sc.get.obs_df(obj, gene_list)
        pred_tmp = model.predict_proba(tmp_data)
        binary_tmp = np.where(pred_tmp[:, 1]>auc_cutoff, 'Positive cells', 'Negative cells')
        obj.obs['predict_'+save_col] = binary_tmp
        obj.obs['predict_prob_'+save_col] = pred_tmp[:, 1]
        print(obj.obs['predict_'+save_col].value_counts())
        sc.pl.embedding(obj, color = 'predict_'+save_col, basis = 'spatial', size=pt_size, save='predict.pdf')
    return obj

def plotPRC_curve(target_test, test_preds, pos_label: int=1, title: str='prolif_gene_default_cutoff', prefix: str='fine_tune_'):
    """
    Plot Precision-Recall Curve for binary classification.

    Parameters
    ----------
    target_test : pandas.DataFrame
        True binary labels of the test set.
    test_preds : array-like
        Predicted probabilities or scores for the positive class.
    pos_label : int
        Label of the positive class.
    title : str, optional
        Title for the plot, by default 'prolif_gene_default_cutoff'.
    prefix : str, optional
        Prefix for saving the plot file, by default 'fine_tune_'.
    """
    precision, recall, thresholds = metrics.precision_recall_curve(target_test, test_preds, pos_label=pos_label)
    auprc = metrics.average_precision_score(target_test, test_preds)

    plt.figure(dpi=300, figsize=(5.5,5))
    plt.title('XGBoostClassifier ('+title+')', fontweight='bold', fontsize = 11)
    plt.plot(recall, precision, 'b', label = 'AUPRC = %0.2f' % auprc)
    plt.legend(loc = 'lower right', fontsize = 8)
    # plt.plot([0, 1], [0, 1],'r--')
    plt.ylabel('Precision', fontsize = 10)
    plt.xlabel('Recall', fontsize = 10)
    plt.xticks(size = 10)
    plt.yticks(size = 10)
    plt.grid(False)
    plt.gcf().savefig(prefix + title+'_PRC.pdf')

def plotROC_curve(target_test, test_preds, pos_label: int=1, title: str='prolif_gene_default_cutoff', prefix: str='fine_tune_'):
    """
    Plot Receiver Operating Characteristic (ROC) Curve for binary classification.

    Parameters
    ----------
    target_test : pandas.DataFrame
        True binary labels of the test set.
    test_preds : array-like
        Predicted probabilities or scores for the positive class.
    pos_label : int
        Label of the positive class.
    title : str, optional
        Title for the plot, by default 'prolif_gene_default_cutoff'.
    prefix : str, optional
        Prefix for saving the plot file, by default 'fine_tune_'.
    """
    fpr, tpr, threshold = metrics.roc_curve(target_test, test_preds, pos_label=pos_label)
    roc_auc = metrics.auc(fpr, tpr)

    plt.figure(dpi=300, figsize=(5.5,5))
    plt.title('XGBoostClassifier ('+title+")", fontweight='bold', fontsize = 11)
    plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
    plt.legend(loc = 'lower right', fontsize = 8)
    plt.plot([0, 1], [0, 1],'r--')
    plt.ylabel('True Positive Rate', fontsize = 10)
    plt.xlabel('False Positive Rate', fontsize = 10)
    plt.xticks(size = 10)
    plt.yticks(size = 10)
    plt.grid(False)
    plt.gcf().savefig(prefix + title+'_roc.pdf')

def confusion_matrix_plot(y_true, y_pred_prob, pos_label: int=1, title: str="prolif_gene_default_cutoff", prefix: str='fine_tune_'):
    """
    Plot Confusion Matrix for binary classification.

    Parameters
    ----------
    y_true : pandas.DataFrame
        True binary labels.
    y_pred_prob : array-like
        Predicted probabilities or scores for the positive class.
    pos_label : int, optional
        Label of the positive class, by default 1.
    title : str, optional
        Title for the plot, by default 'prolif_gene_default_cutoff'.
    prefix : str, optional
        Prefix for saving the plot file, by default 'fine_tune_'.
    """
    fpr, tpr, threshold = metrics.roc_curve(y_true, y_pred_prob, pos_label=pos_label)

    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = threshold[optimal_idx]
    y_pred_binary = np.where(y_pred_prob > optimal_threshold, 1, 0)

    cf_matrix = metrics.confusion_matrix(y_true, y_pred_binary)
    plt.figure(dpi = 120, figsize=(5,5))
    sns.heatmap(cf_matrix, annot=True, fmt='.20g', cmap='Blues')
    plt.xlabel("Prediction")
    plt.ylabel("True")
    plt.title(title+ '/'+ str(optimal_threshold))
    plt.savefig(prefix + title + "_confusion_matrix.pdf")

def shap_genes_plot(model, X_test, feature_perturbation: str="tree_path_dependent", approximate: bool=True, filename: str="shap_genes.pdf"):
    """
    Generate a SHAP summary plot for feature importance.

    Parameters
    ----------
    model : xgb.sklearn.XGBClassifier
        Trained machine learning model.
    X_test : array-like
        Testing data features.
    feature_perturbation : str, optional
        Type of SHAP feature perturbation, by default "tree_path_dependent".
    approximate : bool, optional
        Whether to use the approximate SHAP values, by default True.
    filename : str, optional
        Filename for saving the SHAP summary plot, by default "shap_genes.pdf".
    """
    # Create object that can calculate shap values
    explainer = shap.TreeExplainer(model, feature_perturbation = feature_perturbation, approximate = approximate) # interventional

    # calculate shap values. This is what we will plot.
    # Calculate shap_values for all of val_X rather than a single row, to have more data for plot.
    shap_values_test = explainer.shap_values(X_test)

    shap.summary_plot(shap_values_test, X_test, show=False)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()