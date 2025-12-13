# -*- coding: utf-8 -*-
"""
ORIGINAL AUTHOR: Sanjay Kottapalli
MODIFIED BY: Theo Botella for GoTChAo
OPTIMIZED: Vectorized logic + Analytical Moments + One-Pass KNN
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

from sklearn.neighbors import KernelDensity
from numpy.polynomial.polynomial import Polynomial
from numpy.polynomial.polynomial import polyfit
from scipy.optimize import minimize
from collections import Counter
from sklearn.neighbors import KNeighborsClassifier
from scipy.stats import gmean
import matplotlib as mpl

# SPEED: Turn off interactive plotting
plt.ioff()

mpl.rcParams['figure.dpi'] = 300 
sns.set(font_scale=1.2)
sns.set_style("ticks")
tab20=plt.get_cmap('tab20')
gen_cmap = {'MUT':tab20.colors[0], 'WT':tab20.colors[2],
            'HET':tab20.colors[4], 'NA':tab20.colors[6], '-1':tab20.colors[8]}

def GotchaLabeling(path="", infile="", sample_id=""):
    """
    Main controller function.
    path: The directory where outputs go.
    infile: The CSV filename (must be inside path).
    sample_id: Prefix for output files.
    """
    print("--- Starting Genotyping Analysis ---")
    
    # Ensure path ends with slash for concatenation
    if not path.endswith('/'):
        path += '/'

    full_path = os.path.join(path, infile)
    typing = read_data(full_path)
    
    # UPDATED: Use the provided path directly, do not create a subfolder
    sample_dir = path
    print(f"   Outputs will be saved to: {sample_dir}")
    
    for i in ["WT", "MUT"]:
        print("Noise correcting {} read counts.".format(i))
        if i=="WT":
            typing, wt_min = noise_correct(typing, i, sample_dir, sample_id)
        else:
            typing, mut_min = noise_correct(typing, i, sample_dir, sample_id)
    
    print("Performing quadrant genotyping.")
    typing = quadrant_genotype(typing, wt_min, mut_min)
        
    print("Computing KNN-based clusters.")
    typing = KNN_cluster(typing, wt_min, mut_min, 0.05, sample_dir, sample_id)
    
    out_file = os.path.join(sample_dir, f"{sample_id}_genotype_labels.csv")
    typing.to_csv(out_file)
    print(f"Saved Genotyped Labels to: {out_file}")
    
    print("All analysis complete!")

    return typing

def read_data(infile):
    print(f"Loading: {infile}")
    try:
        genotyping = pd.read_csv(infile)
        if 'Barcode' in genotyping.columns:
            genotyping.set_index('Barcode', inplace=True)
        elif 'Read_Name' not in genotyping.columns: 
            genotyping.set_index(genotyping.columns[0], inplace=True)
            
        genotyping['WTcount'] = genotyping['WT_Count']
        genotyping['MUTcount'] = genotyping['MUT_Count']
        genotyping = genotyping.dropna()
    except Exception as e:
        print(f"Error reading data: {e}")
        raise e
    print("Number of cells: " + str(genotyping.shape[0]))
    return genotyping

def noise_correct(typing, feature="", sample_dir="", sample_id=""):
    np.random.seed(0)
    pseudocount = 1
    X = typing[[feature+'count']]+pseudocount
    logged_counts = np.log(X) 

    typing['transf_{}'.format(feature)] = logged_counts

    # Plot Log Counts
    plt.figure()
    plt.hist(logged_counts, density=True, bins=50)
    plt.title(f"{feature} counts")
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    plt.savefig(os.path.join(sample_dir, f"{sample_id}_log_{feature}_counts.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
        
    bw = 0.1
    kde = KernelDensity(bandwidth=bw)
    kde.fit(typing['transf_{}'.format(feature)].values.reshape(-1, 1))
    
    x_bin = np.histogram(typing['transf_{}'.format(feature)], bins=50)[1]
    kde_x = np.linspace(min(x_bin)-0.5,max(x_bin)+0.5,10000)
    kde_smooth = np.exp(kde.score_samples(kde_x.reshape(-1, 1)))

    # Plot KDE Initial
    plt.figure()
    plt.hist(typing['transf_{}'.format(feature)], density=True, bins=50)
    plt.plot(kde_x, kde_smooth, color='red')
    plt.title('Initial KDE')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    plt.savefig(os.path.join(sample_dir, f"{sample_id}_{feature}_kde_initial.pdf"), dpi=300, bbox_inches="tight")
    plt.close()

    fit = polyfit(kde_x, kde_smooth, 12)
    poly = Polynomial(fit)
    mid = (max(kde_x)-min(kde_x))/2
    
    result = minimize(poly, x0=np.array([mid]), options={'disp':False}, bounds=((min(kde_x), max(kde_x)),), method='Nelder-Mead')
    poly_min = result.x[0]
    
    def kde_scalar_func(x_scalar):
        return np.exp(kde.score_samples([[x_scalar[0]]]))[0]
    
    result = minimize(kde_scalar_func, x0=np.array([poly_min]), options={'disp':True}, bounds=((min(kde_x), max(kde_x)),), method='Nelder-Mead')
    new_min = result.x[0]
    
    if feature=="WT": alt_feature="MUT"
    else: alt_feature="WT"

    if new_min >= max(kde_x) - 0.1:
        zero_counts = typing['{}count'.format(alt_feature)]
        zero_counts = zero_counts[zero_counts==0].index
        if len(zero_counts) > 0:
            new_min = np.percentile(typing.loc[zero_counts,'transf_{}'.format(feature)], 99.99)
        else:
            new_min = mid

    noise_values = logged_counts[logged_counts < new_min].dropna()
    
    # Plot Mixture
    plt.figure()
    plt.hist(typing['transf_{}'.format(feature)], density=True, bins=50)
    plt.plot(kde_x, kde_smooth, color='red')
    plt.axvline(new_min, color='blue', linestyle='--')
    plt.title('Noise vs. Signal Threshold')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    plt.savefig(os.path.join(sample_dir, f"{sample_id}_{feature}_kde_mixture.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
    
    noise_mean = noise_values.values.mean()
    noise_var = noise_values.values.var(ddof=0) + bw**2
    
    added_var = 0.3**2
    new_var = noise_var - added_var
    if new_var < 0: new_var = 0.0001
    new_std = np.sqrt(new_var)

    typing['transf_{}'.format(feature)] = (np.log(X)-noise_mean)/new_std

    # Plot Z-Scores
    plt.figure()
    plt.hist(typing['transf_{}'.format(feature)], bins=50, density=True)
    plt.title(f'{feature} Z-scores')
    plt.ylabel("Probability")
    plt.xlabel("Z-score")
    plt.savefig(os.path.join(sample_dir, f"{sample_id}_{feature}_zscores.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
    
    feat_min = (new_min-noise_mean)/new_std
    return typing, feat_min

def quadrant_genotype(typing, wt_min, mut_min):
    wt_vals = typing['transf_WT'].values
    mut_vals = typing['transf_MUT'].values
    
    cond_na  = (wt_vals < wt_min) & (mut_vals < mut_min)
    cond_mut = (wt_vals < wt_min) & (mut_vals >= mut_min)
    cond_wt  = (wt_vals >= wt_min) & (mut_vals < mut_min)
    cond_het = (wt_vals >= wt_min) & (mut_vals >= mut_min)
    
    conditions = [cond_na, cond_mut, cond_wt, cond_het]
    choices = ['NA', 'MUT', 'WT', 'HET']
    
    typing['quadrant_class'] = np.select(conditions, choices, default='NA')
    return typing

def KNN_cluster(typing, wt_min, mut_min, knn_window, sample_dir, sample_id):
    typing['clusters'] = typing['quadrant_class']
    data = typing[['transf_WT', 'transf_MUT']].values
    
    range_wt = data[:,0].max() - data[:,0].min()
    range_mut = data[:,1].max() - data[:,1].min()
    
    wt_lower = wt_min - knn_window * range_wt
    wt_upper = wt_min + knn_window * range_wt
    mut_lower = mut_min - knn_window * range_mut
    mut_upper = mut_min + knn_window * range_mut
    
    mask_wt = (data[:,0] > wt_lower) & (data[:,0] <= wt_upper)
    mask_mut = (data[:,1] > mut_lower) & (data[:,1] <= mut_upper)
    uncertain_mask = mask_wt | mask_mut
    
    counts = Counter(typing['quadrant_class'].values)
    n = list(counts.values())
    n_neighbors = 5 if not n else round(np.sqrt(gmean(n)))
        
    train_mask = ~uncertain_mask
    X_train = data[train_mask]
    y_train = typing.loc[train_mask, 'quadrant_class']
    
    # Plot Confident Calls
    plt.figure()
    for label in np.unique(y_train):
        subset = data[typing['quadrant_class'] == label]
        subset_train = subset[np.isin(subset, X_train).all(axis=1)]
        color = gen_cmap.get(str(label), 'gray')
        if len(subset) > 0:
            plt.scatter(subset[:, 0], subset[:, 1], label=label, s=7, color=color)
    plt.legend()
    plt.title('Confident Calls (Training Set)')
    plt.xlabel('WT Z-scores')
    plt.ylabel('MUT Z-scores')
    plt.savefig(os.path.join(sample_dir, f"{sample_id}_confident_calls.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
    
    if len(X_train) > n_neighbors:
        knn = KNeighborsClassifier(n_neighbors=n_neighbors, weights='distance', n_jobs=-1)
        knn.fit(X_train, y_train)
        typing['genotype_pred'] = knn.predict(data)
    else:
        print("Warning: Not enough confident data points for KNN. Falling back to quadrants.")
        typing['genotype_pred'] = typing['quadrant_class']
        
    del typing['clusters']
    
    def plot_final(col_name, fname, title):
        plt.figure()
        labels = typing[col_name].unique()
        for label in labels:
            label_str = str(label)
            subset = data[typing[col_name] == label]
            color = gen_cmap.get(label_str, 'gray')
            plt.scatter(subset[:, 0], subset[:, 1], label=label_str, s=7, color=color)
        plt.legend()
        plt.title(title)
        plt.xlabel('WT Z-scores')
        plt.ylabel('MUT Z-scores')
        plt.axhline(y=mut_min, color='r', linestyle='--')
        plt.axvline(x=wt_min, color='r', linestyle='--')
        plt.savefig(os.path.join(sample_dir, f"{sample_id}_{fname}"), dpi=300, bbox_inches="tight")
        plt.close()

    plot_final('genotype_pred', "cluster_genotype.pdf", 'Cluster Genotype')
    plot_final('quadrant_class', "quadrant_genotype.pdf", 'Quadrant Genotype')
    
    return typing