from anndata import AnnData
import numpy as np
import os
import scanpy.api as sc
from scipy.sparse import vstack
from sklearn import preprocessing
import sys 
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics.cluster import adjusted_rand_score

# Run louvain algorithm
def eval_louvain(adata, resol):
    sc.tl.louvain(adata, key_added='louvain', resolution = resol)
    louv_labels = np.array(adata.obs['louvain'].tolist())

    le = preprocessing.LabelEncoder().fit(louv_labels)
    cell_labels_pred  = le.transform(louv_labels)
    return int(np.unique(cell_labels_pred).shape[0])


def eval_bisect_louvain(adata, minresolution, maxresolution, n_clusters):
    M = -1
    k1 = n_clusters + 1
    k2 = 0
    iterations = 0 
    #adata = AnnData(X = data)
    #sc.pp.neighbors(adata, n_neighbors=20, use_rep='X')

    # find minresolution and maxresolution s.t k1 < n_clusters < k2 
    while k1 > n_clusters and iterations < 8:
        minresolution = minresolution/2 
        k1 = eval_louvain(adata, minresolution)
        if k1 == n_clusters:
            return minresolution
        iterations = iterations + 1
    while k2 < n_clusters and iterations < 8:
        maxresolution = 2*maxresolution
        k2 = eval_louvain(adata, maxresolution)
        if k2 == n_clusters:
            return maxresolution
        iterations = iterations + 1

    # bisection main
    while iterations < 40 and abs(maxresolution - minresolution) > 1e-5:
        M = (minresolution + maxresolution)/2 
        iterations = iterations + 1
        k3 = eval_louvain(adata, M)
        if k3 == n_clusters:
            return M 
        elif k3 < n_clusters:
            minresolution = M 
        else:
            maxresolution = M

    if iterations >= 40:
        print("bisection algorithm could not find the right resolution")

def louvain_exact_K(X_dimred, n_clusters):
    adata = AnnData(X=X_dimred)
    sc.pp.neighbors(adata,n_neighbors=10, use_rep='X')
    
    resol = eval_bisect_louvain(adata, 0.5, 1.7, n_clusters)
    
    sc.tl.louvain(adata, key_added='louvain', resolution = resol)
    louv_labels = np.array(adata.obs['louvain'].tolist())
    le = preprocessing.LabelEncoder().fit(louv_labels)
    cell_labels_pred  = le.transform(louv_labels)
    return cell_labels_pred
