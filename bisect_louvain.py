from anndata import AnnData
import numpy as np
import os
import scanpy.api as sc
from scipy.sparse import vstack
from sklearn import preprocessing
import sys 
from scipy import sparse

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

def semi_louvain_exact_K(X_train, y_train, X_test, n_clusters):
    X_new = np.concatenate((X_train, X_test), axis=0)
    adata = AnnData(X=X_new)
    sc.pp.neighbors(adata, n_neighbors=30, use_rep='X')
    
    # Construct similarity matrix
    #for i in range(len(y_train)):
        #for j in range(len(y_train)):
            #if y_train[i] == y_train[j]:
                #adata.uns['neighbors']['connectivities'][i,j] = 1.0
            #else:
                #adata.uns['neighbors']['connectivities'][i,j] = 0.0
                
    # Construct sim matrix: convert to dense matrix
    print('Change similarity matrix')
    matrix_neighbors = adata.uns['neighbors']['connectivities'].toarray()
    for i in range(len(y_train)):
        for j in range(len(y_train)):
            if y_train[i] == y_train[j]:
                matrix_neighbors[i,j] = 1.0
            else:
                matrix_neighbors[i,j] = 0.0
                
    adata.uns['neighbors']['connectivities'] = sparse.csr_matrix(matrix_neighbors)
    

    resol = eval_bisect_louvain(adata, 0.5, 1.7, n_clusters)
    print("resol = ", resol)
    
    sc.tl.louvain(adata, key_added='louvain', resolution = resol)
    louv_labels = np.array(adata.obs['louvain'].tolist())
    le = preprocessing.LabelEncoder().fit(louv_labels)
    cell_labels_pred  = le.transform(louv_labels)
    return cell_labels_pred
            
