import numpy as np
import matplotlib.pyplot as plt
import csv
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy import sparse
import louvain
import igraph as ig
from velocyto.estimation import colDeltaCorSqrt
import bezier
from sklearn.manifold import TSNE, JTSNE

def visualize_protein_markers_tsne(vlm, protein_markers, pc_targets, visualize_clusters=False,
                              colormap='inferno'):
    """
    Visualizes an array of protein markers in protein principal component space. The components to plot are given by
    the list pc_targets. If visualize_clusters is selected, an additional cluster-colored plot is generated.
    Useful for iterative manual procedure to identify clusters based on characteristic markers.
    """
    array_proteins = vlm.adt_names
    pcs = vlm.prot_tsne
    pc_zi = [pc_targets[0]-1, pc_targets[1]-1]
    
    n_addit = int(visualize_clusters)
    
    nrows = int(np.ceil((len(protein_markers)+n_addit)/5))
#     print(nrows)
    
    f, ax = plt.subplots(nrows=nrows,ncols=5,figsize=(12,0.25+2*nrows))
    ax = ax.flatten()


    for j in range(len(protein_markers)):
        prot_name=protein_markers[j]
        sc = ax[j].scatter(pcs[:,pc_zi[0]],pcs[:,pc_zi[1]],s=3,c=np.log(vlm.P[array_proteins == prot_name][0]+1),
                      alpha=0.2,cmap=colormap)
        plt.colorbar(sc)
        ax[j].set_title(prot_name)
    
    if visualize_clusters:
        
        if hasattr(vlm,'cluster_ID') and hasattr(vlm,'COLORS'):
            col=vlm.COLORS[vlm.cluster_ID]
        else:
            COLORS = np.rand(np.amax(cluster_ID)+1,3)
            col=COLORS[vlm.cluster_ID]
            
                    
        ax[-1].scatter(pcs[:,pc_zi[0]],pcs[:,pc_zi[1]],s=3,c=col,alpha=0.9)

        
    for k in range(len(ax)):
        ax[k].axis('off')
        
def fit_jointTSNE(vlm, pc_space_name1, pc_space_name2, n_pcs1, n_pcs2, ts_name, seed=None):
    #bh_tsne = TSNE(random_state=seed)
    data = {'X1': getattr(vlm,pc_space_name1)[:, :n_pcs1], 'X2': getattr(vlm,pc_space_name2)[:, :n_pcs2]}
    bh_tsne = JTSNE(init='pca').fit_transform(data)
    setattr(vlm,ts_name,bh_tsne)
        
