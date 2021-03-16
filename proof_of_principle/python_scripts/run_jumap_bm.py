## Run benchmark j-UMAP
import numpy as np
from numpy import genfromtxt
from sklearn.decomposition import PCA
from Jvis import JUMAP, UMAP
import sys
from timeit import default_timer as timer
import time

data_set = sys.argv[1] #'GeqN1kD1k'

data_name = "output_jumap/" + data_set #output

def shuffle_rows(expression_matrix, prop, random_seed):
    np.random.seed(random_seed)
    n_rows, n_cols = expression_matrix.shape
    n_elem = round(prop*n_cols)
    s = np.arange(n_cols)
    row_id = list(np.random.choice(s, size=n_elem, replace=False))
    v = expression_matrix[:, row_id]
    np.random.shuffle(np.transpose(v))
    expression_matrix[:, row_id] = v
    
def normalize_matrix(A, option = "Frobenius"):
    if option == "Frobenius":
        return (A/np.linalg.norm(A, 'fro'))
    elif option == "max":
        return (A/A.max())
    else:
        "Not implemented yet"

# Read the expr matrix
expr_mat = genfromtxt('../sim_splat/data/'+data_set+'_rna.csv',  delimiter=',')
labels_true = genfromtxt('../sim_splat/data/'+data_set+'_labels.csv', delimiter=',')
# labels_true = labels_true.values.flatten()
expr_mat_log_t = expr_mat
print("RNA_shape", expr_mat.shape)
expr_reduced = PCA(n_components=50).fit_transform(expr_mat_log_t)

# Read ADT
adt_mat = genfromtxt('../sim_splat/data/'+data_set+'_adt.csv', delimiter=',')
print("ADT shape: ", adt_mat.shape)

start = time.time()

for noise_level in [0, 0.1, 0.2, 0.3, 0.4]:
    # Create noise
    print("Dataset: ", data_set)
    print("Noise level: ", noise_level)
    expr_mat_shuffle = np.copy(expr_mat_log_t)
    # expr_mat_shuffle = np.copy(adt_mat)
    shuffle_rows(expr_mat_shuffle.T, prop = noise_level, random_seed=0)

    if expr_mat_shuffle.shape[1] <500:
        expr_reduced_shufle = expr_mat_shuffle
    else:
        expr_reduced_shufle = PCA(n_components=50).fit_transform(expr_mat_shuffle)
    maxIter = 10
    print("Run optimized Jvis")
    for _lambda in [0.2, 0.5, 1, 2]:
        joint_umap_obj = JUMAP(init='random')
        data = {'modal-1': expr_reduced, 'modal-2': adt_mat,  'modal-noise': expr_reduced_shufle}
        joint_umap = joint_umap_obj.fit_transform(X = data, method = 'auto', ld = _lambda, max_iter=maxIter)
        # save data
        
        np.savetxt(data_name +"_OPTjumap_ld" + str(_lambda)+"_noise" + str(noise_level) +".csv", joint_umap, delimiter=",")
        np.savetxt(data_name +"_alphajumap_lambda" + str(_lambda) +"_noise" + str(noise_level)+ ".csv", np.array(joint_umap_obj.alpha), delimiter=",")

    # Run Jvis (uniform)
    print("Run uniform Jvis")
    joint_umap_obj_uni = JUMAP(init='random')
    joint_umap_uni = joint_umap_obj_uni.fit_transform(X = data, method = 'uniform')
    # Save data
    np.savetxt(data_name +"_uniformjumap"+"_noise" + str(noise_level)+ ".csv", joint_umap_uni, delimiter=",")

    # Run concat
    print("Run concat")
    concat = np.concatenate((normalize_matrix(expr_mat_log_t), normalize_matrix(adt_mat), normalize_matrix(expr_mat_shuffle)), axis = 1)
    expr_reduced = PCA(n_components=100).fit_transform(concat)
    concat_umap = UMAP(init='random').fit_transform(expr_reduced)
    # save data
    np.savetxt(data_name +"_concatFroumap"+"_noise" + str(noise_level)+ ".csv", concat_umap, delimiter=",")
    print("Run concat")
    concat = np.concatenate((normalize_matrix(expr_mat_log_t, "max"), normalize_matrix(adt_mat, "max"), normalize_matrix(expr_mat_shuffle, "max")), axis = 1)
    expr_reduced = PCA(n_components=100).fit_transform(concat)
    concat_umap = UMAP(init='random').fit_transform(expr_reduced)
    # save data
    np.savetxt(data_name +"_concatMaxumap"+"_noise" + str(noise_level)+ ".csv", concat_umap, delimiter=",")
    print("Run concat")
    concat = np.concatenate((expr_mat_log_t, adt_mat, expr_mat_shuffle), axis = 1)
    expr_reduced = PCA(n_components=100).fit_transform(concat)
    concat_umap = UMAP(init='random').fit_transform(expr_reduced)
    # save data
    np.savetxt(data_name +"_concatumap"+"_noise" + str(noise_level)+ ".csv", concat_umap, delimiter=",")

print("Total runtime: ", time.time() - start, " s.")
