from skfeature.function.similarity_based import lap_score
import numpy as np
from skfeature.function.similarity_based import SPEC
from skfeature.function.sparse_learning_based import MCFS
from skfeature.function.sparse_learning_based import NDFS
from skfeature.utility import construct_W
from skfeature.utility.sparse_learning import feature_ranking


def fs_lap_score(X,dim):
    kwargs_W = {"metric": "euclidean", "neighbor_mode": "knn", "weight_mode": "heat_kernel", "k": 5, 't': 1}
    W = construct_W.construct_W(X, **kwargs_W)
    score = lap_score.lap_score(X, W=W)
    idx = lap_score.feature_ranking(score)
    num_fea = dim
    selected_features = X[:, idx[0:num_fea]]
    return selected_features

def fs_low_variance(X,dim):
    idx = np.argsort(-np.var(X, axis=0))
    #print(np.var(X, axis=0))
    num_fea = dim
    selected_features = X[:, idx[0:num_fea]]
    return selected_features

def fs_SPEC(X,dim):
    kwargs = {'style': 0}
    score = SPEC.spec(X, **kwargs)
    idx = SPEC.feature_ranking(score, **kwargs)
    num_fea = dim
    selected_features = X[:, idx[0:num_fea]]
    return selected_features

def fs_MCFS(X,dim,n):
    kwargs = {"metric": "euclidean", "neighborMode": "knn", "weightMode": "heatKernel", "k": 5, 't': 1}
    W = construct_W.construct_W(X, **kwargs)
    num_fea = dim
    Weight = MCFS.mcfs(X, n_selected_features=num_fea, W=W, n_clusters=n)
    idx = MCFS.feature_ranking(Weight)
    selected_features = X[:, idx[0:num_fea]]
    return selected_features

def fs_NDFS(X,dim,n):
    kwargs = {"metric": "euclidean", "neighborMode": "knn", "weightMode": "heatKernel", "k": 5, 't': 1}
    W = construct_W.construct_W(X, **kwargs)
    Weight = NDFS.ndfs(X, W=W, n_clusters=n)
    idx = feature_ranking(Weight)
    num_fea = dim
    selected_features = X[:, idx[0:num_fea]]
    return selected_features