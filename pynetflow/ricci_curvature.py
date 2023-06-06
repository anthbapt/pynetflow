# Inspired from GeneralizedFormanCurvature Python library
from scipy.sparse import dok_matrix
from collections import defaultdict
from itertools import combinations
from operator import add
import numpy as np
import xgi


def faces(simplices):
    faceset = set()
    for simplex in simplices :
        numnodes = len(simplex)
        for r in range(numnodes, 0, -1) :
            for face in combinations(simplex, r) :
                faceset.add(tuple(sorted(face)))
                
    return faceset


def n_faces(face_set, n):
    filt = filter(lambda face: len(face)==n+1, face_set)
    
    return filt


def boundary_operator(face_set, i):
    source_simplices = list(n_faces(face_set, i))
    target_simplices = list(n_faces(face_set, i-1))
    if len(target_simplices)==0 :
        S = dok_matrix((1, len(source_simplices)), dtype=np.float64)
        S[0, 0:len(source_simplices)] = 1
    else:
        source_simplices_dict = {source_simplices[j]: j for j in range(len(source_simplices))}
        target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}
        S = dok_matrix((len(target_simplices), len(source_simplices)), dtype=np.float64)
        for source_simplex in source_simplices :
            for a in range(len(source_simplex)) :
                target_simplex = source_simplex[:a]+source_simplex[(a+1):]
                i = target_simplices_dict[target_simplex]
                j = source_simplices_dict[source_simplex]
                S[i, j] = -1 if a % 2==1 else 1
                
    return S


def compute_laplacian(S, p) :
    """
    Compute and Store Hodge Laplacian up to dimension p. 
    """
    laplacian = []
    laplacian.append(np.matmul(boundary_operator(S, 1).toarray(), np.transpose(boundary_operator(S, 1).toarray())))
    for i in range(1, p+2) :
        b1 = boundary_operator(S, i+1).toarray()
        b2 = boundary_operator(S, i).toarray()
        b1_ = np.matmul(b1,np.transpose(b1))
        b2_ = np.matmul(np.transpose(b2), b2)
        laplacian.append(b1_+b2_)
        
    return laplacian


def compute_bochner(laplacian):

    la = laplacian
    boch = []
    for i in range(len(la)) :
        temp = la[i]-np.diag(np.diag(la[i])) #Bochner-Weitzenböck Decomposition
        boch.append(temp + np.diag(np.sum(np.abs(temp), axis=1)))

    return boch


def func_compute_forman(S, p, laplacian, simplex):
    """
    Lookup the values in Hodge Laplacian and Compute the Forman Ricci Curvature
    """
    
    l = len(simplex)
    if l <= p+1 and l >= 2 :
        target_simplices = list(n_faces(S, l-1))
        target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}
        m = target_simplices_dict[simplex]
        mat = np.delete(laplacian[l-1][m], m) 
        val = laplacian[l-1][m][m] - sum(np.abs(mat)) #Bochner-Weitzenböck Decomposition
        
        return val
    

def compute_forman(S, p, laplacian):
    """ Compute Generalised Forman Ricci Curvature for simplices up to dimension p 

    Parameters
    -----------
    p : maximum dimension of simplex to compute curvature.

    Output
    -----------
    forman_dict: dict[simplex, forman curvature value]
        A dictionary of Forman Ricci Curvature values for simplices up to dimension p.
        If p = 0 or p = 1 is entered, it will return the forman ricci curvature for a graph network.
    """

    forman_dict = defaultdict(dict)
    if p < 1 :
        p = 1
    simplices = list(S)
    for i in range(len(simplices)) :
        l = len(simplices[i])
        if l <= p+1 and l >= 2 :
            forman_dict[l-1][simplices[i]] = func_compute_forman(S, p, laplacian, simplices[i])
    target_simplices = list(n_faces(S, 0))
    cnt = np.zeros(len(target_simplices))
    vals = np.zeros(len(target_simplices))
    nbrs = list(n_faces(S, 1))
    for vertex in nbrs :
        for v in vertex :
            vals[v] += forman_dict[1][vertex]
            cnt[v] += 1
    for i in range(len(vals)) :
        if cnt[i] > 0 :
            forman_dict[0][(i,)] = vals[i]/cnt[i]
        else :
            forman_dict[0][(i,)] = 0
    
    return forman_dict


@xgi.nodestat_func
def forman_degree(SC, bunch):
    list_nodes = set([tuple([i]) for i in SC.nodes.keys()])
    list_edges = set(tuple(sorted(list(SC.edges.members()[i]))) for i in range(len(SC.edges.members())))
    S = set(list_nodes).union(list_edges)
    laplacian = compute_laplacian(S, 0)
    dict_forman = compute_forman(S, 0, laplacian)
    node_dict = dict_forman[0]
    new_node_dict = {list(node_dict.keys())[i][0]: node_dict[list(node_dict.keys())[i]] for i in range(len(node_dict.keys()))}
    
    return new_node_dict


@xgi.edgestat_func
def forman_hyperedge(SC, bunch, p):
    list_nodes = set([tuple([i]) for i in SC.nodes.keys()])
    list_edges = set(tuple(sorted(list(SC.edges.members()[i]))) for i in range(len(SC.edges.members())))
    S = set(list_nodes).union(list_edges)
    laplacian = compute_laplacian(S, p)
    dict_forman = compute_forman(S, p, laplacian)
    hyperedge_dict = dict_forman[p]
    dict_to_key = {i:SC.edges.members()[i] for i in SC.edges.filterby('order',p)}
    dict_correspondance = {SC.edges.members()[i]:tuple(sorted(list(SC.edges.members()[i]))) for i in SC.edges.filterby('order',p)}
    new_hyperedge_dict = {k:hyperedge_dict[dict_correspondance[v]] for k,v in dict_to_key.items()}
    
    return new_hyperedge_dict
