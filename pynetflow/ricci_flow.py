import networkx as nx
import pandas as pd
import numpy as np


def drop_weights(G) :
    '''Drop the weights from a networkx weighted graph.'''
    for node, edges in nx.to_dict_of_dicts(G).items():
        for edge, attrs in edges.items():
            attrs.pop('weight', None)
            

def dist_mat(vec) :
    dist = np.zeros((len(vec),len(vec)))
    for i in range(len(vec)) :
        for j in range(len(vec)) :
            if vec[i] != 0 and vec[j] != 0 :
                dist[i,j] = 1/(vec[i]*vec[j])
            else :
                dist[i,j] = 0
    return dist


def RicciFlow_normalised(G, vec, vecf, Niter = 10, eta = 0.01) :   
    Iteration = list()
    Pair = list()
    Distance = list()
    F = list()
    nv = len(G.nodes)
    Dup = np.zeros((nv,nv))
    D = 1/np.outer(vec, vec)
    fr = Forman_ricci(G, D)["ricci"]  
    Df = dist_mat(vecf)
    Frf = Forman_ricci(G, Df)["ricci"]  
    rows = np.where(fr != 0)[0]
    cols = np.where(fr != 0)[1]
    pairs = [(rows[i], cols[i]) for i in range(len(rows))]
    print(Total_ricci_curvature(G, D))
    for it in range(Niter): 
        print (it, end = ",")
        for pair in pairs:
            a = pair[0]
            b = pair[1]
            Dup[a,b] = D[a,b] + eta*(fr[a,b]-Frf[a,b])*D[a,b]
            Iteration.append(it)
            Pair.append((a,b))
            Distance.append(Dup[a,b])
            F.append(fr[a,b])
            if D[a,b] <0:
                print("end process")      
        D = Dup
        D[D == np.inf] = 0
        fr = Forman_ricci(G, D)["ricci"]  
        print(Total_ricci_curvature(G, D))
    df = pd.DataFrame(list(zip(Iteration, Pair, Distance, F)))
    output = {"D": D, "FR": fr, "df": df}
    return output


def RicciFlow(G, vec, Niter = 10, eta = 0.01) :   
    Iteration = list()
    Pair = list()
    Distance = list()
    F = list()
    nv = len(G.nodes)
    Dup = np.zeros((nv,nv))
    D = 1/np.outer(vec, vec)
    fr = Forman_ricci(G, D)["ricci"]  
    rows = np.where(fr != 0)[0]
    cols = np.where(fr != 0)[1]
    pairs = [(rows[i], cols[i]) for i in range(len(rows))]
    print(Total_ricci_curvature(G, D))
    for it in range(Niter): 
        print (it, end = ",")
        for pair in pairs:
            a = pair[0]
            b = pair[1]
            Dup[a,b] = D[a,b] - eta*(fr[a,b])*D[a,b]
            Iteration.append(it)
            Pair.append((a,b))
            Distance.append(Dup[a,b])
            F.append(fr[a,b])
            if D[a,b] <0:
                print("end process")      
        D = Dup
        D[D == np.inf] = 0
        fr = Forman_ricci(G, D)["ricci"]  
        print(Total_ricci_curvature(G, D))
    df = pd.DataFrame(list(zip(Iteration, Pair, Distance, F)))
    output = {"D": D, "FR": fr, "df": df}
    return output


def Forman_ricci(G, dist, node_weight = None) :        
    nv = sorted(list(G.nodes()))
    Adj = nx.adjacency_matrix(G, nodelist = nv)
    nv = np.array(nv)-1
    degree = [np.sum(Adj.todense()[:,i]) for i in range(len(G))]
    ricci_mat = np.zeros((len(nv),len(nv)))
    scalar = [0 for i in nv]
    for c1 in nv :  
        for c2 in nv :
            if (Adj[c1,c2] != 0):
                sum1 = np.array([(dist[c1,i])**(-1/2) for i in list(set(nv) - {c2}) if Adj[c1,i] != 0])
                sum2 = np.array([(dist[i,c2])**(-1/2) for i in list(set(nv) - {c1}) if Adj[i,c2] != 0])
                if node_weight is not None:
                    we1 = (1/node_weight[c1])*(1 - (dist[c1,c2]**(1/2))*np.sum(sum1))
                    we2 = (1/node_weight[c2])*(1 - (dist[c1,c2]**(1/2))*np.sum(sum2))
                else:
                    we1 = (1/degree[c1])*(1 - (dist[c1,c2]**(1/2))*np.sum(sum1))
                    we2 = (1/degree[c2])*(1 - (dist[c1,c2]**(1/2))*np.sum(sum2))
                ricci_mat[c1,c2] =  we1 + we2   
        scalar[c1] = (np.sum(ricci_mat[c1,:]))/degree[c1]
    output = {"ricci": ricci_mat ,"scalar": scalar}
    return output


def Signal_entropy(G, dist) :
    nv = sorted(list(G.nodes()))
    Adj = nx.adjacency_matrix(G, nodelist = nv)
    nv = np.array(nv)-1
    dist = np.reciprocal(dist)*Adj
    summ = np.sum(dist, axis = 0)
    mu = summ/np.sum(summ)
    stoch = (dist@np.diag(1/summ)).T
    ent = [0 for i in nv]
    for c1 in nv :  
        temp_Si = 0
        for c2 in nv :
            temp_proba = stoch[c1,c2]
            if temp_proba != 0:
                temp_Si -= temp_proba*np.log(temp_proba)
        ent[c1] = temp_Si
    SR = np.sum(mu*ent)
    return SR


def Total_ricci_curvature(G, dist):
    nv = sorted(list(G.nodes()))
    Adj = nx.adjacency_matrix(G, nodelist = nv)
    nv = np.array(nv)-1
    curv = Forman_ricci(G, dist)["scalar"]
    dist = np.reciprocal(dist)*Adj
    summ = np.sum(dist, axis = 0)
    mu = summ/np.sum(summ)
    Tot_curv = np.sum(mu*curv)
    return Tot_curv