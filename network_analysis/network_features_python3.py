#!/usr/bin/python
import sys
import networkx as nx
import networkx.algorithms.approximation as nx_approx
import pandas as pd
import numpy as np
from numexpr3 import evaluate as ev
from scipy.spatial.distance import cdist, pdist, squareform
import scipy.integrate as integrate
import time
import math

# giulia's entropy functions
# giulia spatial entropy
def giulia_spatial_entropy(edgelist, nodelist):

    start = time.time()

    g = nx.read_edgelist(edgelist, delimiter=',')
    nodes = list(g.nodes())
    nodes.sort()

    PP = np.array(nx.to_numpy_matrix(g,nodelist=nodes))

    n = len(nodes)
    distance = np.zeros((n,n))

    ens_id_expr_df = pd.read_csv(nodelist,names=['ens_id','expr_val'])
    ens_id_expr_df = ens_id_expr_df.set_index('ens_id')
    ens_id_expr_map = ens_id_expr_df.to_dict()['expr_val']

    for row in range(n):
        for col in range(n):
            x = ens_id_expr_map[nodes[row]]
            y = ens_id_expr_map[nodes[col]]
            distance[row][col] = math.fabs(x-y)

    Nbin = int(math.sqrt(n)+1)

    #linear binning
    mi,ma = np.min(distance[distance>0]),np.max(distance)
    limiti = np.linspace(mi,ma,Nbin+1)
    limiti[-1]+=0.1*ma
    limiti[0]-=0.1*mi

    b = np.searchsorted(limiti,distance)-1

    massimo = np.max(b)+1
   
    # BC gives how many links fall in each bin
    BC = [ np.sum(b[PP>0]==i)/2 for i in range(massimo) ]

    N=PP.shape[0];
    
    connectivity = np.sum(PP,1);
    avg_conn = np.mean(connectivity)

    #Lagrangians
    z = connectivity/(math.sqrt(avg_conn*N))
    w = BC / (avg_conn*N)

    old_z = z
    old_w = w

    loops = 5  # CHANGE to 10000 again
    precision = 1E-5


    for idx in range(loops):
        bigW = w.take(b)

        for i in range(N):  bigW[i,i]=0.

        U = ev("bigW * z")
        UT = U.T    
        D = ev("(UT * z) + 1.")

        UD = ev("U/D")

        del D,U,UT

        for i in range(N):  UD[i,i]=0.

        z = connectivity / np.sum(UD,1)

        zt = z[:].reshape(N,1)
        D2 = ev("(z*zt) / ( bigW *(z*zt) + 1.)")

        B2 = np.array([np.sum(np.where(b==i,D2,0.)) for i in range(len(w))])/2.

        w = np.where( (BC!=0) & (B2!=0),BC/B2,0)
        rz= ev("abs(1.-z/old_z)")
        rw= ev("abs(1.-w/old_w)")
        rz[np.isinf(rz)]=0.
        rw[np.isinf(rw)]=0.
        rz[np.isnan(rz)]=0.
        rw[np.isnan(rw)]=0.

        if max(rz)<precision and max(rw)<precision:
            break

        old_z = z
        old_w = w

    bigW = w.take(b)
    for i in range(N):  bigW[i,i]=0.

    z2 = bigW * np.outer(z,z)
    P = z2 / ( z2 + 1)
    Q=1.-P
    S = -1*(np.sum(np.log(P**P) + np.log(Q**Q) ))/2.

    stop = time.time()

    output = dict()
    output['num_nodes'] = n
    output['num_edges'] = len(g.edges())
    output['giulia_spatial_entropy'] = S
    output['runtime'] = stop-start
    output['nodelist'] = nodelist
    output['edgelist'] = edgelist

    return output

# giulia config entropy
def giulia_config_entropy(edgelist):

    start = time.time()

    g = nx.read_edgelist(edgelist, delimiter=',')
    nodes = list(g.nodes())
    nodes.sort()

    PP = np.array(nx.to_numpy_matrix(g,nodelist=nodes))

    N=PP.shape[0];
    connectivity = np.sum(PP,1);
    avg_conn = np.mean(connectivity)

    #Lagrangians
    z = connectivity/(math.sqrt(avg_conn*N))
    old_z = z
    loops = 10000
    precision = 1

    for idx in range(loops):
        zT = z[:,np.newaxis]
        D = ev("(zT * z) + 1.")
        UD = ev("z/D")
        del D
        for i in range(N):  UD[i,i]=0.
        z = connectivity / np.sum(UD,1)
        rz= ev("abs(1.-z/old_z)")
        rz[np.isinf(rz)]=0.
        rz[np.isnan(rz)]=0.

        if max(rz)<precision:
            break
        old_z = z
    z2 = np.outer(z,z)
    for i in range(N):  z2[i,i]=0.
    P = z2 / ( z2 + 1)
    Q=1.-P
    S = -1*(np.sum(np.log(P**P) + np.log(Q**Q) ))/2.
    stop = time.time()


    output = dict()
    output['num_nodes'] = n
    output['num_edges'] = len(g.edges())
    output['giulia_config_entropy'] = S
    output['runtime'] = stop-start
    output['edgelist'] = edgelist

    return output


# marinka resilience function (and auxiliary functions)
# trapezoidal approximation
def num_integration(func, lb, ub, num_points):
    
    auc = 0.0

    x_coord = np.linspace(lb,ub,num=num_points)
    y_coord = [func(x) for x in x_coord]

    step_size = x_coord[1]-x_coord[0]

    for idx in range(len(y_coord)-1):
        auc += (y_coord[idx] + y_coord[idx+1]) / 2.0

    auc *= step_size
   
    return auc

# marinka's modified shannon
def modified_shannon(g, num_rounds, fail_rate):
    
    entropy_trials = [0]*num_rounds

    for trial in range(num_rounds): 

        graph = g.copy()

        nodes = graph.nodes
        N = len(nodes)

        num_remove = min(int(fail_rate*N),len(nodes))

        failed_nodes = np.random.choice(nodes,num_remove,replace=False)

        for v in failed_nodes:
            graph.remove_edges_from(graph.edges(v))

        entropy = 0

        for cc in nx.connected_components(graph):
            C_i = len(list(cc))
            p_i = float(C_i) / N
            entropy += p_i * math.log(p_i)

        entropy *= (-1. / math.log(N))

        entropy_trials[trial] = entropy

    return np.mean(entropy_trials)

# marinka resilience
def marinka_resilience(edgelist):
    
    start = time.time()

    g = nx.read_edgelist(edgelist, delimiter=',')
    auc_entropy = num_integration(lambda x: modifiedShannon(g,1,x), 0, 1, 50)
    resilience = 1 - auc_entropy

    stop = time.degault_timer()

    output = dict()
    output['num_nodes'] = len(g.nodes())
    output['num_edges'] = len(g.edges())
    output['marinka_resilience'] = resilience
    output['runtime'] = stop-start
    output['edgelist'] = edgelist

    return output

# summary: num_nodes, num_edges, 
#          average_node_degree,
#          num_connected_components
def get_summary(edgelist):

    start = time.time()
    g = nx.read_edgelist(edgelist, delimiter=',')

    output = dict()

    output['num_nodes'] = len(g.nodes())
    output['num_edges'] = len(g.edges())
    output['average_node_degree'] = np.mean([x[1] for x in list(nx.degree(g))])
    output['num_connected_components'] = nx.number_connected_components(g)

    stop = time.time()
    output['runtime'] = stop-start
    output['edgelist'] = edgelist
    
    return output

    

