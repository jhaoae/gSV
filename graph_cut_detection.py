import numpy as np
import gco
import scipy.io as scio
import time
from scipy.sparse import coo_matrix

def readgetAdj2(sizeData):
    numSites = np.prod(sizeData)
    id1 = np.arange(0, numSites-1)
    id2 = np.arange(1, numSites)
    id3 = np.arange(0, numSites - 2)
    id4 = np.arange(2, numSites)
    value1 = np.ones(numSites - 1)
    value2 = np.ones(numSites - 2)
    W = coo_matrix((value1, (id1, id2)), shape=(numSites, numSites))
    W += coo_matrix((value2, (id3, id4)), shape=(numSites, numSites))
    return W

def graph_cut_detection(loss, depth_mean,beta_para,gamma,amplify_para):
    e0 = loss[0]

    beta = depth_mean / beta_para
    #gamma = 1
    amplify = np.ceil(depth_mean / amplify_para)

    # Create a graph
    hMRF = gco.GCO()
    hMRF.create_general_graph(len(e0),2)
   
    # Set smoothness costs
    smooth_cost = np.array([[0, 1],[1, 0]])
    hMRF.set_smooth_cost(smooth_cost)

    # Set neighbors
    AdjMatrix = readgetAdj2(len(e0))
    #print('AdjMartix',AdjMatrix)
    #print(AdjMatrix.tocoo().row,AdjMatrix.tocoo().col,amplify * AdjMatrix.data)
    hMRF.set_all_neighbors(AdjMatrix.tocoo().row,AdjMatrix.tocoo().col,amplify * AdjMatrix.data)

    # Set data costs
    #print(np.column_stack((e0,np.ones(len(e0))*beta)))
    data_cost = (1/gamma)*np.column_stack((e0,np.ones(len(e0))*beta)).astype(np.int32)
    hMRF.set_data_cost(data_cost)
    
    # expansion
    hMRF.expansion()


    # Get the cut results
    #print(hMRF.get_labels())
    label = hMRF.get_labels()
    S = np.reshape(np.where(label==1,1,0),len(e0))

    diff = np.diff(S.astype(int))
    diff1 = np.where(diff == 1)[0]
    diff2 = np.where(diff == -1)[0]
    rrr = list(zip(diff1, diff2))
    hMRF.destroy_graph()
    return rrr,S




