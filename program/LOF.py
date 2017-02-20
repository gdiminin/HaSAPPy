import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors




#knn function gets the dataset and calculates K-Nearest neighbors and distances
def knn(df,k):
    nbrs = NearestNeighbors(n_neighbors=1)
    nbrs.fit(df)
    distances, indices = nbrs.kneighbors(df)
    return distances, indices

#reachDist calculates the reach distance of each point to MinPts around it
def reachDist(df,MinPts,knnDist):
    nbrs = NearestNeighbors(n_neighbors=MinPts)
    nbrs.fit(df)
    distancesMinPts, indicesMinPts = nbrs.kneighbors(df)
    distancesMinPts[:,0] = np.amax(distancesMinPts,axis=1)
    distancesMinPts[:,1] = np.amax(distancesMinPts,axis=1)
    distancesMinPts[:,2] = np.amax(distancesMinPts,axis=1)
    return distancesMinPts, indicesMinPts

#lrd calculates the Local Reachability Density
def lrd(MinPts,knnDistMinPts):
    return (MinPts/np.sum(knnDistMinPts,axis=1))

#Finally lof calculates lot outlier scores
def lof(Ird,MinPts,dsts):
    lof=[]
    for item in dsts:
       tempIrd = np.divide(Ird[item[1:]],Ird[item[0]])
       lof.append(tempIrd.sum()/MinPts)
    return lof

def main (data,Info,m):
    #m = MinPoints
    data = (data-data.mean())/data.std()
    knndist, knnindices = knn(data,1)
    reachdist, reachindices = reachDist(data,m,knndist)
    irdMatrix = lrd(m,reachdist)
    lofScores = lof(irdMatrix,m,reachindices) 
    scores= pd.DataFrame(lofScores,columns=['Score'],index = data.index)
    #mergedData=pd.merge(data,scores,left_index=True,right_index=True)
    #mergedData['flag'] = mergedData.apply(returnFlag,axis=1)
    return scores
    
     
    
    

    


