import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as mPatch
from matplotlib.legend_handler import HandlerLine2D
import itertools
import os



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

#We flag anything with outlier score greater than 1.2 as outlier#This is just for charting purposes
def returnFlag(x):
    if x['Score']>2:
       return 1
    else:
       return 0

#Below section creates the charts
def plot_outliers(mergedData,Info,to_plot):
    ####    
    def plot_3D(outliers,normals,Info,comb_columns):
        line1, = plt.plot([1], marker='o', label='Regular',linestyle='None',color='blue')
        line2, = plt.plot([1], marker='*', label='Outlier',linestyle='None',color='red')
    
        fig=plt.figure(dpi=80, facecolor='w', edgecolor='k')
        fig.legend((line2,line1),('Outliers','Regular'),loc=1,numpoints=1,ncol=2)
    
    
        ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
        ax1.scatter(outliers[comb_columns[0]],outliers[comb_columns[1]],outliers[comb_columns[2]],c='r',marker='*')
        ax1.scatter(normals[comb_columns[0]],normals[comb_columns[1]],normals[comb_columns[2]],c='b',marker='o')
        ax1.set_xlabel(comb_columns[0])
        ax1.set_ylabel(comb_columns[1])
        ax1.set_zlabel(comb_columns[2])
        ax1.set_title('Outliers Vs. Rest\n%s, %s, %s' %(comb_columns[0],comb_columns[1],comb_columns[2]))
    
        plt.tight_layout()
        plt.show()
        fig.savefig(os.path.join(Info.GroupAnalysis.storage_loc,'Outliers_Vs_Rest(%s-%s-%s).svg' %(comb_columns[0],comb_columns[1],comb_columns[2])))
        
    def plot_2D(outliers,normals,Info):
        pass
    ####
    
    if len(mergedData.columns)>2:
        for comb_columns in itertools.combinations(to_plot,3):
            plot_3D(mergedData[(mergedData['flag']==1)],mergedData[(mergedData['flag']==0)],Info,comb_columns)
        else:
            plot_2D([(mergedData['flag']==1)],mergedData[(mergedData['flag']==0)],Info)
            

    #Finally we draw the histogram of scores
    
    fig=plt.figure(dpi=80, facecolor='w', edgecolor='k')
    
    
    ax2 = plt.subplot2grid((1,1), (0,0))
    ax2.hist(mergedData['Score'],bins=100,facecolor='cornflowerblue')
    ax2.set_xlabel('LOF Score')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Outlier Scores')
    fig.savefig(os.path.join(Info.GroupAnalysis.storage_loc,'graph','OutliersScores.svg'))
    
    #plt.tight_layout()
    #plt.show()


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
    
     
    
    

    


