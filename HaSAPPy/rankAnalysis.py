#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 18:08:54 2017

@author: gdiminin
"""

import pandas as pd
import numpy as np
import HaSAPPy.LOF as LOF
import os
import cPickle as pickle


###TEST
#
#with open ('/home/giulio/BIGDATA/Data/Analysis/WNT/Analysis/2017-11-16/raw/#GroupAnalysis.pkl','rb') as load:
#    GroupAnalysis = pickle.load(load)
#with open('/home/giulio/BIGDATA/Data/Analysis/WNT/Analysis/2017-11-16/raw/#RawData.pkl','rb') as load:
#    dataframe = pickle.load(load)
#dataframe.Rank = pd.DataFrame()
###



def calculateRank (GroupAnalysis,DATA):
    ####
    def calculate_ranking(GroupAnalysis,DATA,parameter,DataFrame):
        index_reference = '%s_%s'%(GroupAnalysis.Reference.name,parameter)
        DataFrame['%s_rank'%(index_reference)] = DATA['%s_sum'%(index_reference)].rank(method = 'max',ascending = False)   
        for group in GroupAnalysis.Others.name:
            index_others = '%s_%s'%(group,parameter)
            DataFrame['%s_rank'%(index_others)] = DATA['%s_sum'%(index_others)].rank(method = 'max',ascending = False)
            DataFrame['%s_%s'%(index_others,'rankedFold')] = DataFrame['%s_%s'%(index_reference,'rank')]/DataFrame['%s_%s'%(index_others,'rank')]
        
        return DataFrame
    ####   
    def adjust_Bias_reference(index):
        total_Bias  = index.FW + index.RV
        if total_Bias == 0:
            return 0
        if index.RV == 0:
            index.RV = 1
        ratio = float(index.FW)/index.RV
        if total_Bias < 15:
            if ratio <1:
                ratio = 1
        return ratio
    ####
    def adjust_Bias_selected(index):
        total_Bias  = index.FW + index.RV
        if total_Bias == 0:
            return 0
        elif total_Bias < 15:
            return 1
        else:
            if index.RV <1:
                index.RV = 1
            return float(index.FW) / index.RV
    ####
    
    rank_fold = pd.DataFrame()
    
    if GroupAnalysis.Parameters.II:
        DATA.Rank = calculate_ranking(GroupAnalysis,DATA.II,'II',DATA.Rank)
    if GroupAnalysis.Parameters.KI:
        DATA.Rank = calculate_ranking(GroupAnalysis,DATA.KI,'KI',DATA.Rank)
    if GroupAnalysis.Parameters.Reads:
        DATA.Rank = calculate_ranking(GroupAnalysis,DATA.KI,'Reads',DATA.Rank)
    
    if GroupAnalysis.Parameters.Bias:
            
        reference = pd.DataFrame({'FW':DATA.biasFW['%s_biasFW_sum'%GroupAnalysis.Reference.name],'RV':DATA.biasRV['%s_biasRV_sum'%GroupAnalysis.Reference.name]})
        reference['reference_fold'] = reference.apply(adjust_Bias_reference,axis = 1)
        reference['reference_fold']= reference['reference_fold'].replace(to_replace=0,value = 1)
        DATA.Rank['%s_Bias_rank'%GroupAnalysis.Reference.name] = reference['reference_fold'].rank(method = 'max',ascending = False)
        for group in GroupAnalysis.Others.name:
            selected = pd.DataFrame({'FW':DATA.biasFW['%s_biasFW_sum'%group],'RV':DATA.biasRV['%s_biasRV_sum'%group]})
            selected['selected_fold'] = selected.apply(adjust_Bias_selected,axis = 1)
            DATA.Rank['%s_Bias_rank'%group] = selected['selected_fold'].rank(method = 'max',ascending = False)
            DATA.Rank['%s_Bias_rankedFold'%group] = DATA.Rank['%s_Bias_rank'%GroupAnalysis.Reference.name]/DATA.Rank['%s_Bias_rank'%group]
    return DATA

def calculateOutlierRank(GroupAnalysis,DATA):
    rank_fold = pd.DataFrame()
    for name in DATA.Rank.columns:
        if '_rankedFold' in name:
            rank_fold[name] = DATA.Rank[name]
    for group in GroupAnalysis.Others.name: 
        rank_group = rank_fold[[name for name in rank_fold.columns if '%s_'%group in name]]
        if LOF.testDataForLOF(rank_group):
            outliers =LOF.main(rank_group,len(rank_group.index))
            DATA.Rank =pd.merge(DATA.Rank,outliers,left_index=True,right_index=True)
            DATA.Rank = DATA.Rank.rename(columns ={'Score':'%s_Score_rank' % group})               
    else:
            print '\n!!!Data provided is not compatible with LOF analysis. Calculation of %s_Score_rank will be skipped.\n' % group

    return DATA

