#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 08:49:20 2017

@author: gdiminin
"""

import pandas as pd
import numpy as np
import LOF


def main (GroupAnalysis,DATA):
    ####
    def outlier_InputData (GroupAnalysis,DATA,group):
            ####            
        def prepare_lists (fold_to_concat,group,dataframe,parameter):
                
            temporary_series = dataframe['%s_%s_fold'%(group,parameter)]
            temporary_series.name = '%s_%s_outlier'%(group,parameter)
            fold_to_concat.append(temporary_series)
                
            return fold_to_concat

            ####
                        
        fold_to_concat = []

        if GroupAnalysis.Outlier.Parameters.II:
            fold_to_concat = prepare_lists(fold_to_concat,group,DATA.II,'II')                                             
        if GroupAnalysis.Outlier.Parameters.KI:
            fold_to_concat = prepare_lists (fold_to_concat,group,DATA.KI,'KI')
        if GroupAnalysis.Outlier.Parameters.Bias:
            fold_to_concat = prepare_lists (fold_to_concat,group,DATA.Bias,'Bias')
        if GroupAnalysis.Outlier.Parameters.Reads:
            fold_to_concat = prepare_lists (fold_to_concat,group,DATA.Reads,'Reads')
       
        outlier_fold = pd.concat(fold_to_concat, axis = 1)

        if len(group)> 4:
            outlier_fold = outlier_fold[DATA.II['%s_II_sum' % group] > 4] #!!!Selection of genes with at least 3 II in the selected libraries
        else:
            outlier_fold = outlier_fold[DATA.II['%s_II_sum' % group] > len(group)] 

    
        return outlier_fold

        ####
    
    for group in GroupAnalysis.Others.name:
        outlier_fold = outlier_InputData(GroupAnalysis,DATA,group)

        
    
        if GroupAnalysis.Outlier.fidelity != 0:
            IIadjust = DATA.II['%s_II_sum'%group].copy()
            if len(group)> 4:
                IIadjust = IIadjust[DATA.II['%s_II_sum' % group] > 4] #!!!Selection of genes with at least 3 II in the selected libraries
            else:
                IIadjust = IIadjust[DATA.II['%s_II_sum' % group] > len(group)] 
             
            if GroupAnalysis.Outlier.fidelity == 'n.d.':
                slop_value = DATA.II['%s_II_sum'%group][DATA.II['%s_II_sum'%group]>0].mean()
            else:
                slop_value = GroupAnalysis.Outlier.fidelity
               

            IIadjust = IIadjust.apply(lambda r: 1/(1+np.exp(-(r-slop_value)/10)))
            
            outlier_fold = outlier_fold.mul(IIadjust,axis = 0)

        outliers =LOF.main(outlier_fold,len(outlier_fold.index))
    
        outliers =pd.merge(outlier_fold,outliers,left_index=True,right_index=True)
        outliers = outliers.rename(columns ={'Score':'%s_Score_fold' % group})            
    	DATA.Outlier = pd.concat([DATA.Outlier,outliers],axis = 1)

    return DATA
    
    
    
    
