# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 08:07:37 2016

@author: GDM
"""

#####                  Importing modules               #####
import cPickle as pickle
import pandas as pd
import numpy as np
import os
from scipy.stats  import ttest_ind
from scipy.stats  import fisher_exact
import HaSAPPy.rankAnalysis as rankAnalysis

import re
from HaSAPPy.HaSAPPY_time import *
############################################################

class Analysis:
    
    def __init__(self,classifier,date_today):
        self.name = classifier + date_today
        self.repicates = 0
        self.Reads = pd.DataFrame()
        self.II = pd.DataFrame()
        self.KI = pd.DataFrame()
        self.Bias = pd.DataFrame()
        self.biasFW = pd.DataFrame()
        self.biasRV = pd.DataFrame()
        self.all = pd.DataFrame()
        self.Outlier = pd.DataFrame()
        self.Rank = pd.DataFrame()
        
        
def performing_analysis(Info):
    ####
    def print_save_analysis (string, storage_loc):
        print string
        with open (os.path.join(storage_loc,'analysis_info.txt'), 'a' ) as write:
            print >> write, string
    #### 
    
    """The function recover data from the libraries provided as input and generates for each block a Analysis Class containing different DataFrame for the parameter analysed"""
    ####    
    def create_analysis(Info,group_name,group_experiments,date_today):#group are the differnt libraies in Info.group
        ####
    
        def creating_DataFrame(group_name,experiments,parameter):    
            """Concatenation of the different replicates
            group_name = the name of the original group = Info.GroupAnalysis.Group.name
            experiments_name = list of group_name + replicates name
            parameter= list of Series for a particular parameter
            """
            series = []
            keys_name = []
            for exp in experiments:
                keys_name.append('%s_%s_%s'%(group_name,parameter,exp))
                if parameter == 'II':
                    series.append(experiments[exp].II)
                elif parameter == 'KI':
                    series.append(experiments[exp].KI)
                elif parameter == 'Reads':
                    series.append(experiments[exp].reads)
                elif parameter == 'biasFW':
                    series.append(experiments[exp].bias_FW)
                elif parameter == 'biasRV':
                    series.append(experiments[exp].bias_RV)
                    
            
            #Calculation of sum and mean and stdv (if replicates >1)
            fusion = pd.concat(series, axis = 1, keys= keys_name)#concatantaion of the different experiments
            ongoing = fusion.copy()
            fusion['%s_%s_sum'%(group_name,parameter)] = fusion.sum(axis = 1)
            
            if len(keys_name) > 1:
                fusion['%s_%s_mean'%(group_name,parameter)] = ongoing.mean(axis = 1)
                fusion['%s_%s_stdev'%(group_name,parameter)] = ongoing.std(axis = 1)
                
                
            return fusion #DataFrame containing for a parameter for a group all the replicates, sum, mean and stdev
        ####
            
        processing = Analysis(group_name,date_today) #Generation of Analysis class for each group
        experiments = {}
        replicates_number = 0
        for replicate in group_experiments:#each replicates
            replicates_number +=1
            replicate_input = Info.GroupAnalysis.input_files[Info.GroupAnalysis.lib_names.index(replicate)]#get the input file for the replicate of the group
            with open (replicate_input,'rb') as loading: #recover .pkl file with general informations of the library
                raw_data = pickle.load(loading)
            experiments[replicate] = raw_data #Upload in experiments (connected with library name) general address used in creating_DataFrame

        processing.replicates = replicates_number      
        
        if Info.GroupAnalysis.Parameters.II:
            processing.II = creating_DataFrame(group_name,experiments,'II')#store in a DataFrame concatenation experiment
        if Info.GroupAnalysis.Parameters.KI:
            processing.KI = creating_DataFrame(group_name,experiments,'KI')
        if Info.GroupAnalysis.Parameters.Bias:
            processing.biasFW = creating_DataFrame(group_name,experiments,'biasFW')
            processing.biasRV = creating_DataFrame(group_name,experiments,'biasRV')
        if Info.GroupAnalysis.Parameters.Reads:
            processing.Reads = creating_DataFrame(group_name,experiments,'Reads')
        
        return processing
        ####
    ####
    def comparing_experiments(Info,categories,date_today):
        ####    
        def fold_ttest (Info,experiments,replicates,parameter):
            
            list_to_concat = []
            
            control = experiments[Info.GroupAnalysis.Reference.name]
            list_to_concat.append(control)

            for exp in Info.GroupAnalysis.Others.name:
                list_to_concat.append(experiments[exp])
                control_series = control['%s_%s_sum' %(Info.GroupAnalysis.Reference.name,parameter)]
                exp_series = experiments[exp]['%s_%s_sum' %(exp,parameter)]                    
                temporary_dataframe = pd.concat([control_series,exp_series],axis = 1)
                temporary_dataframe.columns = ['control','exp']

                temporary_dataframe['control'][temporary_dataframe['control']==0] = 1

                temporary_dataframe['fold'] =  temporary_dataframe['exp'] / temporary_dataframe['control']
          
                fold_series = pd.DataFrame(temporary_dataframe['fold'])
                fold_series.columns = ["%s_%s_fold" %(exp,parameter)]
                list_to_concat.append(fold_series)

                #Calculate pvalue
                if replicates[Info.GroupAnalysis.Reference.name] >2 and replicates[exp] >2:
                    replicates_control = []
                    replicates_exp = []
                    for column in control.columns:
                        if column.find('_stdev') == -1 and column.find('_mean') == -1:
                            replicates_control.append(control[column])
                    for column in experiments[exp].columns:
                        if column.find('_stdev') == -1 and column.find('_mean') == -1:
                            replicates_exp.append(experiments[exp][column])
                    x,p = ttest_ind (replicates_control, replicates_exp)
                    ttest_series = pd.DataFrame(p,index = control.index,columns =["%s_%s_ttest" %(exp,parameter)])
                    list_to_concat.append(ttest_series)
                       
            return pd.concat(list_to_concat,axis =1)                
        ####
            
        summary = Analysis('Summary',date_today)
        if Info.GroupAnalysis.Parameters.II:
            on_going_experiments = {}
            on_going_replicates = {}
            for group in categories:
                on_going_replicates[group] = categories[group].replicates
                on_going_experiments[group] = categories[group].II
            summary.II = fold_ttest(Info,on_going_experiments,on_going_replicates,'II')
            
        if Info.GroupAnalysis.Parameters.KI:    
            on_going_experiments = {}
            on_going_replicates = {}
            for group in categories:
                on_going_replicates[group] = categories[group].replicates
                on_going_experiments[group] = categories[group].KI
            summary.KI = fold_ttest(Info,on_going_experiments,on_going_replicates,'KI')
            
        if Info.GroupAnalysis.Parameters.Bias:
            on_going_experiments = {}
            on_going_replicates = {}
            for group in categories:
                on_going_replicates[group] = categories[group].replicates
                on_going_experiments[group] = categories[group].biasFW
            summary.biasFW = fold_ttest(Info,on_going_experiments,on_going_replicates,'biasFW')
            on_going_experiments = {}
            on_going_replicates = {}
            for group in categories:
                on_going_replicates[group] = categories[group].replicates
                on_going_experiments[group] = categories[group].biasRV
            summary.biasRV = fold_ttest(Info,on_going_experiments,on_going_replicates,'biasRV')
        
        if Info.GroupAnalysis.Parameters.Reads:    
            on_going_experiments = {}
            on_going_replicates = {}
            for group in categories:
                on_going_replicates[group] = categories[group].replicates
                on_going_experiments[group] = categories[group].Reads
            summary.Reads = fold_ttest(Info,on_going_experiments,on_going_replicates,'Reads') 
        summary.replicates = on_going_replicates
        
        return summary
    
    ####
    def define_Bias(DATA,Info):
        ####
        def analyase_Bias_categories(DATA,Bias,Info,group_name, group_experiments):
            for exp in group_experiments:
                ongoing = pd.DataFrame([DATA.biasFW['%s_biasFW_%s'%(group_name,exp)],DATA.biasRV['%s_biasRV_%s'%(group_name,exp)]])
                Bias['%s_Bias_%s'%(group_name,exp)] = ongoing.apply(lambda r: '%d/%d'%(r['%s_biasFW_%s'%(group_name,exp)],r['%s_biasRV_%s'%(group_name,exp)]))
            if len(group_experiments) > 1:
                for factor in ['sum','mean']:
                    ongoing = pd.DataFrame([DATA.biasFW['%s_biasFW_%s'% (group_name,factor)],DATA.biasRV['%s_biasRV_%s'% (group_name,factor)]])
                    Bias['%s_Bias_%s'%(group_name,factor)] = ongoing.apply(lambda r: '%d/%d'%(r['%s_biasFW_%s'%(group_name,factor)],r['%s_biasRV_%s'%(group_name,factor)]))
            else:
                for factor in ['sum']:
                    ongoing = pd.DataFrame([DATA.biasFW['%s_biasFW_%s'% (group_name,factor)],DATA.biasRV['%s_biasRV_%s'% (group_name,factor)]])
                    Bias['%s_Bias_%s'%(group_name,factor)] = ongoing.apply(lambda r: '%d/%d'%(r['%s_biasFW_%s'%(group_name,factor)],r['%s_biasRV_%s'%(group_name,factor)]))
            return Bias
        ####
        def calculate_Bias_fold (DATA,reference,group):
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
#            ongoing = pd.DataFrame([DATA.biasFW['%s_biasFW_sum'%group],DATA.biasRV['%s_biasRV_sum'%group],DATA.biasFW['%s_biasFW_sum'%reference],DATA.biasRV['%s_biasRV_sum'%reference]])
            reference = pd.DataFrame({'FW':DATA.biasFW['%s_biasFW_sum'%reference],'RV':DATA.biasRV['%s_biasRV_sum'%reference]})

            reference['reference_fold'] = reference.apply(adjust_Bias_reference,axis = 1)
            reference['reference_fold']= reference['reference_fold'].replace (0,1)
            
            group = pd.DataFrame({'FW':DATA.biasFW['%s_biasFW_sum'%group],'RV':DATA.biasRV['%s_biasRV_sum'%group]})
            
            group['group_fold'] = group.apply(adjust_Bias_selected,axis = 1)
            return group['group_fold']/reference['reference_fold']
        
        ####
                                                                           
        Bias = pd.DataFrame()
        
        Bias = analyase_Bias_categories(DATA,Bias,Info,Info.GroupAnalysis.Reference.name, Info.GroupAnalysis.Reference.experiments)
        for group in Info.GroupAnalysis.Others.name:
            Bias = analyase_Bias_categories(DATA,Bias,Info,group, Info.GroupAnalysis.Others.experiments[Info.GroupAnalysis.Others.name.index(group)])
            Bias['%s_Bias_fold'%group] = calculate_Bias_fold(DATA,Info.GroupAnalysis.Reference.name,group)
        return Bias

    #### Running commands ####

    date_today = getDay()   
    
    ###Printing statements
    
    strings = []    
    strings.append('\n***\tPerform Group Analysis\t***\t\tDate: %s' % date_today)
    
    string = '\n\t{:25s}\t'.format('Reference group:') +  Info.GroupAnalysis.Reference.name
    strings.append(string)  
    string = '\t\t{:25s}\t'.format('Numbers of replicates:') + str(len(Info.GroupAnalysis.Reference.experiments))
    strings.append(string) 

    for replicate in Info.GroupAnalysis.Reference.experiments:
        string = '\t\t- %s  (%s)'%(replicate, Info.GroupAnalysis.input_files[Info.GroupAnalysis.lib_names.index(replicate)])
        strings.append(string)
    string = '\n\t{:25s}\t'.format('Analysed groups number:') +  str(len(Info.GroupAnalysis.Others.name))
    strings.append(string)
    for exp in Info.GroupAnalysis.Others.name:
        pos = Info.GroupAnalysis.Others.name.index(exp)
        string = '\t\t{:20s}\t'.format(str(pos+1)+')'+ exp) 
        strings.append(string)
        string  = '\t\t\t{:20s}\t'.format('Numbers of replicates:') + str(len(Info.GroupAnalysis.Others.experiments[pos]))
        strings.append(string)
        for  replicate in Info.GroupAnalysis.Others.experiments[pos]:
            string = '\t\t\t- %s  (%s)'%(replicate, Info.GroupAnalysis.input_files[Info.GroupAnalysis.lib_names.index(replicate)])
            strings.append(string)
        
    string = '\n\t{:25s}\t'.format('Analysis Parameters')
    strings.append(string)
    string = '\t\t{:8s}:\t'.format('II') + '%s' % Info.GroupAnalysis.Parameters.II
    strings.append(string)
    string = '\t\t{:8s}:\t'.format('KI') + '%s' % Info.GroupAnalysis.Parameters.KI
    strings.append(string)
    string = '\t\t{:8s}:\t'.format('Bias') + '%s' % Info.GroupAnalysis.Parameters.Bias
    strings.append(string)
    string = '\t\t{:8s}:\t'.format('Reads') + '%s' % Info.GroupAnalysis.Parameters.Reads
    strings.append(string)
    
    string = '\t\t{:8s}:\t'.format('Outlier') + '%s' % Info.GroupAnalysis.Outlier.perform
    strings.append(string)
    if Info.GroupAnalysis.Outlier.perform:
        string = '\t\t\t{:8s}:\t'.format('II') + '%s' % Info.GroupAnalysis.Outlier.Parameters.II
        strings.append(string)
        string = '\t\t\t{:8s}:\t'.format('KI') + '%s' % Info.GroupAnalysis.Outlier.Parameters.KI
        strings.append(string)
        string = '\t\t\t{:8s}:\t'.format('Bias') + '%s' % Info.GroupAnalysis.Outlier.Parameters.Bias
        strings.append(string)
        string = '\t\t\t{:8s}:\t'.format('Reads') + '%s' % Info.GroupAnalysis.Outlier.Parameters.Reads
        strings.append(string)
        string = '\t\t\tFidelity correction: %s' % str(Info.GroupAnalysis.Outlier.fidelity)
        strings.append(string)

    for string in strings:
        print_save_analysis (string, Info.GroupAnalysis.storage_loc)

    ###


    ### Run Commands:
    
    categories = {}
    
    startTime = getCurrTime()      
    string = '\nGenearte group analysis\n\tStarted: %s' % startTime
    print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    categories[Info.GroupAnalysis.Reference.name] = create_analysis (Info,Info.GroupAnalysis.Reference.name,Info.GroupAnalysis.Reference.experiments,date_today)
    
    for group in Info.GroupAnalysis.Others.name:
        categories[group] = create_analysis (Info,group,Info.GroupAnalysis.Others.experiments[Info.GroupAnalysis.Others.name.index(group)],date_today)
    
    string = '\tRunTime: %s' % computeRunTime(startTime, getCurrTime())
    print_save_analysis (string, Info.GroupAnalysis.storage_loc)    
    

    
    
    startTime = getCurrTime()
    string = '\nStatistycal analysis of the groups\n\tStarted: %s' % startTime
    print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    summary = comparing_experiments(Info,categories,date_today)

    if Info.GroupAnalysis.Parameters.Bias:
        summary.Bias = define_Bias(summary,Info)

    summary = rankAnalysis.calculateRank(Info.GroupAnalysis,summary) #calculate rank values

    if Info.GroupAnalysis.Parameters.KI:
        import HaSAPPy.Fisher as Fisher
        summary.KI = Fisher.main(Info.GroupAnalysis,summary.KI)
        string = '\t-Fisher\'s exact test on KI insertion'
        print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    if Info.GroupAnalysis.Outlier.perform:
	if Info.GroupAnalysis.Outlier.Approach.fold:
            import HaSAPPy.Outlier_fold as Outlier_fold
	    summary = Outlier_fold.main(Info.GroupAnalysis,summary)
            string = '\t-Outlier_fold analysis'
            print_save_analysis (string, Info.GroupAnalysis.storage_loc)
	if Info.GroupAnalysis.Outlier.Approach.rank:
            summary = rankAnalysis.calculateOutlierRank(Info.GroupAnalysis,summary)
            string = '\t-Outlier_rank analysis'
            print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    string = '\tRunTime: %s' % computeRunTime(startTime, getCurrTime())
    print_save_analysis (string, Info.GroupAnalysis.storage_loc)

    summary.all = pd.concat([ parameter for parameter in [summary.II,summary.KI,summary.Bias,summary.biasFW,summary.biasRV,summary.Reads,summary.Outlier,summary.Rank] if not parameter.empty],axis = 1)

    string = '\nColumns: \n\t%s' %(' | ').join(summary.all.columns)
    print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    
    with open (os.path.join(Info.GroupAnalysis.storage_loc,'raw', 'RawData.pkl'),'wb') as saving:
        pickle.dump(summary,saving)
        
    string = '\nSaved RawData analysis in : %s' %(os.path.join(Info.GroupAnalysis.storage_loc,'raw', 'RawData.pkl'))
    print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    
    with open (os.path.join(Info.GroupAnalysis.storage_loc,'raw', 'GroupAnalysis.pkl'),'wb') as saving:
        pickle.dump(Info.GroupAnalysis,saving)
    
    string = 'Saved GroupAnalysis file (necessary for table generation) in : %s' %(os.path.join(Info.GroupAnalysis.storage_loc,'raw', 'GroupAnalysis.pkl'))
    print_save_analysis (string, Info.GroupAnalysis.storage_loc)       
    
    for exp in Info.GroupAnalysis.input_files:
        location = re.findall('^(.+)/raw/(.+)_GenesData.pkl',exp)[0]
	if os.path.isfile(os.path.join(location[0],location[1] + '_info.txt')):
            with open (os.path.join(location[0],location[1] + '_info.txt'),'a') as write:
                string = '\t%s :\t %s' % (date_today,Info.GroupAnalysis.storage_loc)
                print >> write,string
    string = 'Writing in Input Info files their usage for this analysis'
    print string
            
    string = '***\tEND Perform Group Analysis\t***' 
    print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    

    
        
    

    
    
    
    



####

        

    


