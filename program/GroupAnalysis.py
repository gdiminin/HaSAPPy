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
import re
from HaSAPPY_time import *
############################################################


class Analysis:
    
    def __init__(self,classifier,date_today):
        self.name = classifier + date_today
        self.repicates = 0
        self.Reads = pd.DataFrame()
        self.II = pd.DataFrame()
        self.KI = pd.DataFrame()
        self.Bias = pd.DataFrame()
        self.all = pd.DataFrame()
        self.Outlier = pd.DataFrame()
        
        
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
                elif parameter == 'Bias':
                    series.append(experiments[exp].bias)
                    
                
            fusion = pd.concat(series, axis = 1, keys= keys_name)#concatantaion of the different experiments
            
            if len(keys_name) > 1:
                
                fusion['%s_%s_mean'%(group_name,parameter)] = fusion.mean(axis = 1)
                fusion['%s_%s_stdev'%(group_name,parameter)] = fusion.std(axis = 1)      
            return fusion
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
            processing.bias = creating_DataFrame(group_name,experiments,'Bias')
        if Info.GroupAnalysis.Parameters.Reads:
            processing.reads = creating_DataFrame(group_name,experiments,'Reads')
        
        return processing
        ####
    ####
    def comparing_experiments(Info,categories,date_today):
        ####    
        def fold_ttest (Info,experiments,replicates,parameter):
            
            list_to_concat = []
            
            control = experiments[Info.GroupAnalysis.Reference.name]
            list_to_concat.append(control)
            control_replicates = replicates[Info.GroupAnalysis.Reference.name]
            
            for exp in Info.GroupAnalysis.Others.name:
                list_to_concat.append(experiments[exp])
                exp_replicates = replicates[exp]
                temporary_dataframe = pd.DataFrame()
                substitution = 1/float(control_replicates)
                if control_replicates > 1:
                    control_series = control['%s_%s_mean' %(Info.GroupAnalysis.Reference.name,parameter)]
                else:
                    control_series = control
                if exp_replicates > 1:
                    exp_series = experiments[exp]['%s_%s_mean' %(exp,parameter)]
                else:
                    exp_series = experiments[exp]
                                       
                temporary_dataframe = pd.concat([control_series,exp_series],axis = 1)
                temporary_dataframe.columns = ['control','exp']  
                grouped = temporary_dataframe.groupby(by= [temporary_dataframe.control == 0, temporary_dataframe.exp == 0])
                if (True,True)in grouped.groups:
                    for index in grouped.groups[(True,True)]:
                        temporary_dataframe.ix[index,'control'] = 1
                if (False,True)in grouped.groups:
                    for index in grouped.groups[(False,True)]:
                        temporary_dataframe.ix[index,'exp'] = substitution
                if (True,False)in grouped.groups:
                    for index in grouped.groups[(True,False)]:
                        temporary_dataframe.ix[index,'control'] = substitution

                temporary_dataframe['fold'] =  temporary_dataframe['exp'] / temporary_dataframe['control']

                
                fold_series = pd.DataFrame(temporary_dataframe['fold'])
                fold_series.columns = ["%s_%s_fold" %(exp,parameter)]
                list_to_concat.append(fold_series)
                
                if control_replicates >2 and exp_replicates >2:
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
                on_going_experiments[group] = categories[group].bias
            
            summary.Bias = fold_ttest(Info,on_going_experiments,on_going_replicates,'Bias') 
            
        if Info.GroupAnalysis.Parameters.Reads:    
            on_going_experiments = {}
            on_going_replicates = {}
            for group in categories:
                on_going_replicates[group] = categories[group].replicates
                on_going_experiments[group] = categories[group].reads
            summary.Reads = fold_ttest(Info,on_going_experiments,on_going_replicates,'Reads') 
        summary.replicates = on_going_replicates
        
        return summary
    ####
    
    ####    
    def outlier_analysis (Info,summary):
        import LOF
        ####
        def outlier_InputData (Info,summary,group):
            ####            
            def prepare_lists (mean_to_concat,fold_to_concat,summary,group,dataframe,parameter):
                
                if summary.replicates[group] == 1:
                    value = Info.GroupAnalysis.Others.experiments[0][0]
                    mean_to_concat.append(dataframe['%s_%s_%s'%(group,parameter,value)]-dataframe['%s_%s_%s'%(Info.GroupAnalysis.Reference.name,parameter,Info.GroupAnalysis.Reference.experiments[0])])       
                else:
                    mean_to_concat.append(dataframe['%s_%s_mean'%(group,parameter)]-dataframe['%s_%s_mean'%(Info.GroupAnalysis.Reference.name,parameter)])
                fold_to_concat.append(dataframe['%s_%s_fold'%(group,parameter)])
                    
                
                return mean_to_concat,fold_to_concat
            ####
                        
            mean_to_concat = []
            fold_to_concat = []
            II_selection = pd.DataFrame()
            if Info.GroupAnalysis.Outlier.Parameters.II:
                mean_to_concat,fold_to_concat = prepare_lists (mean_to_concat,fold_to_concat,summary,group,summary.II,'II')
                if summary.replicates[group] == 1:
                    value = Info.GroupAnalysis.Others.experiments[0][0]
                    name = '%s_%s_%s'%(group,'II',value)
                    II_selection[name] = summary.II[name]
                else:
                    name = '%s_%s_mean'%(group,'II')
                    II_selection[name] = summary.II[name]
                                              
            if Info.GroupAnalysis.Outlier.Parameters.KI:
                mean_to_concat,fold_to_concat = prepare_lists (mean_to_concat,fold_to_concat,summary,group,summary.KI,'KI')
            if Info.GroupAnalysis.Outlier.Parameters.Bias:
                mean_to_concat,fold_to_concat = prepare_lists (mean_to_concat,fold_to_concat,summary,group,summary.Bias,'Bias')
            if Info.GroupAnalysis.Outlier.Parameters.Reads:
               mean_to_concat,fold_to_concat = prepare_lists (mean_to_concat,fold_to_concat,summary,group,summary.Reads,'Reads')
                       
            outlier_mean = pd.concat(mean_to_concat, axis = 1)
            outlier_mean[outlier_mean < 0] = 0
            outlier_fold = pd.concat(fold_to_concat, axis = 1)
            if summary.replicates[group] == 1:
                value = Info.GroupAnalysis.Others.experiments[0][0]
            else:
                value = 'mean'
            
                   
            outlier_fold = outlier_fold[II_selection['%s_II_%s' % (group,value)] >5]
    
            outlier_mean = outlier_mean[II_selection['%s_II_%s' % (group,value)] >5]
            
            
            return outlier_mean, outlier_fold

            
        ####
        for group in Info.GroupAnalysis.Others.name:        
            outlier_mean,outlier_fold = outlier_InputData (Info,summary,group)
            
            value_of_trustability = 0.4
            
            
            trustability = LOF.main(outlier_mean,Info,len(outlier_mean.index))
            trustability = trustability * value_of_trustability
            trustability = pd.merge(outlier_mean,trustability,left_index=True,right_index=True)
            

            outliers =LOF.main(outlier_fold,Info,len(outlier_fold.index))
            outliers=pd.merge(outlier_fold,outliers,left_index=True,right_index=True)
            
            summary.Outlier['%s_Outliers'%group] = outliers ['Score']
            
            summary.Outlier['%s_Score'%group] = outliers ['Score']*trustability['Score']         


            
            return summary
            

                
            


    
        ####
        
    
    
    ####

    #### Running commands ####

    date_today = getDay()   
    
    #Printing statements    
    
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

    for string in strings:
        print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    
            
    
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
    string = '\tRunTime: %s' % computeRunTime(startTime, getCurrTime())
    print_save_analysis (string, Info.GroupAnalysis.storage_loc) 
    
    startTime = getCurrTime()    
    if Info.GroupAnalysis.Outlier.perform:
        string = '\nOutlier analysis\n\tStarted: %s' %startTime
        print_save_analysis (string, Info.GroupAnalysis.storage_loc)

        summary = outlier_analysis(Info,summary)
        
        
        string = '\tRunTime: %s' % computeRunTime(startTime,getCurrTime())
        print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    
    
    
    summary.all = pd.concat([ parameter for parameter in [summary.II,summary.KI,summary.Bias,summary.Reads,summary.Outlier] if not parameter.empty],axis = 1)
    
    
    
    
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
        with open (os.path.join(location[0],location[1] + '_info.txt'),'a') as write:
            string = '\t%s :\t %s' % (date_today,Info.GroupAnalysis.storage_loc)
            print >> write,string
    string = 'Writing in Input Info files their usage for this analysis'
    print string
            
    string = '***\tEND Perform Group Analysis\t***' 
    print_save_analysis (string, Info.GroupAnalysis.storage_loc)
    

    
        
    

    
    
    
    



####

        

    


