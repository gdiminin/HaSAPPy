# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 08:47:13 2016

@author: GDM
"""
#####                  Importing modules               #####
import HTSeq
import cPickle as pickle
import pandas as pd
import os
from HaSAPPY_time import *
############################################################



class Genes_library():
    """...explanation..."""
    def __init__(self,exp,input_): 
        self.name = exp
        self.input = input_
        self.reads = pd.Series()
        self.II = pd.Series()
        self.KI = pd.Series()
        self.bias_FW = pd.Series()
        self.bias_RV = pd.Series()
        self.bias =pd.Series()
        
        
def library_analysis(Info):
    ######      Internal Functions Definitions     ######
    def identify_genes_insertions (exp,Info,genes_ref,genes,genes_sense,exons,introns_sense,introns_antisense):
        #the function devides in two branches according if also Reads parameter should be analysed
        #if type_of_analysis['Reads']: #Perform this counting if are also present Reads informations
        startTime = getCurrTime()        
        string = 'Localize insertions in gene categories:\n\tStarted: %s' %startTime
        Info.print_save(exp,string)
        
        library = Genes_library(exp,Info.GeneDefinition.input_files[Info.GeneDefinition.lib_names.index(exp)])
                    
        with open (library.input,'rb') as loading:
            raw_data = pickle.load(loading)
            
        library.II = pd.Series(data = 0, index = genes_ref.index)        
        library.II['_no_feature'] = 0
                    
        library.KI = pd.Series(data = 0, index = genes_ref.index)
        library.KI['_no_feature'] = 0
            
        library.bias_FW = pd.Series(data = 0, index = genes_ref.index)
        library.bias_FW['_no_feature'] = 0
            
        library.bias_RV = pd.Series(data = 0, index = genes_ref.index)
        library.bias_RV['_no_feature'] = 0
            
        library.bias = pd.Series(data = 0, index = genes_ref.index)
        library.bias['_no_feature'] = 0
            
        library.reads = pd.Series(data = 0, index = genes_ref.index)
        library.reads['_no_feature'] = 0
    
        for i in raw_data.raw.keys():
            i_nosense = HTSeq.GenomicPosition(i.chrom,i.start,strand =".")
            if Info.GeneDefinition.Parameters.Reads: #Counting all Genes Insertions to be added to library.reads if paired-end
                gene_ids = set(genes[i_nosense])
                if gene_ids == set([]):
                    library.reads["_no_feature"] += raw_data.raw[i]
                else:
                    for gene_id in gene_ids: 
                        library.reads[gene_id] += raw_data.raw[i]
            if Info.GeneDefinition.Parameters.II:
                update_library(genes,i_nosense,library.II)
            if Info.GeneDefinition.Parameters.KI:    
                update_library(genes_sense,i,library.KI) #Counting all Sense Insertions to be added to library.KI
                update_library(exons,i,library.KI) #Counting all Exon Insertions NOT sense to be added to library.KI
            if Info.GeneDefinition.Parameters.Bias:
                update_library(introns_sense,i,library.bias_FW) #Counting all Intron Insertions sense to be added to library.bias_FW
                update_library(introns_antisense,i,library.bias_RV) #Counting all Intron Insertions NOT sense to be added to library.bias_RV
                    
       
        temporary_dataframe = pd.concat([library.bias_FW,library.bias_RV],axis = 1,keys = ['FW','RV'])
        grouped = temporary_dataframe.groupby(by= [temporary_dataframe.FW == 0, temporary_dataframe.RV == 0])
        if (True,True)in grouped.groups:
            for index in grouped.groups[(True,True)]:
                temporary_dataframe.ix[index,'RV'] = 1
        if (False,True)in grouped.groups:
            for index in grouped.groups[(False,True)]:
                temporary_dataframe.ix[index,'RV'] = 1
        if (True,False)in grouped.groups:
            for index in grouped.groups[(True,False)]:
                temporary_dataframe.ix[index,'FW'] = 1

        library.bias = temporary_dataframe['FW'] / temporary_dataframe['RV']
            
        #Store the class
        location = os.path.join(Info.General.storing_loc,exp + '_' +Info.General.date,'raw',exp + '_GenesData.pkl')   
        with open (location,'wb') as saving:
            pickle.dump(library,saving)
            
        library_evaluation(Info,library,len(raw_data.raw.keys()))
                    
        
        string = '\tRunTime: %s\n' %computeRunTime(startTime, getCurrTime())  
        
        string = 'Informations stored in %s\n***\tEND of Localizing insertions into gene models\t***'
        Info.print_save(exp,string)
        
        string = '\nThe file has been used for the following analysis:'        
        Info.print_save(exp,string)
        return library
        
    ###
    def update_library(genes_location, i, library_field):
        """The function upload insertions (i) according to gene structures provided by genes_location in the
        specific library_field"""
        gene_ids = set(genes_location[i])
        if gene_ids == set([]):
            library_field["_no_feature"] += 1
        else:
            for gene_id in gene_ids:
                library_field[gene_id] += 1 
                
    def library_evaluation(Info,library,count_II):
        strings_to_print = []
        strings_to_print.append('\t-Total number of I.I.: %i (%.2f%%)'%(count_II, 100*(float(count_II)/count_II)))
        strings_to_print.append('\t-I.I. mapped in genes: %i (%.2f %%)'%((library.II.sum()-library.II["_no_feature"]),100*(float(library.II.sum()-library.II["_no_feature"])/count_II)))
        coverage = float(len(library.II[library.II > 0])) / len(library.II)
        strings_to_print.append('\t-Coverage(genes with at least 1 I.I.): %.2f %%'%(coverage*100))
        strings_to_print.append('\t-Killing Insertions: %i (%.2f %%)'%((library.KI.sum()-library.KI["_no_feature"]),100*(float(library.KI.sum()-library.KI["_no_feature"])/count_II)))
        strings_to_print.append('\t-I.I. in introns FW: %i (%.2f %%)'%((library.bias_FW.sum()-library.bias_FW["_no_feature"]),100*(float(library.bias_FW.sum()-library.bias_FW["_no_feature"])/count_II)))
        strings_to_print.append('\t-I.I. in introns RV: %i (%.2f %%)'%((library.bias_RV.sum()-library.bias_RV["_no_feature"]),100*(float(library.bias_RV.sum()-library.bias_RV["_no_feature"])/count_II)))
        
        
        for string in strings_to_print:
            Info.print_save(exp,string)
        
    
    
    ######      Internal Functions Definitions_END!!!  ######    
    
    type_analysis = []
    if Info.GeneDefinition.Parameters.II:
        type_analysis.append('II')
    if Info.GeneDefinition.Parameters.KI:
        type_analysis.append('KI')
    if Info.GeneDefinition.Parameters.Bias:
        type_analysis.append('Bias')
    if Info.GeneDefinition.Parameters.Reads:
        type_analysis.append('Reads')
        
    
    for exp in Info.GeneDefinition.lib_names:
        string = '***\tLocalizing insertions into gene models\t***\n\n- Input file: %s\n- Gene reference location: %s\n- Parameter analyzed: %s' %(Info.GeneDefinition.input_files[Info.GeneDefinition.lib_names.index(exp)],Info.GeneDefinition.reference,' '.join(type_analysis))
        Info.print_save(exp,string)    
    
    
    # Restoring and preparation gene categories for subsequent analysis
    startTime = getCurrTime()    
    string = 'Restoring gene reference library:\n\tStarted: %s' % startTime    
    for exp in Info.GeneDefinition.lib_names:
        Info.print_save(exp,string)
     
    with open(Info.GeneDefinition.reference, 'rb') as handle: # To open and restore the Genes Dictionary (!!!CHANGE "file_name!!!)
        genes_ref = pickle.load(handle)
    
    #Creation of GenomicArrays where to store information and structure of the genes
    genes = HTSeq.GenomicArrayOfSets("auto", stranded = False) #Upload genes location, without taking into account gene strand

    genes_sense = HTSeq.GenomicArrayOfSets("auto", stranded = True)#Upload genes locations acccording to oriantaion: used to define KI insertions
    
    exons = HTSeq.GenomicArrayOfSets("auto", stranded = True)#Upload exons locations opposite oriantation: used to define KI insertions
    
    introns_sense = HTSeq.GenomicArrayOfSets("auto", stranded = True)#Upload intron locations sense oriantation: used to define bias_FW insertions
    
    introns_antisense = HTSeq.GenomicArrayOfSets("auto", stranded = True)#Upload intron locations anti-sense oriantation: used to define bias_RV insertions    
    
    
    for gene in genes_ref.index: #Uploading genes
        genes[genes_ref.ix[gene,'genomic_interval']] += gene
        
        genes_sense[genes_ref.ix[gene,'genomic_interval']] += gene
        
        for interval in genes_ref.ix[gene,'exon_specific']:
            if interval.strand == "+":
                interval_rev = HTSeq.GenomicInterval(interval.chrom,interval.start,interval.end,strand = "-")
            else:
                interval_rev = HTSeq.GenomicInterval(interval.chrom,interval.start,interval.end,strand = "+")
            exons[interval_rev] += gene
            
        for interval in genes_ref.ix[gene,'introns_all']:
            introns_sense[interval] +=gene
            if interval.strand == "+":
                interval_rev = HTSeq.GenomicInterval(interval.chrom,interval.start,interval.end,strand = "-")
            else:
                interval_rev = HTSeq.GenomicInterval(interval.chrom,interval.start,interval.end,strand = "+")
            introns_antisense[interval_rev] +=gene
    
    string = '\tRunTime: %s\n' %computeRunTime(startTime, getCurrTime())   
    for exp in Info.GeneDefinition.lib_names:
        Info.print_save(exp,string)
    
    for exp in Info.GeneDefinition.lib_names:
        library = identify_genes_insertions (exp,Info,genes_ref,genes,genes_sense,exons,introns_sense,introns_antisense)

    return library
