# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 08:33:32 2016

@author: GDM
"""

#####                  Importing modules               #####
import HTSeq
import cPickle as pickle
import pandas as pd
from collections import Counter 
import os
from HaSAPPy.HaSAPPY_time import *
import itertools
############################################################
#1) Defining the Class Library
class Library():
    """...defintion..."""
    def __init__(self,exp,input_): 
        self.name = exp
        self.input = input_
        self.informations = {'Total':'nd','Aligned':'nd','Unique_reads':'nd','Insertions':'nd','I.I.':'nd'}        
        self.raw = pd.Series()
        

#2) Function library_gneration
def library_generation (exp, Info):
    #Generation of the Class Library specific for "exp"
    library = Library(exp,Info.IIDefinition.input_files[Info.IIDefinition.lib_names.index(exp)])
    print library.name
    string ='\n\n***\tInedependent Insertions (I.I.) definition\t***\n\n- Input file: %s\n- Pair ends: %s\n- Alignment cutoff: %s\n- Remove duplicates: %s\n- Insertion cutoff: %i' %(library.input, Info.General.pair_ends,Info.IIDefinition.fidelity_limit,Info.IIDefinition.reads_duplicate,Info.IIDefinition.ins_iv)        
    Info.print_save(exp,string)
    
    startTime = getCurrTime()
    string = '\tSelection of Insertions (I.): %s' %startTime
    Info.print_save(exp,string)
    aligned_file = HTSeq.SAM_Reader(library.input)

    #aligned_file = [seq for seq in itertools.islice(aligned_file,100000)]
    
    insertions_counts = Counter()

    count_aligned = 0
    count_GoodQualityAlignment = 0
    count_total = 0
    for algnt in aligned_file:
        if algnt.aligned:
            if algnt.iv.chrom.startswith('chr'):
                chromosome_style = ''
            else:
                chromosome_style = 'chr'
                break
	
    if Info.General.pair_ends: #Pair ends library
        for bundle in HTSeq.pair_SAM_alignments(aligned_file, bundle=True):
	    if len(bundle) != 1:
	        continue # Skip multiple alignments
            first_almnt, second_almnt = bundle[0] # extract pair 
	    if first_almnt.aligned and second_almnt.aligned:
		if first_almnt.aQual >= Info.IIDefinition.fidelity_limit:
		    ins = HTSeq.GenomicPosition('%s%s' %(chromosome_style,str(first_almnt.iv.chrom)),first_almnt.iv.start_d,first_almnt.iv.strand)
		    insertions_counts[ins] +=1
		    count_GoodQualityAlignment +=1
                count_aligned +=1    
	    count_total +=1
    
    else:	#Single ends library
        for algnt in aligned_file:
            if algnt.aligned:
                if algnt.aQual >= Info.IIDefinition.fidelity_limit:
                    ins = HTSeq.GenomicPosition('%s%s' %(chromosome_style,str(algnt.iv.chrom)),algnt.iv.start_d,algnt.iv.strand)
                    insertions_counts[ins] +=1
                    count_GoodQualityAlignment +=1                    
                count_aligned +=1                    
            count_total +=1
	

    del aligned_file
            
    string = '\t-Total reads: %i\n\t-Aligned reads: %i\n\t-Aligned Reads trusted: %i\n\t-Insertions identified: %i' %(count_total,count_aligned,count_GoodQualityAlignment, len(insertions_counts.keys()))
    Info.print_save(exp,string)
          
    string = '\tRunTime: %s' % computeRunTime(startTime, getCurrTime())
    Info.print_save(exp,string)
                    
    ### To collapse insertions in insertion array that are in the same interval (4bps)
    
    
    startTime = getCurrTime()
    string = 'Define Independent Insertions\n\tStarted: %s' %startTime
    Info.print_save(exp,string)
    insertions_series = pd.Series(insertions_counts, index = insertions_counts.keys())
    del insertions_counts
    insertions_order = insertions_series.copy()
    insertions_order.sort_values(ascending = False)
    insertions_genomicarray = HTSeq.GenomicArray("auto",stranded = True)
    
    count_indipendent_insertions = 0
    count_indipendent_insertions_aborted = 0

    insertions_tuple = zip(insertions_order.index,insertions_order.values)
    del insertions_order
    del insertions_series

    for ins in insertions_tuple:
        insertions_genomicarray[ins[0]] = ins[1]	    

        
    insertions_collapsed = {}
    
    for n in insertions_tuple:
        i = n[0]
        if insertions_genomicarray[i]>0:
            counted = 0
            iv_i = HTSeq.GenomicInterval(i.chrom,i.start-2,i.start+2,i.strand)
            for i_2 in iv_i.xrange(step=1):
                try:
                    counted += insertions_genomicarray[i_2]
                    insertions_genomicarray[i_2] = 0
                except IndexError:
                    string =  "\t!!!Skipped from analysis: %s" % i_2
                    Info.print_save(exp,string)
                    continue
                    
            if counted >= Info.IIDefinition.ins_iv:
                if insertions_collapsed.has_key(i):
                    insertions_collapsed[i] += counted
                else:
                    insertions_collapsed[i] = counted
                count_indipendent_insertions +=1
            else:
                count_indipendent_insertions_aborted +=1
    
    string = '\t-Total insertions: %i\n\t-Independent Insertions (I.I.): %i' % ((count_indipendent_insertions + count_indipendent_insertions_aborted), count_indipendent_insertions)
    Info.print_save(exp,string)
            
    string = '\tRunTime: %s' % computeRunTime(startTime, getCurrTime())
    Info.print_save(exp,string)

    ###Storing data in library class that will be returned modifed as result of the function  
          
    library.informations['Total'] = count_total
    library.informations['Aligned'] = count_aligned
    library.informations['Insertions'] = count_indipendent_insertions
    library.informations['II'] = count_indipendent_insertions
    if Info.IIDefinition.reads_duplicate:
        library.informations['Unique_reads'] = count_reads

    
            
    library.raw = pd.Series(insertions_collapsed, index = insertions_collapsed.keys())
    
    #####Store the class!!!!!#####
    location = os.path.join(Info.General.storing_loc,exp + '_' +Info.General.date,'raw',exp + '_IIRawdata.pkl')   
    with open (location,'wb') as saving:
            pickle.dump(library,saving)
    #####END the program#####
    string = 'Informations stored in %s\n***\tEND of Inedependent Insertions (I.I.) definition\t***' % location 
    Info.print_save(exp,string)

    return library


def load(Info):
    for exp in Info.IIDefinition.lib_names:
        library_generation(exp,Info)
