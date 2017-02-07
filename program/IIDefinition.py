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
from HaSAPPY_time import *
############################################################
#1) Defining the Class Library
class Library():
    """...defintion..."""
    def __init__(self,exp,input_): 
        self.name = exp
        self.input = input_
        self.informations = {'Total':'nd','Aligned':'nd','Unique_reads':'nd','Insertions':'nd','I.I.':'nd'}        
        self.row = pd.Series()
        

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
    insertions_counts = Counter()

    if not Info.IIDefinition.reads_duplicate: # single_end sequencing!!!!
        count_aligned = 0
        count_GoodQualityAllignment = 0
        count_total = 0
        for algnt in aligned_file:
            if algnt.aligned:
                if algnt.aQual >= Info.IIDefinition.fidelity_limit:
		    ins = HTSeq.GenomicPosition(algnt.iv.chrom,algnt.iv.start_d,algnt.iv.strand)
#                    ins = HTSeq.GenomicPosition('chr%s' %str(algnt.iv.chrom),algnt.iv.start_d,algnt.iv.strand)
                    insertions_counts[ins] +=1
                    count_GoodQualityAllignment +=1                    
                count_aligned +=1                    
            count_total +=1
            
        string = '\t-Total reads: %i\n\t-Alligned reads: %i\n\t-Alligned Reads trusted: %i\n\t-Insertions identified: %i' %(count_total,count_aligned,count_GoodQualityAllignment, len(insertions_counts.keys()))
        Info.print_save(exp,string)
      
    else: # pair_end sequencing = Allow analysis of number of reads       
        insertions = {}
        count_total = 0
        count_aligned = 0
        count_position = 0
        count_duplicates = 0
        count_GoodQualityAllignment = 0

        for algnt in aligned_file:
            count_total +=1 #all the read in file
            
            if algnt.pe_which == "first" and algnt.proper_pair:#just consider sequenced 1)alligned, 2)for which a paired_end has been identified, 3) coming from the "first" file (in theory p5 strand)
                if algnt.aQual >= Info.IIDefinition.fidelity_limit:
                    count_GoodQualityAllignment +=1
                    i_p5 = HTSeq.GenomicPosition('chr%s'%str(algnt.iv.chrom),algnt.iv.start_d,algnt.iv.strand) # Viral insertion genomic position
                    if insertions.has_key(i_p5):
                        if algnt.inferred_insert_size in insertions[i_p5]:#where algnt.mate_start is the starting genomic position of the reverse paired read (that for us..is the LAM_PCR end)
                            count_duplicates +=1 # Discard sequence (PCR artifact)      
                        else:
                            insertions[i_p5].add(algnt.inferred_insert_size) #Add a new element to the insertion's list (=no redoundant read)
                    else: #Identification of a new insertion
                        count_position +=1
                        insertions[i_p5] = set() #Creation of a list were to store length of different reads
                        insertions[i_p5].add(algnt.inferred_insert_size)
                        
                count_aligned +=1 #all the reads with pair_end alignment and coming from the p5 file
                        
        
        count_reads = 0        
        for ins in insertions:
            length = len(insertions[ins])
            insertions_counts[ins] = length
            count_reads += length
            
        string = '\t-Total reads: %i\n\t-Alligned reads: %i\n\t-Alligned Reads trusted: %i\n\t-Detected duplicates: %i\n\tInsertions identified: %i' %(count_total,count_aligned,count_GoodQualityAllignment,count_duplicates,len(insertions_counts.keys()))
        Info.print_save(exp,string)
        
    string = '\tRunTime: %s' % computeRunTime(startTime, getCurrTime())
    Info.print_save(exp,string)
                    
    ### To collapse insertions in insertion array that are in the same interval (4bps)
    
    
    startTime = getCurrTime()
    string = 'Define Independent Insertions\n\tStarted: %s' %startTime
    Info.print_save(exp,string)
    insertions_series = pd.Series(insertions_counts, index = insertions_counts.keys())
    insertions_series = insertions_series.order(ascending = False)
    insertions_genomicarray = HTSeq.GenomicArray("auto",stranded = True)
    
    count_indipendent_insertions = 0
    count_indipendent_insertions_aborted = 0
    
    for ins in insertions_series.index:
        insertions_genomicarray[ins] = insertions_series[ins]
        
    insertions_collapsed = {}
    
    for i in insertions_series.index:
        
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

    
            
    library.row = pd.Series(insertions_collapsed, index = insertions_collapsed.keys())
    
    #####Store the class!!!!!#####
    location = os.path.join(Info.General.storing_loc,exp + '_' +Info.General.date,'row',exp + '_IIRowdata.pkl')   
    with open (location,'wb') as saving:
            pickle.dump(library,saving)
    #####END the program#####
    string = 'Informations stored in %s\n***\tEND of Inedependent Insertions (I.I.) definition\t***' % location 
    Info.print_save(exp,string)

    return library


def load(Info):
    for exp in Info.IIDefinition.lib_names:
        library_generation(exp,Info)
