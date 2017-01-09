# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 07:28:14 2016

@author: GDM
"""

import subprocess
import os
import HTSeq
from HaSAPPY_time import *

#############################    
####Programs involved in align library with Bowtie2 cheack__BEGIN

####
def load(Info):

    ####
    def Align_library(Info,Info_selected,exp,storing_loc,storing_suffix,program,pair_end = False):
        #Initializing
    
        input_file_I = Info_selected.input_files_I[Info_selected.lib_names.index(exp)]
        if pair_end:
            input_file_II = Info_selected.input_files_II[Info_selected.lib_names.index(exp)]
        
        output_file_I = os.path.join(storing_loc,exp + storing_suffix)  
        if pair_end:
            output_file_II = os.path.join(storing_loc,exp + '_pairend' + storing_suffix)
        
        if Info.Type_analysis.AllignPhix:
            if Info_selected == Info.AllignPhix:
                command = "bowtie2 --time --threads 8 -x %s %s -S %s" %(Info_selected.reference,input_file_I,output_file_I) #bowtie2 command
                if pair_end:
                    command_pair_end = "bowtie2 --time --threads 8 -x %s %s -S %s" %(Info_selected.reference,input_file_II,output_file_II) #bowtie2 command

        if Info.Type_analysis.AllignGenome:
            if Info_selected == Info.AllignGenome : 
                if program == 'bowtie2':
                    if pair_end:
                        command = "bowtie2 --time --local --threads 8 -x %s -1 %s -2 %s -S %s" %(Info_selected.reference,input_file_I,input_file_II,output_file_I)
                    else:
                        command = "bowtie2 --time --local --threads 8 -x %s %s -S %s" %(Info_selected.reference,input_file_I,output_file_I)
                if program == 'NextGenMap':
                    if pair_end:
                        pass #understand the command to allow pair-end allignment
                    else:
                        command = 'ngm -q %s -r %s -o %s -t 8 -g' %(input_file_I,Info_selected.reference,output_file_I)
                if program == 'nvBowtie':
                    if pair_end:
                        pass #understand the command to allow pair-end allignment
                    else:
                        command = 'nvBowtie --local --device 0 -U %s -x %s -S %s' % (input_file_I,Info_selected.reference,output_file_I)
        
        startTime = getCurrTime()
        run_output = subprocess.check_output(command,shell= True,stderr=subprocess.STDOUT)
        string = "\t\tCommand used: %s\n\t\t%s\n\tStarted: %s\tRunTime: %s\n\t\t\tStored in: %s\n" %(command,run_output,startTime,computeRunTime(startTime, getCurrTime()),output_file_I)
        Info.print_save(exp,string)
        
        if Info.Type_analysis.AllignPhix:
            if Info_selected == Info.AllignPhix:
                if pair_end:
                    startTime = getCurrTime()
                    run_output = subprocess.check_output(command_pair_end,shell= True,stderr=subprocess.STDOUT)
                    string = "\t\tCommand used: %s\n\t\t%s\n\t\tStarted: %s\tRunTime: %s\n\t\tStored in: %s\n" %(command_pair_end,run_output,startTime,computeRunTime(startTime, getCurrTime()),output_file_II)
                    Info.print_save(exp,string)
    ####
                
    def Trim_Phix(Info,exp,storing_loc,pair_end):
        to_clean = []
        to_clean.append(os.path.join(storing_loc,exp + '_Alligned.sam'))
        if pair_end:
            to_clean.append(os.path.join(storing_loc,exp + '_pairend_Alligned.sam'))
            
        for library in to_clean:                                                                               
            Phix_blast = HTSeq.SAM_Reader(library)
            if to_clean.index(library) == 0:
                storage_file = os.path.join(storing_loc,exp + '_PhixCleaned.fastq')
            elif to_clean.index(library) == 1:
                storage_file = os.path.join(storing_loc,exp + '_pairend_PhixCleaned.fastq')
                
            with open(storage_file,"w") as Library_Trimmed_Phix:
                count_total = 0
                count_selected = 0
                for read in Phix_blast:
                    count_total +=1
                    if not read.aligned:
                        selected_fastq = read.read
                        selected_fastq.write_to_fastq_file(Library_Trimmed_Phix)
                        count_selected +=1
                    if count_total%1000000 == 0:
                        print "Analyzed ", count_total, "sequences."
    
            string = "\tSaved in: %s\n\t\tSelected(NO Phix) sequences %i of %i."% (storage_file,count_selected,count_total)
            Info.print_save(exp,string)
    ####
            
            
    if Info.Type_analysis.AllignPhix:
        for exp in Info.AllignPhix.lib_names:      
            string = "\n\n***\tAllignment\t***\n"
            Info.print_save(exp,string)
    else:
        if Info.Type_analysis.AllignGenome:
            for exp in Info.AllignGenome.lib_names:
                string = "\n\n***\tAllignment\t***\n"
                Info.print_save(exp,string)
                
    
            
    if Info.Type_analysis.AllignPhix:        
        for exp in Info.AllignPhix.lib_names:
            startTime = getCurrTime()
            string = 'Allign to Phix genome\n\tStarted: %s' %(startTime)
            Info.print_save(exp,string)
            storing_loc = os.path.join(Info.General.storing_loc,exp + '_' +Info.General.date,'row')
            storing_suffix = '_Alligned.sam'
            
            Align_library(Info,Info.AllignPhix,exp,storing_loc,storing_suffix,program = 'bowtie2',pair_end = Info.General.pair_ends)
            
            string = '\tRunTime: %s' % computeRunTime(startTime, getCurrTime())
            Info.print_save(exp,string)
            
            
            startTime = getCurrTime()
            string = "\tCleaning library from Phix genes.\n\t\t\tStarted: %s"% startTime
            Info.print_save(exp,string)
            
            Trim_Phix(Info,exp,storing_loc,pair_end = Info.General.pair_ends)
            
            string = '\t\tRunTime: %s' % computeRunTime(startTime, getCurrTime())
            Info.print_save(exp,string)
            
            if Info.Type_analysis.Trim:
                if not Info.Trim.store:
                    os.remove(Info.AllignPhix.input_files_I[Info.AllignPhix.lib_names.index(exp)])
                    if not os.path.exists(Info.AllignPhix.input_files_I[Info.AllignPhix.lib_names.index(exp)]):
                        string = "\tCancelled file:\n\t\t-%s" % Info.AllignPhix.input_files_I[Info.AllignPhix.lib_names.index(exp)]
                        Info.print_save(exp,string)
                    if Info.General.pair_ends:
                        os.remove(Info.AllignPhix.input_files_II[Info.AllignPhix.lib_names.index(exp)])
                        if not os.path.exists(Info.AllignPhix.input_files_I[Info.AllignPhix.lib_names.index(exp)]):
                            string = "\tCancelled file:\n\t\t-%s" % Info.AllignPhix.input_files_I[Info.AllignPhix.lib_names.index(exp)]
                            Info.print_save(exp,string)
                        
                
    
    if Info.Type_analysis.AllignGenome:
        for exp in Info.AllignGenome.lib_names:
            print exp
            startTime = getCurrTime()
            string = 'Allign to Reference genome:\n\tStarted: %s' %startTime
            Info.print_save(exp,string)
            storing_loc = os.path.join(Info.General.storing_loc,exp + '_' +Info.General.date,'row')
            storing_suffix = '_Alligned.sam'
            
            Align_library(Info,Info.AllignGenome,exp,storing_loc,storing_suffix,program = Info.AllignGenome.program_type,pair_end = Info.General.pair_ends)
            
            string = '\tRunTime: %s' % computeRunTime(startTime, getCurrTime())
            Info.print_save(exp,string)

    if Info.Type_analysis.AllignPhix:    
        if not Info.AllignPhix.store:
            for exp in Info.AllignGenome.lib_names:
                 os.remove(os.path.join(storing_loc,exp + '_PhixCleaned.fastq'))
                 if not os.path.exists(os.path.join(storing_loc,exp + '_PhixCleaned.fastq')):
                     string = "\tCancelled files:\n\t\t-%s" % os.path.join(storing_loc,exp + '_PhixCleaned.fastq')
                     Info.print_save(exp,string)
                 if Info.General.pair_ends:
                     os.remove(os.path.join(storing_loc,exp + '_pairend_PhixCleaned.fastq'))
                     if not os.path.exists(os.path.join(storing_loc,exp + '_pairend_PhixCleaned.fastq')):
                         string = "\t\t-%s" % os.path.join(storing_loc,exp + '_pairend_PhixCleaned.fastq')
                         Info.print_save(exp,string)
                     
                     
    if Info.Type_analysis.AllignGenome:
        for exp in Info.AllignGenome.lib_names:      
            string = "***\tEND Allignment\t***"
            Info.print_save(exp,string)
    else:
        if Info.Type_analysis.AllignPhix:
            for exp in Info.AllignPhix.lib_names:   
                string = "***\tEND Allignment\t***"
                Info.print_save(exp,string)
                 
    

#############################    
