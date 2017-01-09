# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 19:26:59 2016

@author: GDM
"""

import re
import os
import subprocess
from HaSAPPY_time import *


def load(Info):
    for exp in Info.Trim.lib_names:
        startTime = getCurrTime()
        string = '\n\n***\tTrim adaptor and select for Quality\t***\n\tStarted: %s\n' %(startTime)
        Info.print_save(exp,string)
        input_file_I = Info.Trim.input_files_I[Info.Trim.lib_names.index(exp)]
        output_file_I = os.path.join(Info.General.storing_loc,exp + '_' + Info.General.date,'row',exp + '_QS.fastq')
        command = './PreprocessReads -i %s -o %s -len 25 -a %s -f' %(input_file_I,output_file_I,Info.Trim.p7adaptor)
        run_output = subprocess.check_output(command,shell= True,stderr=subprocess.STDOUT)
        total_reads = int(re.findall('\\nReads in input file \S* : (\d*)\\n',run_output)[0])
        selected_reads= int(re.findall('\\nReads written to output file \S* : (\d*)\\n',run_output)[0])
        string = "\tCommand used: %s\n\tAnalysed reads: %d\n\tReads passed selection: %d\n\tRunTime: %s\n\tStored in: %s\n" %(command,total_reads,selected_reads,computeRunTime(startTime,getCurrTime()),output_file_I)
        Info.print_save(exp,string)
        if Info.General.pair_ends:
            startTime = getCurrTime()
            string = '\tII strain: %s' % startTime
            Info.print_save(exp,string)
            input_file_II = Info.Trim.input_files_II[Info.Trim.lib_names.index(exp)]
            output_file_II = os.path.join(Info.General.storing_loc,exp + '_' + Info.General.date,'row',exp + '_pairend_QS.fastq')
            
            command = './PreprocessReads -i %s -o %s -len 25 -a %s -f' %(input_file_II,output_file_II,Info.Trim.p5adaptor)
            print command
            
            run_output = subprocess.check_output(command,shell= True,stderr=subprocess.STDOUT)
            total_reads = int(re.findall('\\nReads in input file \S* : (\d*)\\n',run_output)[0])
            selected_reads= int(re.findall('\\nReads written to output file \S* : (\d*)\\n',run_output)[0])
            
            string = "\tCommand used: %s\n\tAnalysed reads: %d\n\tReads passed selection: %d\n\tRunTime: %s\n\tStored in: %s\n" %(command,total_reads,selected_reads,computeRunTime(startTime,getCurrTime()),output_file_II)
            Info.print_save(exp,string)
        
        string = '***\tEND of Trim adaptor and select for Quality\t***'
        Info.print_save(exp,string)

            
        
        
    
    
