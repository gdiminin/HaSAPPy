# -*- coding: utf-8 -*-
"""
Created on Tue May 24 08:20:07 2016

@author: GDM
"""
#### Importing modules ####
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import HTSeq
import cPickle as pickle
import os
import re
import pandas as pd
import itertools
from mpl_toolkits.mplot3d import Axes3D
mpl.interactive(False)
from HaSAPPy.HaSAPPY_time import *
####



#### Class definition ####
class GeneII():
    def __init__(self,group_name):
        self.name = group_name
        self.FW =[]
        self.RV =[]
####
        
#### Functions definition ####
        
def start(Info):
    ####
    def print_save_analysis (string, storage_loc):
        print string
        with open (os.path.join(storage_loc,'analysis_info.txt'), 'a' ) as write:
            print >> write, string
    ####
    ####
    def design_gene_insertions(Info,GroupAnalysis):
        ####
        def group_generation(Info,genes,GroupAnalysis,GroupAnalysis_group,rank):
            ####
            def upload_informations(Info,block):
                for exp in block['input']: #each replicate experiment
                    with open (exp, 'rb') as handle: #open in raw Series
                        insertions = pickle.load(handle).raw
                    for gene in genes:#iteration for each gene of interest 
                        if type(gene) == HTSeq._HTSeq.GenomicInterval:
                            iv = gene
                        else:
                            iv = genes_ref.loc[gene,'genomic_interval'].copy()#genomic interval of the gene of intersts
                            iv.strand = '.' #remove the strand parameter (for the moment)
			selected_ins  = []
			insertion_list = zip(insertions.index,insertions.values)                     
			for n in insertion_list:
			    i = n[0]
			    if iv.contains(i):
				selected_ins.append((i,n[1]))
				

#                        selected_ins= [(i,insertions[i]) for i in insertions.index if iv.contains(i)] #list of insertions in the replicate experiment containined in gene
                        for i in selected_ins: #divide in sense and anti-sense insertions
                            if type(gene) == HTSeq._HTSeq.GenomicInterval:
                                if i[0].strand == '+':
                                    block[gene].FW.append(i)
                                else:
                                    block[gene].RV.append(i)
                            else:        
                                if i[0].strand == genes[gene]['genomic_interval'].strand:
                                    block[gene].FW.append(i)
                                else:
                                    block[gene].RV.append(i)
                return block
            ####
            
            block = {}
            block['input'] = []
            if rank == 'NaN':            
                block['name'] = GroupAnalysis_group.name
            else:
                block['name'] = GroupAnalysis_group.name[rank]
            for gene in genes:
                block[gene] = GeneII(block['name'])
            if rank == 'NaN':
                experiments = GroupAnalysis_group.experiments
            else:
                experiments = GroupAnalysis_group.experiments[rank]
            for exp in experiments:
                pos = GroupAnalysis.lib_names.index(exp)
                address = GroupAnalysis.input_files[pos]
                with open (address,'rb') as handle:
                    block['input'].append(pickle.load(handle).input)     
            block = upload_informations(Info,block)
            return block
            ####
            
            
        def draw_gene(Info,name,gene,group_reference,group_other,storage_loc):
            
            
            if type(name) == HTSeq._HTSeq.GenomicInterval:
                genomic_interval = True
            else:
                genomic_interval = False
                
            fig1 = plt.figure()
            if genomic_interval:
                fig1.suptitle('%s'%name)
            else:
                fig1.suptitle('%s (%s)'%(name,gene['genomic_interval']))
            ax1 = fig1.add_subplot(311)    
            ax2 = fig1.add_subplot(312)
            ax3 = fig1.add_subplot(313)    
            
            if genomic_interval:
                start_pos = name.start_d
                end_pos = name.end_d
            else:
                start_pos = gene['genomic_interval'].start_d
                end_pos = gene['genomic_interval'].end_d
            
            #To define max read value in the two libraries for get ymax position
            try:
                x_all,y_all = zip(*(group_reference.FW + group_reference.RV + group_other.FW +group_other.RV))
                ymax = max(y_all)
            except ValueError:
                ymax = 10
                
            
            #Eventually add here log schale if ymax > ...
            
            #Design Scatter for 1st experiment    
            ax1.axis(xmin =start_pos, xmax = end_pos, ymin = 0, ymax = ymax) #First scatter plot
            ax1.set_ylabel(group_other.name) #Name of the experiment to show in ylabel
            ax1.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom='off',      # ticks along the bottom edge are off
                top='off',         # ticks along the top edge are off
                labelbottom='off') # labels along the bottom edge are off
            if ymax > 200:
                ax1.set_yscale('log')
                ax1.axis(xmin =start_pos, xmax = end_pos, ymin = 1, ymax = ymax)
                
            for i in group_other.FW: # Plotting points according thier value. Devided in two colors: red if sense, green if anti-sense
                ax1.scatter(i[0].pos,i[1],s = 10, color = 'r')
            for i in group_other.RV:
                ax1.scatter(i[0].pos,i[1],s = 10, color = 'g')
                
            
            
            #Design Scatter for 2nd experiment 
            ax3.axis(xmin =start_pos, xmax = end_pos, ymin = 0, ymax =ymax)
            ax3.invert_yaxis()
        
            ax3.set_ylabel(group_reference.name) #Name of the experiment to show in ylabel
            ax3.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom='off',      # ticks along the bottom edge are off
                top='off',         # ticks along the top edge are off
                labelbottom='off') # labels along the bottom edge are off
            if ymax > 200:
                ax3.set_yscale('log')
                ax3.axis(xmin =start_pos, xmax = end_pos, ymin = 1, ymax = ymax)
                ax3.invert_yaxis()
                
            for i in group_reference.FW: # Plotting points according thier value. Devided in two colors: red if sense, green if anti-sense
                ax3.scatter(i[0].pos,i[1],s = 10, color = 'r')
            for i in group_reference.RV:
                ax3.scatter(i[0].pos,i[1],s = 10, color = 'g')
            
            #Design gene models
            
            if genomic_interval:
                transcripts = gene
            else:
                transcripts = gene['variants']
            
            ax2.axis([start_pos,end_pos,0,len(transcripts)+1])
            ax2.axis('off')
            
            y_value = 0 #location of gene in y axes, acccording to transcripts number
            for transcript in transcripts:
                y_value +=1 #move 1 up
                
                ax2.text((end_pos), # Transcript name starting position x-axis (end of gene)
                        (y_value-0.2), # Transcript name starting position y-axis
                        ('   ' + transcript),fontsize = 10)
                        
                #line rapresenting all the transcript length        
                ax2.plot([min([exon.start_d for exon in transcripts[transcript]]), max([exon.end_d for exon in transcripts[transcript]])],[y_value,y_value],'k',linewidth = 2./len(transcripts))
                
                for exon in transcripts[transcript]:
                    ax2.add_patch(patches.Rectangle(
                    (exon.start,(y_value-0.2)), exon.length,0.4,linewidth = 0.1,facecolor='black'))
                
        
            fig1.show()
            if genomic_interval:
                fig1.savefig(os.path.join(storage_loc,'graph','%s_%svs%s.svg'%('interval',group_other.name,group_reference.name)),dpi=300, bbox_inches='tight')
                string = '\t\t- %s' %os.path.join(storage_loc,'graph','%s_%svs%s.svg'%('interval',group_other.name,group_reference.name))
            else:
                fig1.savefig(os.path.join(storage_loc,'graph','%s_%svs%s.svg'%(name,group_other.name,group_reference.name)),dpi=300, bbox_inches='tight')
                string = '\t\t- %s' %os.path.join(storage_loc,'graph','%s_%svs%s.svg'%(name,group_other.name,group_reference.name))
            print_save_analysis (string, GroupAnalysis.storage_loc)
            ####
        
        strings = []    
        strings.append('\nPlot I.I. in gene models')
        strings.append('\tSelected genes (or intervals):') 
        for gene in Info.Design.Gene_model.genes:
            strings.append('\t\t-%s' %gene) 
        
        for string in strings:
            print_save_analysis (string, GroupAnalysis.storage_loc)   
            
        startTime = getCurrTime()      
        string= '\tRestoring Gene Reference file\n\t\tStarted: %s' % startTime
        print_save_analysis (string, GroupAnalysis.storage_loc) 
                
        genome = HTSeq.GenomicArrayOfSets("auto", stranded = False)
        with open(Info.Design.Gene_model.reference, 'rb') as handle:
            genes_ref = pickle.load(handle)
        for gene in genes_ref.index:
            genome[genes_ref.ix[gene,'genomic_interval']] += gene
        
        string = '\t\tRunTime: %s' % computeRunTime(startTime,getCurrTime()) 
        print_save_analysis (string, GroupAnalysis.storage_loc)           
        
        
        genes ={} 
        for gene in Info.Design.Gene_model.genes:
            if re.search ('\((.+)_(.+)_(.+)\)',gene):
                value = re.findall ('\((.+)_(.+)_(.+)\)',gene)[0]
                iv = HTSeq.GenomicInterval(value[0],int(value[1]),int(value[2]),'.')
                genes_in_iv = [x[1] for x in genome[iv].steps()]
                genes_selected = set()
                for group in genes_in_iv:
                    for gene in group:
                        genes_selected.add(gene)
                
                genes[iv] = {}                
                for gene in genes_selected:
                    genes[iv][gene] = genes_ref.loc[gene,'variants'].values()[0]  
            else:
                genes[gene] = genes_ref.loc[gene]
        
        reference = group_generation(Info,genes,GroupAnalysis,GroupAnalysis.Reference,'NaN')
    
        others = []
        for group in GroupAnalysis.Others.name:
            others.append(group_generation(Info,genes,GroupAnalysis,GroupAnalysis.Others,GroupAnalysis.Others.name.index(group)))
        
        string = '\tGenes saved in:'
        print_save_analysis (string, GroupAnalysis.storage_loc)
        for name in genes:
            for group in others:
                draw_gene(Info,name,genes[name],reference[name],group[name],GroupAnalysis.storage_loc)
    
    
    ####
    def design_insertions_distribution(Info,GroupAnalysis):
        ####
        def plot_distribution(GroupaAnalysis,data,to_plot,group):
            ####    
            def scatter3D(dataframe,comb_columns,color):
                ax.scatter(dataframe[comb_columns[0]],dataframe[comb_columns[1]],dataframe[comb_columns[2]],c =color)
            def scatter2D(dataframe,comb_columns,color):
                ax.scatter(dataframe[comb_columns[0]],dataframe[comb_columns[1]],c =color)
            ####
            
            
            
            
            if len(to_plot)>2:
                string = '\t3D plots:'
                print_save_analysis (string, GroupAnalysis.storage_loc)
                
                files_saved = []

                for comb_columns in itertools.combinations(to_plot,3):
                    fig=plt.figure(dpi=600, facecolor='w', edgecolor='k')
                    ax = Axes3D(fig)
                    ax = plt.subplot2grid((1,1), (0,0),projection='3d')
                    if Info.Design.Distribution.outlier and GroupAnalysis.Outlier.perform:
                        #Outlier design
                        scatter3D(data[data['Score'] >= Info.Design.Distribution.outlier_value],comb_columns,'r')
                        #Normal design
                        scatter3D(data[data['Score'] < Info.Design.Distribution.outlier_value],comb_columns,'b')
                    else:
                        #All design
                        scatter3D(data,comb_columns,'b')
                    #annotate
                    if Info.Design.Distribution.annotate:
                        if Info.Design.Distribution.annotate_value.isdigit(): 
                            if Info.Design.Distribution.outlier and GroupAnalysis.Outlier.perform:
                                index_to_annotate = [index for index in data[data['Score']>= int(Info.Design.Distribution.annotate_value)].index]
                        else:
                            index_to_annotate = [index.lstrip(' ').rstrip(' ') for index in Info.Design.Distribution.annotate_value.split(',')]
                        
                        for index in index_to_annotate:
                            try:
                                ax.text(data.ix[index,comb_columns[0]],data.ix[index,comb_columns[1]],data.ix[index,comb_columns[2]],index, size=10, zorder=40,color='k')
                            except KeyError:
                                string =  'Warning: %s was not found' %index
                                print_save_analysis (string, GroupAnalysis.storage_loc)

                    ax.set_xlabel(comb_columns[0])
                    ax.set_ylabel(comb_columns[1])
                    ax.set_zlabel(comb_columns[2])
                    ax.set_title('Gene Distribution\n%s, %s, %s' %(comb_columns[0],comb_columns[1],comb_columns[2]))
            
                    fig.savefig(os.path.join(GroupAnalysis.storage_loc,'graph','Outliers3D_%s(%s-%s-%s).jpg' %(group,comb_columns[0],comb_columns[1],comb_columns[2])),dpi= 200)
                    files_saved.append(os.path.join(GroupAnalysis.storage_loc,'graph','Outliers3D_%s(%s-%s-%s).jpg' %(group,comb_columns[0],comb_columns[1],comb_columns[2])))
                
                string =  '\t\tFile saved:'
                print_save_analysis (string, GroupAnalysis.storage_loc)
                    
                for file_ in files_saved:
                    string =  '\t\t- %s'% file_
                    print_save_analysis (string, GroupAnalysis.storage_loc)
                    
            string = '\t2D plots:'
            print_save_analysis (string, GroupAnalysis.storage_loc)
            
            files_saved = []
            for comb_columns in itertools.combinations(to_plot,2):
                fig=plt.figure(facecolor='w', edgecolor='k')
                ax = fig.add_subplot(111) 
                if Info.Design.Distribution.outlier and GroupAnalysis.Outlier.perform:
                    #Outlier design
                    scatter2D(data[data['Score'] >= Info.Design.Distribution.outlier_value],comb_columns,'r')
                    #Normal design
                    scatter2D(data[data['Score'] < Info.Design.Distribution.outlier_value],comb_columns,'b')
                else:
                    #All design
                    scatter2D(data,comb_columns,'b')
                #annotate
                if Info.Design.Distribution.annotate:
                    if Info.Design.Distribution.annotate_value.isdigit(): 
                        if Info.Design.Distribution.outlier and GroupAnalysis.Outlier.perform:
                            index_to_annotate = [index for index in data[data['Score']>= int(Info.Design.Distribution.annotate_value)].index]
                    else:
                        index_to_annotate = [index.lstrip(' ').rstrip(' ') for index in Info.Design.Distribution.annotate_value.split(',')]
                    
                    for index in index_to_annotate:
                        try:
                            ax.text(data.ix[index,comb_columns[0]],data.ix[index,comb_columns[1]],index, size=10, zorder=40,color='k')
                        except KeyError:
                            string =  'Warning: %s was not found' %index
                            print_save_analysis (string, GroupAnalysis.storage_loc) 

                ax.set_xlabel(comb_columns[0])
                ax.set_ylabel(comb_columns[1])
                ax.set_title('Gene Distribution\n%s, %s' %(comb_columns[0],comb_columns[1]))
                         
                fig.savefig(os.path.join(GroupAnalysis.storage_loc,'graph','Outliers2D_%s(%s-%s).jpg' %(group,comb_columns[0],comb_columns[1])),dpi= 200)
                files_saved.append(os.path.join(GroupAnalysis.storage_loc,'graph','Outliers2D_%s(%s-%s).jpg' %(group,comb_columns[0],comb_columns[1])))
                
            string =  '\t\tFile saved:'
            print_save_analysis (string, GroupAnalysis.storage_loc)
                    
            for file_ in files_saved:
                string =  '\t\t- %s'% file_
                print_save_analysis (string, GroupAnalysis.storage_loc)
                
        ####
        strings = []    
        
        strings.append('\nPlot gene distributions')
        if Info.Design.Distribution.outlier:
            strings.append('\tOutlier value: %s' % Info.Design.Distribution.outlier_value)
        if Info.Design.Distribution.annotate:
            strings.append('\tAnnotated genes or value: %s' % Info.Design.Distribution.annotate_value)
        for string in strings:
            print_save_analysis (string, GroupAnalysis.storage_loc)   
            
        startTime = getCurrTime()      
        string= '\nRestoring RawData files\n\tStarted: %s' % startTime
        print_save_analysis (string, GroupAnalysis.storage_loc)
            
        with open (os.path.join(GroupAnalysis.storage_loc,'raw','RawData.pkl'), 'rb') as loading:
            GroupAnalysis_table = (pickle.load(loading)).all
    
        string= '\tRowData.pkl location: %s' % os.path.join(GroupAnalysis.storage_loc,'raw', 'RawData.pkl')
        print_save_analysis (string, GroupAnalysis.storage_loc)
        string = '\tRunTime: %s' % computeRunTime(startTime,getCurrTime())
        print_save_analysis (string, GroupAnalysis.storage_loc)
        
        startTime = getCurrTime()      
        string= '\nPlot Generation\n\tStarted: %s' % startTime
        print_save_analysis (string, GroupAnalysis.storage_loc)    
        
        for group in GroupAnalysis.Others.name:
            data = pd.DataFrame()
            to_plot = []
            if GroupAnalysis.Parameters.II:
                to_plot.append('II')
                if GroupAnalysis.Outlier.perform and GroupAnalysis.Outlier.Parameters.II:
                    data['II'] = GroupAnalysis_table['%s_%s_outlier' %(group,'II')]
                else:
                    data['II'] = GroupAnalysis_table['%s_%s_fold' %(group,'II')]
            if GroupAnalysis.Parameters.KI:
                to_plot.append('KI')
                if GroupAnalysis.Outlier.perform and GroupAnalysis.Outlier.Parameters.KI:
                    data['KI'] = GroupAnalysis_table['%s_%s_outlier' %(group,'KI')]
                else:
                    data['KI'] = GroupAnalysis_table['%s_%s_fold' %(group,'KI')]
            if GroupAnalysis.Parameters.Bias:
                to_plot.append('Bias')
                if GroupAnalysis.Outlier.perform and GroupAnalysis.Outlier.Parameters.Bias:
                    data['Bias'] = GroupAnalysis_table['%s_%s_outlier' %(group,'Bias')]
                else:
                    data['Bias'] = GroupAnalysis_table['%s_%s_fold' %(group,'Bias')]
            if GroupAnalysis.Parameters.Reads:
                to_plot.append('Reads')
                if GroupAnalysis.Outlier.perform and GroupAnalysis.Outlier.Parameters.Reads:
                    data['Reads'] = GroupAnalysis_table['%s_%s_outlier' %(group,'Reads')]
                else:
                    data['Reads'] = GroupAnalysis_table['%s_%s_fold' %(group,'Reads')]
            if GroupAnalysis.Outlier.perform:
                data['Score'] = GroupAnalysis_table['%s_%s' %(group,'Score_fold')]

            data = data[data['II']>0]
            plot_distribution(GroupAnalysis,data,to_plot,group)
        
        string = '\tRunTime: %s' % computeRunTime(startTime,getCurrTime())
        print_save_analysis (string, GroupAnalysis.storage_loc)
        
    ####
    date_today = getDay()
    strings = []
    strings.append('\n***\tPlot generation\t***\t\tDate: %s' % date_today)
    startTime = getCurrTime()      
    strings.append('\nRestoring GroupAnalysis information\n\tStarted: %s' % startTime)

    
    with open (Info.Design.input_files,'rb') as handle: #loading GroupAnlaysis class exential to get all information on the samples
        GroupAnalysis = pickle.load(handle)
    strings.append('\tGroupAnalysis location: %s' % Info.Design.input_files)
    strings.append('\tRunTime: %s' % computeRunTime(startTime,getCurrTime()))
    for string in strings:
        print_save_analysis (string, GroupAnalysis.storage_loc)
        
    if Info.Design.Gene_model.perform:
        design_gene_insertions(Info,GroupAnalysis)
    if Info.Design.Distribution.perform:
        design_insertions_distribution(Info,GroupAnalysis)
        
    string = '\n***\tEND Plot Generation\t***' 
    print_save_analysis (string, GroupAnalysis.storage_loc)
        
        
            
    
            
            
            
            






