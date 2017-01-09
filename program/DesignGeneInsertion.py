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
mpl.interactive(False)
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
    def group_generation(Info,genes,GroupAnalysis,GroupAnalysis_group,rank):
        ####
        def upload_informations(Info,block):
            for exp in block['input']: #each replicate experiment
                with open (exp, 'rb') as handle: #openin Row Series
                    insertions = pickle.load(handle).row
                for gene in genes:#iteration for each gene of interest 
                    if type(gene) == HTSeq._HTSeq.GenomicInterval:
                        iv = gene
                    else:
                        iv = genes_ref.loc[gene,'genomic_interval'].copy()#genomic interval of the gene of intersts
                        iv.strand = '.' #remove the strand parameter (for the moment)
                    
                    selected_ins= [(i,insertions[i]) for i in insertions.index if iv.contains(i)] #list of insertions in the replicate experiment containined in gene
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
        x_all,y_all = zip(*(group_reference.FW + group_reference.RV + group_other.FW +group_other.RV))
        ymax = max(y_all)
            
        
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
            ax1.scatter(i[0].pos,i[1],s = 1, color = 'r')
        for i in group_other.RV:
            ax1.scatter(i[0].pos,i[1],s = 1, color = 'g')
            
        
        
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
            ax3.scatter(i[0].pos,i[1],s = 1, color = 'r')
        for i in group_reference.RV:
            ax3.scatter(i[0].pos,i[1],s = 1, color = 'g')
        
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
                    ('   ' + transcript),fontsize = 4)
                    
            #line rapresenting all the transcript length        
            ax2.plot([min([exon.start_d for exon in transcripts[transcript]]), max([exon.end_d for exon in transcripts[transcript]])],[y_value,y_value],'b',linewidth = 2./len(transcripts))
            
            for exon in transcripts[transcript]:
                ax2.add_patch(patches.Rectangle(
                (exon.start,(y_value-0.2)), exon.length,0.4,linewidth = 0.1))
            
    
        fig1.show()
        if genomic_interval:
            fig1.savefig(os.path.join(storage_loc,'%s_%svs%s.svg'%('interval',group_other.name,group_reference.name)))
        else:
            fig1.savefig(os.path.join(storage_loc,'%s_%svs%s.svg'%(name,group_other.name,group_reference.name)))
        ####
        
    genome = HTSeq.GenomicArrayOfSets("auto", stranded = False)
    with open(Info.Design.reference, 'rb') as handle:
        genes_ref = pickle.load(handle)
    for gene in genes_ref.index:
        genome[genes_ref.ix[gene,'genomic_interval']] += gene
               
    
    
    genes ={} 
    for gene in Info.Design.genes:
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
        
    
    with open (Info.Design.input_files,'rb') as handle: #loading GroupAnlaysis class exential for get all information on the samples
        GroupAnalysis = pickle.load(handle)
    
    reference = group_generation(Info,genes,GroupAnalysis,GroupAnalysis.Reference,'NaN')

    others = []
    for group in GroupAnalysis.Others.name:
        others.append(group_generation(Info,genes,GroupAnalysis,GroupAnalysis.Others,GroupAnalysis.Others.name.index(group)))
    
    for name in genes:
        for group in others:
            draw_gene(Info,name,genes[name],reference[name],group[name],GroupAnalysis.storage_loc)
            
            






