# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 08:04:13 2016

@author: GDM
"""
import re
import os
import datetime
import time


class Upperlevel:
    def __repr__(self):        
        return self.__class__.__name__ + ': ' + '(' +', '.join([n for n in dir(self) if not n.startswith('__')]) +')'
   
class Type_analysis(Upperlevel):
    def __init__(self,dictionary):
        self.Trim = dictionary[0]['A']
        self.AlignPhix = dictionary[0]['B']
        self.AlignGenome = dictionary[0]['C']
        self.IIDefinition = dictionary[0]['D']
        self.GeneDefinition = dictionary[0]['E']
        self.GroupAnalysis = dictionary[0]['F']
        self.Tables = dictionary[0]['G']
        self.Design = dictionary[0]['H']
                       
class General(Upperlevel):
    def __init__(self,dictionary):                      
        self.operator = dictionary[1]['A']
        self.storing_loc = dictionary[1]['B']
        self.pair_ends = True if (dictionary[2]['B']+dictionary[3]['D']+dictionary[4]['E']+dictionary[5]['F']+dictionary[6]['G']) >= 1 else False
        self.date = (datetime.date.today()).isoformat()

class QS (Upperlevel):
    def __init__(self,first_base, sequence_3end):
        self.first_base = first_base
        self.sequence_3end = sequence_3end

class Trim(Upperlevel):
    def __init__(self,dictionary,pair_ends):         
        self.lib_numb = dictionary[2]['A']
        self.lib_names = dictionary[2]['C']
        self.input_files_I = dictionary[2]['D']
        self.p7adaptor = dictionary[2]['F']
        if pair_ends:
            self.input_files_II = dictionary[2]['E']
            self.p5adaptor = dictionary[2]['G']    
        self.QS = QS(dictionary[2]['H'],dictionary[2]['J'])
        self.store = dictionary[2]['K']       
        
class AlignPhix(Upperlevel):
    def __init__(self,dictionary,pair_ends,starting):                 
        self.reference = dictionary[3]['A']
        self.store = dictionary[3]['B']
        if starting:
            self.lib_numb = dictionary[3]['C']
            self.lib_names = dictionary[3]['E']
            self.input_files_I = dictionary[3]['F']
            if pair_ends:
                self.input_files_II = dictionary[3]['G']
        else:
            self.lib_numb = int()
            self.lib_names = []
            self.input_files_I = []
            if pair_ends:
                self.input_files_II = []
                
class AlignGenome(Upperlevel):
    def __init__(self,dictionary,pair_ends,starting):
        self.program_type = dictionary[4]['A']
        self.reference = dictionary[4]['B']
        self.store = dictionary[4]['C']
        if starting:
            self.lib_numb = dictionary[4]['D']
            self.lib_names = dictionary[4]['F']
            self.input_files_I = dictionary[4]['G']
            if pair_ends:
                self.input_files_II = dictionary[4]['H']
        else:
            self.lib_numb = int()
            self.lib_names = []
            self.input_files_I = []
            if pair_ends:
                self.input_files_II = []
                
class IIDefinition(Upperlevel):
    def __init__(self,dictionary,pair_ends,starting):  
        if starting:
            self.lib_numb = dictionary[5]['E']
            self.lib_names = dictionary[5]['G']
            self.input_files = dictionary[5]['H'] 
        else:
            self.lib_numb = int()
            self.lib_names = []
            self.input_files = []    
        self.ins_iv = dictionary[5]['A']
        if pair_ends:
            self.reads_duplicate = dictionary[5]['B']
        else:
            self.reads_duplicate = False
        self.fidelity = dictionary[5]['C']
        if self.fidelity:
            self.fidelity_limit = dictionary[5]['D']
        else:
            self.fidelity_limit = 0

class Parameters (Upperlevel):
    def __init__(self,II,KI,Bias,Reads):
        self.II = II
        self.KI = KI
        self.Bias = Bias
        self.Reads = Reads

class GeneDefinition (Upperlevel):
    def __init__(self,dictionary,starting):
        self.reference = dictionary[6]['A']
        self.Parameters = Parameters(II = dictionary[6]['B'],KI = dictionary[6]['C'],Bias = dictionary[6]['D'],Reads =dictionary[6]['E'])
        if starting:
            self.lib_numb = dictionary[6]['F']
            self.lib_names = dictionary[6]['H']
            self.input_files = dictionary[6]['I']
        else:
            self.lib_numb = int()
            self.reads_norm = False
            self.lib_names = []
            self.input_files = []
            
class Outlier (Upperlevel):
    def __init__(self,dictionary):
        self.perform = dictionary[7]['J']
        if self.perform:
            self.Parameters = Parameters(II = dictionary[7]['K'],KI = dictionary[7]['L'],Bias = dictionary[7]['M'],Reads =dictionary[7]['N'])
            self.Limits = Parameters(II = dictionary[7]['O'],KI = dictionary[7]['P'],Bias = dictionary[7]['Q'],Reads =dictionary[7]['R'])
                
        
class Experiments (Upperlevel):
    def __init__(self,name,experiments):
        self.name = name
        self.experiments = experiments
            
class GroupAnalysis (Upperlevel):
    def __init__(self,dictionary,starting):  
        self.group_numb = dictionary[7]['A']
        self.storage_loc = ''
        self.Reference = Experiments(dictionary[7]['B'],dictionary[7]['C'])
        self.Others = Experiments(dictionary[7]['D'],dictionary[7]['E'])
        self.Parameters = Parameters(II = dictionary[7]['F'],KI = dictionary[7]['G'],Bias = dictionary[7]['H'],Reads =dictionary[7]['I'])
        self.Outlier = Outlier(dictionary)
        if starting: 
            self.lib_numb = dictionary[7]['S']
            self.lib_names = dictionary[7]['T']
            self.input_files = dictionary[7]['U']
        else:
            self.lib_numb = int()
            self.lib_names = []
            self.input_files = []           

class Tables(Upperlevel):    
    def __init__(self,dictionary,starting):
        self.names = dictionary[8]['A']
        self.keys = dictionary[8]['B']
        self.filters = dictionary[8]['C']
        if starting:
            self.input_files = dictionary[8]['D']
        else:
            self.input_files = ''

            
class Gene_model(Upperlevel):
    def __init__ (self,perform,reference,genes):
        self.perform = perform
        self.reference = reference
        self.genes = []
        for n in genes.split(','):
            self.genes.append(n.rstrip(' ').lstrip(' '))

class Insertions_distribution(Upperlevel):
    def __init__ (self,perform,outlier,outlier_value,annotate,annotate_value):
        self.perform = perform
        self.outlier = outlier
        self.outlier_value = outlier_value
        self.annotate = annotate
        self.annotate_value = annotate_value
        
            
class Design (Upperlevel):
    def __init__ (self,dictionary,starting):
        self.Gene_model = Gene_model(dictionary[9]['A'],dictionary[9]['C'],dictionary[9]['D'])
        self.Distribution = Insertions_distribution(dictionary[9]['B'],dictionary[9]['E'],dictionary[9]['F'],dictionary[9]['G'],dictionary[9]['H']) 
        if starting:
            self.input_files = dictionary[9]['I']
        else:
            self.input_files = ''
        
class Info(Upperlevel):
    def __init__(self, dictionary):
        self.Type_analysis = Type_analysis(dictionary)
        self.General = General(dictionary)
            
    def upload_informations(self,dictionary):
        if self.Type_analysis.Trim:
            self.Trim = Trim(dictionary,self.General.pair_ends)
        if self.Type_analysis.AlignPhix:
            self.AlignPhix = AlignPhix(dictionary,self.General.pair_ends,starting = not(self.Type_analysis.Trim))
        if self.Type_analysis.AlignGenome:
            self.AlignGenome = AlignGenome (dictionary,self.General.pair_ends,starting = False if (self.Type_analysis.Trim + self.Type_analysis.AlignPhix) == 1 else True)
        if self.Type_analysis.IIDefinition:  
            self.IIDefinition = IIDefinition (dictionary,self.General.pair_ends, starting = not(self.Type_analysis.AlignGenome))
        if self.Type_analysis.GeneDefinition:  
            self.GeneDefinition = GeneDefinition(dictionary, starting = not(self.Type_analysis.IIDefinition))
        if self.Type_analysis.GroupAnalysis:
            self.GroupAnalysis = GroupAnalysis(dictionary, starting = not(self.Type_analysis.GeneDefinition))
        if self.Type_analysis.Tables:
            self.Tables = Tables (dictionary, starting = not(self.Type_analysis.GroupAnalysis))
        if self.Type_analysis.Design:
            self.Design = Design (dictionary, starting = not(self.Type_analysis.GroupAnalysis))
            
    def starting (self):
        if self.Type_analysis.Trim:
            self.Type_analysis.starting = self.Trim
        elif self.Type_analysis.AlignPhix:
            self.Type_analysis.starting = self.AlignPhix
        elif self.Type_analysis.AlignGenome:
            self.Type_analysis.starting = self.AlignGenome
        elif self.Type_analysis.IIDefinition:
            self.Type_analysis.starting = self.IIDefinition
        elif self.Type_analysis.GeneDefinition:
            self.Type_analysis.starting = self.GeneDefinition
        elif self.Type_analysis.GroupAnalysis:
            self.Type_analysis.starting = self.GroupAnalysis
        elif self.Type_analysis.Tables:
            self.Type_analysis.starting = self.Tables
        elif self.Type_analysis.Design:
            self.Type_analysis.starting = self.Design
            
            
    def to_do  (self):
        to_do = []
        if self.Type_analysis.Trim:
            if self.Type_analysis.starting != self.Trim:
                to_do.append(self.Trim)
        if self.Type_analysis.AlignPhix:
            if self.Type_analysis.starting != self.AlignPhix:
                to_do.append(self.AlignPhix)
        if self.Type_analysis.AlignGenome:
            if self.Type_analysis.starting != self.AlignGenome:
                to_do.append(self.AlignGenome)
        if self.Type_analysis.IIDefinition:
            if self.Type_analysis.starting != self.IIDefinition:
                to_do.append(self.IIDefinition)
        if self.Type_analysis.GeneDefinition:
            if self.Type_analysis.starting != self.GeneDefinition:
                to_do.append(self.GeneDefinition)
        if self.Type_analysis.GroupAnalysis:
            if self.Type_analysis.starting != self.GroupAnalysis:
                to_do.append(self.GroupAnalysis)
        if self.Type_analysis.Tables:
            if self.Type_analysis.starting != self.Tables:
                to_do.append(self.Tables)
        if self.Type_analysis.Design:
            if self.Type_analysis.starting != self.Design:
                to_do.append(self.Design)
                      
        self.Type_analysis.to_do = to_do
        
                
        
        
            
    def make_folders(self):
        """Function to create folders were to store subsequent data:
        After checking if the storing_location does not exist (and eventually asking for a new name), creates:
        the main folder: self.General.storing_loc + '/' exp_name + '_' + self.General.date
        the graph folder: main folder +'/graph/'
        the raw folder: main folder +'/raw/' """
        
        do_folder_experiments = True
        if self.Type_analysis.GroupAnalysis:
            if self.Type_analysis.starting == self.GroupAnalysis:
                do_folder_experiments = False
        if self.Type_analysis.Tables:
            if self.Type_analysis.starting == self.Tables:
                do_folder_experiments = False
        if self.Type_analysis.Design:
            if self.Type_analysis.starting == self.Design:
                do_folder_experiments = False
        
        if do_folder_experiments: 
            for exp in self.Type_analysis.starting.lib_names:
                location = os.path.join(self.General.storing_loc,exp + '_' +self.General.date)   
                while os.path.isdir(location):
                    print 'The directory %s already exists: '% location
                    new_name = raw_input('Change the experiment name (the previous one was %s): ' % exp )
                    new_loc = os.path.join(self.General.storing_loc,new_name + '_' +self.General.date)
                    while os.path.isdir(new_loc):
                        new_name = raw_input('The directory %s already exists. Change the experiment name (the previous one was %s): '%(new_loc,new_name))
                        new_loc = os.path.join(self.General.storing_loc,new_name + '_' +self.General.date)
                    location = new_loc
                    self.Type_analysis.starting.lib_names = [new_name if x == exp else x for x in self.Type_analysis.starting.lib_names]
                                        
                os.mkdir(location)
                graph_loc = os.path.join(location, "graph")
                os.mkdir(graph_loc)
                
                raw = os.path.join(location, "raw")
                os.mkdir(raw)
                        
        if self.Type_analysis.GroupAnalysis:            
            analysis_main = os.path.join(self.General.storing_loc,'Analysis')  
            if not os.path.isdir(analysis_main):
               os.mkdir(analysis_main)
            analysis_loc = os.path.join(analysis_main, self.General.date) 
            while os.path.isdir(analysis_loc):
                print 'The directory %s already exists: '% analysis_loc
                new_name = raw_input('Change the experiment name (the previous one was %s): ' %  self.General.date)
                analysis_loc = os.path.join(analysis_main,new_name)
                while os.path.isdir(analysis_loc):
                    new_name = raw_input('The directory %s already exists. Change the experiment name (the previous one was %s): '%(analysis_loc,new_name))
                    analysis_loc = os.path.join(analysis_main,new_name)
            self.GroupAnalysis.storage_loc = analysis_loc
            os.mkdir(analysis_loc)
            raw_loc = os.path.join(analysis_loc, "raw")
            os.mkdir(raw_loc)
            graph_loc = os.path.join(analysis_loc, "graph")
            os.mkdir(graph_loc)
            
    def print_info (self):
        print "\n\n***\tInformations on experiments' analysis\t***"
        print 'Generals'
        print '\t{:20s}:\t'.format('Operator') + '%s'%self.General.operator
        print '\t{:20s}:\t'.format('Date') + '%s'%self.General.date
        print '\t{:20s}:\t'.format('Storing location') + '%s'%self.General.storing_loc
        print '\t{:20s}:\t'.format('Pair ends') + '%s'% self.General.pair_ends
        print 'Type of analysis'
        print '\t{:20s}:\t'.format('Trim') + '%s' % self.Type_analysis.Trim
        print '\t{:20s}:\t'.format('AlignPhix') + '%s'%self.Type_analysis.AlignPhix
        print '\t{:20s}:\t'.format('AlignGenome') + '%s'%self.Type_analysis.AlignGenome
        print '\t{:20s}:\t'.format('IIDefinition') + '%s'%self.Type_analysis.IIDefinition
        print '\t{:20s}:\t'.format('GeneDefinition') + '%s'%self.Type_analysis.GeneDefinition
        print '\t{:20s}:\t'.format('GroupAnalysis') + '%s'%self.Type_analysis.GroupAnalysis
        print '\t{:20s}:\t'.format('Tables') + '%s'%self.Type_analysis.Tables
        if self.Type_analysis.Trim:
            print 'Trim'
            print '\t{:20s}:\t'.format('Libraries number') + '%s' % self.Trim.lib_numb
            print '\t{:20s}:\t'.format('Libraries name') + '%s' % '\n\t{:20s}\t'.format('').join(self.Trim.lib_names)
            print '\t{:20s}:\t'.format('Input') + '%s' % '\n\t{:20s}\t'.format('').join(self.Trim.input_files_I)
            if self.General.pair_ends:
                print '\t{:20s}:\t'.format('Input_pairends') + '%s' % '\n\t{:20s}\t'.format('').join(self.Trim.input_files_II)
            print '\t{:20s}:\t'.format('P7 adaptor') + '%s' % self.Trim.p7adaptor
            if self.General.pair_ends:
                print '\t{:20s}:\t'.format('P5 adaptor') + '%s' % self.Trim.p5adaptor
            print '\tQS parmeters'
            print '\t\t{:20s}:\t'.format('First base Q') + '%s' % self.Trim.QS.first_base
            print '\t\t{:20s}:\t'.format('Extension limit') + '%s' % self.Trim.QS.sequence_3end
            print '\t{:20s}:\t'.format('Permanently store') + '%s' % self.Trim.store
            
        if self.Type_analysis.AlignPhix:
            print 'AlignPhix'
            print '\t{:20s}:\t'.format('Phix reference') + '%s' % self.AlignPhix.reference
            print '\t{:20s}:\t'.format('Permanently store') + '%s' % self.AlignPhix.store
            print '\t{:20s}:\t'.format('Libraries number') + '%s' % self.AlignPhix.lib_numb
            print '\t{:20s}:\t'.format('Libraries name') + '%s' % '\n\t{:20s}\t'.format('').join(self.AlignPhix.lib_names)
            print '\t{:20s}:\t'.format('Lib input') + '%s' % '\n\t{:20s}\t'.format('').join(self.AlignPhix.input_files_I)
            if self.General.pair_ends:
                print '\t{:20s}:\t'.format('Lib input pair_end') + '%s' % '\n\t{:20s}\t'.format('').join(self.AlignPhix.input_files_II)
            
        if self.Type_analysis.AlignGenome:
            print 'AlignGenome'
            print '\t{:20s}:\t'.format('Alignment Program') + '%s' % self.AlignGenome.program_type
            print '\t{:20s}:\t'.format('Reference') + '%s' % self.AlignGenome.reference
            print '\t{:20s}:\t'.format('Permanently store') + '%s' % self.AlignGenome.store
            print '\t{:20s}:\t'.format('Libraries number') + '%s' % self.AlignGenome.lib_numb
            print '\t{:20s}:\t'.format('Libraries name') + '%s' % '\n\t{:20s}\t'.format('').join(self.AlignGenome.lib_names)
            print '\t{:20s}:\t'.format('Lib input') + '%s' % '\n\t{:20s}\t'.format('').join(self.AlignGenome.input_files_I)
            if self.General.pair_ends:
                print '\t{:20s}:\t'.format('Lib input pair_end') + '%s' % '\n\t{:20s}\t'.format('').join(self.AlignGenome.input_files_II)
        
        if self.Type_analysis.IIDefinition:
            print 'I.I. Definition'
            print '\t{:20s}:\t'.format('Libraries number') + '%s' % self.IIDefinition.lib_numb
            print '\t{:20s}:\t'.format('Libraries name') + '%s' % '\n\t{:20s}\t'.format('').join(self.IIDefinition.lib_names)
            print '\t{:20s}:\t'.format('Lib input') + '%s' % '\n\t{:20s}\t'.format('').join(self.IIDefinition.input_files)
            print '\t{:20s}:\t'.format('Bases to define I.I.') + '%s' % self.IIDefinition.ins_iv
            print '\t{:20s}:\t'.format('Alignment fidelity limit') + '%s' % self.IIDefinition.fidelity
            if self.IIDefinition.fidelity:
                print '\t{:20s}:\t'.format('Fidelity limit value') + '%s' % self.IIDefinition.fidelity_limit
            
            if self.General.pair_ends:
                print '\t{:20s}:\t'.format('Remove duplicates') + '%s' % self.IIDefinition.reads_duplicate

        if self.Type_analysis.GeneDefinition:
            print 'Genes Definition'
            print '\t{:20s}:\t'.format('Genes annotation') + '%s' % self.GeneDefinition.reference
            print '\t{:20s}'.format('Parameters')
            print '\t\t{:8s}:\t'.format('II') + '%s' % self.GeneDefinition.Parameters.II
            print '\t\t{:8s}:\t'.format('KI') + '%s' % self.GeneDefinition.Parameters.KI
            print '\t\t{:8s}:\t'.format('Bias') + '%s' % self.GeneDefinition.Parameters.Bias
            print '\t\t{:8s}:\t'.format('Reads') + '%s' % self.GeneDefinition.Parameters.Reads
            print '\t{:20s}:\t'.format('Libraries number') + '%s' % self.GeneDefinition.lib_numb
            print '\t{:20s}:\t'.format('Libraries name') + '%s' % '\n\t{:20s}\t'.format('').join(self.GeneDefinition.lib_names)
            print '\t{:20s}:\t'.format('Lib input') + '%s' % '\n\t{:20s}\t'.format('').join(self.GeneDefinition.input_files)
            
        if self.Type_analysis.GroupAnalysis:
            print 'Group analysis'
            print '\t{:20s}:\t'.format('Storage location') + '%s' % self.GroupAnalysis.storage_loc
            print '\t{:20s}:\t'.format('Libraries number') + '%s' % self.GroupAnalysis.lib_numb
            print '\t{:20s}:\t'.format('Libraries name') + '%s' % '\n\t{:20s}\t'.format('').join(self.GroupAnalysis.lib_names)
            print '\t{:20s}:\t'.format('Lib input') + '%s' % '\n\t{:20s}\t'.format('').join(self.GroupAnalysis.input_files)
            print '\t{:20s}:\t'.format('Groups number') + '%s' % self.GroupAnalysis.group_numb
            print '\tReference exp.'
            print '\t\t{:20s}:\t'.format('Name') + '%s' % self.GroupAnalysis.Reference.name
            print '\t\t{:20s}:\t'.format('Experiments') + '%s' % '\n\t\t{:20s}\t'.format('').join(self.GroupAnalysis.Reference.experiments)
            print '\tOther groups'
            for name in self.GroupAnalysis.Others.name:
                print '\t\t{:20s}:\t'.format('Name') + '%s' % name 
                print '\t\t{:20s}:\t'.format('Experiments') + '%s' % '\n\t\t{:20s}\t'.format('').join(self.GroupAnalysis.Others.experiments[self.GroupAnalysis.Others.name.index(name)])
            print '\t{:20s}'.format('Parameters')
            print '\t\t{:8s}:\t'.format('II') + '%s' % self.GroupAnalysis.Parameters.II
            print '\t\t{:8s}:\t'.format('KI') + '%s' % self.GroupAnalysis.Parameters.KI
            print '\t\t{:8s}:\t'.format('Bias') + '%s' % self.GroupAnalysis.Parameters.Bias
            print '\t\t{:8s}:\t'.format('Reads') + '%s' % self.GroupAnalysis.Parameters.Reads
            print '\t{:20s}:\t'.format('Outlier analysis') + '%s' % self.GroupAnalysis.Outlier.perform
            if self.GroupAnalysis.Outlier.perform:
                print '\t\t{:20s}'.format('Parameters')
                print '\t\t\t{:8s}:\t'.format('II') + '%s' % self.GroupAnalysis.Outlier.Parameters.II,
                if self.GroupAnalysis.Outlier.Parameters.II:
                    print '=> Value limit: %s' % self.GroupAnalysis.Outlier.Limits.II
                else:
                    print ''
                print '\t\t\t{:8s}:\t'.format('KI') + '%s' % self.GroupAnalysis.Outlier.Parameters.KI,
                if self.GroupAnalysis.Outlier.Parameters.KI:
                    print '=> Value limit: %s' % self.GroupAnalysis.Outlier.Limits.KI
                else:
                    print ''
                print '\t\t\t{:8s}:\t'.format('Bias') + '%s' % self.GroupAnalysis.Outlier.Parameters.Bias,
                if self.GroupAnalysis.Outlier.Parameters.Bias:
                    print '=> Value limit: %s' % self.GroupAnalysis.Outlier.Limits.Bias
                else:
                    print ''
                print '\t\t\t{:8s}:\t'.format('Reads') + '%s' % self.GroupAnalysis.Outlier.Parameters.Reads,
                if self.GroupAnalysis.Outlier.Parameters.Reads:
                    print '=> Value limit: %s' % self.GroupAnalysis.Outlier.Limits.Reads
                else:
                    print ''
        if self.Type_analysis.Tables:
            print 'Tables'
            print '\t{:20s}:\t'.format('Data_reference') + self.Tables.input_files
            for table in self.Tables.names:
                position = self.Tables.names.index(table)
                print '\t{:20s}:\t'.format('Name') + table
                print '\t\t{:20s}:\t\t'.format('Keys') + str(self.Tables.keys[position]).strip('[]').replace("'","").replace(" ","")
                filter_list = []
                if position < len(self.Tables.filters):       
                    for n in self.Tables.filters[position]:
                        parameter = str(n['parameter']).strip('[]').replace("'","").replace(" ","")
                        if n.has_key('number'):
                            filter_list.append('%s%s%i' % (parameter,n['operation'],n['number']))
                        else:
                            filter_list.append('%s%s' % (parameter,n['operation']))
                print '\t\t{:20s}:\t\t'.format('Filters') + ", ".join(filter_list)
        if self.Type_analysis.Design:
            print 'Gene rappresentation'
            print '\t{:20s}:\t'.format('Data_reference') + self.Design.input_files 
            print '\t{:20s}:\t'.format('Type of plot')
            print '\t\t{:20s}:\t'.format('I.I. in gene models') + '%s' % self.Design.Gene_model.perform
            print '\t\t{:20s}:\t'.format('Genes distribution') + '%s' % self.Design.Distribution.perform
            if self.Design.Gene_model.perform:
                print '\t{:20s}:\t'.format('I.I. in gene models')
                print '\t{:20s}:\t'.format('Genes annotation') + '%s' % self.Design.Gene_model.reference
                print '\t{:20s}:\t'.format('Genes (or intervals)')
                for name in self.Design.Gene_model.genes:
                    print '\t\t%s' % name
            if self.Design.Distribution.perform:
                print '\t{:20s}:\t'.format('Genes distribution')
                print '\t\t{:20s}:\t'.format('Mark Outliers') + '%s' % self.Design.Distribution.outlier
                if self.Design.Distribution.outlier:
                    print '\t\t\t{:20s}:\t'.format('Outlier value') + '%s' % self.Design.Distribution.outlier_value
                print '\t\t{:20s}:\t'.format('Annotate genes') + '%s' % self.Design.Distribution.annotate
                if self.Design.Distribution.annotate:
                    print '\t\t\t{:20s}:\t'.format('Value or Gene list') + '%s' % self.Design.Distribution.annotate_value
                    
            
    def fill_up(self):
        def fill_up_categories(self,category,terminal_string,pair_end = True):
            category.lib_numb = self.Type_analysis.starting.lib_numb
            category.lib_names = self.Type_analysis.starting.lib_names 
            if pair_end:
                input_ =[]
                for exp in self.Type_analysis.starting.lib_names:
                    file_input = os.path.join(self.General.storing_loc,exp + '_' +self.General.date,'raw',exp + terminal_string)
                    input_.append(file_input)
                category.input_files_I = input_             
                if self.General.pair_ends:
                    input_ = []
                    for exp in self.Type_analysis.starting.lib_names:
                        file_input = os.path.join(self.General.storing_loc,exp + '_' +self.General.date,'raw',exp + '_pairend' + terminal_string)
                        input_.append(file_input)
                    category.input_files_II = input_
            else:
                input_ =[]
                for exp in self.Type_analysis.starting.lib_names:
                    file_input = os.path.join(self.General.storing_loc,exp + '_' +self.General.date,'raw',exp + terminal_string)
                    input_.append(file_input)
                category.input_files = input_             
                
        
        for category in self.Type_analysis.to_do:
            if self.Type_analysis.AlignPhix:
                if category == self.AlignPhix:
                    fill_up_categories(self,category,'_QS.fastq',pair_end = True)
            if self.Type_analysis.AlignGenome:                    
                if category == self.AlignGenome:
                    if self.Type_analysis.AlignPhix:
                        fill_up_categories(self,category,'_PhixCleaned.fastq',pair_end = True)
                    else:
                        fill_up_categories(self,category,'_QS.fastq',pair_end = True)
            if self.Type_analysis.IIDefinition:
                if category == self.IIDefinition:
                    fill_up_categories(self,category,'_Aligned.sam',pair_end = False)
            if self.Type_analysis.GeneDefinition:
                if category == self.GeneDefinition:
                    fill_up_categories(self,category,'_IIRawdata.pkl',pair_end = False)
            if self.Type_analysis.GroupAnalysis:
                if category == self.GroupAnalysis:
                    fill_up_categories(self,category,'_GenesData.pkl',pair_end = False)
            if self.Type_analysis.Tables:
                if category == self.Tables:
                    file_input = os.path.join(self.GroupAnalysis.storage_loc,'raw','GroupAnalysis.pkl')
                    self.Tables.input_files = file_input
            if self.Type_analysis.Design:
                if category == self.Design:
                    file_input = os.path.join(self.GroupAnalysis.storage_loc,'raw','GroupAnalysis.pkl')
                    self.Design.input_files = file_input
                    
    def print_save(self,exp,string):
        """Print the selected string on the screen and store the information in txt file in the main folder"""
        print string
        with open (os.path.join(self.General.storing_loc,exp + '_' + self.General.date,exp +'_info.txt'), 'a' ) as write:
                print >> write, string
         
def read_txt(informations,text):
    
    def get_TRUE(section,task,line,informations):    
        if not informations.has_key(section):
            informations[section] = {}
        if re.search('^\s*@\d[A-Z]\).*([YyNn]).*',line):
            value = re.findall ('^\s*@\d[A-Z]\).*([YyNn]).*',line)[0]
            if value.lower() == 'y':
                informations[section][task] = True
            else:
                informations[section][task] = False
        else:
            informations[section][task] = False
        return informations
        
    def get_NUMBER(section,task,line,informations):
        if not informations.has_key(section):
            informations[section] = {}
        if re.search ('^\s*@\d[A-Z]\)\D*([\d.]+)\D*',line):
            value = re.findall ('^\s*@\d[A-Z]\)\D*([\d.]+)\D*',line)[0]
            if value:
                informations[section][task] = float(value)
            else:
                informations[section][task] = 'n.d'
        else:
             informations[section][task] = 'n.d'
        return informations
    
    def get_STRING(section,task,line,informations):
        if not informations.has_key(section):
            informations[section] = {}
    
        string = re.findall ('^\s*@\d[A-Z]\)\s*(.*)',line)[0]
        if string:
            informations[section][task] = string
        else:
            informations[section][task] = 'n.d.'
        return informations
        
    def get_LIST(section,task,line,informations):
        if not informations.has_key(section):
            informations[section] = {}
        if not informations[section].has_key(task):
            informations[section][task] = []
        string = re.findall ('^\s*@\d[A-Z]\)\s*(.*)',line)[0]
        if string:
            informations[section][task].append(string)
        return informations
        
    def extract_line (section,task,line,informations,true=[],number=[],list_=[],string=[]):
        if task in number:
            informations = get_NUMBER(section,task,line,informations)    
        elif task in true:
            informations = get_TRUE(section,task,line,informations)
        elif task in list_:
            informations = get_LIST(section,task,line,informations)
        elif task in string:  
            informations = get_STRING(section,task,line,informations)
        return informations    
    
    hand = open(text,'rb')
    line_count = 0
    for line in hand:
        
        line_count +=1
        line = line.rstrip()
        find = re.findall('^\s*@(\d)([A-Z])\)',line)
        
        if find:
            section, task = find[0]
            section = int(section)


            if section == 0:
                informations = get_TRUE(section,task,line,informations)
                
            elif section == 1:
               informations = get_STRING(section,task,line,informations) 
            
            elif section == 2:
                informations = extract_line (section,task,line,informations,true =['B','K'],number =['A','H','J'],list_=['C','D','E'],string=['F','G'])
            
            elif section == 3:        
                informations = extract_line (section,task,line,informations,true =['B','D'],number =['C'],list_=['E','F','G'],string=['A'])
                            
            elif section == 4:        
                informations = extract_line (section,task,line,informations,true =['C','E'],number =['D'],list_=['F','G','H'],string=['A','B'])  
             
            elif section == 5:        
                informations = extract_line (section,task,line,informations,true =['B','C','F'],number =['A','D','E'],list_=['G','H'],string=[])  
        
            elif section == 6:        
                informations = extract_line (section,task,line,informations,true =['B','C','D','E','G'],number =['F'],list_=['H','I'],string=['A'])  
            
            elif section == 7:
                informations = extract_line (section,task,line,informations,true =['F','G','H','I','J','K','L','M','N'],number =['A','P','O','P','Q','R','S'],list_=['D','T','U'],string =['B'])  
                if task in ['C']:
                    if re.search('^\s*@\d[A-Z]\)\s*(\S+.*\S*)', line):
                        string = [n.lstrip().rstrip() for n in re.findall('^\s*@\d[A-Z]\)\s*(\S+.*\S*)', line)[0].split(',')]
                        informations[section][task] = string
                
                if task in ['E']:
                    if re.search('^\s*@\d[A-Z]\)\s*(\S+.*\S*)', line):
                        if not informations[section].has_key(task):
                            informations[section][task] =[]
                        string = [n.lstrip().rstrip() for n in re.findall('^\s*@\d[A-Z]\)\s*(\S+.*\S*)', line)[0].split(',')]
                        informations[section][task].append(string)
            elif section == 8:
                informations = extract_line (section,task,line,informations,list_=['A'],string =['D'])
                if task == 'B':
                    if not informations[section].has_key(task):
                        informations[section][task] =[]    
                    key = []
                    string = re.findall('\((.*?)\)',line)
                    for param in string:
                        variable = tuple(n.lstrip().rstrip() for n in param.split(',') if n != '')
                        if len(variable) == 3:
                            key.append(variable)
                        elif variable[1] == 'Score' or variable[1] == 'Outliers':
                            key.append(variable)
                        else:
                            print "ERROR! Line %d Section %s Task %d: The table parameter '%s' doesn't contain the correct structure" %(line_count,section,task,param)
                    informations[section][task].append(key)
                elif task == 'C':
                    if not informations[section].has_key(task):
                        informations[section][task] =[]
                    filters = []
                    items = re.findall('[fF]ilter\s*?\[(.*?)\]',line)
                    for q in items:
                        filter_ = {}
                        if re.search('\((.*?)\)\s*?,\s*?(ascending|descending)',q):
                            parameters = re.findall('\((.*?)\)\s*?,\s*?(ascending|descending)',q)[0]
                            filter_['operation'] = parameters[1]
                            variable = tuple(n.lstrip().rstrip() for n in parameters[0].split(',') if n != '')
                            if len(variable) == 3:
                                filter_['parameter']= variable
                            elif variable[1] == 'Score':
                                filter_['parameter']= variable
                            else:
                                print "ERROR! Line %i Section %i Task %s: The table parameter '%s' doesn't contain the correct structure" %(line_count,section,task,q)
                                continue
                        elif re.search('\((.*?)\)\s*,\s*(\S+)\s*,\s*(\S+)',q):
                            parameters = re.findall('\((.*?)\)\s*,\s*(\S+)\s*,\s*(\S+)',q)[0]
                            filter_['operation'] = parameters[1]
                            filter_['number'] = float(parameters[2])
                            variable = tuple(n.lstrip().rstrip() for n in parameters[0].split(',') if n != '')
                            if len(variable) == 3:
                                filter_['parameter']= variable
                            else:
                                print "ERROR! Line %d Section %d Task %s: The table parameter '%s' doesn't contain the correct structure" %(line_count,section,task,q)
                                continue                        
                        else: 
                            print "ERROR! Line %d Section %d Task %s: The table parameter '%s' doesn't contain the correct structure" %(line_count,section,task,q)
                            continue
                        filters.append(filter_)
                    informations[section][task].append(filters)
        
            elif section == 9:
                informations = extract_line (section,task,line,informations,true =['A','B','E','G'],number =['F'],string=['C','D','H','I'])           
    
    return informations

        
       
       

                
                
                
            
            
            
