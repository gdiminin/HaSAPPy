#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 19:26:59 2016

@author: GDM
"""

import HaSAPPy.INFOloads as INFOloads
import argparse

parser = argparse.ArgumentParser(description='Launching command of HaSAPPy program')
print
print "***********************************************************"
print "***HaSAPPy: Haploid Screening Analysis Package in Python***"
print '***********************************************************\n\n'


parser.add_argument('path', help = 'Provide PATH of LoadModule file to start analysis. For more details visit HaSAPPy webpage on "https://github.com/gdiminin/HaSAPPy" ', action="store")


text = parser.parse_args()

if text.path == None:
    print '\nWARNING: informations provided are not sufficent.\nCheck -h option to have more details on requested parameters'


informations =  {}
informations = INFOloads.read_txt(informations,text.path)
analysis = INFOloads.Info(informations)
analysis.upload_informations(informations)
analysis.starting()
analysis.to_do()
analysis.make_folders()
analysis.fill_up()
analysis.print_info()

if analysis.Type_analysis.Trim:
    import HaSAPPy.Trim as Trim
    Trim.load(analysis)

if analysis.Type_analysis.AlignPhix or analysis.Type_analysis.AlignGenome:
    import HaSAPPy.Align as Align
    Align.load(analysis)

if analysis.Type_analysis.IIDefinition:
    import HaSAPPy.IIDefinition as IIDefinition
    IIDefinition.load(analysis)

if analysis.Type_analysis.GeneDefinition:
    import HaSAPPy.GeneDefinition as GeneDefinition
    GeneDefinition.library_analysis(analysis)

if analysis.Type_analysis.GroupAnalysis:
    import HaSAPPy.GroupAnalysis as GroupAnalysis
    GroupAnalysis.performing_analysis(analysis)

if analysis.Type_analysis.Tables:
    import HaSAPPy.Tables as Tables
    Tables.main(analysis)

if analysis.Type_analysis.Design:
    import HaSAPPy.DesignGeneInsertion as Design
    Design.start(analysis)





    
    

