#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 19:26:59 2016

@author: GDM
"""

import INFOloads
import argparse

parser = argparse.ArgumentParser(description='Launching command of HaSAPPy program')
print
print "***********************************************************"
print "***HaSAPPy: Haploid Screening Analysis Package in Python***"
print '***********************************************************\n\n'


parser.add_argument('path', help = 'Provide PATH of LoadModule file to start analysis. For more details visit HaSAPPy webpage on "https://github.com/gdiminin/HaSAPPy" ', action="store")

text = parser.parse_args()
if text == None:
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
    import Trim
    Trim.load(analysis)

if analysis.Type_analysis.AlignPhix or analysis.Type_analysis.AlignGenome:
    import Align
    Align.load(analysis)

if analysis.Type_analysis.IIDefinition:
    import IIDefinition
    IIDefinition.load(analysis)

if analysis.Type_analysis.GeneDefinition:
    import GeneDefinition
    GeneDefinition.library_analysis(analysis)

if analysis.Type_analysis.GroupAnalysis:
    import GroupAnalysis
    GroupAnalysis.performing_analysis(analysis)

if analysis.Type_analysis.Tables:
    import Tables
    Tables.main(analysis)

if analysis.Type_analysis.Design:
    import DesignGeneInsertion as Design
    Design.start(analysis)





    
    

