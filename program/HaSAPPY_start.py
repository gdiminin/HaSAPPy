# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 19:26:59 2016

@author: GDM
"""

import INFOloads
import argparse

parser = argparse.ArgumentParser(description='Lunching command to HaSA program')

parser.add_argument('path', action="store")
text = (parser.parse_args())

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

if analysis.Type_analysis.Allign:
    import Allign
    Allign.load(analysis)

if analysis.Type_analysis.IIDefinition:
    import IIDefinition
    IIDefinition.main(analysis)

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





    
    

