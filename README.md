#HaSAPPy (Haploid Screening Analysis Python Package)



##INTRODUCTION

Forward genetic screens represent one of the most important approaches to unveil the complexity of cellular mechanisms. The derivation of different mammalian haploid cell models (from tumors and from embryos) opened the possibility to perform forward genetic screens using insertional mutagenesis in mammals.
HaSAPPy software can be used to identify insertion locations in the whole genome, map them at the level of genes, and classify insertions according to their effects on gene expression. The code conforms to current Python programming guidelines and can be freely adapted and extended according to user needs.

**Modules:**
 - Trim adaptor and sequence by quality
 - Align to reference genome
 - Identify Independent Insertions (I.I.)
 - Classify insertions at the genes level
 - Enrichment analysis
 - Data presentation


##REQUISITES

For the correct functionality of HaSAPPy program the following Python Packages are necessary and requested to download/update using pip -instal or pip install --upgrade command:
 - numpy
 - HTSeq
 - matplotlib
 - pandas
 - scipy
 - sklearn
 - xlsxwriter
 
Moreover, alignment module requires installation in PATH of the following programs:
- Bowtie2: to remove Phix sequences from libraries
- Bowtie2, nvBowtie and Nextgenemap: to align libraries against a reference genome

Packages and installation details can be found following the links:

| Program     | Source                                                |
| ----------- | ------------------------------------------------------|
| Bowtie2     | http://bowtie-bio.sourceforge.net/bowtie2/index.shtml |
| nvBowtie    | http://nvlabs.github.io/nvbio/                        |
| Nextgenemap | http://cibiv.github.io/NextGenMap/                    |

File containing genome sequence of the organism of interest should be provided to build the genome reference used for the alignment. Fasta files can be found in Illumina browser. Use the UCSC source (http://support.illumina.com/sequencing/sequencing_software/igenome.html). 


##INSTALLATION

To use HaSAPPy program download the package from GitHub repository. Programs can be started from the 'HaSAPPy/program' subfolder using command prompt


##GENERATE GENES REFERENCE LIBRARY FOR HaSAPPy SOFTWARE

After installation of HaSAPPY program, Genes Reference Library must be generated using GeneReference_built.py
The program requires two variables:

**-i** (INPUT) 	location of .txt file containing gene annotations according to UCSC browser. In the download folder, users can find the mm10_REFSEQgenes.txt file built for the mouse genome according to the assembly of Dec. 2011 (GCRm38/mm10). Alternatively, annotations can be obtained from UCSC browser (http://genome.ucsc.edu/cgi-bin/hgTables?command=start). Provide the following informations:	

| Task | Selection |
| --- | --- |
| *clade*	| Mammal |
| *genome* | Mouse or Human |
| *assembly* | (according to the last version) |
| *group*	| Genes and Gene Predictions |
| *track*	| RefSeq Genes |
| *table*	| refGene |
| *region* | genome |
| *output format*	| all fields from selected table |
| *output file* |	… |
| *file type returned* | gzip compressed |

**-o** (OUTPUT)	location where to store pandas library containing gene annotations necessary for HaSAPPy program. Informations are stored in a .pkl file. Add .pkl extension to the PATH provided

Ex.
```User$ python GeneReference_built.py -i User/HaSAPPy/Doc/mm10_REFSEQgenes.txt -o User/HaSAPPy/mm10_HaSAPPY_refernce.pkl```

Write output PATH in LoadModule.txt where requested (see next)


##RUN HaSAPPy SOFTWARE

Compile LoadModule.txt file according to instructions and save the modified file (Don’t overwrite it). Run the program using HaSAPPy_start.py file. The program requests as parameter the location of LoadModule file.

Ex.
```User$ python HaSAPPy_start.py User/HaSAPPy/Commands/LoadModule_170101.txt```




