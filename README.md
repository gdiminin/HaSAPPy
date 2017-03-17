# HaSAPPy (Haploid Screening Analysis Package in Python)



## INTRODUCTION

Forward genetic screens are becoming increasingly useful as approaches to unveil the complexity of cellular mechanisms. Recent introduction of mammalian haploid cells opened the possibility to perform forward genetic screens using insertional mutagenesis in mammals. Haploid cells have originally been obtained from human tumors and subsequently have been established from haploid embryos of several mammals.

Mutagenesis of haploid cells leads to hemizygous mutantions which ensures phenotypic exposure. This is in contrast to diploid cells, where in general two copies of each chromosome are present, and effects of mutations in specific genes can be masked by the functional copy on the second chromosome.

For screening viral gene trap vectors are commonly used as mutagens for their large mutational effect and utility for providing a sequence tag for the genomic insertion site. Large haploid mutant pools are established by viral transduction following by selection of a phenotype of interest. Subsequent identification of viral insertions that become enriched after selection can define gene mutations and pathways that contribute to the phenotype under investigation. Readout of many thousands of insertion sites has become possible by deep sequencing of DNA extracted from large cell pools using specialized Next Generation Sequencing (NGS) protocols.

HaSAPPy performs analysis of NGS datasets form pooled haploid mammalian cell screens and is used to identify insertion locations in the whole genome, map them at the level of genes, and classify insertions according to their effects on gene function. Customizable output of all calculated parameters and ranking of candidate genes is performed. The implementation conforms to current Python programming guidelines and is open to be adapted and extended according to experimental setup. The data processing and analysis is performed by the following modules:

**Modules:**
 - Read trimming of adaptor and low quality sequences in **Trim.py** (requires PreprocessReads)
 - Read alignment to PhiX control and reference genomes in **Align.py** (bowtie2, NextGenMap, or nvBowtie are required)
 - Identification of Independent Insertions (I.I.) in **IIDefinition.py**
 - Counting insertions at the genes level in **GenesDefinition.py**
 - Enrichment analysis in **GroupAnalysis.py**
 - Data presentation **Tables.py** and **DesignGeneInsertions.py**

**Output**
Output of the program is written in standard .xlsx and .svg file formats and can be customized in a command script.


## PREREQUISITES

To use HaSAPPy the Python 2.7 version has to be installed in the local machine (Python 3 is not compatible with some used pakages).

HaSAPPy program requests the following Python Packages that can be downloaded/updated using the `pip -install` or `pip install --upgrade` command:
 - numpy
 - HTSeq
 - matplotlib
 - pandas
 - scipy
 - sklearn
 - xlsxwriter
 
For using the alignment module NGS readmappers are required and should be in the PATH:
- Bowtie2: to remove Phix sequences from libraries
- Bowtie2, nvBowtie and Nextgenemap: to align libraries against a reference genome

Instructions how to install these read mappers can be found at:

| Program     | Source                                                |
| ----------- | ------------------------------------------------------|
| Bowtie2     | http://bowtie-bio.sourceforge.net/bowtie2/index.shtml |
| nvBowtie    | http://nvlabs.github.io/nvbio/                        |
| Nextgenemap | http://cibiv.github.io/NextGenMap/                    |

Fasta files for the PhiX and reference genomes can be obtained from the Illumina browser, for use with HaSAPPy select UCSC as the source:

http://support.illumina.com/sequencing/sequencing_software/igenome.html 


## INSTALLATION OF HaSAPPy

To use HaSAPPy program download the package from GitHub repository. If you have git installed you can use the following command:
```
cd
git clone https://github.com/gdiminin/HaSAPPy.git
```
This will install all Python files in the 'HaSAPPy' folder under your home folder. For read preprocessing a precompiled executable PreprocessReads is supplied. This is an optional step in data processing and if desired requires to set the execute property on this file:
```
cd HaSAPPy/program
chmod +x PreprocessReads
```

This step is not required if adaptor and quality trimming of the reads is not desired.


## GENERATE THE GENE ANNOTATION REFERENCE FOR HaSAPPy

After installation of HaSAPPY program, Genes Annotation Reference must be generated using GeneReference_built.py
The program requires two variables:

`python GeneReference_built.py -i <path-to-input-file.txt> -o <path-to-output-anotation.pkl>`

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

An example command line for generating an annotation for the mouse genome:
```
python GeneReference_built.py -i User/HaSAPPy/docs/mm10_REFSEQgenes.txt -o User/HaSAPPy/docs/GeneReference_Mouse-MM10.pkl
```

The [Generate Gene Reference file for Human genome](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/CreateHumanGeneAnnotationReference.md) tutorial explains how to download annotations for the human genome and generate Genes Annotation Reference.

To use a Gene Annotation Reference .pkl file with HaSAPPy its file path is entered in the command script (see the [Performing a test run of your HaSAPPy installation](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/TestRunHaSAPPY.md) tutorial).


## Running HaSAPPy

Running a data analysis requires a command script that is pecified on the command line as input for HaSAPPy. An empty script is provided in docs/LoadModule.txt. Copy this file and enter the experimental parameters in the form (see more detail on how to generate a command script in the [Controlling the HaSAPPy workflow throught the command script](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/RunningHaSAPPyWorkflowsScript.md) tutorial).

The analysis workflow can then be strated using HaSAPPY_start.py with the appropriate command script, see the following example:

```
python HaSAPPY_start.py User/HaSAPPy/Commands/LoadModule.txt
```

After installation it is adviced to run HaSAPPy using the LoadModule_Test.txt file to confirm program and dependencies installation success. For instructions on how to perform the test run refer to the tutorial [Performing a test run of your HaSAPPy installation](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/TestRunHaSAPPY.md).

To explore HaSAPPy potentialities and to get experience in pipline managment the tutorial [Time for a real case: Using HaSAPPy for analysing a screen for identifying candidates involved in X chromosome inactivation](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/XCI_ScreenAnalysis.md) guides the user on how performing the analysis on true datasets.


## REPORTING OF ERRORS, FEATURE REQUEST, AND REQUESTS FOR HELP

Any issues with using HaSAPPy should be raised at https://github.com/gdiminin/HaSAPPy/issues

Feature requests and pull request are appreciated https://github.com/gdiminin/HaSAPPy/pulls
and be attended to as quickly as resources allow.

Send additional enquiries to `giulio.diminin@biol.ethz.ch`



