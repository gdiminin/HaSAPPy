#Generate Gene Reference file for Human genome
In this tutorial, I will show you how to download gene annotations from the USCS browser and how to generate the GeneReference.pkl file used by HaSAPPy software.
##Preliminary
HaSAPPy software should be already installed on your computer and all requested Python modules already updated (see README.md file).

HaSAPPy software was installed in:

```
/Users/User/Analysis
```

The HaSAPPy directory should have this structure:

```
|- HaSAPPy
    |- LICENSE
    |- README.md
    |- program
        |- …
    |- docs
        |- test
            |- …
        |- LoadModule.txt
        |- mm10REFSEQgenes.txt
```

#Download the gene reference file form the UCSC browser
In a Web browser, go to the UCSC link: http://genome.ucsc.edu/cgi-bin/hgTables?command=start
Compile the page request using the following parameters:

| Task | Selection |
| --- | --- |
| *clade*	| Mammal |
| *genome* | Human |
| *assembly* | (according to the last version) |
| *group*	| Genes and Gene Predictions |
| *track*	| RefSeq Genes |
| *table*	| refGene |
| *region* | genome |
| *output format*	| all fields from selected table |
| *output file* | RefSeq_human_GRCh38-hg38 |
| *file type returned* | gzip compressed |

The Web browser window should look like that:

![alt-text](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/Figures/Generate_human_genome_reference_1.png)

Press ‘get output’

Move downloaded file in:

```
/Users/User/Analysis/HaSAPPy/docs
```

##Generate the GeneReference.pkl file using GeneReference_built.py module
Open the terminal and move to HaSAPPy program directory

```
cd /Users/User/Analysis/HaSAPPy/program
```

Start the GeneReference_built.py module using as input (-i) the downloaded file. Save the output (-o) in the same directory creating a file with the name GeneReference_Homo.pkl

```
python GeneReference_built.py -i  /Users/User/Analysis/HaSAPPy/docs/RefSeq_human_GRCh38-hg38.txt -o /Users/User/Analysis/HaSAPPy/docs/GeneReference_Homo.pkl
```

##Inspect the GeneReference.pkl file generatedG
The GeneReference_Homo.pkl should have been saved in the selected folder. To verify the integrity of the file, open the python interface

```
python
```

The Python console started

```
Python 2.7.12 |Anaconda 4.2.0 (x86_64)| (default, Jul  2 2016, 17:43:17) 
[GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
Anaconda is brought to you by Continuum Analytics.
Please check out: http://continuum.io/thanks and https://anaconda.org
>>> 
```

import the module that we are going to use

```python
import pickle
import pandas as pd
```

Load the data staored in the GeneReference_Homo.pkl file

```python
with open ('/Users/User/Analysis/HaSAPPy/docs/GeneReference_Homo.pkl', 'rb') as load:
 	gene_reference = pickle.load(load)
```

Inspect 
  *the file structure

```python
gene_reference.head(10)
```

![alt-text](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/Figures/Generate_human_genome_reference_2.png)

  *the columns name

```python
for column in gene_reference.columns:
	print column```

  *the number of genes

```python
print len(gene_reference.index)
```

The expected output should be

```
>>> for column in gene_reference.columns:
...     print column
... 
reference
genomic_interval
variants
exon_specific
introns_all
>>> 
>>> 
>>> print len(gene_reference.index)
32471
```

The GeneReference.pkl file was correctly generated and now it is ready to be used in HaSAPPy program
