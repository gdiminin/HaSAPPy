# Generate a Gene Reference file for the Human genome

This tutorial provides a step-by-step guide how to download the genome annotation using the UCSC genome browser and how to generate the **GeneReference.pkl** file used by HaSAPPy.

The HaSAPPy software should be already installed on your computer with all dependencies (see the [README.md](https://github.com/gdiminin/HaSAPPy/blob/master/README.md) file for how to install) and tested ([Performing a test of your HaSAPPy installation](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/Test_HaSAPPy_installation.md)).

This tutorial assumes that HaSAPPy was installed in:

```
/Users/User
```

The HaSAPPy directory should have this structure:

```
/Users/User/
│
└── HaSAPPy
├── LICENSE
├── README.md
├── program
│   └── ...
└── docs
├── test
│   └── ...
├── LoadModule.txt
└── mm10REFSEQgenes.txt
```

## Download the gene reference file form the UCSC browser

In a Web browser, go to the UCSC link: http://genome.ucsc.edu/cgi-bin/hgTables?command=start
The web page will show a number of input fields requesting information on the genome assembly, type and extend of annotation to download. Enter the following parameters:

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

Move the downloaded annotation file into the HaSAPPy/docs folder:

```
/Users/User/HaSAPPy/docs
```

## Generate the GeneReference.pkl file using GeneReference_built.py
Open the terminal and move to the HaSAPPy program folder

```
cd /Users/User/HaSAPPy/program
```

Start the **GeneReference_built.py** module using as input (option `-i`) the downloaded file. Save the output (option `-o`) in the same folder creating a file with the name GeneReference_Homo.pkl

```
python GeneReference_built.py -i  /Users/User/HaSAPPy/docs/RefSeq_human_GRCh38-hg38.txt -o /Users/User/HaSAPPy/docs/GeneReference_Homo.pkl
```

## Inspect the GeneReference.pkl file generated
The **GeneReference_Homo.pkl** should have been saved in the selected folder. To verify the integrity of the file, open the python console:

```
python
```

The Python console started

```
Python 2.7.12 |Anaconda 4.2.0 (x86_64)| (default, Jul  2 2016, 17:43:17) 
[GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
Anaconda is brought to you by Continuum Analytics.
Please check out: http://continuum.io/thanks and https://anaconda.org
>>> 
```

Import the **pickle** and **pandas** modules necessary to use for working with the gene annotation reference file:

```python
import pickle
import pandas as pd
```

Load the data stored in the **GeneReference_Homo.pkl** file with the following command:

```python
with open ('/Users/User/HaSAPPy/docs/GeneReference_Homo.pkl', 'rb') as load:
gene_reference = pickle.load(load)
```

Inspect 
* the file structure

```python
gene_reference.head(10)
```

![alt-text](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/Figures/Generate_human_genome_reference_2.png)

* the columns name

```python
for column in gene_reference.columns:
print column
```

* the number of genes

```python
print len(gene_reference.index)
```

The expected output should be

```python
>>> for column in gene_reference.columns:
...     print column
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

The GeneReference.pkl file was correctly generated and can now be used for working with the human genome in HaSAPPy. To specify the human genome annotation for a workflow enter the file PATH into a HaSAPPy command script. **NOTE:** read mapping must be performed with to the same assembly of the human genome as was used to download the annotation!


[**RETURN TO THE MAIN PAGE**](https://github.com/gdiminin/HaSAPPy/blob/master/README.md)
