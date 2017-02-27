# Performing a test run of your HaSAPPy installation

This tutorial shows how to quickly verify that HaSAPPy is correctly installed on your system and all dependencies are met. For testing the core functionality of HaSAPPy a pre-filled **LoadModule.txt** command script and a **Sequence_test** file are included in the docs/test folder in the installation.
It is recommended to run this text for a new installation before starting to use HaSAPPy for data analysis workflows.

The working directory for this tutorial is:

```
/Users/User
```

## Prerequisites

1. Installation of HaSAPPy package
2. Necessary Python packages correctly updated
3. Installation of bowtie2 

### Create a folder where the experiments and run information will be stored
Enter into the HaSAPPy folder and create a new **experiments/** folder:
```
cd HaSAPPy
mkdir experiments
cd experiments
mkdir test
```
The folder structure should look like this:
```
/Users/User/
|
└── HaSAPPy
    ├── ...
    └── experiments
        └── test
            └── ...
```

### Download and unpack sequence archives for PhiX and Mus Musculus genomes
Add a **reference/** folder where the genome sequences and indices will be stored:
```
mkdir reference
```
The folder structure should look like this:
```
/Users/User/
|
└── HaSAPPy
    ├── ...
    └── experiments
    |   └── test
    |       └── ...
    └── reference
	    └── ...
```

Genome sequences and indices for read mapping for common species are available at:

http://support.illumina.com/sequencing/sequencing_software/igenome.html

Using a web browser download the PhiX genome (NCBI 1993-04-28) and the mouse genome (GCRm38/mm10) to your ~/HaSAPPy/reference folder. Then unpack the archives in the terminal:

```
cd reference
tar -xzvf PhiX_NCBI_1993-04-28.tar.gz
tar -xzvf Mus_musculus_UCSC_mm10.tar.gz
```


### Set execute permission on PreprocessReads

PreprocessReads is a precompiled executable file for read trimming (the source code is available at [github](https://github.com/zanton123/PreprocessReads). You need to set the execute permission on Linux based file systems to be able to run the program. Enter the /HaSAPPy/program/** folder and set the execute property with the following commands:

```
cd ..
cd program
chmod +x PreprocessReads 
```

### Build the gene annotation database using GeneReference_built.py

A gene annotation file for the mouse genome (GCRm38/mm10) is supplied with the source.

>**NOTE:** If you need other genomes suitable annotation files can be obtained from the UCSC genome browser (http://genome.ucsc.edu/cgi-bin/hgTables?command=start). 
Refer to the [Generate a human gene reference](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/CreateHumanGeneAnnotationReference.md) tutorial for a detailed example:

```
python GeneReference_built.py -i /Users/User/HaSAPPy/docs/mm10_REFSEQgenes.txt -o /Users/User/HaSAPPy/docs/GeneReference_mouse_mm10.pkl
```

A new file GeneReference_mouse_mm10.pkl should appear in the current folder.

Now we are ready to start. Your working directory (**/Users/User/HaSAPPy**) should have this structure:

* the HaSAPPy directory

```
/Users/User/
|
└── HaSAPPy
    ├── ...
    └── experiments
    |   └── test
    |       └── ...
    └── reference
    |   └── ...
    ├── ...
    ├── program
    |   └── ...
    └── docs
        ├── LoadModule.txt
        ├── mm10REFSEQgenes.txt
        ├── GeneReference_mouse_mm10.pkl
        ├── Tutorials
        |   └── ...
        └── test
             ├── Aligned.sam
             ├── Sequence.fastq
             ├── LoadModule_test_from_Trim.txt
             └── LoadModule_test_from_IIdefinition.txt
```

* the reference directory

```
└── reference
    ├── Phix
    |   └── NCBI
    |       └── 1993-04-28
    |           ├── Annotation
    |           |   └── ...
    |           ├── Sequence
    |           ├── ...
    |	        └── Bowtie2Index
    |	            ├── genome.1.ebwt
    |	            └── ...
    └── Mus_musculus
        └── UCSC
            └── mm10
                 ├── Annotation
                 |   └── ...
                 ├── Sequence
                 ├── ...
     	         └── Bowtie2Index
     	             ├── genome.1.ebwt
     	             └── ...      
```

* the experiments directory

```
└── experiments
    └── test
         └── ...
```

## Edit the HaSAPPy command script - the LoadModule file

In this tutorial, we are going to test the functionality of the HaSAPPy packages (starting from the **Trim.py** module) using the **LoadModule_test_from_Trim.txt** and the **Sequence.fastq** files in the `docs/test` folder.

```
/Users/User/HaSAPPy/docs/test
```

In addition, you can test HaSAPPy function from the **IIdefinition.py** module using the **LoadModule_test_from_IIdefinition.txt** and the **Aligned.sam** files.

Open the LoadModule_test_from_Trim.txt file in a text editor. As you can see most of the tasks have already been filled out. You need to provide the missing PATH information where files are stored in your file system (absolute PATH from root `/` is required). Complete the input fields of the form in the following way:

* Section 1
* Task 1A and 1B

```
Operator Name: 
@1A) User
Storing location (provide a correct path):
@1B) /Users/User/HaSAPPy/experiments/test
```
* Section 2
* 2D

```
Location of input file 1 (add additional lines if necessary):
@2D) /Users/User/HaSAPPy/docs/test/Sequence.fastq
@2D) /Users/User/HaSAPPy/docs/test/Sequence.fastq
```

* Section 3
*  3A

```
Location of Phix reference genome:
@3A)/Users/User/HaSAPPy/reference/PhiX/NCBI/1993-04-28/Sequence/Bowtie2Index/genome
```

> **NOTE:** For Bowtie2 the path and filename for the genome index should be provided without the file extension!

* Section 4
* 4B

```
Location of reference genome:
@4B) /Users/User/HaSAPPy/reference/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome
```

* Section 6
* 6A

```
Location of gene reference:
@6A)/Users/User/HaSAPPy/docs/GeneReference_Mouse-MM10.pkl
```

* Section 9
* 9C

```
FOR Plot I.I. in gene models:
Location of gene reference:
@9C)/Users/User/HaSAPPy/docs/GeneReference_Mouse-MM10.pkl
```

Save the file and copy the file name to the clipboard.

## Run the test

Move into the HaSAPPy/program folder:

```
cd /Users/User/HaSAPPy/program
```

Start the analysis using **HaSAPPY_start.py**  with the **LoadModule_test_from_Trim.txt** file path as command line argument:

```
python HaSAPPY_start.py /Users/User/HaSAPPy/docs/test/LoadModule_test_from_Trim.txt
```

If everything is correctly installed the analysis should finish without any errors.

## Inspect the output

In the `/Users/User/HaSAPPy/experiments/test` folder the following files and folders should have been created

```
/Users/User/HaSAPPy/experiments/test
    ├── test_1_yyyy-mm-dd
    |   ├── test_1_info.txt
    |   ├── graph
    |   └── raw
    |       ├── test_1_Aligned.sam
    |       ├── ...
    |       └── test_1_IIRawdata.pkl
    ├── test_2_yyyy-mm-dd
    |   ├── test_2_info.txt
    |   ├── graph
    |   └── raw
    |       ├── test_2_Aligned.sam
    |       ├── ...
    |       └── test_2_IIRawdata.pkl
    └── Analysis
        └── yyyy-mm-dd
            ├── analysis_info.txt
            ├── graph
            |   ├── Xist_SelectedvsControl.svg
            |   └── ...
            ├── raw
            |   ├── GroupAnalysis.pkl
            |   └── RawData.pkl
            └── Table_yyyy-mm-dd.xlsx
```

Open the **Table_yyyy-mm-dd.xlsx** file in a spreadsheet program (eg. LibreOffice Calc, WPS Office Spreadsheet, or Microsoft Excel) and the **Xist_SelectedvsControl.svg** file in a vector graphics program (eg. Adobe Illustratore, or Inkscape). You should see the output in table format and a graphic view of the insertions in the _Xist_ gene similar as following output:

![alt text] (https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/Figures/Test_HaSAPPy_installation.png)

If this is what you see, you have now a working HaSAPPy installation. Congratulation!!!!


[**RETURN TO THE MAIN PAGE**](https://github.com/gdiminin/HaSAPPy/blob/master/README.md)
