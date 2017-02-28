# Time for a real case: Using HaSAPPy for analysing a screen for identifying candidates involved in X chromosome inactivation

In this tutorial, we are going to perform an anlysis with HaSAPPy software of a published datasets from a haploid embryonic stem cell screen to identify modulators of X chromosome inactivation. For the details of the screening refer to the paper ([Monfort et al., 2015](https://www.ncbi.nlm.nih.gov/pubmed/26190100)).

## Prerequisites
1. Installation of HaSAPPy package
2. Necessary Python packages correctly updated
3. Installation of Bowtie2
4. Download of reference genomes for the alignment programs
5. Maked executable the PreprocessReads
6. Built the gene annotation database using GeneReference_built.py
6. Performed the Test_HaSAPPy_installation tutorial

The working directory is
```
Users/User/HaSAPPy
```

Your working directory should look like that:
```
/Users/User/HaSAPPy/experiments/
│
├── *library*_yyyy-mm-dd
│   │
│   ├── graph
│   │   └── <graphic output of run>
│   └── raw
│       └── <intermediary files of data processing>
│
└── Analysis
└── HaSAPPy
    ├── reference
    │   ├── Phix
    │   │   └── ...
    │   └── Mus_musculus
    │       └── ...
    ├── experiments
    │   └── ...
    ├── program
    │   └── ...
    └── docs
        ├── LoadModule.txt
        ├── mm10REFSEQgenes.txt
        ├── GeneReference_Mouse-MM10.pkl
        ├── Tutorials
        │    └── ...
        └── test
             └── ...
```         

## Collecting FASTQ file from the SRA repository
We will use data from the SRA database associated with the haploid ES cell screen and store the read FASTQ files in a data folder in a selected/ and unselected/ subfolder.

```
cd ..
mkdir data
cd data
mkdir unselected
mkdir selected
```

The SRA database can be accessed at:
http://www.ncbi.nlm.nih.gov/sra

For the dataset of the unselected cell population enter **SRX1060416** into the search box. The selected dataset has the accession number **SRX1060407**. Download all read files in FASTQ format using the SRA toolkit. Install the SRA toolkit from: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
Choose Ubuntu Linux 64 bit architecture and spcify the ~/Downloads folder as destination. Open a new terminal window, enter the ~/Downloads folder and unzip the toolkit.

```
cd
cd Downloads
tar -xzvf sratoolkit.2.8.1-2-ubuntu64.tar.gz
exit
```

This closes this terminal window and returns you to the first terminal window. Now run the SRA toolkit **fastq-dump** command for downloading 7 unselected read files:

```
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064917
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064920
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064924
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064927
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064928
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064929
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064930
```

Then download the 7 selected read files to the data/selected/ folder:

```
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064886
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064896
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064898
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064899
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064902
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064905
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064908
```

## Prepare the LoadModule.txt file
This script will specify the reference and sequence files, as well as the processing steps and analysis that HaSAPPy will perform. 
User can utilize a pre-filled LoadMoule (find it in ```/Users/User/HaSAPPy/docs/test/LoadModule_XCI-tutorial.txt```) or try to use and fill up an empty template. The empty template and an example script are supplied with a source (see the tutorial [Controlling the Workflow of HaSAPPy with a command script](https://github.com/gdiminin/HaSAPPy/blob/master/docs/Tutorials/RunningHaSAPPyWorkflowsScript.md)).

### Create in experiments a folder where to store analysis data
```
cd /Users/User/HaSAPPy/experiments/
mkdir XCI_tutorial
```

### Use the pre-filled LoadModule
To use the pre-filled LoadModule open the LoadModule_XCI-tutorial.txt file in a text editor. As you can see most of the tasks have already been filled out. You need to provide the missing PATH information where files are stored in your file system (absolute PATH from root `/` is required). Complete the input fields of the form in the following sections:

* Section 1
* Task 1A and 1B

```
Operator Name: 
@1A) User
Storing location (provide a correct path):
@1B) /Users/User/HaSAPPy/experiments/XCI_tutorial
```
* Section 2
* 2D

```
Location of input file 1 (add additional lines if necessary):
@2D) /Users/User/HaSAPPy/data/unselected/SRR2064917.fastq.gz
@2D) ...
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

### Generate a new LoadModule for the analysis

Copy the **LoadModule.txt** command script template into a new command folder:
```
cd ..
mkdir command
cd command
cp ../docs/LoadModule.txt ./LoadModule_XCI-tutorial.txt
```

Use a text editor (eg. gedit or vi) to edit the LoadModule_XCI-tutorial.txt file in the ~/HaSAPPy/command/ folder. For read trimming before alignment follow the instructions for PreprocessReads, a useful adaptor sequence to trim from the 3'-end is the Ilumina P7 adaptor reverse complement (ATCTCGTATGCCGTCTTCTGCTT). 

## Start the analysis
Start the analysis with the following command:
> **NOTE:** Provide to HaSAPPY_start.py the correct LoadMoudle_XCI-tutorial.txt path. It is differnt if you are using the pre-filled (```/Users/User/HaSAPPy/docs/test/LoadModule_XCI-tutorial.txt```) or a new LoadModule file (```/Users/User/HaSAPPy/command/LoadModule_XCI-tutorial.txt```)

```
cd ..
cd program
python HaSAPPY_start.py /Users/User/HaSAPPy/../LoadModule_XCI-tutorial.txt
```

This will start the HaSAPPy program. The output will be stored under the experiments folder. You can compare it with the one present in HaSAPPy manual.


## Using nvBowtie with HaSAPPy

GPU accelerated alignment can speedup the alignment of the reads to the target genome. For using nvBowtie from NVIDIA specify nvBowtie in the command file as the aligner for the target reference genome. Both the nvBowtie and nvBWT executables will need to be installed on your system and in the path. You can generate the genome index for nvBowtie using nvBWT as follows:
```
cd
cd HaSAPPy/reference/Mus_musculus/UCSC/mm10/Sequence
mkdir nvBowtieIndex
nvBWT WholeGenomeFasta/genome.fa nvBowtieIndex/mm10
```

In the command script provide for the Mus musculus genome reference the following path: 
```
/home/anton/HaSAPPy/rereference/Mus_musculus/UCSC/mm10/Sequence/nvBowtieIndex/mm10
```

That is it! You can now run the analysis as before and might experience a reduction of the overall runtime of the analysis depending on the prowess of your GPU.

Note: If you have more than one GPU you can edit the Allign.py source file to specify multiple devices for the nvBowtie launch command specifying multiple --device options (--device 0 --device 1 --device 3).

## Using NextGenMap with HaSAPPy

NextGenMap is an aligner that supports the OpenCL API for using graphic processors and is designed to allow mapping of reads containing higher numbers of mismatches. The latter characteristic might make this aligner especially suitable for mapping insertion site libraries, when sequencing or PCR errors are suspected. Acceleration on AMD and NVIDIA GPUs as well as multi-GPU systems is supported. Use the --gpu 0,1,2 parameter for specifying the graphic processor device identifiers and -t 8 for specifying the number of threads (on hyperthreaded CPUs use double the number of CPU cores).

```
ngm --local --gpu 0,1 -t 8 -r <reference> -q <read file> -o <output file>
```

Specify ngm as the aligner in the comand script and /home/ … /HaSAPPy/reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa as the path for the target genome.



[**RETURN TO THE MAIN PAGE**](https://github.com/gdiminin/HaSAPPy/blob/master/README.md)
