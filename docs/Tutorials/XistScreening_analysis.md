We will use data from the SRA database associated with the haploid ES cell screen in Monfort et al. 2015 and store the read FASTQ files in a data folder/ in a selected/ and unselected/ subfolder.

cd ..
mkdir data
cd data
mkdir unselected
mkdir selected

The SRA database can be accessed at:

http://www.ncbi.nlm.nih.gov/sra

For the dataset of the unselected cell population enter SRX1060416 into the search box. The selected dataset has the accession number SRX1060407. Download all read files in FASTQ format using the SRA toolkit. Install the SRA toolkit from: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
Choose Ubuntu Linux 64 bit architecture and spcify the ~/Downloads folder as destination. Open a new terminal window, enter the ~/Downloads folder and unzip the toolkit.

cd
cd Downloads
tar -xzvf sratoolkit.2.8.1-2-ubuntu64.tar.gz
exit

This closes this terminal window and returns you to the first terminal window. Now run the SRA toolkit fastq-dump command for downloading 7 unselected read files:

~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064917
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064920
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064924
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064927
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064928
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064929
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O unselected/ --gzip SRR2064930

Then download the 7 selected read files to the data/selected/ folder:

~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064886
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064896
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064898
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064899
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064902
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064905
~/Downloads/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump -O selected/ --gzip SRR2064908

Next make a HaSAPPy command script for your project. This script will specify the reference and sequence files, as well as the processing steps and analysis that HaSAPPy will perform. An empty template and an example script are supplied with a source. Copy the template into a new command folder:

cd ..
mkdir command
cd command
cp ../docs/LoadModule.txt ./myProj.txt

Use a text editor (eg. gedit or vi) to edit the myProj.txt file in the ~/HaSAPPy/command/ folder. For read trimming before alignment follow the instructions for PreprocessReads, a useful adaptor sequence to trim from the 3'-end is the Ilumina P7 adaptor reverse complement (ATCTCGTATGCCGTCTTCTGCTT). Start the analysis with the following command:

cd ..
mkdir analysis
cd program
python HaSAPPY_start.py ../command/myProj.txt

This will start the HaSAPPy program the output will be stored under the analysis/ folder.

Update the HaSAPPy version from github in the HaSAPPy folder using:

git pull origin


Using nvBowtie with HaSAPPy

GPU accelerated alignment can speedup the alignment of the reads to the target genome. For using nvBowtie from NVIDIA specify nvBowtie in the command file as the aligner for the target reference genome. Both the nvBowtie and nvBWT executables will need to be installed on your system and in the path. You can generate the genome index for nvBowtie using nvBWT as follows:

cd
cd HaSAPPy/reference/Mus_musculus/UCSC/mm10/Sequence
mkdir nvBowtieIndex
nvBWT WholeGenomeFasta/genome.fa nvBowtieIndex/mm10

In the command script provide for the Mus musculus genome reference the following path: /home/anton/HaSAPPy/rereference/Mus_musculus/UCSC/mm10/Sequence/nvBowtieIndex/mm10

That is it! You can now run the analysis as before and might experience a reduction of the overall runtime of the analysis depending on the prowess of your GPU.

Note: If you have more than one GPU you can edit the Allign.py source file to specify multiple devices for the nvBowtie launch command specifying multiple --device options (--device 0 --device 1 --device 3).

Using NextGenMap with HaSAPPy

NextGenMap is an aligner that supports the OpenCL API for using graphic processors and is designed to allow mapping of reads containing higher numbers of mismatches. The latter characteristic might make this aligner especially suitable for mapping insertion site libraries, when sequencing or PCR errors are suspected. Acceleration on AMD and NVIDIA GPUs as well as multi-GPU systems is supported. Use the --gpu 0,1,2 parameter for specifying the graphic processor device identifiers and -t 8 for specifying the number of threads (on hyperthreaded CPUs use double the number of CPU cores).

ngm --local --gpu 0,1 -t 8 -r <reference> -q <read file> -o <output file>

Specify ngm as the aligner in the comand script and /home/ … /HaSAPPy/reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa as the path for the target genome.
