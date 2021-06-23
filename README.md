# TTV-primer-analysis
A  Nextflow DSL 1 pipeline for analysis of primer sequences in Torque teno virus.

This pipeline focuses on finding three primers in the untranslated regions (UTR) that are commonly used for PCR analysis of Torque teno virus (TTV). 
It finds primers AMTS UTR, AMTAS UTR and AMTPTU UTR. Primers found in sequences are based on exact matching, one mismatch and two mismatches. 
Outputs are calculations of sequences and tables with summaries of sequences that contain at least one primer. 
Different bar charts are also created based on the results from the tables.

There are four tables created. One table only for extact primer results, one table for 0-1 mismatch, 
one table that includes 0-2 mismatches and a table with the combinations of all three tables displaying which sequnces has a exact match with the label "M0", 
"M1" for sequences with one mismatch, "M2" with two mismatches and "-" for no primer match. 

![alt text](/img/ttv_primer_comb_table.png)

This image displayes the table with the combination of all tables

In img folder you can find different bar chart plots that are produced by the pipeline. 
There are two versions of each bar chart. First bar chart is UTR primers in relation to each other and the second is the relation to all primer and total sequences used in the pipeline. 

#### Universal TM-PCR - Primer sequences used in pipeline

| Oligonucleotide | Sequence (5'-3')   | Genome localization | Position (nt) |
| --------------  | ------------------ | ------------------  |  -----------  |
| AMTS            | GTGCCGIAGGTGAGTTTA | UTR                 | 177–194       |
| AMTAS           | AGCCCGGCCAGTCC     | UTR                 | 226–239       |
| AMTPTU          | TCAAGGGGCAATTCGGGCT| UTR                 | 205–223       |

## Please note
The pipeline does not reads sequences from both the positive and negative strand. It only reads sequences based on positive. If you want to read the reverse complement of your sequences I would advise that you create reverse sequences and keep the same fasta header. All calculations are based on the header and if you get a primer match in any of the sequences with the same header it counts as one match. 

 ## Software requirements 
 All versions of the softwares should be compatible with the pipeline. Currently I don't have any conflicts between softwares. 
 - [Nextflow DSL1](https://www.nextflow.io/)
 - [Python3](https://www.python.org/downloads/)
    - numpy (pyhton 3 library)
    - matplotlib (pyhton 3 library)
 - [seqtk](https://github.com/lh3/seqtk)
 - [seqkit](https://bioinf.shenwei.me/seqkit/)

#### Software versions used when developing pipeline
Softwares
| Software | Version |
| -------- | ------- |
| Nextflow | 19.07.0 |
| Python   | 3.7.8   |
| Seqtk    | 1.3-r106|
| Seqkit   | 0.13.2  |

Python3 modules
| Package   | Version |
| --------- | ------- |
| numpy     | 1.17.2  |
| matplotlib| 3.3.0   |

## Prerequisites
To ensure that resources are available or to increase the number of threads, you can change the configuration for profile 'amanj' contained in the conf folder. It is good to remember that all three processes contained in the configuration file are the only ones capable of multithreading and runs simultaneously. Remaining processes are not listed in the configuration file because they run in one thread only.

In the following code, the number can be changed after 'cpus' (threads):
```
process {
    withName: grep_misMatch0{
        cpus = 13
    }
    withName: grep_misMatch1{
        cpus = 13
    }
    withName: grep_misMatch2{
        cpus = 13
    }
}
```

## Installing TTV-primer-analysis using Conda environment
You could also use the "ttv-primer-env.yml" file to install all software requirements in the following terminal commands below.
```
conda env create -f ttv-primer-env.yml
conda activate ttv-primer
```
Once you have activated the Conda environment "ttv-primer", you can type the following command in the terminal.
```
which python3

example output would be: /home/user/anaconda3/envs/ttv-primer/bin/python3
```
Now you should replace the shebang line ``` "#!/bin/python3" ``` to ```"#!/home/user/anaconda3/envs/ttv-primer/bin/python3"``` in this file "ttv_utr_primer_analysis.nf". There should be 8 python3 shebang lines in the file. You can use any text editor to find and replace all lines. 

## Running TTV-primer-analysis
The user should create a directory called 'input_data' and store all samples must be in fasta format with the file extension ".fasta". 
It is also possible to change the following code line to the preferred filename and file extension, but the content of the files should be in fasta format:
```
/* input files */
//contig sequences
contig_files = Channel.fromFilePairs("*.fasta",size:1) 
```
It is preferred to use longer sequences as input data but works fine with short reads.
You should take into account that the execution time is linear.

To run the pipeline in command line:
```
nextflow -C ttv_utr_primer_analysis.config run ttv_utr_primer_analysis.nf -profile amanj
```
To run the pipeline in command line and resume from cache memory:
```
nextflow -C ttv_utr_primer_analysis.config run ttv_utr_primer_analysis.nf -profile amanj -resume
```
