# TTV-primer-analysis
A  Nextflow DSL 1 pipeline for analysis of primer sequences in Torque teno virus.

This pipeline focuses on finding three primers in the untranslated regions (UTR) that are commonly used for PCR analysis of Torque teno virus (TTV). 
It finds primers AMTS UTR, AMTAS UTR and AMTPTU UTR. Primers found in sequences are based on exact matching, one mismatch and two mismatches. 
Outputs are calculations of sequences and tables with summaries of sequences that contain at least one primer. 
Different bar charts are also created based on the results from the tables.

There are four tables created. One table only for extact primer results, one table for both exact match and 0-1 mismatch, 
one table that includes 0-2 mismatches and a table with the combinations of all three tables displaying which sequnces has a exact match with the label "M0", 
"M1" for sequences with one mismatch, "M2" with two mismatches and "-" for no primer match. 

![alt text](/img/ttv_primer_comb_table.png)

This image displayes the table with the combination of all tables

In img folder you can find different bar chart plots that are produced by the pipeline. 
There are two versions of each bar chart. First bar chart is UTR primers in relation to each other and the second is the relation to all primer and total sequences used in the pipeline. 

## Please note
The pipeline does not reads sequences from both the positive and negative strand. It only reads sequences based on positive. If you want to read the reverse complement of your sequences I would advise that you create reverse sequences and keep the same fasta header. All calculations are based on the header and if you get a primer match in any of the sequences with the same header it counts as one match. 

 ## Software requirements 
 - Nextflow
 - python3
 - numpy (pyhton 3 library)
 - matplotlib (pyhton 3 library)
 - seqtk
 - seqkit
