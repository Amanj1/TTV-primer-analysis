# TTV-primer-analysis
A  Nextflow DSL 1 pipeline for analysis of primer sequences in Torque teno virus

This pipeline focuses on finding three primers in the untranslated regions (UTR) that are commonly used for PCR analysis of Torque teno virus (TTV). 
It finds primers AMTS UTR, AMTAS UTR and AMTPTU UTR. Primers found in sequences are based on exact matching, one mismatch and two mismatches. 
Outputs are calculations of sequences and tables with summaries of sequences that contain at least one primer. 
Different bar charts are also created based on the results from the tables.

There are four tables created. One table only for extact primer results, one table for both exact match and 0-1 mismatch, 
one table that includes 0-2 mismatches and a table with the combinations of all three tables displaying which sequnces has a exact match with the label "M0", 
"M1" for sequences with one mismatch, "M2" with two mismatches and "-" for no primer match. 

![alt text](primar_analysis_ttvUTR.png)
