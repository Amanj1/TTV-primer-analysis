#!/usr/bin/env nextflow

/*
How to run:
nextflow -C ttv_utr_primer_analysis.config run ttv_utr_primer_analysis.nf
*/

params.input='input_data'

/* input files */
//contig sequences
contig_files = Channel.fromFilePairs("*.fasta",size:1) 

process filter_contig_size{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}", mode:'link'

  input:
  set sample_id, seq from contig_files
  
  output:
  set sample_id, "${sample_id}_contigs_filt.fasta"  into filter_contigs_out
  
  script:
""" 
seqtk seq -L 2 ${seq[0]} > ${sample_id}_contigs_filt.fasta
"""
}

filter_contigs_out.into{
misMatch0_in;
misMatch1_in;
misMatch2_in;
filt_contig_size_count_in;
}

process number_of_sequences{
 tag {"${sample_id}"}
 
 publishDir "${params.publish_base_dir}/${sample_id}", mode:'link'
  input:

  set sample_id, seq from filt_contig_size_count_in

  output:
  set sample_id, "${sample_id}_nrOfSeq.txt"  into nr_of_seq_out

  script:
"""
cat ${seq} | grep ">" | awk '!seen[\$0]++' | wc -l > "${sample_id}_nrOfSeq.txt"

"""
}

process collect_nr_of_seq{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All", mode:'link'

input:
file nrOfSeq from nr_of_seq_out.map{it[1]}.collect()

output:

file "all_samples_collected_nrOfSeq.txt" into collected_nrOfSeq_out

script:
"""
IFS=' ' read -r -a arr <<< \$(echo ${nrOfSeq})
total=0
for i in "\${arr[@]}"
do
	total=\$(expr \$total + \$(cat \$i | bc))
done

echo \$total > "all_samples_collected_nrOfSeq.txt"

"""
}

collected_nrOfSeq_out.into{
total_nrOfSeq_M0_in;
total_nrOfSeq_M1_in;
total_nrOfSeq_M2_in;
toatl_nrOfSeq_all_in
}

process grep_misMatch0{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/misMatch0", mode:'link'

  input:
  set sample_id, seq from misMatch0_in 
  
  output:
  set sample_id, "${sample_id}_list_all_headers_M0.txt", "${sample_id}_AMTS_UTR_M0.fasta", "${sample_id}_AMTAS_UTR_M0.fasta", "${sample_id}_AMTPTU_UTR_M0.fasta" into misMatch0_out

  script:
""" 
#!/bin/bash
#Universal TM-PCR
#AMTS
seqkit grep --by-seq --max-mismatch 0 --threads 13 --pattern "GTGCCGAAGGTGAGTTTA" ${seq} > AMTS_UTR.fasta
sleep 2
seqkit grep --by-seq --max-mismatch 0 --threads 13 --pattern "GTGCCGTAGGTGAGTTTA" ${seq} >> AMTS_UTR.fasta
sleep 2
seqkit grep --by-seq --max-mismatch 0 --threads 13 --pattern "GTGCCGGAGGTGAGTTTA" ${seq} >> AMTS_UTR.fasta
sleep 2
seqkit grep --by-seq --max-mismatch 0 --threads 13 --pattern "GTGCCGCAGGTGAGTTTA" ${seq} >> AMTS_UTR.fasta
sleep 2
cat AMTS_UTR.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > tmp_list.txt
sleep 2
seqtk subseq ${seq} tmp_list.txt > ${sample_id}_AMTS_UTR_M0.fasta
sleep 2
rm tmp_list.txt AMTS_UTR.fasta
#AMTAS
seqkit grep --by-seq --max-mismatch 0 --threads 13 --pattern "AGCCCGGCCAGTCC" ${seq} > AMTAS_UTR.fasta
sleep 2
cat AMTAS_UTR.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > tmp_list.txt
sleep 2
seqtk subseq ${seq} tmp_list.txt > ${sample_id}_AMTAS_UTR_M0.fasta
sleep 2
rm tmp_list.txt AMTAS_UTR.fasta
#AMTPTU
seqkit grep --by-seq --max-mismatch 0 --threads 13 --pattern "TCAAGGGGCAATTCGGGCT" ${seq} > AMTPTU_UTR.fasta
sleep 2
cat AMTPTU_UTR.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > tmp_list.txt
seqtk subseq ${seq} tmp_list.txt > ${sample_id}_AMTPTU_UTR_M0.fasta
sleep 2
rm tmp_list.txt AMTPTU_UTR.fasta

cat *M0.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > ${sample_id}_list_all_headers_M0.txt
"""
}

process grep_misMatch1{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/misMatch1", mode:'link'

  input:
  set sample_id, seq from misMatch1_in 
  
  output:
  set sample_id, "${sample_id}_list_all_headers_M1.txt", "${sample_id}_AMTS_UTR_M1.fasta", "${sample_id}_AMTAS_UTR_M1.fasta", "${sample_id}_AMTPTU_UTR_M1.fasta" into misMatch1_out

  script:
""" 
#!/bin/bash
#Universal TM-PCR
#AMTS
seqkit grep --by-seq --max-mismatch 1 --threads 13 --pattern "GTGCCGAAGGTGAGTTTA" ${seq} > AMTS_UTR.fasta
sleep 2
seqkit grep --by-seq --max-mismatch 1 --threads 13 --pattern "GTGCCGTAGGTGAGTTTA" ${seq} >> AMTS_UTR.fasta
sleep 2
seqkit grep --by-seq --max-mismatch 1 --threads 13 --pattern "GTGCCGGAGGTGAGTTTA" ${seq} >> AMTS_UTR.fasta
sleep 2
seqkit grep --by-seq --max-mismatch 1 --threads 13 --pattern "GTGCCGCAGGTGAGTTTA" ${seq} >> AMTS_UTR.fasta
sleep 2
cat AMTS_UTR.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > tmp_list.txt
seqtk subseq ${seq} tmp_list.txt > ${sample_id}_AMTS_UTR_M1.fasta
sleep 2
rm tmp_list.txt AMTS_UTR.fasta
#AMTAS
seqkit grep --by-seq --max-mismatch 1 --threads 13 --pattern "AGCCCGGCCAGTCC" ${seq} > AMTAS_UTR.fasta
sleep 2
cat AMTAS_UTR.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > tmp_list.txt
seqtk subseq ${seq} tmp_list.txt > ${sample_id}_AMTAS_UTR_M1.fasta
sleep 2
rm tmp_list.txt AMTAS_UTR.fasta
#AMTPTU
seqkit grep --by-seq --max-mismatch 1 --threads 13 --pattern "TCAAGGGGCAATTCGGGCT" ${seq} > AMTPTU_UTR.fasta
sleep 2
cat AMTPTU_UTR.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > tmp_list.txt
seqtk subseq ${seq} tmp_list.txt > ${sample_id}_AMTPTU_UTR_M1.fasta
sleep 2
rm tmp_list.txt AMTPTU_UTR.fasta
cat *M1.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > ${sample_id}_list_all_headers_M1.txt
"""
}

process grep_misMatch2{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/misMatch2", mode:'link'

  input:
  set sample_id, seq from misMatch2_in 
  
  output:
  set sample_id, "${sample_id}_list_all_headers_M2.txt", "${sample_id}_AMTS_UTR_M2.fasta", "${sample_id}_AMTAS_UTR_M2.fasta", "${sample_id}_AMTPTU_UTR_M2.fasta" into misMatch2_out

  script:
""" 
#!/bin/bash
#Universal TM-PCR
#AMTS
seqkit grep --by-seq --max-mismatch 2 --threads 13 --pattern "GTGCCGAAGGTGAGTTTA" ${seq} > AMTS_UTR.fasta
sleep 2
seqkit grep --by-seq --max-mismatch 2 --threads 13 --pattern "GTGCCGTAGGTGAGTTTA" ${seq} >> AMTS_UTR.fasta
sleep 2
seqkit grep --by-seq --max-mismatch 2 --threads 13 --pattern "GTGCCGGAGGTGAGTTTA" ${seq} >> AMTS_UTR.fasta
sleep 2
seqkit grep --by-seq --max-mismatch 2 --threads 13 --pattern "GTGCCGCAGGTGAGTTTA" ${seq} >> AMTS_UTR.fasta
sleep 2
cat AMTS_UTR.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > tmp_list.txt
seqtk subseq ${seq} tmp_list.txt > ${sample_id}_AMTS_UTR_M2.fasta
sleep 2
rm tmp_list.txt AMTS_UTR.fasta
#AMTAS
seqkit grep --by-seq --max-mismatch 2 --threads 13 --pattern "AGCCCGGCCAGTCC" ${seq} > AMTAS_UTR.fasta
sleep 2
cat AMTAS_UTR.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > tmp_list.txt
seqtk subseq ${seq} tmp_list.txt > ${sample_id}_AMTAS_UTR_M2.fasta
sleep 2
rm tmp_list.txt AMTAS_UTR.fasta
#AMTPTU
seqkit grep --by-seq --max-mismatch 2 --threads 13 --pattern "TCAAGGGGCAATTCGGGCT" ${seq} > AMTPTU_UTR.fasta
sleep 2
cat AMTPTU_UTR.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > tmp_list.txt
seqtk subseq ${seq} tmp_list.txt > ${sample_id}_AMTPTU_UTR_M2.fasta
sleep 2
rm tmp_list.txt AMTPTU_UTR.fasta
cat *M2.fasta | grep ">" | cut -d ">" -f 2 | awk '!seen[\$0]++' > ${sample_id}_list_all_headers_M2.txt
"""
}

process create_table_misMatch0{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/tables", mode:'link'

  input:
    set sample_id, list, AMTS_UTR, AMTAS_UTR, AMTPTU_UTR from misMatch0_out
  
  output:
  set sample_id, "${sample_id}_table_misMatch0.tsv", list into misMatch0_table_out

  script:
""" 
#!/bin/python3

fw = open("${sample_id}_table_misMatch0.tsv", "a")
header = "sample_id" + "\\t" + "seq_id" + "\\t" + "AMTS_UTR" + "\\t" + "AMTAS_UTR" + "\\t" + "AMTPTU_UTR" + "\\n"
fw.write(header)
sampleid = "${sample_id}"
f = open("${list}", "r")
f_list = f.readlines()
f.close()
f = open("${AMTS_UTR}", "r")
AMTS = f.readlines()
f.close()
f = open("${AMTAS_UTR}", "r")
AMTAS = f.readlines()
f.close()
f = open("${AMTPTU_UTR}", "r")
AMTPTU = f.readlines()
f.close()
for seq in f_list:
    seq_id = seq.split()[0].replace('\\n', '')
    AMTS_prime = 0
    AMTAS_prime = 0
    AMTPTU_prime = 0
    for prim_seq in AMTS:
        if seq_id.replace('\\n', '') == prim_seq.replace('\\n', '').replace('>', ''):
            AMTS_prime = 1
    for prim_seq in AMTAS:
        if seq_id.replace('\\n', '') == prim_seq.replace('\\n', '').replace('>', ''):
            AMTAS_prime = 1
    for prim_seq in AMTPTU:
        if seq_id.replace('\\n', '') == prim_seq.replace('\\n', '').replace('>', ''):
            AMTPTU_prime = 1

    str_line = str(sampleid) + "\\t" + str(seq_id) + "\\t" + str(AMTS_prime) + "\\t" + str(AMTAS_prime) + "\\t" + str(AMTPTU_prime) + "\\n"
    fw.write(str_line)
fw.close()
"""
}

process create_table_misMatch1{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/tables", mode:'link'

  input:
    set sample_id, list, AMTS_UTR, AMTAS_UTR, AMTPTU_UTR from misMatch1_out
  
  output:
  set sample_id, "${sample_id}_table_misMatch1.tsv", list into misMatch1_table_out

  script:
""" 
#!/bin/python3

fw = open("${sample_id}_table_misMatch1.tsv", "a")
header = "sample_id" + "\\t" + "seq_id" + "\\t" + "AMTS_UTR" + "\\t" + "AMTAS_UTR" + "\\t" + "AMTPTU_UTR" + "\\n"
fw.write(header)
sampleid = "${sample_id}"
f = open("${list}", "r")
f_list = f.readlines()
f.close()
f = open("${AMTS_UTR}", "r")
AMTS = f.readlines()
f.close()
f = open("${AMTAS_UTR}", "r")
AMTAS = f.readlines()
f.close()
f = open("${AMTPTU_UTR}", "r")
AMTPTU = f.readlines()
f.close()

for seq in f_list:
    seq_id = seq.split()[0].replace('\\n', '')
    AMTS_prime = 0
    AMTAS_prime = 0
    AMTPTU_prime = 0
    for prim_seq in AMTS:
        if seq_id.replace('\\n', '') == prim_seq.replace('\\n', '').replace('>', ''):
            AMTS_prime = 1
    for prim_seq in AMTAS:
        if seq_id.replace('\\n', '') == prim_seq.replace('\\n', '').replace('>', ''):
            AMTAS_prime = 1
    for prim_seq in AMTPTU:
        if seq_id.replace('\\n', '') == prim_seq.replace('\\n', '').replace('>', ''):
            AMTPTU_prime = 1

    str_line = str(sampleid) + "\\t" + str(seq_id) + "\\t" + str(AMTS_prime) + "\\t" + str(AMTAS_prime) + "\\t" + str(AMTPTU_prime) + "\\n"
    fw.write(str_line)
fw.close()
"""
}

process create_table_misMatch2{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/tables", mode:'link'

  input:
    set sample_id, list, AMTS_UTR, AMTAS_UTR, AMTPTU_UTR from misMatch2_out
  
  output:
  set sample_id, "${sample_id}_table_misMatch2.tsv", list into misMatch2_table_out

  script:
""" 
#!/bin/python3

fw = open("${sample_id}_table_misMatch2.tsv", "a")
header = "sample_id" + "\\t" + "seq_id" + "\\t" + "AMTS_UTR" + "\\t" + "AMTAS_UTR" + "\\t" + "AMTPTU_UTR" + "\\n"
fw.write(header)
sampleid = "${sample_id}"
f = open("${list}", "r")
f_list = f.readlines()
f.close()
f = open("${AMTS_UTR}", "r")
AMTS = f.readlines()
f.close()
f = open("${AMTAS_UTR}", "r")
AMTAS = f.readlines()
f.close()
f = open("${AMTPTU_UTR}", "r")
AMTPTU = f.readlines()
f.close()

for seq in f_list:
    seq_id = seq.split()[0].replace('\\n', '')
    AMTS_prime = 0
    AMTAS_prime = 0
    AMTPTU_prime = 0
    for prim_seq in AMTS:
        if seq_id.replace('\\n', '') == prim_seq.replace('\\n', '').replace('>', ''):
            AMTS_prime = 1
    for prim_seq in AMTAS:
        if seq_id.replace('\\n', '') == prim_seq.replace('\\n', '').replace('>', ''):
            AMTAS_prime = 1
    for prim_seq in AMTPTU:
        if seq_id.replace('\\n', '') == prim_seq.replace('\\n', '').replace('>', ''):
            AMTPTU_prime = 1

    str_line = str(sampleid) + "\\t" + str(seq_id) + "\\t" + str(AMTS_prime) + "\\t" + str(AMTAS_prime) + "\\t" + str(AMTPTU_prime) + "\\n"
    fw.write(str_line)
fw.close()
"""
}

misMatch0_table_out.into{
misMatch0_table_res_plot_in;
misMatch0_combine_all_in;
misMatch0_mix_misMatches_in;
}

misMatch1_table_out.into{
misMatch1_table_res_plot_in;
misMatch1_combine_all_in;
misMatch1_mix_misMatches_in;
}

misMatch2_table_out.into{
misMatch2_table_res_plot_in;
misMatch2_combine_all_in;
misMatch2_mix_misMatches_in;
}

tmp = misMatch0_mix_misMatches_in.combine(misMatch1_mix_misMatches_in, by: 0)
all_misMatches_in = tmp.combine(misMatch2_mix_misMatches_in, by: 0)

process create_table_all_misMatches{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/tables", mode:'link'

  input:
  set sample_id, table_M0, list_M0, table_M1, list_M1, table_M2, list_M2 from all_misMatches_in
  
  output:
  set sample_id, "${sample_id}_table_misMatches_combined.tsv", list_M2 into all_misMatches_table
  set sample_id, list_M0, list_M1, list_M2 into nr_of_matched_seq

  script:
""" 
#!/bin/python3

def isEmpty(myString):
    if myString and myString.strip():
        #Is not empty
        return False
    #Is empty
    return True

fw = open("${sample_id}_table_misMatches_combined.tsv", "a")
header = "sample_id" + "\\t" + "seq_id" + "\\t" + "AMTS_UTR" + "\\t" + "AMTAS_UTR" + "\\t" + "AMTPTU_UTR" + "\\n"
fw.write(header)
f = open("${list_M2}", "r")
f_list = f.readlines()
f.close()
f = open("${table_M0}", "r")
table0 = f.readlines()
f.close()
f = open("${table_M1}", "r")
table1 = f.readlines()
f.close()
f = open("${table_M2}", "r")
table2 = f.readlines()
f.close()

sampleid = "${sample_id}"

for seq in f_list:
    seq_id = seq.split()[0].replace('\\n', '')
    AMTS_prime = "-"
    AMTAS_prime = "-"
    AMTPTU_prime = "-"
    for t0 in table0:
        if isEmpty(t0) != True:
            if str(seq_id) == str(t0.split()[1]):
                if int(t0.split()[2]) != 0:
                    AMTS_prime = "M0"
                if int(t0.split()[3]) != 0:
                    AMTAS_prime = "M0"
                if int(t0.split()[4]) != 0:
                    AMTPTU_prime = "M0"
    for t0 in table1:
        if isEmpty(t0) != True:
            if str(seq_id) == str(t0.split()[1]):
                if int(t0.split()[2]) != 0 and AMTS_prime != "M0":
                    AMTS_prime = "M1"
                if int(t0.split()[3]) != 0 and AMTAS_prime != "M0":
                    AMTAS_prime = "M1"
                if int(t0.split()[4]) != 0 and AMTPTU_prime != "M0":
                    AMTPTU_prime = "M1"
    for t0 in table2:
        if isEmpty(t0) != True:
            if str(seq_id) == str(t0.split()[1]):
                if int(t0.split()[2]) != 0 and AMTS_prime != "M0" and AMTS_prime != "M1":
                    AMTS_prime = "M2"
                if int(t0.split()[3]) != 0 and AMTAS_prime != "M0" and AMTAS_prime != "M1":
                    AMTAS_prime = "M2"
                if int(t0.split()[4]) != 0 and AMTPTU_prime != "M0" and AMTPTU_prime != "M1":
                    AMTPTU_prime = "M2"
    str_line = str(sampleid) + "\\t" + str(seq_id) + "\\t" + str(AMTS_prime) + "\\t" + str(AMTAS_prime) + "\\t" + str(AMTPTU_prime) + "\\n"
    fw.write(str_line)

fw.close()
"""
}

process collect_misMatch_tables_M0{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tables", mode:'link'

input:
file tables from misMatch0_combine_all_in.map{it[1]}.collect()

output:

file "all_samples_collected_table_misMatches_M0.tsv" into collected_misMatch_table_M0

script:
"""
IFS=' ' read -r -a arr <<< \$(echo ${tables})
readarray -t sortedArr < <(for i in "\${arr[@]}"; do echo "\$i"; done | sort)
first_elem=true
for i in "\${sortedArr[@]}"
do
	if [ "\$first_elem" = true ] ; then
		cat \$i >> "all_samples_collected_table_misMatches_M0.tsv"
		first_elem=false
	else
		cat \$i | tail -n +2 >> "all_samples_collected_table_misMatches_M0.tsv"
	fi
   
done
"""
}

process collect_misMatch_tables_M1{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tables", mode:'link'

input:
file tables from misMatch1_combine_all_in.map{it[1]}.collect()

output:

file "all_samples_collected_table_misMatches_M1.tsv" into collected_misMatch_table_M1

script:
"""
IFS=' ' read -r -a arr <<< \$(echo ${tables})
readarray -t sortedArr < <(for i in "\${arr[@]}"; do echo "\$i"; done | sort)
first_elem=true
for i in "\${sortedArr[@]}"
do
	if [ "\$first_elem" = true ] ; then
		cat \$i >> "all_samples_collected_table_misMatches_M1.tsv"
		first_elem=false
	else
		cat \$i | tail -n +2 >> "all_samples_collected_table_misMatches_M1.tsv"
	fi
   
done
"""
}

process collect_misMatch_tables_M2{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tables", mode:'link'

input:
file tables from misMatch2_combine_all_in.map{it[1]}.collect()

output:

file "all_samples_collected_table_misMatches_M2.tsv" into collected_misMatch_table_M2

script:
"""
IFS=' ' read -r -a arr <<< \$(echo ${tables})
readarray -t sortedArr < <(for i in "\${arr[@]}"; do echo "\$i"; done | sort)
first_elem=true
for i in "\${sortedArr[@]}"
do
	if [ "\$first_elem" = true ] ; then
		cat \$i >> "all_samples_collected_table_misMatches_M2.tsv"
		first_elem=false
	else
		cat \$i | tail -n +2 >> "all_samples_collected_table_misMatches_M2.tsv"
	fi
   
done
"""
}

process collect_misMatch_tables_comb{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tables", mode:'link'

input:
file tables from all_misMatches_table.map{it[1]}.collect()

output:

file "all_samples_collected_table_misMatches_combined.tsv" into collected_misMatch_table_comb

script:
"""
IFS=' ' read -r -a arr <<< \$(echo ${tables})
readarray -t sortedArr < <(for i in "\${arr[@]}"; do echo "\$i"; done | sort)
first_elem=true
for i in "\${sortedArr[@]}"
do
	if [ "\$first_elem" = true ] ; then
		cat \$i >> "all_samples_collected_table_misMatches_combined.tsv"
		first_elem=false
	else
		cat \$i | tail -n +2 >> "all_samples_collected_table_misMatches_combined.tsv"
	fi
   
done
"""
}


nr_of_matched_seq.into{
nr_of_matched_seq_M0;
nr_of_matched_seq_M1;
nr_of_matched_seq_M2
}

process collect_nr_matched_seq_M0{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/lists", mode:'link'

input:
file lists from nr_of_matched_seq_M0.map{it[1]}.collect()

output:

set "all_samples_collected_listsM0.txt", "nr_of_matched_seq_M0.txt" into collected_lists_and_nr_of_matched0

script:
"""
IFS=' ' read -r -a arr <<< \$(echo ${lists})
readarray -t sortedArr < <(for i in "\${arr[@]}"; do echo "\$i"; done | sort)
first_elem=true
for i in "\${sortedArr[@]}"
do
	cat \$i >> "all_samples_collected_listsM0.txt"
   
done

cat "all_samples_collected_listsM0.txt" | awk '!seen[\$0]++' | wc -l > "nr_of_matched_seq_M0.txt"
"""
}

process collect_nr_matched_seq_M1{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/lists", mode:'link'

input:
file lists from nr_of_matched_seq_M1.map{it[2]}.collect()

output:

set "all_samples_collected_listsM1.txt", "nr_of_matched_seq_M1.txt" into collected_lists_and_nr_of_matched1

script:
"""
IFS=' ' read -r -a arr <<< \$(echo ${lists})
readarray -t sortedArr < <(for i in "\${arr[@]}"; do echo "\$i"; done | sort)
for i in "\${sortedArr[@]}"
do
	cat \$i >> "all_samples_collected_listsM1.txt"
   
done

cat "all_samples_collected_listsM1.txt" | awk '!seen[\$0]++' | wc -l > "nr_of_matched_seq_M1.txt"
"""
}

process collect_nr_matched_seq_M2{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/lists", mode:'link'

input:
file lists from nr_of_matched_seq_M2.map{it[3]}.collect()

output:

set "all_samples_collected_listsM2.txt", "nr_of_matched_seq_M2.txt" into collected_lists_and_nr_of_matched2

script:
"""
IFS=' ' read -r -a arr <<< \$(echo ${lists})
readarray -t sortedArr < <(for i in "\${arr[@]}"; do echo "\$i"; done | sort)
for i in "\${sortedArr[@]}"
do
	cat \$i >> "all_samples_collected_listsM2.txt"
   
done

cat "all_samples_collected_listsM2.txt" | awk '!seen[\$0]++' | wc -l > "nr_of_matched_seq_M2.txt"
"""
}

collected_lists_and_nr_of_matched0.into{
collected_lists_and_nr_of_matched0_comb;
collected_lists_and_nr_of_matched0_plot
}

collected_lists_and_nr_of_matched1.into{
collected_lists_and_nr_of_matched1_comb;
collected_lists_and_nr_of_matched1_plot
}

collected_lists_and_nr_of_matched2.into{
collected_lists_and_nr_of_matched2_comb;
collected_lists_and_nr_of_matched2_plot
}

collected_lists_and_nr_of_matched_comb = collected_lists_and_nr_of_matched0_comb.combine(collected_lists_and_nr_of_matched1_comb).combine(collected_lists_and_nr_of_matched2_comb)


plot_M0_in = collected_misMatch_table_M0.combine(total_nrOfSeq_M0_in).combine(collected_lists_and_nr_of_matched0_plot)
plot_M1_in = collected_misMatch_table_M1.combine(total_nrOfSeq_M1_in).combine(collected_lists_and_nr_of_matched1_plot)
plot_M2_in = collected_misMatch_table_M2.combine(total_nrOfSeq_M2_in).combine(collected_lists_and_nr_of_matched2_plot)
plot_comb_in = collected_misMatch_table_comb.combine(toatl_nrOfSeq_all_in).combine(collected_lists_and_nr_of_matched_comb)

process plot_comb{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/plots", mode:'link'

input:
set table, nrOfSeq, list0, nrOfUniqueSeq0, list1, nrOfUniqueSeq1, list2, nrOfUniqueSeq2 from plot_comb_in

output:

set "plot_comb.png", "plot_comb_with_total_seq.png" into plot_comb_out

script:
"""
#!/bin/python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

f = open("${nrOfUniqueSeq0}", "r")
uniqueSeq0 = int(f.readline().replace("\\n", ""))
f.close()
f = open("${nrOfUniqueSeq1}", "r")
uniqueSeq1 = int(f.readline().replace("\\n", ""))
f.close()
f = open("${nrOfUniqueSeq2}", "r")
uniqueSeq2 = int(f.readline().replace("\\n", ""))
f.close()

uniqueSeq2 = uniqueSeq2 - uniqueSeq1
uniqueSeq1 = uniqueSeq1 - uniqueSeq0

f = open("${table}", "r")
header = f.readline().replace("\\n", "").split('\t')[2:]
f_table = f.readlines()
f.close()
f = open("${nrOfSeq}", "r")
total_seq = int(f.readline())
f.close()
primer_count_arr = [0]*len(header)
primer_count_arr0 = [0]*len(header)
primer_count_arr1 = [0]*len(header)
primer_count_arr2 = [0]*len(header)

for line in f_table:
    tmp = line.replace("\\n", "").split('\t')[2:]
    for i in range(len(tmp)):
        if "M0" in tmp[i]:
            primer_count_arr0[i] = primer_count_arr0[i] + 1
            primer_count_arr[i] = primer_count_arr[i] + 1
        if "M1" in tmp[i]:
            primer_count_arr1[i] = primer_count_arr1[i] + 1
            primer_count_arr[i] = primer_count_arr[i] + 1
        if "M2" in tmp[i]:
            primer_count_arr2[i] = primer_count_arr2[i] + 1
            primer_count_arr[i] = primer_count_arr[i] + 1

for i in range(len(primer_count_arr)): 
    for j in range(0, len(primer_count_arr)-i-1): 
        if primer_count_arr[j] < primer_count_arr[j+1]: 
            primer_count_arr[j], primer_count_arr[j+1] = primer_count_arr[j+1], primer_count_arr[j]
            primer_count_arr0[j], primer_count_arr0[j+1] = primer_count_arr0[j+1], primer_count_arr0[j]
            primer_count_arr1[j], primer_count_arr1[j+1] = primer_count_arr1[j+1], primer_count_arr1[j]  
            primer_count_arr2[j], primer_count_arr2[j+1] = primer_count_arr2[j+1], primer_count_arr2[j]  
            header[j], header[j+1] = header[j+1], header[j] 

dataset = np.array(primer_count_arr)
dataset0 = np.array(primer_count_arr0)
dataset1 = np.array(primer_count_arr1)
dataset2 = np.array(primer_count_arr2)

fig, ax = plt.subplots(figsize=(20,10))
width = 0.35 
ax.bar(header, dataset0, width, label='M0')
ax.bar(header, dataset1, width, bottom=dataset0, label='M1')
ax.bar(header, dataset2, width, bottom=dataset0+dataset1, label='M2')
ax.set_xticklabels(header,rotation='vertical')
ax.set_ylabel('Number of occurrences')
ax.set_title('Number of occurrences of the most common primer sequences')
ax.legend()
plt.subplots_adjust(bottom=0.24)
plt.savefig("plot_comb.png")


header.append("total_seq")
primer_count_arr.append(total_seq)
primer_count_arr0.append(total_seq)
primer_count_arr1.append(0)
primer_count_arr2.append(0)
header.append("total_nr_unique_seq")
all_unique = uniqueSeq0 + uniqueSeq1 + uniqueSeq2
primer_count_arr.append(all_unique)
primer_count_arr0.append(uniqueSeq0)
primer_count_arr1.append(uniqueSeq1)
primer_count_arr2.append(uniqueSeq2)

for i in range(len(primer_count_arr)): 
    for j in range(0, len(primer_count_arr)-i-1): 
        if primer_count_arr[j] < primer_count_arr[j+1]: 
            primer_count_arr[j], primer_count_arr[j+1] = primer_count_arr[j+1], primer_count_arr[j]
            primer_count_arr0[j], primer_count_arr0[j+1] = primer_count_arr0[j+1], primer_count_arr0[j]
            primer_count_arr1[j], primer_count_arr1[j+1] = primer_count_arr1[j+1], primer_count_arr1[j]  
            primer_count_arr2[j], primer_count_arr2[j+1] = primer_count_arr2[j+1], primer_count_arr2[j]  
            header[j], header[j+1] = header[j+1], header[j] 

_dataset = np.array(primer_count_arr)
_dataset0 = np.array(primer_count_arr0)
_dataset1 = np.array(primer_count_arr1)
_dataset2 = np.array(primer_count_arr2)

fig, ax = plt.subplots(figsize=(20,10))
width = 0.35 
ax.bar(header, _dataset0, width, label='M0')
ax.bar(header, _dataset1, width, bottom=_dataset0, label='M1')
ax.bar(header, _dataset2, width, bottom=_dataset0+_dataset1, label='M2')
ax.set_xticklabels(header,rotation='vertical')
ax.set_ylabel('Number of occurrences')
ax.set_title('Number of occurrences of the most common primer sequences')
ax.legend()
plt.subplots_adjust(bottom=0.24)
plt.savefig("plot_comb_with_total_seq.png")
"""
}

process plot_M0{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/plots", mode:'link'

input:
set table, nrOfSeq, list, nrOfUniqueSeq from plot_M0_in

output:

set "plot_M0.png", "plot_M0_with_total_seq.png" into plot_M0_out

script:
"""
#!/bin/python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

f = open("${nrOfUniqueSeq}", "r")
uniqueSeq = int(f.readline().replace("\\n", ""))
f.close()
f = open("${table}", "r")
header = f.readline().replace("\\n", "").split('\t')[2:]
f_table = f.readlines()
f.close()
f = open("${nrOfSeq}", "r")
total_seq = int(f.readline())
f.close()
primer_count_arr = [0]*len(header)
for line in f_table:
    tmp = line.replace("\\n", "").split('\t')[2:]
    for i in range(len(tmp)):
        if int(tmp[i]) != 0:
            primer_count_arr[i] = primer_count_arr[i] + 1

for i in range(len(primer_count_arr)): 
    for j in range(0, len(primer_count_arr)-i-1): 
        if primer_count_arr[j] < primer_count_arr[j+1]: 
            primer_count_arr[j], primer_count_arr[j+1] = primer_count_arr[j+1], primer_count_arr[j] 
            header[j], header[j+1] = header[j+1], header[j] 

x = np.arange(len(header))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(20,10))
rects1 = ax.bar(x, primer_count_arr, width)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of occurrences')
ax.set_title('Number of occurrences of the most common primer sequences')
ax.set_xticks(x)
ax.set_xticklabels(header,rotation='vertical')
ax.legend().remove()

def autolabel(rects):
    #Attach a text label above each bar in *rects*, displaying its height.
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
fig.tight_layout()
plt.savefig("plot_M0.png")

header.append("total_seq")
primer_count_arr.append(total_seq)
header.append("total_nr_unique_seq")
primer_count_arr.append(uniqueSeq)

for i in range(len(primer_count_arr)): 
    for j in range(0, len(primer_count_arr)-i-1): 
        if primer_count_arr[j] < primer_count_arr[j+1]: 
            primer_count_arr[j], primer_count_arr[j+1] = primer_count_arr[j+1], primer_count_arr[j] 
            header[j], header[j+1] = header[j+1], header[j] 

x = np.arange(len(header))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(20,10))
rects1 = ax.bar(x, primer_count_arr, width)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of occurrences')
ax.set_title('Number of occurrences of the most common primer sequences')
ax.set_xticks(x)
ax.set_xticklabels(header,rotation='vertical')
ax.legend().remove()

def autolabel(rects):
    #Attach a text label above each bar in *rects*, displaying its height.
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
fig.tight_layout()
plt.savefig("plot_M0_with_total_seq.png")
"""
}

process plot_M1{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/plots", mode:'link'

input:
set table, nrOfSeq, list, nrOfUniqueSeq from plot_M1_in

output:

set "plot_M1.png", "plot_M1_with_total_seq.png" into plot_M1_out

script:
"""
#!/bin/python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

f = open("${nrOfUniqueSeq}", "r")
uniqueSeq = int(f.readline().replace("\\n", ""))
f.close()
f = open("${table}", "r")
header = f.readline().replace("\\n", "").split('\t')[2:]
f_table = f.readlines()
print(f_table)
print(header)
f.close()
f = open("${nrOfSeq}", "r")
total_seq = int(f.readline())
f.close()
primer_count_arr = [0]*len(header)
for line in f_table:
    tmp = line.replace("\\n", "").split('\t')[2:]
    for i in range(len(tmp)):
        if int(tmp[i]) != 0:
            primer_count_arr[i] = primer_count_arr[i] + 1

for i in range(len(primer_count_arr)): 
    for j in range(0, len(primer_count_arr)-i-1): 
        if primer_count_arr[j] < primer_count_arr[j+1]: 
            primer_count_arr[j], primer_count_arr[j+1] = primer_count_arr[j+1], primer_count_arr[j] 
            header[j], header[j+1] = header[j+1], header[j] 

x = np.arange(len(header))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(20,10))
rects1 = ax.bar(x, primer_count_arr, width)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of occurrences')
ax.set_title('Number of occurrences of the most common primer sequences')
ax.set_xticks(x)
ax.set_xticklabels(header,rotation='vertical')
ax.legend().remove()

def autolabel(rects):
    #Attach a text label above each bar in *rects*, displaying its height.
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
fig.tight_layout()
plt.savefig("plot_M1.png")

header.append("total_seq")
primer_count_arr.append(total_seq)
header.append("total_nr_unique_seq")
primer_count_arr.append(uniqueSeq)

for i in range(len(primer_count_arr)): 
    for j in range(0, len(primer_count_arr)-i-1): 
        if primer_count_arr[j] < primer_count_arr[j+1]: 
            primer_count_arr[j], primer_count_arr[j+1] = primer_count_arr[j+1], primer_count_arr[j] 
            header[j], header[j+1] = header[j+1], header[j] 

x = np.arange(len(header))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(20,10))
rects1 = ax.bar(x, primer_count_arr, width)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of occurrences')
ax.set_title('Number of occurrences of the most common primer sequences')
ax.set_xticks(x)
ax.set_xticklabels(header,rotation='vertical')
ax.legend().remove()

def autolabel(rects):
    #Attach a text label above each bar in *rects*, displaying its height.
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
fig.tight_layout()
plt.savefig("plot_M1_with_total_seq.png")
"""
}

process plot_M2{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/plots", mode:'link'

input:
set table, nrOfSeq, list, nrOfUniqueSeq from plot_M2_in

output:

set "plot_M2.png", "plot_M2_with_total_seq.png" into plot_M2_out

script:
"""
#!/bin/python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

f = open("${nrOfUniqueSeq}", "r")
uniqueSeq = int(f.readline().replace("\\n", ""))
f.close()
f = open("${table}", "r")
header = f.readline().replace("\\n", "").split('\t')[2:]
f_table = f.readlines()
print(f_table)
print(header)
f.close()
f = open("${nrOfSeq}", "r")
total_seq = int(f.readline())
f.close()
primer_count_arr = [0]*len(header)
for line in f_table:
    tmp = line.replace("\\n", "").split('\t')[2:]
    for i in range(len(tmp)):
        if int(tmp[i]) != 0:
            primer_count_arr[i] = primer_count_arr[i] + 1

for i in range(len(primer_count_arr)): 
    for j in range(0, len(primer_count_arr)-i-1): 
        if primer_count_arr[j] < primer_count_arr[j+1]: 
            primer_count_arr[j], primer_count_arr[j+1] = primer_count_arr[j+1], primer_count_arr[j] 
            header[j], header[j+1] = header[j+1], header[j] 

x = np.arange(len(header))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(20,10))
rects1 = ax.bar(x, primer_count_arr, width)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of occurrences')
ax.set_title('Number of occurrences of the most common primer sequences')
ax.set_xticks(x)
ax.set_xticklabels(header,rotation='vertical')
ax.legend().remove()

def autolabel(rects):
    #Attach a text label above each bar in *rects*, displaying its height.
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
fig.tight_layout()
plt.savefig("plot_M2.png")

header.append("total_seq")
primer_count_arr.append(total_seq)
header.append("total_nr_unique_seq")
primer_count_arr.append(uniqueSeq)

for i in range(len(primer_count_arr)): 
    for j in range(0, len(primer_count_arr)-i-1): 
        if primer_count_arr[j] < primer_count_arr[j+1]: 
            primer_count_arr[j], primer_count_arr[j+1] = primer_count_arr[j+1], primer_count_arr[j] 
            header[j], header[j+1] = header[j+1], header[j] 

x = np.arange(len(header))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(20,10))
rects1 = ax.bar(x, primer_count_arr, width)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of occurrences')
ax.set_title('Number of occurrences of the most common primer sequences')
ax.set_xticks(x)
ax.set_xticklabels(header,rotation='vertical')
ax.legend().remove()

def autolabel(rects):
    #Attach a text label above each bar in *rects*, displaying its height.
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
fig.tight_layout()
plt.savefig("plot_M2_with_total_seq.png")
"""
}




