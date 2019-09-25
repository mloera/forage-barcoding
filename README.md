 # forage-barcoding
Forage DNA barcoding: data analysis workflow

This workflow was used to analyze the DNA barcoding data from 16 forage plant species  (SWFRG project at BOLD v4). The workflow involves the following steps:
1. Downloading sequences from BOLD
2. Re-formatting the headers of BOLD sequences
3. blastn
4. Match calling with a "Leave-One-Out" approach

## Dependencies:

### python 2.7
BioPython 1.68

### R version 3.5.1 (2018-07-02)
tidyverse 1.2.1,
stringr 1.4.0

### blast 2.2.30

## Step 1: Downloading sequences from BOLD
Poaceae and Fabaceae sequences for barcodes <i>matK</i>, <i>rbcLa</i>, and <i>trnH-psbA</i> were downloaded from BOLD v4. Only sequences without reported contaminants and with sequence length > 200 bp were downloaded.

The downloaded fasta files are located in the docs folder:

1. docs/matK-may2019.fas
2. docs/rbcLa-may2019.fas
3. docs/trnH-psbA-may2019.fas

## Step 1: Re-formatting the headers of BOLD sequences
This step is meant to reformat the headers of the downloaded BOLD sequences in order to have them in a blastn-friendly configuration (i.e., no spaces). The fasta_name_reformat.py script found in this repository was used for this purpose:

<code>python fasta_name_reformat.py barcode.fasta > barcode.reformatted.fasta</code>

The reformatted fasta files are also located in the docs folder:

1. docs/matK-shname.SWFRG.fas
2. docs/rbcLa-shname.SWFRG.fas
3. docs/trnH-psbA-shname.SWFRG.fas


## Step 2: blastn
A blast database was built for each reformatted fasta file.

<code>for BCODE in trnH-psbA matK rbcLa; do  \
makeblastdb -in ${BCODE}.reformatted.fasta \
  -input_type fasta \
  -dbtype nucl \
  -title ${BCODE} \
  -out ${BCODE}.SWFRG; 
done;</code>

Then, a fasta file was created for the SWFRG sequences of each barcode. 

<code>for BCODE in trnH-psbA matK rbcLa; do \
 
 grep -A1 'SWFRG' ${BCODE}-shname.fas > ${BCODE}.SWFRG.only.fasta</code>

Each SWFRG fasta file was then blasted against its corresponding blast database with the flag max_target_seqs = 5. The blast output files are in a tabular format (outfmt = 6).

<code>for BCODE in trnH-psbA matK rbcLa; do \
blastn -query ${BCODE}.SWFRG.only.fasta \
  -db ${BCODE}.SWFRG  \
  -max_target_seqs 5 \
  -outfmt 6 \
  -out ${BCODE}.blastn.SWFRG; 
done;</code>

The blastn output tables are also located in the docs folder:

1. docs/matK.blastn.SWFRG
2. docs/rbcLa.blastn.SWFRG
3. docs/trnH-psbA.blastn.SWFRG

## Step 3: Match calling on the SWFRG sequences
The blast output tables were parsed with the blastn_matcher.R script, ran on Rstudio. The script removes self-hits (i.e., blast hits between the same sequence ID), and also corrects some misspellings in the species name of queries and database hits. The script then compares the taxonomy of the query and the database hits at the species- and genus-level. A match is called when the taxonomy of a query sequence matches the taxonomy of the highest scoring hits.
