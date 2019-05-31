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

## Step 1: Re-formatting the headers of BOLD sequences
This step is meant to reformat the headers of the downloaded BOLD sequences in order to have them in a blastn-friendly configuration (i.e., no spaces). The fasta_name_reformat.py script found in this repository was used for this purpose:

<code>python fasta_name_reformat.py barcode.fasta > barcode.reformatted.fasta</code>

## Step 2: blastn
A blast database was built for each reformatted fasta file.

<code>for BCODE in trnH-psbA matK rbcLa; do  \\
makeblastdb -in ${BCODE}.reformatted.fasta \
  -input_type fasta \
  -dbtype nucl \
  -title ${BCODE} \
  -out ${BCODE}.SWFRG; 
done;</code>

Then, each reformatted fasta file was blasted against its own blast database. The blast output files are in a tabular format (outfmt = 6).

<code>for BCODE in trnH-psbA matK rbcLa; do \
blastn -query ${BCODE}-shname-may2019.SWFRG.fas \
  -db ${BCODE}.SWFRG  \
  -max_target_seqs 5 \
  -outfmt 6 \
  -out ${BCODE}.blastn.SWFRG; 
done;</code>

## Step 3: Match calling
