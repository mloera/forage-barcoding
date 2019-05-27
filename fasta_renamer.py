########################################################################################################################################
# This script is meant to modify the name of each entry of the original BOLD-downloaded sequence headers.
#
# The BOLD-downloaded sequence headers have this format:
# >SEQUENCE ID | TAXON | SAMPLE ID | FAMILY | SPECIES (with spaces)
#
# e.g., :
# >UHURU053-14|Rhynchosia malacophylla|E7_K1207_Unknown_sp|Fabaceae|Rhynchosia malacophylla
#
# This script reduces the headers to:
#>SEQUENCE ID | FAMILY | SPECIES (no spaces)
#
# Use as follows:
# python fasta_renamer.py input.fasta > output.fasta
########################################################################################################################################

#!/usr/bin/python

import sys
from Bio import SeqIO

infasta = sys.argv[1]

for r in SeqIO.parse(infasta, 'fasta'):

        s_id = r.description.split('|')[0]

        family = r.description.split('|')[-2]
        
        species = '_'.join(r.description.split('|')[-1].split('_')[:2])

        

        print('>%s|%s|%s\n%s' % (s_id, family, species, str(r.seq).replace('-','')))  
