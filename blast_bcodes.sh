

for BCODE in trnH-psbA matK rbcLa ITS2 trnL; do  \
makeblastdb -in ${BCODE}-shname-may2019.SWFRG.fas \
  -input_type fasta -dbtype nucl -title ${BCODE} \
  -out ${BCODE}.SWFRG; 
done;

for BCODE in trnH-psbA matK rbcLa ITS2 trnL; do \
blastn -query ${BCODE}-shname-may2019.SWFRG.fas \
  -db ${BCODE}.SWFRG  \
  -max_target_seqs 5 \
  -outfmt 6 \
  -out ${BCODE}.blastn.SWFRG; 
done;
