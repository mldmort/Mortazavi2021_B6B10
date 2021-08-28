#!/bin/bash

echo "+++++++++++++++++++ SV ++++++++++++++++"
FDR_FILE=fdr_sv_inDEG_b6b10.txt
awk 'BEGIN{FS="\t";OFS="\t"}$15<0.05 && $5=="LUMPY_CALL"{print $8,"SV_"$4,$9,$1,$2,$3,$15}' $FDR_FILE | sort -k4,4 -V

echo "+++++++++++++++++++ SNP ++++++++++++++++"
FDR_FILE=fdr_vep_inDEG_b6b10.txt
awk 'BEGIN{FS="\t";OFS="\t"}$19<0.05{print $11,$6,$9,$1,$2,$3,$19}' $FDR_FILE | sort -k5,5 -V

echo "+++++++++++++++++++ CNV ++++++++++++++++"
FDR_FILE=fdr_cnv_inDEG_b6b10.txt
awk 'BEGIN{FS="\t";OFS="\t"}$37<0.05{print $4,$1,$2,$3,$37}' $FDR_FILE | sort -k2,2 -k1,1 -u -V
