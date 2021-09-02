#!/bin/bash

### Use this script if you want to compute the total length of all exons of the genes in a gtf file
### you should provide the ANNOT_GTF file bellow

ANNOT_GTF=Mus_musculus.GRCm38.84.chr.gtf  ### your favorate gtf file
FILE_OUT=Mus_musculus.GRCm38.84.all_exons.merged.sizes.txt
if [ -f $FILE_OUT ]; then
	rm -v $FILE_OUT
fi

ANNOT_GTF_exons=Mus_musculus.GRCm38.84.all_exons.txt
awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="exon"){split($9,fld,";");split(fld[1],gene," ");print "chr"$1,$4,$5,gene[2]}}' $ANNOT_GTF | sed 's/"//g' | sort -k 1,1 -k 2,2n -k 3,3n > $ANNOT_GTF_exons

ANNOT_GTF_exons_merged=Mus_musculus.GRCm38.84.all_exons.merged.txt
if [ -f $ANNOT_GTF_exons_merged ]; then
	rm -v $ANNOT_GTF_exons_merged
fi
for gene in `awk '{print $4}' $ANNOT_GTF_exons | sort -k 4 | uniq`; do grep $gene $ANNOT_GTF_exons | bedtools merge -c 4 -o distinct >> $ANNOT_GTF_exons_merged; done


for gene in `awk '{print $4}' $ANNOT_GTF_exons_merged | sort -u`; do grep $gene $ANNOT_GTF_exons_merged | awk 'BEGIN{sum=0;OFS="\t"}{sum+=($3-$2)}END{print $4,sum}' >> $FILE_OUT; done
