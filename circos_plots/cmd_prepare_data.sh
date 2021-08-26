#!/bin/bash

# for this commend you need to have vcftools

DIR_DATA=../data
DIR_DATA_CIR=files

if [ ! -d $DIR_DATA_CIR ]; then
	mkdir -v $DIR_DATA_CIR
fi
 
#====================== SNPs ======================
VCF_SNP_NO_HET=$DIR_DATA/speedseqsnp30x_merged_realigned_region_snp_no_Het.vcf
vcftools --vcf $VCF_SNP_NO_HET --SNPdensity 1000000 --out $DIR_DATA_CIR/SNP_density_NO_HET

SNP_IN=$DIR_DATA_CIR/SNP_density_NO_HET.snpden
SNP_OUT=$DIR_DATA_CIR/snp_density_no_het.txt
echo "snp file: $SNP_OUT"
awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+1e6,$3}}' $SNP_IN > $SNP_OUT

# strain specific SNPs
STRAINS="R2 R52 R72 R84 R104 R114 R124 R134 R144 R22 R32 R42 R62 R94"
for strain in $STRAINS;
do
	sample=$strain\_merged_realigned.bam
	echo $sample
	strain_vcf=$DIR_DATA_CIR/$strain\_SNP_no_Het.vcf
	bcftools view -s $sample $VCF_SNP_NO_HET | bcftools view -i 'GT=="1/1"' > $strain_vcf
	vcftools --vcf $strain_vcf --SNPdensity 1000000 --out $DIR_DATA_CIR/$strain\_SNP_density_NO_HET
	SNP_IN=$DIR_DATA_CIR/$strain\_SNP_density_NO_HET.snpden
	SNP_OUT=$DIR_DATA_CIR/$strain\_snp_density_no_het.txt
	echo "snp file: $SNP_OUT"
	awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+1e6,$3}}' $SNP_IN > $SNP_OUT
done

#====================== INDELs ======================
VCF_INDEL_NO_HET=$DIR_DATA/speedseqsnp30x_merged_realigned_region_indel_no_Het.vcf
vcftools --vcf $VCF_INDEL_NO_HET --SNPdensity 1000000 --out $DIR_DATA_CIR/INDEL_density_NO_HET

INDEL_IN=$DIR_DATA_CIR/INDEL_density_NO_HET.snpden
INDEL_OUT=$DIR_DATA_CIR/indel_density_no_het.txt
echo "INDEL file: $INDEL_OUT"
awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+1e6,$3}}' $INDEL_IN > $INDEL_OUT

# strain specific INDELs
STRAINS="R2 R52 R72 R84 R104 R114 R124 R134 R144 R22 R32 R42 R62 R94"
for strain in $STRAINS;
do
	sample=$strain\_merged_realigned.bam
	echo $sample
	strain_vcf=$DIR_DATA_CIR/$strain\_INDEL_no_Het.vcf
	bcftools view -s $sample $VCF_INDEL_NO_HET | bcftools view -i 'GT=="1/1"' > $strain_vcf
	vcftools --vcf $strain_vcf --SNPdensity 1000000 --out $DIR_DATA_CIR/$strain\_INDEL_density_NO_HET
	INDEL_IN=$DIR_DATA_CIR/$strain\_INDEL_density_NO_HET.snpden
	INDEL_OUT=$DIR_DATA_CIR/$strain\_indel_density_no_het.txt
	echo "INDEL file: $INDEL_OUT"
	awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+1e6,$3}}' $INDEL_IN > $INDEL_OUT
done

#====================== STRs ======================
VCF_STR=$DIR_DATA/ALL_chrom_hipstr.vcf
vcftools --vcf $VCF_STR --SNPdensity 1000000 --out $DIR_DATA_CIR/STR_density
STR_IN=$DIR_DATA_CIR/STR_density.snpden
STR_OUT=$DIR_DATA_CIR/str_density.txt
echo "str file: $STR_OUT"
awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+1e6,$3}}' $STR_IN > $STR_OUT

# strain specific STRs
# NOTE: there is a typo in the STR file: instead of the correct R84 it has R82
STRAINS="R2 R52 R72 R82 R104 R114 R124 R134 R144 R22 R32 R42 R62 R94"
VCF_IN=$DIR_DATA/ALL_chrom_hipstr.vcf
for strain in $STRAINS;
do
	sample=$strain
	echo $sample
	strain_vcf=$DIR_DATA_CIR/STR_$strain.vcf
	bcftools view -s $sample $VCF_IN | bcftools view -e 'GT=="0|0"' > $strain_vcf
	vcftools --vcf $strain_vcf --SNPdensity 1000000 --out $DIR_DATA_CIR/STR_$strain

	STR_IN=$DIR_DATA_CIR/STR_$strain.snpden
	STR_OUT=$DIR_DATA_CIR/STR_$strain.txt
	echo "STR file: $STR_OUT"
	awk 'BEGIN{OFS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+1e6,$3}}' $STR_IN > $STR_OUT
done

#====================== SVs ======================
VCF_SV_IN=$DIR_DATA/merged_cnv_lumpy_SVs_sorted.vcf.gz
VCF_SV_OUT=$DIR_DATA_CIR/merged_cnv_lumpy_SVs_sorted_NO_BND.vcf
bcftools view -e 'SVTYPE=="BND"' $VCF_SV_IN > $VCF_SV_OUT
vcftools --gzvcf $VCF_SV_OUT --SNPdensity 1000000 --out $DIR_DATA_CIR/SV_density

SV_IN=$DIR_DATA_CIR/SV_density.snpden
SV_OUT=$DIR_DATA_CIR/sv_density.txt
echo "sv file: $SV_OUT"
awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+1e6,$3}}' $SV_IN > $SV_OUT

# strain specific SVs
STRAINS="R2 R52 R72 R84 R104 R114 R124 R134 R144 R22 R32 R42 R62 R94"
VCF_IN=$DIR_DATA_CIR/merged_cnv_lumpy_SVs_sorted_NO_BND.vcf
for strain in $STRAINS;
do
	sample=$strain
	echo $sample
	strain_vcf=$DIR_DATA_CIR/SV_$strain.vcf
	bcftools view -s $sample $VCF_IN | bcftools view -i 'GT=="1/1"' > $strain_vcf
	vcftools --vcf $strain_vcf --SNPdensity 1000000 --out $DIR_DATA_CIR/SV_$strain

	SV_IN=$DIR_DATA_CIR/SV_$strain.snpden
	SV_OUT=$DIR_DATA_CIR/SV_$strain.txt
	echo "SV file: $SV_OUT"
	awk 'BEGIN{OFS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+1e6,$3}}' $SV_IN > $SV_OUT
done

#====================== DEGs ======================
ANOVA_FILE=$DIR_DATA/anova_b6b10_normalized_CPM1.txt
GENES_DEG_VCF=$DIR_DATA_CIR/Genes_DEG.vcf

cat vcf.header > $GENES_DEG_VCF
awk 'BEGIN{OFS="\t"}{if($1!="Gene" && $11<0.05){print "chr"$3,$4,".","A","T","999","PASS","PV="$11,"GT","1/1"}}' $ANOVA_FILE | sort -k1,1V -k2,2n >> $GENES_DEG_VCF
vcftools --vcf $GENES_DEG_VCF --SNPdensity 1000000 --out $DIR_DATA_CIR/GENES_DEG_density
GENES_DEG_IN=$DIR_DATA_CIR/GENES_DEG_density.snpden
GENES_DEG_OUT=$DIR_DATA_CIR/Gene_DEG_density.txt
echo "Genes DEG file: $GENES_DEG_OUT"
awk 'BEGIN{OFS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+1e6,$3}}' $GENES_DEG_IN > $GENES_DEG_OUT

