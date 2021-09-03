#!/bin/bash

### to run this script you need plink 1.9

###===================== SNPs =========================
VCF_MAIN=../data/speedseqsnp30x_merged_realigned_region_snp_no_Het.vcf

plink --vcf $VCF_MAIN --const-fid --allow-extra-chr --make-bed --out snp_plink_out

cat snp_plink_out.fam | awk '{gsub("_merged_realigned.bam", "", $2); print $0}' > snp_plink_out.fam2
mv snp_plink_out.fam2 snp_plink_out.fam

awk 'BEGIN{OFS="\t"}{print $1, "snp_"$1"_"$4, $3, $4, $5, $6}' snp_plink_out.bim > snp_plink_out.bim2
mv snp_plink_out.bim2 snp_plink_out.bim

plink --bfile snp_plink_out --indep-pairwise 50 5 0.5 --out snp_plink_out_LD

plink --bfile snp_plink_out --extract snp_plink_out_LD.prune.in --make-bed --out snp_plink_out_LD

###===================== STRs =========================
VCF_MAIN=../data/ALL_chrom_hipstr.vcf

plink --vcf $VCF_MAIN --const-fid --allow-extra-chr --make-bed --out str_plink_out

# correct the fam file R82->R84
awk 'BEGIN{OFS="\t"}{gsub("R82", "R84"); print $0}' str_plink_out.fam > str_plink_out.fam2
mv str_plink_out.fam2 str_plink_out.fam

###===================== SVs =========================
VCF_ORIG=../data/merged_cnv_lumpy_SVs_sorted.vcf.gz
VCF_MOD=merged_cnv_lumpy_SVs_MOD.vcf

bcftools view $VCF_ORIG | bcftools view -i 'LUMPY_CALL==1' | bcftools view -e 'SVTYPE=="BND"' | bcftools view -g ^miss > $VCF_MOD
bcftools view $VCF_ORIG | bcftools view -i 'CNV_CALL==1' | bcftools +missing2ref | bcftools view -H >> $VCF_MOD

VCF_MOD_SORTED=merged_cnv_lumpy_SVs_MOD_sorted.vcf
bcftools sort $VCF_MOD > $VCF_MOD_SORTED

plink --vcf $VCF_MOD_SORTED --const-fid --allow-extra-chr --make-bed --out sv_plink_out
awk 'BEGIN{OFS="\t"}{if($2=="."){print $1, "cnv_"$1"_"$4, $3, $4, $5, $6}else{print $1, "lumpy_"$1"_"$4, $3, $4, $5, $6}}' sv_plink_out.bim > sv_plink_out.bim2
mv sv_plink_out.bim2 sv_plink_out.bim

###===================== merge =========================
echo -e "snp_plink_out_LD\nstr_plink_out\nsv_plink_out" > files_merge.txt
plink --merge-list files_merge.txt --make-bed --out snp_str_sv_plink_out

###===================== IBS matrix =========================
plink --bfile snp_str_sv_plink_out --cluster --matrix --out snp_str_sv_matrix

