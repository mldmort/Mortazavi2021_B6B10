# Mortazavi2021_B6B10
**Polymorphic SNPs, short tandem repeats and structural variants cause differential gene expression among inbred C57BL/6 and C57BL/10 substrains**
***

1. Clone the repository
2. `cd Mortazavi2021_B6B10`
3. `mkdir data`
4. Download the data needed from [here](https://data.mendeley.com/drafts/39sw8xcrmv) into the folder `data`

## Circos plots
```
Programs needed:
bcftools
vcftools
circos
```
Files needed:
```
speedseqsnp30x_merged_realigned_region_snp_no_Het.vcf
speedseqsnp30x_merged_realigned_region_indel_no_Het.vcf
ALL_chrom_hipstr.vcf
merged_cnv_lumpy_SVs_sorted.vcf.gz
anova_b6b10_normalized_CPM1.txt
```
How to run:
1. Download the files needed above from [here](https://data.mendeley.com/drafts/39sw8xcrmv) to a directory named `Mortazavi2021_B6B10/data`
2. `cd circos_plots`
3. `bash cmd_prepare_data.sh`
4. `bash cmd_run_circos.sh`
***

## Association tests
```
python modules needed:
pandas
numpy
matplotlib
statsmodels
limix
```
Files needed:
```
K_chr*.csv
anova_b6b10_normalized_CPM1.txt
rpkm_b6b10_normalized_CPM1.txt
gene_cnv_pair.bed
SV_feature_pair.bed
STR_feature_pair.bed
vep_lof.bed
vep_mis.bed
NoInters_SVs_closestGene.bed
NoInters_STRs_closestGene.bed
```
How to run:
1. Download the files needed above from [here](https://data.mendeley.com/drafts/39sw8xcrmv) to a directory named `Mortazavi2021_B6B10/data`
2. `cd associations`
3. `python run_tests.py`
4. `bash cmd_get_significant.sh`
***

## RNA data analysis
```
python modules needed:
pandas

R packages needed:
edgeR
```
Files needed:
```
b6b10_rnaseq_raw.txt
b6b10_rnaseq_ids.txt
Mus_musculus.GRCm38.84.all_exons.merged.sizes.txt
ensembl_symbol.txt
```
How to run:
1. Download the files needed above from [here](https://data.mendeley.com/drafts/39sw8xcrmv) to a directory named `Mortazavi2021_B6B10/data`
2. `cd rna_analysis`
3. `bash cmd_run_analysis.sh`

## Dendrogram
```
programs needed:
plink 1.9
```
Files needed:
```
speedseqsnp30x_merged_realigned_region_snp_no_Het.vcf
ALL_chrom_hipstr.vcf
merged_cnv_lumpy_SVs_sorted.vcf.gz
```
How to prepare files:
1. Download the files needed above from [here](https://data.mendeley.com/drafts/39sw8xcrmv) to a directory named `Mortazavi2021_B6B10/data`
2. `cd dendrogram`
3. `bash cmd_prepare_files.sh`

To generate dendrogram plot refer to `dendrogram.ipynb` in `Colab_Notebooks`
