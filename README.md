# Mortazavi2021_B6B10
**Polymorphic SNPs, short tandem repeats and structural variants cause differential gene expression among inbred C57BL/6 and C57BL/10 substrains**
***

1. Clone the repository
2. `cd Mortazavi2021_B6B10`
3. `mkdir data`
4. Download the data [here](https://drive.google.com/drive/folders/1g6WIabQRq3H0IpUBDZSswbRIbpRDYjY6?usp=sharing) into the folder `data`

## Circos plots
Programs needed:
```
bcftools
vcftools
circos
```
How to run:
1. `cd circos_plots`
2. `bash cmd_prepare_data.sh`
3. `bash cmd_run_circos.sh`
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
How to run:
1. `cd associations`
2. `python run_tests.py`
3. `bash cmd_get_significant.sh`
***

## RNA data analysis
```
python modules needed:
pandas

R packages needed:
edgeR
```
How to run:
1. `cd rna_analysis`
2. `bash cmd_run_analysis.sh`

## Dendrogram
```
programs needed:
plink 1.9
```
How to prepare files:
1. `cd dendrogram`
2. `bash cmd_prepare_files.sh`

To generate dendrogram plot refer to `dendrogram.ipynb` in `Colab_Notebooks`
