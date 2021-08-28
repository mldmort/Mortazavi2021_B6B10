import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from limix.qtl import scan

Data_Dir = '../data/'

def run_test_sv():

	print('++++++++++++++++++++ run SV tests ++++++++++++++++++')
	### general purpose ###
	gt_order = ['R2', 'R52', 'R72', 'R84', 'R104', 'R114', 'R124', 'R134', 'R144', 'R22', 'R32', 'R42', 'R62', 'R94']
	strain_dict = {'B6N-TyrC_BrdCrlCrl':'R134', 'C57BL_6ByJ':'R72', 'C57BL_6J':'R2', 'C57BL_6JBomTac':'R144', 'C57BL_6JEiJ':'R84', 'C57BL_6NCrl':'R114', 'C57BL_6NHsd':'R124', 'C57BL_6NJ':'R52', 'C57BL_6NTac':'R104', 'C57BL_10ScNHsd':'R94', 'C57BL_10ScSnJ':'R22', 'C57BL_10SnJ':'R32', 'C57BL_10ScCr':'R42', 'C57BL_10J':'R62'}

	#############################################################################
	sv_data_file = Data_Dir+'SV_feature_pair.bed'
	sv_data = pd.read_table(sv_data_file, header=0, index_col=0)
	print('sv data:')
	print(sv_data)

	sv_data['genotype_b6'] = sv_data['genotype'].apply(lambda x: ';'.join(x.split(';')[:9]))
	sv_data['genotype_b10'] = sv_data['genotype'].apply(lambda x: ';'.join(x.split(';')[9:14]))
	#print('sv data:.loc[:, ["genotype", "genotype_b6", "genotype_b10"]]')
	#print(sv_data.loc[:, ['genotype', 'genotype_b6', 'genotype_b10']])

	sv_iGene_data_file = Data_Dir+'NoInters_SVs_closestGene.bed'
	sv_iGene_data = pd.read_table(sv_iGene_data_file, header=None, index_col=None)
	sv_iGene_data['temp'] = sv_iGene_data.iloc[:,5:19].agg(';'.join, axis=1)
	def f_cnv_lum(s):
		return ';'.join(['1' if (ss=='1/1') else '0' for ss in s.split(';')])
	sv_iGene_data['genotype'] = sv_iGene_data['temp'].apply(lambda s : f_cnv_lum(s))
	sv_iGene_data.drop([5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,'temp'], axis=1, inplace=True)
	sv_iGene_data.columns = ['chrom', 'start', 'end', 'SVTYPE', 'caller', 'feature', 'Gene_Sym', 'Gene_ID', 'genotype']
	sv_iGene_data = sv_iGene_data[['chrom', 'start', 'end', 'SVTYPE', 'caller', 'genotype', 'Gene_ID', 'Gene_Sym', 'feature']]
	#print('sv iGene data:')
	#print(sv_iGene_data)

	expr_file = Data_Dir+'rpkm_b6b10_normalized_CPM1.txt'
	expr_data = pd.read_table(expr_file, header=0, index_col=None)
	#print('expr data:')
	#print(expr_data)

	anova_file = Data_Dir+'anova_b6b10_normalized_CPM1.txt'
	anova_data = pd.read_table(anova_file, header=0, index_col=None)
	#print('anova data:')
	#print(anova_data)
	anova_data_b6 = anova_data.loc[anova_data.anova_fdr_b6<0.05]
	anova_data_b10 = anova_data.loc[anova_data.anova_fdr_b10<0.05]
	anova_data_b6b10 = anova_data.loc[anova_data.anova_fdr_b6b10<0.05]
	#print('anova data_b6:')
	#print(anova_data_b6)
	#print('anova data_b10:')
	#print(anova_data_b10)
	#print('anova data_b6b10:')
	#print(anova_data_b6b10)

	#############################################################################
	### read LOCO kinship matrix
	K_dict = {}
	for i in list(range(1,20,1))+[23,24]:
		k_file = Data_Dir + 'K_chr'+str(i)+'.csv'
		name = 'chr'+str(i)
		if i==23:
			name = 'chrX'
		if i==24:
			name = 'chrY'
		K_dict[name] = np.loadtxt(k_file)
	#############################################################################

	###################################### SV in DEGenes #######################################
	sv_data_inDEG_b6b10 = sv_data.loc[sv_data.Gene_ID.isin(anova_data_b6b10.Gene.tolist())]
	print('SVs in DEGenes:')
	print(sv_data_inDEG_b6b10)

	print('test SVs in DEGenes...')
	temp = pd.DataFrame()
	for index in sv_data_inDEG_b6b10.index.tolist():
		#print(str(index)+'/'+str(sv_data_inDEG_b6b10.shape[0]))
		gene_id = sv_data_inDEG_b6b10.loc[index, 'Gene_ID']
		gene_sym = sv_data_inDEG_b6b10.loc[index, 'Gene_Sym']
		chrom = sv_data_inDEG_b6b10.loc[index, 'chrom']
		gt_dict = {gt_order[ix]:float(x) for ix, x in enumerate(str(sv_data_inDEG_b6b10.loc[index, 'genotype']).split(';'))}
		sv_expr_slice = pd.DataFrame(expr_data[['Strain', 'Group', gene_id]])
		sv_expr_slice['genotype'] = [gt_dict[strain_dict[x]] for x in sv_expr_slice['Strain'].tolist()]
		sv_expr_slice['Group_bin'] = sv_expr_slice['Group'].apply(lambda x: 0. if x=='B6' else 1.)
		#sv_expr_slice['PC1'] = [pca_data.loc[x, 'PC1'] for x in sv_expr_slice['Strain'].tolist()]
		#sv_expr_slice['PC2'] = [pca_data.loc[x, 'PC2'] for x in sv_expr_slice['Strain'].tolist()]
		#print('sv_expr_slice:')
		#print(sv_expr_slice)

		y = pd.DataFrame(sv_expr_slice[gene_id])
		X = pd.DataFrame(sv_expr_slice['genotype'])

		### without Group_bin as covariate
		M = pd.DataFrame(np.ones((X.shape)), columns=['offset'])
		res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
		temp.loc[index, 'pval_noCV'] = float(res.stats.pv20)

		### with Group_bin as covariate
		M = pd.DataFrame(sv_expr_slice['Group_bin']).assign(offset=1.)
		res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
		temp.loc[index, 'pval_B6B10CV'] = float(res.stats.pv20)
		#print('res:')
		#print(res)
		#print('res.stats:')
		#print(res.stats)
		#print('res.stats.pv20:')
		#print(float(res.stats.pv20))
	sv_data_inDEG_b6b10 = pd.merge(sv_data_inDEG_b6b10, temp, left_index=True, right_index=True)
	#print('sv_data_inDEG_b6b10:')
	#print(sv_data_inDEG_b6b10)

	file_out = 'pval_sv_inDEG_b6b10.txt'
	print('file out: ', file_out)
	sv_data_inDEG_b6b10.to_csv(file_out, sep='\t', header=True, index=False)
	######################################################################################################

	###################################### SV iGene in DEGenes #######################################
	sv_iGene_data_inDEG_b6b10 = sv_iGene_data.loc[sv_iGene_data.Gene_ID.isin(anova_data_b6b10.Gene.tolist())]
	print('intergenic SVs near DEGenes:')
	print(sv_iGene_data_inDEG_b6b10)

	print('test intergenic SVs near DEGenes...')
	temp = pd.DataFrame()
	for index in sv_iGene_data_inDEG_b6b10.index.tolist():
		#print(str(index)+'/'+str(sv_iGene_data_inDEG_b6b10.shape[0]))
		gene_id = sv_iGene_data_inDEG_b6b10.loc[index, 'Gene_ID']
		gene_sym = sv_iGene_data_inDEG_b6b10.loc[index, 'Gene_Sym']
		chrom = sv_iGene_data_inDEG_b6b10.loc[index, 'chrom']
		gt_dict = {gt_order[ix]:float(x) for ix, x in enumerate(str(sv_iGene_data_inDEG_b6b10.loc[index, 'genotype']).split(';'))}
		sv_expr_slice = pd.DataFrame(expr_data[['Strain', 'Group', gene_id]])
		sv_expr_slice['genotype'] = [gt_dict[strain_dict[x]] for x in sv_expr_slice['Strain'].tolist()]
		sv_expr_slice['Group_bin'] = sv_expr_slice['Group'].apply(lambda x: 0. if x=='B6' else 1.)
		#sv_expr_slice['PC1'] = [pca_data.loc[x, 'PC1'] for x in sv_expr_slice['Strain'].tolist()]
		#sv_expr_slice['PC2'] = [pca_data.loc[x, 'PC2'] for x in sv_expr_slice['Strain'].tolist()]
		#print('sv_expr_slice:')
		#print(sv_expr_slice)

		y = pd.DataFrame(sv_expr_slice[gene_id])
		X = pd.DataFrame(sv_expr_slice['genotype'])

		### without Group_bin as covariate
		M = pd.DataFrame(np.ones((X.shape)), columns=['offset'])
		res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
		temp.loc[index, 'pval_noCV'] = float(res.stats.pv20)

		### with Group_bin as covariate
		M = pd.DataFrame(sv_expr_slice['Group_bin']).assign(offset=1.)
		res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
		temp.loc[index, 'pval_B6B10CV'] = float(res.stats.pv20)
		#print('res:')
		#print(res)
		#print('res.stats:')
		#print(res.stats)
		#print('res.stats.pv20:')
		#print(float(res.stats.pv20))
	sv_iGene_data_inDEG_b6b10 = pd.merge(sv_iGene_data_inDEG_b6b10, temp, left_index=True, right_index=True)
	#print('sv_iGene_data_inDEG_b6b10:')
	#print(sv_iGene_data_inDEG_b6b10)

	file_out = 'pval_sv_iGene_inDEG_b6b10.txt'
	print('file out: ', file_out)
	sv_iGene_data_inDEG_b6b10.to_csv(file_out, sep='\t', header=True, index=False)
	######################################################################################################

def run_test_cnv():

	print('++++++++++++++++++++ run CNV tests ++++++++++++++++++')
	### general purpose ###
	gt_order = ['R2', 'R52', 'R72', 'R84', 'R104', 'R114', 'R124', 'R134', 'R144', 'R22', 'R32', 'R42', 'R62', 'R94']
	strain_dict = {'B6N-TyrC_BrdCrlCrl':'R134', 'C57BL_6ByJ':'R72', 'C57BL_6J':'R2', 'C57BL_6JBomTac':'R144', 'C57BL_6JEiJ':'R84', 'C57BL_6NCrl':'R114', 'C57BL_6NHsd':'R124', 'C57BL_6NJ':'R52', 'C57BL_6NTac':'R104', 'C57BL_10ScNHsd':'R94', 'C57BL_10ScSnJ':'R22', 'C57BL_10SnJ':'R32', 'C57BL_10ScCr':'R42', 'C57BL_10J':'R62'}

	cnv_cols = ['cnv_R2', 'cnv_R52', 'cnv_R72', 'cnv_R84', 'cnv_R104', 'cnv_R114', 'cnv_R124', 'cnv_R134', 'cnv_R144', 'cnv_R22', 'cnv_R32', 'cnv_R42', 'cnv_R62', 'cnv_R94']
	std_cols = ['std_R2', 'std_R52', 'std_R72', 'std_R84', 'std_R104', 'std_R114', 'std_R124', 'std_R134', 'std_R144', 'std_R22', 'std_R32', 'std_R42', 'std_R62', 'std_R94']

	#############################################################################
	cnv_file = Data_Dir+'gene_cnv_pair.bed'
	cnv_data = pd.read_table(cnv_file, sep='\t', header=0, index_col=None)
	print('cnv_data:')
	print(cnv_data)

	expr_file = Data_Dir+'rpkm_b6b10_normalized_CPM1.txt'
	expr_data = pd.read_table(expr_file, header=0, index_col=None)
	#print('expr data:')
	#print(expr_data)
	
	anova_file = Data_Dir+'anova_b6b10_normalized_CPM1.txt'
	anova_data = pd.read_table(anova_file, header=0, index_col=None)
	#print('anova data:')
	#print(anova_data)
	anova_data_b6 = anova_data.loc[anova_data.anova_fdr_b6<0.05]
	anova_data_b10 = anova_data.loc[anova_data.anova_fdr_b10<0.05]
	anova_data_b6b10 = anova_data.loc[anova_data.anova_fdr_b6b10<0.05]
	#print('anova data_b6:')
	#print(anova_data_b6)
	#print('anova data_b10:')
	#print(anova_data_b10)
	#print('anova data_b6b10:')
	#print(anova_data_b6b10)

	#############################################################################
	### read LOCO kinship matrix
	K_dict = {}
	for i in list(range(1,20,1))+[23,24]:
		k_file = Data_Dir + 'K_chr'+str(i)+'.csv'
		name = 'chr'+str(i)
		if i==23:
			name = 'chrX'
		if i==24:
			name = 'chrY'
		K_dict[name] = np.loadtxt(k_file)
	#############################################################################

	###################################### CNV in DEGenes #######################################
	cnv_data_inDEG_b6b10 = cnv_data.loc[cnv_data.Gene_ID.isin(anova_data_b6b10.Gene.tolist())]
	print('CNVs in DEGenes:')
	print(cnv_data_inDEG_b6b10)

	print('test CNV in DEGenes...')
	temp = pd.DataFrame()
	for index in cnv_data_inDEG_b6b10.index.tolist():
		#print(str(index)+'/'+str(cnv_data_inDEG_b6b10.shape[0]))
		gene_id = cnv_data_inDEG_b6b10.loc[index, 'Gene_ID']
		gene_sym = cnv_data_inDEG_b6b10.loc[index, 'Gene_Sym']
		chrom = cnv_data_inDEG_b6b10.loc[index, 'chrom']
		gt_dict = {gt_order[ix]:float(x) for ix, x in enumerate(cnv_data_inDEG_b6b10.loc[index, cnv_cols].tolist())}
		std_dict = {gt_order[ix]:float(x) for ix, x in enumerate(cnv_data_inDEG_b6b10.loc[index, std_cols].tolist())}
		cnv_expr_slice = pd.DataFrame(expr_data[['Strain', 'Group', gene_id]])
		cnv_expr_slice['genotype'] = [gt_dict[strain_dict[x]] for x in cnv_expr_slice['Strain'].tolist()]
		cnv_expr_slice['gt_std'] = [std_dict[strain_dict[x]] for x in cnv_expr_slice['Strain'].tolist()]
		cnv_expr_slice['Group_bin'] = cnv_expr_slice['Group'].apply(lambda x: 0. if x=='B6' else 1.)
		#cnv_expr_slice['PC1'] = [pca_data.loc[x, 'PC1'] for x in cnv_expr_slice['Strain'].tolist()]
		#cnv_expr_slice['PC2'] = [pca_data.loc[x, 'PC2'] for x in cnv_expr_slice['Strain'].tolist()]
		#print('cnv_expr_slice:')
		#print(cnv_expr_slice)

		y = pd.DataFrame(cnv_expr_slice[gene_id])
		X = pd.DataFrame(cnv_expr_slice['genotype'])

		### without Group_bin as covariate
		M = pd.DataFrame(np.ones((X.shape)), columns=['offset'])
		res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
		temp.loc[index, 'pval_noCV'] = float(res.stats.pv20)

		### with Group_bin as covariate
		M = pd.DataFrame(cnv_expr_slice['Group_bin']).assign(offset=1.)
		res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
		temp.loc[index, 'pval_B6B10CV'] = float(res.stats.pv20)
		#print('res:')
		#print(res)
		#print('res.stats:')
		#print(res.stats)
		#print('res.stats.pv20:')
		#print(float(res.stats.pv20))
	cnv_data_inDEG_b6b10 = pd.merge(cnv_data_inDEG_b6b10, temp, left_index=True, right_index=True)
	#print('cnv_data_inDEG_b6b10:')
	#print(cnv_data_inDEG_b6b10)

	file_out = 'pval_cnv_inDEG_b6b10.txt'
	print('file out: ', file_out)
	cnv_data_inDEG_b6b10.to_csv(file_out, sep='\t', header=True, index=False)
	######################################################################################################

def run_test_snp():

	print('++++++++++++++++++++ run SNP tests ++++++++++++++++++')
	### general purpose ###
	gt_order = ['R2', 'R52', 'R72', 'R84', 'R104', 'R114', 'R124', 'R134', 'R144', 'R22', 'R32', 'R42', 'R62', 'R94']
	strain_dict = {'B6N-TyrC_BrdCrlCrl':'R134', 'C57BL_6ByJ':'R72', 'C57BL_6J':'R2', 'C57BL_6JBomTac':'R144', 'C57BL_6JEiJ':'R84', 'C57BL_6NCrl':'R114', 'C57BL_6NHsd':'R124', 'C57BL_6NJ':'R52', 'C57BL_6NTac':'R104', 'C57BL_10ScNHsd':'R94', 'C57BL_10ScSnJ':'R22', 'C57BL_10SnJ':'R32', 'C57BL_10ScCr':'R42', 'C57BL_10J':'R62'}

	#############################################################################
	expr_file = Data_Dir+'rpkm_b6b10_normalized_CPM1.txt'
	expr_data = pd.read_table(expr_file, header=0, index_col=None)
	#print('expr data:')
	#print(expr_data)
	
	anova_file = Data_Dir+'anova_b6b10_normalized_CPM1.txt'
	anova_data = pd.read_table(anova_file, header=0, index_col=None)
	#print('anova data:')
	#print(anova_data)
	anova_data_b6 = anova_data.loc[anova_data.anova_fdr_b6<0.05]
	anova_data_b10 = anova_data.loc[anova_data.anova_fdr_b10<0.05]
	anova_data_b6b10 = anova_data.loc[anova_data.anova_fdr_b6b10<0.05]
	#print('anova data_b6:')
	#print(anova_data_b6)
	#print('anova data_b10:')
	#print(anova_data_b10)
	#print('anova data_b6b10:')
	#print(anova_data_b6b10)

	#############################################################################
	### read LOCO kinship matrix
	K_dict = {}
	for i in list(range(1,20,1))+[23,24]:
		k_file = Data_Dir + 'K_chr'+str(i)+'.csv'
		name = 'chr'+str(i)
		if i==23:
			name = 'chrX'
		if i==24:
			name = 'chrY'
		K_dict[name] = np.loadtxt(k_file)
	#############################################################################

	###################################### SNP-LoF in DEGenes #######################################
	vep_file = Data_Dir+'vep_lof.bed'
	vep_data = pd.read_table(vep_file, sep='\t', header=None, index_col=None)
	vep_data.columns = ['chrom', 'start', 'end', 'ref', 'alt', 'consequence', 'lof', 'mis', 'syn', 'non', 'Gene_Sym', 'Gene_ID', 'genotype']
	print('LOF SNPs:')
	print(vep_data)

	## order of genotypes in the vep vcf file:
	##order = ['R2', 'R22', 'R32', 'R42', 'R52', 'R62', 'R72', 'R84', 'R94', 'R104', 'R114', 'R124', 'R134', 'R144']
	#########   0		1	  2		3		4		5	6		7		8	9		10		11		12		13
	genotype_b6_ind = [0, 4, 6, 7, 9, 10, 11, 12, 13]
	genotype_b10_ind = [1, 2, 3, 5, 8]

	temp = pd.DataFrame()
	gt_set = set()
	for index in vep_data.index.tolist():
		gt_list = str(vep_data.loc[index, 'genotype']).split(',')
		if {'./.', '0/2', '1/2', '2/2', '0/3', '1/3', '2/3', '3/3'}.isdisjoint( set(gt_list) ):
			temp = pd.concat([temp, pd.DataFrame([vep_data.loc[index,]])], axis=0, ignore_index=True)
			for x in gt_list: gt_set.add(x)
		#if '0/1' in gt_list:
		#	any_het_df = pd.concat([any_het_df, pd.DataFrame([vep_data.loc[index,['chrom', 'start', 'end', 'lof', 'mis', 'syn', 'non', 'Gene_Sym', 'Gene_ID', 'genotype']]])], axis=0, ignore_index=True)
	vep_data = temp
	#print('vep_data: only 0/0, 0/1, 1/1 genotypes included')
	#print(vep_data)
	#print('gt_set:')
	#print(gt_set)

	def f_gt(gt_str):
		gts = gt_str.split(',')
		gt_list = []
		for gt in gts:
			if gt == '1/1':
				gt_list.append('2')
			elif gt == '0/1':
				gt_list.append('1')
			else:
				gt_list.append('0')
		gt_ret = ','.join(gt_list)
		return gt_ret

	temp = pd.DataFrame(vep_data)
	temp['genotype'] = vep_data.genotype.apply(lambda x: f_gt(x))
	vep_data = temp
	
	def f_gt_b6(gt_str):
		gts = gt_str.split(',')
		ret = ','.join([gts[i] for i in genotype_b6_ind])
		return ret

	def f_gt_b10(gt_str):
		gts = gt_str.split(',')
		ret = ','.join([gts[i] for i in genotype_b10_ind])
		return ret

	vep_data['genotype_b6'] = vep_data.genotype.apply(lambda x: f_gt_b6(x))
	vep_data['genotype_b10'] = vep_data.genotype.apply(lambda x: f_gt_b10(x))
	#print('vep_data: 1/1->2,  0/1->1,  0/0->0')
	#print(vep_data)

	####### eliminate duplications in genes and unfold multiple genes
	df_add = pd.DataFrame()
	for index in vep_data.index.tolist():
		vep_data.loc[index, 'Gene_Sym'] = ','.join(list(set(str(vep_data.loc[index, 'Gene_Sym']).split(','))))
		vep_data.loc[index, 'Gene_ID'] = ','.join(list(set(str(vep_data.loc[index, 'Gene_ID']).split(','))))
		if len(str(vep_data.loc[index, 'Gene_Sym']).split(',')) > 1:
			genes_sym_id = zip(str(vep_data.loc[index, 'Gene_Sym']).split(','), str(vep_data.loc[index, 'Gene_ID']).split(','))
			for gene_sym_id in genes_sym_id:
				gene_sym = gene_sym_id[0]
				gene_id = gene_sym_id[1]
				row = pd.DataFrame([vep_data.loc[index,]])
				row.Gene_Sym = gene_sym
				row.Gene_ID = gene_id
				df_add = pd.concat([df_add, row], axis=0, ignore_index=True)
	vep_data = pd.concat([vep_data, df_add], axis=0, ignore_index=True)
	#print('vep_data: one entry for each (SNP, gene) pair')
	#print(vep_data)

	################### SNP: fit linear regression models ######################
	vep_data_inDEG_b6b10 = vep_data.loc[vep_data.Gene_ID.isin(anova_data_b6b10.Gene.tolist())]
	print('LOF SNPs in DEGenes:')
	print(vep_data_inDEG_b6b10)
	
	###################### LMM-GRM ####################
	## order of genotypes in the vep vcf file:
	##order = ['R2', 'R22', 'R32', 'R42', 'R52', 'R62', 'R72', 'R84', 'R94', 'R104', 'R114', 'R124', 'R134', 'R144']
	#########   0		1	  2		3		4		5	6		7		8	9		10		11		12		13
	print('test LOF SNPs in DEGenes...')
	temp = pd.DataFrame()
	for index in vep_data_inDEG_b6b10.index.tolist():
		#print(str(index)+'/'+str(vep_data_inDEG_b6b10.shape[0]))
		gene_id = str(vep_data_inDEG_b6b10.loc[index, 'Gene_ID'])
		gt_list_raw = str(vep_data_inDEG_b6b10.loc[index, 'genotype']).split(',')
		gt_list = [gt_list_raw[0], gt_list_raw[4], gt_list_raw[6], gt_list_raw[7], gt_list_raw[9], gt_list_raw[10], gt_list_raw[11], gt_list_raw[12], gt_list_raw[13], gt_list_raw[1], gt_list_raw[2], gt_list_raw[3], gt_list_raw[5], gt_list_raw[8]]
		gt_dict = {gt_order[ix]:float(x) for ix, x in enumerate(gt_list)}
		chrom = vep_data_inDEG_b6b10.loc[index, 'chrom']
		vep_expr_slice = pd.DataFrame(expr_data[['Strain', 'Group', gene_id]])
		vep_expr_slice['Group_bin'] = vep_expr_slice['Group'].apply(lambda x: 0. if x=='B6' else 1.)
		vep_expr_slice['genotype'] = [gt_dict[strain_dict[x]] for x in vep_expr_slice['Strain'].tolist()]
		#vep_expr_slice['PC1'] = [pca_data.loc[x, 'PC1'] for x in vep_expr_slice['Strain'].tolist()]
		#vep_expr_slice['PC2'] = [pca_data.loc[x, 'PC2'] for x in vep_expr_slice['Strain'].tolist()]
		#print('vep_expr_slice:')
		#print(vep_expr_slice)
		
		if (len(set(vep_expr_slice.genotype))>1):
			y = pd.DataFrame(vep_expr_slice[gene_id])
			X = pd.DataFrame(vep_expr_slice['genotype'])

			### without Group_bin as covariate
			M = pd.DataFrame(np.ones((X.shape)), columns=['offset'])
			res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
			temp.loc[index, 'pval_noCV'] = float(res.stats.pv20)

			### with Group_bin as covariate
			M = pd.DataFrame(vep_expr_slice['Group_bin']).assign(offset=1.)
			res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
			temp.loc[index, 'pval_B6B10CV'] = float(res.stats.pv20)
			#print('res:')
			#print(res)
			#print('res.stats:')
			#print(res.stats)
			#print('res.stats.pv20:')
			#print(float(res.stats.pv20))

	vep_data_inDEG_b6b10 = pd.merge(vep_data_inDEG_b6b10, temp, left_index=True, right_index=True)
	#print('vep_data_inDEG_b6b10:')
	#print(vep_data_inDEG_b6b10)
	
	file_out = 'pval_vep_inDEG_b6b10.txt'
	print('file out: ', file_out)
	vep_data_inDEG_b6b10.to_csv(file_out, sep='\t', header=True, index=False)

	###################################### SNP-mis in DEGenes #######################################
	vep_file = Data_Dir+'vep_mis.bed'
	vep_data = pd.read_table(vep_file, sep='\t', header=None, index_col=None)
	vep_data.columns = ['chrom', 'start', 'end', 'ref', 'alt', 'consequence', 'lof', 'mis', 'syn', 'non', 'Gene_Sym', 'Gene_ID', 'genotype']
	print('MIS SNPs:')
	print(vep_data)

	## order of genotypes in the vep vcf file:
	##order = ['R2', 'R22', 'R32', 'R42', 'R52', 'R62', 'R72', 'R84', 'R94', 'R104', 'R114', 'R124', 'R134', 'R144']
	#########   0		1	  2		3		4		5	6		7		8	9		10		11		12		13

	temp = pd.DataFrame()
	gt_set = set()
	for index in vep_data.index.tolist():
		gt_list = str(vep_data.loc[index, 'genotype']).split(',')
		if {'./.', '0/2', '1/2', '2/2', '0/3', '1/3', '2/3', '3/3'}.isdisjoint( set(gt_list) ):
			temp = pd.concat([temp, pd.DataFrame([vep_data.loc[index,]])], axis=0, ignore_index=True)
			for x in gt_list: gt_set.add(x)
	vep_data = temp
	#print('vep_data: only 0/0, 0/1, 1/1 genotypes included')
	#print(vep_data)
	#print('gt_set:')
	#print(gt_set)

	def f_gt(gt_str):
		gts = gt_str.split(',')
		gt_list = []
		for gt in gts:
			if gt == '1/1':
				gt_list.append('2')
			elif gt == '0/1':
				gt_list.append('1')
			else:
				gt_list.append('0')
		gt_ret = ','.join(gt_list)
		return gt_ret

	temp = pd.DataFrame(vep_data)
	temp['genotype'] = vep_data.genotype.apply(lambda x: f_gt(x))
	vep_data = temp
	
	#print('vep_data: 1/1->2,  0/1->1,  0/0->0')
	#print(vep_data)

	####### eliminate duplications in genes and unfold multiple genes
	df_add = pd.DataFrame()
	for index in vep_data.index.tolist():
		vep_data.loc[index, 'Gene_Sym'] = ','.join(list(set(str(vep_data.loc[index, 'Gene_Sym']).split(','))))
		vep_data.loc[index, 'Gene_ID'] = ','.join(list(set(str(vep_data.loc[index, 'Gene_ID']).split(','))))
		if len(str(vep_data.loc[index, 'Gene_Sym']).split(',')) > 1:
			genes_sym_id = zip(str(vep_data.loc[index, 'Gene_Sym']).split(','), str(vep_data.loc[index, 'Gene_ID']).split(','))
			for gene_sym_id in genes_sym_id:
				gene_sym = gene_sym_id[0]
				gene_id = gene_sym_id[1]
				row = pd.DataFrame([vep_data.loc[index,]])
				row.Gene_Sym = gene_sym
				row.Gene_ID = gene_id
				df_add = pd.concat([df_add, row], axis=0, ignore_index=True)
	vep_data = pd.concat([vep_data, df_add], axis=0, ignore_index=True)
	#print('vep_data: one entry for each (SNP, gene) pair')
	#print(vep_data)

	################### SNP: fit linear regression models ######################
	vep_data_inDEG_b6b10 = vep_data.loc[vep_data.Gene_ID.isin(anova_data_b6b10.Gene.tolist())]
	print('MIS SNPs in DEGenes:')
	print(vep_data_inDEG_b6b10)
	
	###################### LMM-GRM ####################
	## order of genotypes in the vep vcf file:
	##order = ['R2', 'R22', 'R32', 'R42', 'R52', 'R62', 'R72', 'R84', 'R94', 'R104', 'R114', 'R124', 'R134', 'R144']
	#########   0		1	  2		3		4		5	6		7		8	9		10		11		12		13
	print('test MIS SNPs in DEGenes...')
	temp = pd.DataFrame()
	for index in vep_data_inDEG_b6b10.index.tolist():
		#print(str(index)+'/'+str(vep_data_inDEG_b6b10.shape[0]))
		gene_id = str(vep_data_inDEG_b6b10.loc[index, 'Gene_ID'])
		gt_list_raw = str(vep_data_inDEG_b6b10.loc[index, 'genotype']).split(',')
		gt_list = [gt_list_raw[0], gt_list_raw[4], gt_list_raw[6], gt_list_raw[7], gt_list_raw[9], gt_list_raw[10], gt_list_raw[11], gt_list_raw[12], gt_list_raw[13], gt_list_raw[1], gt_list_raw[2], gt_list_raw[3], gt_list_raw[5], gt_list_raw[8]]
		gt_dict = {gt_order[ix]:float(x) for ix, x in enumerate(gt_list)}
		chrom = vep_data_inDEG_b6b10.loc[index, 'chrom']
		vep_expr_slice = pd.DataFrame(expr_data[['Strain', 'Group', gene_id]])
		vep_expr_slice['Group_bin'] = vep_expr_slice['Group'].apply(lambda x: 0. if x=='B6' else 1.)
		vep_expr_slice['genotype'] = [gt_dict[strain_dict[x]] for x in vep_expr_slice['Strain'].tolist()]
		#vep_expr_slice['PC1'] = [pca_data.loc[x, 'PC1'] for x in vep_expr_slice['Strain'].tolist()]
		#vep_expr_slice['PC2'] = [pca_data.loc[x, 'PC2'] for x in vep_expr_slice['Strain'].tolist()]
		#print('vep_expr_slice:')
		#print(vep_expr_slice)
		
		if (len(set(vep_expr_slice.genotype))>1):
			y = pd.DataFrame(vep_expr_slice[gene_id])
			X = pd.DataFrame(vep_expr_slice['genotype'])

			### without Group_bin as covariate
			M = pd.DataFrame(np.ones((X.shape)), columns=['offset'])
			res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
			temp.loc[index, 'pval_noCV'] = float(res.stats.pv20)

			### with Group_bin as covariate
			M = pd.DataFrame(vep_expr_slice['Group_bin']).assign(offset=1.)
			res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
			temp.loc[index, 'pval_B6B10CV'] = float(res.stats.pv20)
			#print('res:')
			#print(res)
			#print('res.stats:')
			#print(res.stats)
			#print('res.stats.pv20:')
			#print(float(res.stats.pv20))

	vep_data_inDEG_b6b10 = pd.merge(vep_data_inDEG_b6b10, temp, left_index=True, right_index=True)
	#print('vep_data_inDEG_b6b10:')
	#print(vep_data_inDEG_b6b10)
	
	file_out = 'pval_vep_mis_inDEG_b6b10.txt'
	print('file out: ', file_out)
	vep_data_inDEG_b6b10.to_csv(file_out, sep='\t', header=True, index=False)
	
def run_test_str():

	print('++++++++++++++++++++ run STR tests ++++++++++++++++++')
	### general purpose ###
	gt_order = ['R2', 'R52', 'R72', 'R84', 'R104', 'R114', 'R124', 'R134', 'R144', 'R22', 'R32', 'R42', 'R62', 'R94']
	strain_dict = {'B6N-TyrC_BrdCrlCrl':'R134', 'C57BL_6ByJ':'R72', 'C57BL_6J':'R2', 'C57BL_6JBomTac':'R144', 'C57BL_6JEiJ':'R84', 'C57BL_6NCrl':'R114', 'C57BL_6NHsd':'R124', 'C57BL_6NJ':'R52', 'C57BL_6NTac':'R104', 'C57BL_10ScNHsd':'R94', 'C57BL_10ScSnJ':'R22', 'C57BL_10SnJ':'R32', 'C57BL_10ScCr':'R42', 'C57BL_10J':'R62'}

	#############################################################################
	expr_file = Data_Dir+'rpkm_b6b10_normalized_CPM1.txt'
	expr_data = pd.read_table(expr_file, header=0, index_col=None)
	#print('expr data:')
	#print(expr_data)
	
	anova_file = Data_Dir+'anova_b6b10_normalized_CPM1.txt'
	anova_data = pd.read_table(anova_file, header=0, index_col=None)
	#print('anova data:')
	#print(anova_data)
	anova_data_b6 = anova_data.loc[anova_data.anova_fdr_b6<0.05]
	anova_data_b10 = anova_data.loc[anova_data.anova_fdr_b10<0.05]
	anova_data_b6b10 = anova_data.loc[anova_data.anova_fdr_b6b10<0.05]
	#print('anova data_b6:')
	#print(anova_data_b6)
	#print('anova data_b10:')
	#print(anova_data_b10)
	#print('anova data_b6b10:')
	#print(anova_data_b6b10)

	#############################################################################
	### read LOCO kinship matrix
	K_dict = {}
	for i in list(range(1,20,1))+[23,24]:
		k_file = Data_Dir + 'K_chr'+str(i)+'.csv'
		name = 'chr'+str(i)
		if i==23:
			name = 'chrX'
		if i==24:
			name = 'chrY'
		K_dict[name] = np.loadtxt(k_file)
	#############################################################################

	################### Process STR data, large effect  ######################
	str_file = Data_Dir+'STR_feature_pair.bed'
	str_data = pd.read_table(str_file, sep='\t', header=0, index_col=0)
	str_data = str_data.loc[str_data.feature.isin(['exon', 'TSS', 'UTR5', 'promoter'])]
	print('str_data:')
	print(str_data)

	BPDevs = []
	for index in str_data.index.tolist():
		if str_data.loc[index, 'chrom'] in ['chrX', 'chrY']:
			bp = str_data.loc[index, 'GB'][:-1]
			BPDevs.append(bp)
		else:
			bp = ';'.join([ str(int(xx.split('|')[0])+int(xx.split('|')[1])) for xx in str_data.loc[index, 'GB'].split(';')[:-1] ])
			BPDevs.append(bp)
	str_data['BPDev'] = BPDevs
	#print('str_data:')
	#print(str_data)

	str_data_inDEG_b6b10 = str_data.loc[str_data.Gene_ID.isin(anova_data_b6b10.Gene.tolist())]
	#print('STRs in DEGenes:')
	#print(str_data_inDEG_b6b10)

	str_data_inDEG_b6b10 = str_data_inDEG_b6b10.loc[str_data_inDEG_b6b10.Period>1]
	print('STRs in DEGenes:')
	print(str_data_inDEG_b6b10)

	#**##order: R104    R114    R124    R134    R144    R2      R22     R32     R42     R52     R62     R72     R84     R94
	#**#####     0       1       2       3       4       5       6       7       8       9       10      11      12      13

	###################### LMM-GRM ####################
	## order of genotypes in the str vcf file (both HipSTR and GangSTR):
	##order: R104    R114    R124    R134    R144    R2      R22     R32     R42     R52     R62     R72     R84     R94
	#####     0       1       2       3       4       5       6       7       8       9       10      11      12      13
	print('test STRs in DEGenes...')
	temp = pd.DataFrame()
	for index in str_data_inDEG_b6b10.index.tolist():
		gene_id = str(str_data_inDEG_b6b10.loc[index, 'Gene_ID'])
		gt_list_raw = str(str_data_inDEG_b6b10.loc[index, 'BPDev']).split(';')
		gt_list = [gt_list_raw[5], gt_list_raw[9], gt_list_raw[11], gt_list_raw[12], gt_list_raw[0], gt_list_raw[1], gt_list_raw[2], gt_list_raw[3], gt_list_raw[4], gt_list_raw[6], gt_list_raw[7], gt_list_raw[8], gt_list_raw[10], gt_list_raw[13]]
		gt_dict = {gt_order[ix]:float(x) for ix, x in enumerate(gt_list)}
		str_expr_slice = pd.DataFrame(expr_data[['Strain', 'Group', gene_id]])
		str_expr_slice['Group_bin'] = str_expr_slice['Group'].apply(lambda x: 0. if x=='B6' else 1.)
		str_expr_slice['genotype'] = [gt_dict[strain_dict[x]] for x in str_expr_slice['Strain'].tolist()]
		chrom = str_data_inDEG_b6b10.loc[index, 'chrom']
		
		if (len(set(str_expr_slice.genotype))>1):
			y = pd.DataFrame(str_expr_slice[gene_id])
			X = pd.DataFrame(str_expr_slice['genotype'])
			M = pd.DataFrame(str_expr_slice['Group_bin']).assign(offset=1.)
			### without Group_bin as covariate
			#M = pd.DataFrame(np.ones((X.shape)))
			res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
			temp.loc[index, 'pval'] = float(res.stats.pv20)
			#print('M:')
			#print(M)
			#print('X:')
			#print(X['genotype'].tolist())
			#print('res:')
			#print(res)
			#print('res.stats:')
			#print(res.stats)
			#print('res.stats.pv20:')
			#print(float(res.stats.pv20))

	str_data_inDEG_b6b10 = pd.merge(str_data_inDEG_b6b10, temp, left_index=True, right_index=True)
	#print('str_data_inDEG_b6b10:')
	#print(str_data_inDEG_b6b10)

	file_out = 'pval_str_inDEG_b6b10_exon_TSS_UTR5_promoter_PeriodLT1.txt'
	print('file out: ', file_out)
	str_data_inDEG_b6b10.to_csv(file_out, sep='\t', header=True, index=False)

	################### Process STR data, intergenic  ######################
	str_file = Data_Dir+'NoInters_STRs_closestGene.bed'
	str_data = pd.read_table(str_file, sep='\t', header=None, index_col=None)
	str_data.columns = ['chrom', 'start', 'end', 'GT', 'GB', 'Q', 'REF', 'ALT', 'Period', 'ft_chr', 'ft_s', 'ft_e', 'feature', 'Gene_Sym', 'Gene_ID']
	#print('str_data:')
	#print(str_data)

	BPDevs = []
	for index in str_data.index.tolist():
		if str_data.loc[index, 'chrom'] in ['chrX', 'chrY']:
			bp = str_data.loc[index, 'GB'][:-1]
			BPDevs.append(bp)
		else:
			bp = ';'.join([ str(int(xx.split('|')[0])+int(xx.split('|')[1])) for xx in str_data.loc[index, 'GB'].split(';')[:-1] ])
			BPDevs.append(bp)
	str_data['BPDev'] = BPDevs
	#print('str_data:')
	#print(str_data)

	str_data_inDEG_b6b10 = str_data.loc[str_data.Gene_ID.isin(anova_data_b6b10.Gene.tolist())]
	print('intergenic STRs near DEGenes:')
	print(str_data_inDEG_b6b10)

	#**##order: R104    R114    R124    R134    R144    R2      R22     R32     R42     R52     R62     R72     R84     R94
	#**#####     0       1       2       3       4       5       6       7       8       9       10      11      12      13

	###################### LMM-GRM ####################
	## order of genotypes in the str vcf file (both HipSTR and GangSTR):
	##order: R104    R114    R124    R134    R144    R2      R22     R32     R42     R52     R62     R72     R84     R94
	#####     0       1       2       3       4       5       6       7       8       9       10      11      12      13
	print('test intergenic STRs near DEGenes...')
	temp = pd.DataFrame()
	#for index in str_data_inDEG_b6b10.loc[str_data_inDEG_b6b10.Gene_Sym == "Aim2"].index.tolist():
	for index in str_data_inDEG_b6b10.index.tolist():
		gene_id = str(str_data_inDEG_b6b10.loc[index, 'Gene_ID'])
		gt_list_raw = str(str_data_inDEG_b6b10.loc[index, 'BPDev']).split(';')
		gt_list = [gt_list_raw[5], gt_list_raw[9], gt_list_raw[11], gt_list_raw[12], gt_list_raw[0], gt_list_raw[1], gt_list_raw[2], gt_list_raw[3], gt_list_raw[4], gt_list_raw[6], gt_list_raw[7], gt_list_raw[8], gt_list_raw[10], gt_list_raw[13]]
		gt_dict = {gt_order[ix]:float(x) for ix, x in enumerate(gt_list)}
		str_expr_slice = pd.DataFrame(expr_data[['Strain', 'Group', gene_id]])
		str_expr_slice['Group_bin'] = str_expr_slice['Group'].apply(lambda x: 0. if x=='B6' else 1.)
		str_expr_slice['genotype'] = [gt_dict[strain_dict[x]] for x in str_expr_slice['Strain'].tolist()]
		chrom = str_data_inDEG_b6b10.loc[index, 'chrom']
		
		if (len(set(str_expr_slice.genotype))>1):
			y = pd.DataFrame(str_expr_slice[gene_id])
			X = pd.DataFrame(str_expr_slice['genotype'])
			M = pd.DataFrame(str_expr_slice['Group_bin']).assign(offset=1.)
			### without Group_bin as covariate
			#M = pd.DataFrame(np.ones((X.shape)))
			res = scan(X, y, "normal", K=K_dict[chrom], M=M, verbose=False)
			temp.loc[index, 'pval'] = float(res.stats.pv20)
			#print('M:')
			#print(M)
			#print('X:')
			#print(X['genotype'].tolist())
			#print('res:')
			#print(res)
			#print('res.stats:')
			#print(res.stats)
			#print('res.stats.pv20:')
			#print(float(res.stats.pv20))

	str_data_inDEG_b6b10 = pd.merge(str_data_inDEG_b6b10, temp, left_index=True, right_index=True)
	#print('str_data_inDEG_b6b10:')
	#print(str_data_inDEG_b6b10)

	file_out = 'pval_str_inDEG_b6b10_iGene.txt'
	print('file out: ', file_out)
	str_data_inDEG_b6b10.to_csv(file_out, sep='\t', header=True, index=False)

def plot_QQ():

	print('++++++++++++++++ QQplot +++++++++++++++++')
	file_sv_inDEG = 'pval_sv_inDEG_b6b10.txt'
	file_sv_iGene_inDEG = 'pval_sv_iGene_inDEG_b6b10.txt'
	file_cnv_inDEG = 'pval_cnv_inDEG_b6b10.txt'
	file_vep_inDEG = 'pval_vep_inDEG_b6b10.txt'
	file_vep_mis_inDEG = 'pval_vep_mis_inDEG_b6b10.txt'
	file_str_inDEG = 'pval_str_inDEG_b6b10_exon_TSS_UTR5_promoter_PeriodLT1.txt'
	file_str_iGene_inDEG = 'pval_str_inDEG_b6b10_iGene.txt'

	pval_sv_inDEG_B6B10CV = np.array(pd.read_table(file_sv_inDEG, header=0, index_col=None)['pval_B6B10CV'])
	pval_sv_inDEG_B6B10CV.sort()

	pval_sv_iGene_inDEG_B6B10CV = np.array(pd.read_table(file_sv_iGene_inDEG, header=0, index_col=None)['pval_B6B10CV'])
	pval_sv_iGene_inDEG_B6B10CV.sort()

	pval_cnv_inDEG_B6B10CV = np.array(pd.read_table(file_cnv_inDEG, header=0, index_col=None)['pval_B6B10CV'])
	pval_cnv_inDEG_B6B10CV.sort()

	pval_vep_inDEG_B6B10CV = np.array(pd.read_table(file_vep_inDEG, header=0, index_col=None)['pval_B6B10CV'])
	pval_vep_inDEG_B6B10CV.sort()

	pval_vep_mis_inDEG_B6B10CV = np.array(pd.read_table(file_vep_mis_inDEG, header=0, index_col=None)['pval_B6B10CV'])
	pval_vep_mis_inDEG_B6B10CV.sort()

	pval_str_inDEG_B6B10CV = np.array(pd.read_table(file_str_inDEG, header=0, index_col=None)['pval'])
	pval_str_inDEG_B6B10CV.sort()

	pval_str_iGene_inDEG_B6B10CV = np.array(pd.read_table(file_str_iGene_inDEG, header=0, index_col=None)['pval'])
	pval_str_iGene_inDEG_B6B10CV.sort()

	### B6B10CV
	fig = plt.figure()

	unif_n = pval_cnv_inDEG_B6B10CV.shape[0]
	unif_x = np.linspace(0,1,unif_n+2)[1:-1]
	plt.plot(-np.log10(unif_x), -np.log10(pval_cnv_inDEG_B6B10CV), 'o', fillstyle='none', label='SegDups')
	qt = np.linspace(0.1,0.9,9)
	data_qt = np.quantile(pval_cnv_inDEG_B6B10CV, qt)
	unif_qt = np.quantile(unif_x, qt)
	plt.plot(-np.log10(unif_qt), -np.log10(data_qt), 'ok')

	unif_n = pval_vep_inDEG_B6B10CV.shape[0]
	unif_x = np.linspace(0,1,unif_n+2)[1:-1]
	plt.plot(-np.log10(unif_x), -np.log10(pval_vep_inDEG_B6B10CV), 'o', fillstyle='none', label='SNP/INDEL, LoF')
	qt = np.linspace(0.1,0.9,9)
	data_qt = np.quantile(pval_vep_inDEG_B6B10CV, qt)
	unif_qt = np.quantile(unif_x, qt)
	plt.plot(-np.log10(unif_qt), -np.log10(data_qt), 'ok')

	unif_n = pval_sv_inDEG_B6B10CV.shape[0]
	unif_x = np.linspace(0,1,unif_n+2)[1:-1]
	plt.plot(-np.log10(unif_x), -np.log10(pval_sv_inDEG_B6B10CV), 'o', fillstyle='none', label='SV, Genic')
	qt = np.linspace(0.1,0.9,9)
	data_qt = np.quantile(pval_sv_inDEG_B6B10CV, qt)
	unif_qt = np.quantile(unif_x, qt)
	plt.plot(-np.log10(unif_qt), -np.log10(data_qt), 'ok')

	unif_n = pval_str_inDEG_B6B10CV.shape[0]
	unif_x = np.linspace(0,1,unif_n+2)[1:-1]
	plt.plot(-np.log10(unif_x), -np.log10(pval_str_inDEG_B6B10CV), 'o', fillstyle='none', label='STR, Genic')
	qt = np.linspace(0.1,0.9,9)
	data_qt = np.quantile(pval_str_inDEG_B6B10CV, qt)
	unif_qt = np.quantile(unif_x, qt)
	plt.plot(-np.log10(unif_qt), -np.log10(data_qt), 'ok')

	unif_n = pval_vep_mis_inDEG_B6B10CV.shape[0]
	unif_x = np.linspace(0,1,unif_n+2)[1:-1]
	plt.plot(-np.log10(unif_x), -np.log10(pval_vep_mis_inDEG_B6B10CV), 'o', fillstyle='none', label='SNP/INDEL, mis')
	qt = np.linspace(0.1,0.9,9)
	data_qt = np.quantile(pval_vep_mis_inDEG_B6B10CV, qt)
	unif_qt = np.quantile(unif_x, qt)
	plt.plot(-np.log10(unif_qt), -np.log10(data_qt), 'ok')

	unif_n = pval_sv_iGene_inDEG_B6B10CV.shape[0]
	unif_x = np.linspace(0,1,unif_n+2)[1:-1]
	plt.plot(-np.log10(unif_x), -np.log10(pval_sv_iGene_inDEG_B6B10CV), 'o', fillstyle='none', label='SV, Intergenic')
	qt = np.linspace(0.1,0.9,9)
	data_qt = np.quantile(pval_sv_iGene_inDEG_B6B10CV, qt)
	unif_qt = np.quantile(unif_x, qt)
	plt.plot(-np.log10(unif_qt), -np.log10(data_qt), 'ok')

	unif_n = pval_str_iGene_inDEG_B6B10CV.shape[0]
	unif_x = np.linspace(0,1,unif_n+2)[1:-1]
	plt.plot(-np.log10(unif_x), -np.log10(pval_str_iGene_inDEG_B6B10CV), 'o', fillstyle='none', label='STR, Intergenic')
	qt = np.linspace(0.1,0.9,9)
	data_qt = np.quantile(pval_str_iGene_inDEG_B6B10CV, qt)
	unif_qt = np.quantile(unif_x, qt)
	plt.plot(-np.log10(unif_qt), -np.log10(data_qt), 'ok')

	unif_n = max([pval_sv_iGene_inDEG_B6B10CV.shape[0], pval_sv_inDEG_B6B10CV.shape[0], pval_cnv_inDEG_B6B10CV.shape[0], pval_vep_inDEG_B6B10CV.shape[0], pval_str_inDEG_B6B10CV.shape[0], pval_str_iGene_inDEG_B6B10CV.shape[0]])
	unif_x = np.linspace(0,1,unif_n+2)[1:-1]
	plt.plot([-np.log10(unif_x[0]), -np.log10(unif_x[-1])], [-np.log10(unif_x[0]), -np.log10(unif_x[-1])], '-r')

	plt.xlabel('-log10(Expected)')
	plt.ylabel('-log10(Observed)')
	plt.legend()
	file_names = ['QQplot_LMM_B6B10CV.png', 'QQplot_LMM_B6B10CV.pdf']
	for file_name in file_names:
		print(file_name)
		fig.savefig(file_name, bbox_inches='tight')
	plt.close()

def compute_fdr():
	print('++++++++++++++++ Compute FDR +++++++++++++++++')
	file_sv_inDEG = 'pval_sv_inDEG_b6b10.txt'
	data_sv_inDEG = pd.read_table(file_sv_inDEG, header=0, index_col=None)
	data_sv_inDEG['fdr_noCV'] = sm.stats.multipletests(data_sv_inDEG.pval_noCV, alpha=0.1, method='fdr_bh')[1]
	data_sv_inDEG['fdr_B6B10CV'] = sm.stats.multipletests(data_sv_inDEG.pval_B6B10CV, alpha=0.1, method='fdr_bh')[1]
	#print('data_sv_inDEG:')
	#print(data_sv_inDEG)
	file_out = 'fdr_sv_inDEG_b6b10.txt'
	print('file out: ', file_out)
	data_sv_inDEG.to_csv(file_out, sep='\t', header=True, index=False)

	file_vep_inDEG = 'pval_vep_inDEG_b6b10.txt'
	data_vep_inDEG = pd.read_table(file_vep_inDEG, header=0, index_col=None)
	data_vep_inDEG['fdr_noCV'] = sm.stats.multipletests(data_vep_inDEG.pval_noCV, alpha=0.1, method='fdr_bh')[1]
	data_vep_inDEG['fdr_B6B10CV'] = sm.stats.multipletests(data_vep_inDEG.pval_B6B10CV, alpha=0.1, method='fdr_bh')[1]
	#print('data_vep_inDEG:')
	#print(data_vep_inDEG)
	file_out = 'fdr_vep_inDEG_b6b10.txt'
	print('file out: ', file_out)
	data_vep_inDEG.to_csv(file_out, sep='\t', header=True, index=False)

	file_cnv_inDEG = 'pval_cnv_inDEG_b6b10.txt'
	data_cnv_inDEG = pd.read_table(file_cnv_inDEG, header=0, index_col=None)
	data_cnv_inDEG['fdr_noCV'] = sm.stats.multipletests(data_cnv_inDEG.pval_noCV, alpha=0.1, method='fdr_bh')[1]
	data_cnv_inDEG['fdr_B6B10CV'] = sm.stats.multipletests(data_cnv_inDEG.pval_B6B10CV, alpha=0.1, method='fdr_bh')[1]
	#print('data_cnv_inDEG:')
	#print(data_cnv_inDEG)
	file_out = 'fdr_cnv_inDEG_b6b10.txt'
	print('file out: ', file_out)
	data_cnv_inDEG.to_csv(file_out, sep='\t', header=True, index=False)

	file_str_inDEG = 'pval_str_inDEG_b6b10_exon_TSS_UTR5_promoter_PeriodLT1.txt'
	data_str_inDEG = pd.read_table(file_str_inDEG, header=0, index_col=None)
	data_str_inDEG['fdr'] = sm.stats.multipletests(data_str_inDEG.pval, alpha=0.1, method='fdr_bh')[1]
	#print('data_str_inDEG:')
	#print(data_str_inDEG)
	file_out = 'fdr_str_inDEG_b6b10.txt'
	print('file out: ', file_out)
	data_str_inDEG.to_csv(file_out, sep='\t', header=True, index=False)

if __name__ == '__main__':

	run_test_sv()

	run_test_cnv()

	run_test_snp()

	run_test_str()

	plot_QQ()

	compute_fdr()
