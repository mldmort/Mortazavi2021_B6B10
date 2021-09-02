import pandas as pd

Data_Dir = '../data/'

def process_rnaseq_raw():
	expr_raw_file = Data_Dir+'b6b10_rnaseq_raw.txt'
	expr_raw_data = pd.read_table(expr_raw_file, header=0)
	print("expr_raw_data:")
	print(expr_raw_data)

	id_file = Data_Dir+'b6b10_rnaseq_ids.txt'
	id_data = pd.read_table(id_file, header=0)
	shape = id_data.shape
	id_data = id_data.iloc[range(0,shape[0],4),:]
	print('id_data:')
	print(id_data)

	expr_data = pd.DataFrame(expr_raw_data['gene'])
	for ic in range(1,expr_raw_data.shape[1],4):
		expr_data = pd.concat([expr_data, expr_raw_data.iloc[:,ic:ic+4].sum(axis=1)], axis=1)
	expr_data.columns = ['gene',*list(id_data['ID'])]
	#print('expr_data:')
	#print(expr_data)

	# remove 2 samples with too many reads (samples 58826 & 58827)
	id_data.set_index('ID', inplace=True)
	id_data.drop([58826, 58827], axis=0, inplace=True)
	#print('id_data:')
	#print(id_data)
	
	expr_data.set_index('gene', inplace=True)
	expr_data.drop([58826, 58827], axis=1, inplace=True)
	expr_data = expr_data.T
	expr_data = pd.concat([id_data['Strain'], expr_data],axis=1)
	#print('expr_data:')
	#print(expr_data)

	strain_dict_b6 = {'B6N-TyrC_BrdCrlCrl':'R134', 'C57BL_6ByJ':'R72', 'C57BL_6J':'R2', 'C57BL_6JBomTac':'R144', 'C57BL_6JEiJ':'R84', 'C57BL_6NCrl':'R114', 'C57BL_6NHsd':'R124', 'C57BL_6NJ':'R52', 'C57BL_6NTac':'R104'}
	strain_dict_b10 = {'C57BL_10ScNHsd':'R94', 'C57BL_10ScSnJ':'R22', 'C57BL_10SnJ':'R32', 'C57BL_10ScCr':'R42', 'C57BL_10J':'R62'}

	expr_data.loc[expr_data['Strain'].isin(strain_dict_b6.keys()), 'Group'] = 'B6'
	expr_data.loc[expr_data['Strain'].isin(strain_dict_b10.keys()), 'Group'] = 'B10'
	columns = ['Strain', 'Group', *expr_data.columns.tolist()[1:-1]]
	expr_data = expr_data[columns]
	#print('expr_data:')
	#print(expr_data)

	file_out = 'b6b10_counts_raw.txt'
	print(file_out)
	expr_data.to_csv(file_out, sep='\t')

if __name__ == '__main__':

	process_rnaseq_raw()
	
