import pandas as pd
import numpy as np

alleles= ['A0101', 'A0201', 'A0203', 'A0204', 'A0207','A0301',
          'A2402', 'A2902', 'A3101', 'A6802', 'B3501','B4402',
          'B4403', 'B5101', 'B5401']
threshold_netmhc = [150, 50, 50, 250, 500, 150, 250, 50, 50, 50,
					100, 300, 350, 500, 100]
threshold_np = [0.5, 0.7, 0.55, 0.7, 0.7, 0.5, 0.65, 0.65, 0.5,
                0.55, 0.55, 0.5, 0.65, 0.65, 0.6]
netmhc = {}
np = {}

for idx, allele in enumerate(alleles):
	netmhc[allele] = threshold_netmhc[idx]
	np[allele] = threshold_np[idx]

for allele in alleles:
    df_1 = pd.read_table(allele+".seq.test.1col_netMHCpan3_eluted", header=0, index_col=False, sep="\t")
    df_2 = pd.read_table(allele+"_elute_monotrained_MSscorei_5xNeg", header=0, index_col=False, sep="\t")

    thresh1 = netmhc[allele]
    thresh2 = np[allele]
    df_combine = pd.DataFrame()
    peptide = []
    netmhc_pos = []
    np_pos = []
    real_pos = []

    for index, row in df_2.iterrows():
    	peptide.append(row[1])
    	if row[3] >= thresh2:
    		np_pos.append(1)
    	else:
    		np_pos.append(0)

    	real_pos.append(row[2])

    	netmhc_row = df_1[df_1['Peptide'] == row[1]]
    	netmhc_pred = netmhc_row.iloc[0]['nM']
    	if netmhc_pred < thresh1:
    		netmhc_pos.append(1)
    	else:
    		netmhc_pos.append(0)

    df_combine.loc[:,'peptide'] = peptide
    df_combine.loc[:,'np'] = np_pos
    df_combine.loc[:,'netmhc'] = netmhc_pos
    df_combine.loc[:,'real'] = real_pos

    np_out = len(df_combine[df_combine['np'] == 1])
    netmhc_out = len(df_combine[df_combine['netmhc'] == 1])
    co1 = len(df_combine[(df_combine['np']==1) & (df_combine['netmhc']==1)])
    co2 = len(df_combine[(df_combine['np']==1) & (df_combine['netmhc']==1) & (df_combine['real']==1)])
    
    print "%s\t%d\t%d\t%d\t%d" % (allele, np_out, netmhc_out, co1, co2)
