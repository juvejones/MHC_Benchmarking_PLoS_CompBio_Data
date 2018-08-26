import pandas as pd
import argparse
import re
import numpy as np

parser = argparse.ArgumentParser(description="Read input file")
parser.add_argument('-i', dest='file', help='assign filename to read')

args = parser.parse_args()
hlaName = []

with open(args.file, 'r') as f:
	for line in f.readlines():
		hlaName.append(line.strip())

sf = open("./Rdata/summary.txt",'a')

exclude = ['B3801','B1501','A3201','A0206','A0205','A0203',
		'B4403','B2705','A3201','A1101','A0301','B0702']
for item in hlaName:
	print(item)

	refData = pd.read_csv("../input_seq/HLA-"+item+".txt",header=None,sep='\t')
	refData.columns = ["HLA","Peptide","meas_contin"]
	meas_bi = np.zeros(len(refData['Peptide']))
	meas_bi[refData['meas_contin'].astype(float) < 500.0] = 1
	meas_bi[refData['meas_contin'].astype(float) >= 500.0] = 0
	refData['meas_bi'] = meas_bi

	data = pd.read_csv("./output/HLA-"+item+".pred",header=0,skiprows=11,sep="\t")

	dataMerge = pd.concat([refData,data[['Peptide','Max_score','Rank']]],axis=1)

	assert len(dataMerge) == len(refData), "%d\t%d" % (len(dataMerge), len(refData))
	
	item = re.sub('[:*]','',item)

	dataMerge[['Peptide','meas_bi','meas_contin','Max_score','Rank']].to_csv("./Rdata/HLA-"+item+".pred",header=True,index=False,sep=',')

	if not item in exclude:
		sf.writelines("##"+item+"\n")
		dataMerge[['Peptide','meas_bi','meas_contin','Max_score','Rank']].to_csv(sf, header=True, index=False, sep=",")

sf.close()
