import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser(description="Read input file")
parser.add_argument('-i', dest='file', help='assign filename to read')

args = parser.parse_args()

with open(args.file, 'r') as f:
	headline = f.readline()
	hlaName = re.sub('[:]','',headline.strip().split()[0])

data = pd.read_csv(args.file, header=0, skiprows=1, sep='\t')
print(hlaName)

refData = pd.read_csv("../mhcflurry_pan/Rdata/"+hlaName+".pan.self.txt",header=0,index_col=0,sep=',')
refData['peptide'] = refData['peptide'].str.strip()

dataMerge = refData.merge(data[['Peptide', '1-log50k']], left_index=True,right_index=True,how='inner')

assert len(dataMerge) == len(refData), "%d\t%d" % (len(dataMerge), len(refData))

dataMerge[['peptide','meas','1-log50k']].to_csv("./Rdata/"+hlaName+".NetMHCpan4.txt",header=True,index=True,sep=',')
