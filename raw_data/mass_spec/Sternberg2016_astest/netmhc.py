import pandas as pd
import numpy as np
import random
import subprocess

dbDir = "/SFS/user/ctc/zhaoweil/ElutionDataTest/"
inputDir = dbDir+"MS_curated_data/Sternberg2016_astest/"
alleles=["A0301", 
	 "B0702",
	 "A2402",
	 "B1501",
	 "B4001",
	 "B5801",
	 "B0801",
	 "B3901",
	 "A0101",
	 "A2601",
	 "A0201",
	 "B2705"]

'''
Perform predictions by netMHCpan3
'''
f = open(inputDir+"netMHC_precision", 'w')
for item in alleles:
    alleleinput = item[:3] + ':' + item[-2:]
    print alleleinput
    filename1 = inputDir+item+"_Sternberg2016_test"
    filename2 = inputDir+item+"_Sternberg2016.9mer"
    table = pd.read_table(filename1, header=None, names = ["peptide", "eluted", "ic50"], index_col=False)
    table = table.loc[table['peptide'].str.len() == 9]
    eluted = table['eluted'].tolist()
    with open(filename2,'w') as ninemerfile:
        ninemerfile.writelines(["%s\n" % line for line in table['peptide'].tolist()])
    output = inputDir+item+"_Sternberg2016_test.out"
    subprocess.call(['netMHCpan', '-a', "HLA-"+str(alleleinput), '-l 9', '-p', '-f', filename2, '-inptype 1', '-xls', '-xlsfile', output])
    print "Finishing %s\n" % (alleleinput)
    
    df = pd.read_csv(output, sep="\s+", skiprows=1, header=0)
    df.loc[:,'eluted'] = eluted
    pred500 = len(df.loc[df['nM']<500])
    pred50 = len(df.loc[df['nM']<50])
    real500 = len(df.loc[(df['eluted'] == 1) & (df['nM']<500)])
    real50 = len(df.loc[(df['eluted'] == 1) & (df['nM']<50)])
    if not pred500 == 0:
        precision500 = float(real500)/float(pred500)
    else:
        precision500 = 0
    if not pred50 == 0:
        precision50 = float(real50)/float(pred50)
    else:
        precision50 = 0
    f.writelines("%s\t%d\t%d\t%.5f\t%d\t%d\t%.5f\n" % (item, pred500, real500, precision500, pred50, real50, precision50))
f.close()
