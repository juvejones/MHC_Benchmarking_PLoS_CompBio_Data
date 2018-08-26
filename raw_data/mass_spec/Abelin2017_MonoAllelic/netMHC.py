import subprocess
import re
import pandas as pd

alleles = ['A01:01', 'A02:01', 'A02:03', 'A02:04', 'A02:07','A03:01',
           'A24:02', 'A29:02', 'A31:01', 'A68:02', 'B35:01','B44:02',
           'B44:03', 'B51:01', 'B54:01']

for allele in alleles:
    allelename = re.sub('[:]','',allele)
    print "Calculating %s\n" % (allele)
    filename = allelename+".seq.test.1col"
    seq = {}
    with open(allelename+'.seq.test', 'r') as f:
        for line in f:
            seq[line.strip().split()[0].upper()] = line.strip().split()[1]
    with open(filename, 'w') as f:
        f.writelines(["%s\n" % key for key in seq.keys()])
    print len(seq.keys())

    output = filename+"_netMHCpan3"
    subprocess.call(['netMHCpan', '-a', "HLA-"+str(allele), '-l 9', '-p', '-f', filename, '-inptype 1', '-xls', '-xlsfile', output])
    print "Finishing %s\n" % (allele)
    
    df = pd.read_csv(output, sep="\s+", skiprows=1, header=0)
    eluted = []
    for x in seq.keys():
        eluted.append(seq[x])
    df.loc[:,'eluted'] = eluted
    #df.to_csv(output+"_eluted", sep="\t", header=True)    
