#!/usr/bin/python

######################################
### Summarize Allele-Specific Data ###
######################################

import subprocess, os, re

def WriteMat(f, key, mat_auc, mat_srcc, i):
	found = False
	for line in f:
		cols = line.split()
		cols[1] = re.sub('HLA-','',cols[1])
		if cols[1] == key:
			found = True
			try:
				mat_auc[i].append(float(cols[2]))
			except:
				mat_auc[i].append(0)
			try:
				mat_srcc[i].append(float(cols[3]))
			except:
				mat_srcc[i].append(0)
	if not found:
		mat_auc[i].append(0)
		mat_srcc[i].append(0)

__location__ = os.path.realpath(
    os.path.join(os.getcwd(),os.path.dirname(__file__)))
methods = ["smm","smmpmbec","ann","NetMHC4","PickPocket","consensus","NetMHCpan2.8","NetMHCpan3","NetMHCpan4","NetMHCcons", "MixMHCpred"]
mat_auc = []
mat_srcc = []
allele_dict1 = {}
allele_dict2 = {}
allele_dict3 = {}

with open("../CountLogStrip",'r') as ref_file:
	reflines = ref_file.readlines()
	for refline in reflines:
		refcols = refline.strip().split(" ")
		if int(refcols[1]) > 2500:
			key = refcols[0]
			allele_dict1[key] = refcols[1]
		elif int(refcols[1]) > 500:
			key = refcols[0]
			allele_dict2[key] = refcols[1]
		else:
			key = refcols[0]
			allele_dict3[key] = refcols[1]	
print(allele_dict1)
print(allele_dict2)
print(allele_dict3)

for i,method in enumerate(methods):
	mat_auc.append([])
	mat_srcc.append([])
	with open(os.path.join(__location__, "../", method, "Rdata/allele_stat.txt"),'r') as f:
		for key, value in allele_dict1.iteritems():
			WriteMat(f, key, mat_auc, mat_srcc, i)
			f.seek(0)
		for key, value in allele_dict2.iteritems():
			WriteMat(f, key, mat_auc, mat_srcc, i)
			f.seek(0)
		for key, value in allele_dict3.iteritems():
			WriteMat(f, key, mat_auc, mat_srcc, i)
			f.seek(0)							

#with open("table_all_alleles_auc_MixMHCpred.txt","w") as f:
#	f.write("method")
#	f.writelines(["\t%s" % key for key, value in allele_dict1.iteritems()])
#	f.writelines(["\t%s" % key for key, value in allele_dict2.iteritems()])
#	f.writelines(["\t%s" % key for key, value in allele_dict3.iteritems()])
#	f.write("\n")
#	for i,method in enumerate(methods):
#		f.write("%s" % method)
#		f.writelines(["\t%.4f" % float(item) for item in mat_auc[i]])
#		f.write("\n")

with open("table_all_alleles_srcc_MixMHCpred.txt","w") as f:
	f.write("method")
	f.writelines(["\t%s" % key for key, value in allele_dict1.iteritems()])
	f.writelines(["\t%s" % key for key, value in allele_dict2.iteritems()])
	f.writelines(["\t%s" % key for key, value in allele_dict3.iteritems()])
	f.write("\n")
	for i,method in enumerate(methods):
		f.write("%s" % method)
		f.writelines(["\t%.4f" % float(item) for item in mat_srcc[i]])
		f.write("\n")
