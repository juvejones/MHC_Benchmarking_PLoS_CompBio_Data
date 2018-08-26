#!/usr/bin/python

######################################
### Summarize Allele-Specific Data ###
######################################

import subprocess, os, re

def WriteMat(f, refline.cols, found, mat_auc, mat_srcc, i):
	for line in f:
		cols = line.split()
		cols[1] = re.sub('HLA-','',cols[1])
		if cols[1] == refline.cols[0]:
			found = True
			try:
				mat_auc[i].append(cols[2])
			except:
				mat_auc[i].append(0)
			try:
				mat_srcc[i].append(cols[3])
			except:
				mat_srcc[i].append(0)
	return found

__location__ = os.path.realpath(
    os.path.join(os.getcwd(),os.path.dirname(__file__)))
methods = ["smm","smmpmbec","ann","NetMHC4","PickPocket","consensus","NetMHCpan2.8","NetMHCpan3","NetMHCcons", "mhcflurry"]
mat_auc = []
mat_srcc = []

for i,method in enumerate(methods):
	mat_auc.append([])
	mat_srcc.append([])
	with open(os.path.join(__location__, "../", method, "Rdata/allele_stat.txt"),'r') as f:
		with open("../CountLogStrip",'r') as ref_file:
			reflines = ref_file.readlines()
			for refline in reflines:
				refline.cols = refline.strip().split(" ")
				found = False
				if int(refline.cols[1]) > 2500:
					found = WriteMat(f, refline.cols, found, mat_auc, mat_srcc, i)							
				if not found:
					mat_auc[i].append(0)
					mat_srcc[i].append(0)	
				f.seek(0)
			reflines = ref_file.readlines()
			for refline in reflines:
				refline.cols = refline.strip().split(" ")
				found = False
				if int(refline.cols[1]) > 500 and int(refline.cols[1]) < 2500:
					found = WriteMat(f, refline.cols, found, mat_auc, mat_srcc, i)
				if not found:
					mat_auc[i].append(0)
					mat_srcc[i].append(0)
				f.seek(0)
			reflines = ref_file.readlines()
			for refline in reflines:
				refline.cols = refline.strip().split(" ")
				found = False
				if int(refline.cols[1]) < 500:
					found = WriteMat(f, refline.cols, found, mat_auc, mat_srcc, i)
				if not found:
					mat_auc[i].append(0)
					mat_srcc[i].append(0)
				f.seek(0)

with open("table_all_alleles_auc.txt","w") as f:
	f.write("method")
	f.writelines(["\t%s" % refline.strip() for refline in reflines])
	f.write("\n")
	for i,method in enumerate(methods):
		f.write("%s" % method)
		f.writelines(["\t%.4f" % float(item) for item in mat_auc[i]])
		f.write("\n")

with open("table_all_alleles_srcc.txt","w") as f:
        f.write("method")
        f.writelines(["\t%s" % refline.strip() for refline in reflines])
        f.write("\n")
        for i,method in enumerate(methods):
                f.write("%s" % method)
                f.writelines(["\t%.4f" % float(item) for item in mat_srcc[i]])
                f.write("\n")
