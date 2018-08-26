import pandas as pd
import numpy as np
import matplotlib 
from matplotlib import pyplot as plt
import seaborn as sns
import sklearn
from sklearn import metrics
import sys

sys.path.append('/SFS/user/ctc/zhaoweil/bin/pyvenn/')
import venn

alleles= ['A0101', 'A0201', 'A0203', 'A0204', 'A0207','A0301',
          'A2402', 'A2902', 'A3101', 'A6802', 'B3501','B4402',
          'B4403', 'B5101', 'B5401']
# alleles = ['A2902', 'B5101', 'B5401']
thresholds_ic50 = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
thresholds_np = np.arange(0.2,1,0.05)

allele_result_netmhc = {}
allele_result_np = {}
cutoff_np = []
cutoff_netmhc = []
pred_labels = {}

for allele in alleles:
    # NetMHCpan results
    df = pd.read_table(allele+".seq.test.1col_netMHCpan3_eluted", header=0, sep="\t")
    df_report = pd.DataFrame()
    precision_list = []
    recall_list = []
    f1 = []
    auc = []
    for thresh in thresholds_ic50:

        try:
        	precision=float(len(df[(df['eluted']==1) & (df['nM']<=thresh)]))/float(len(df[df['nM']<= thresh]))
        	precision_list.append(precision)
        except:
        	precision = 0
        	precision_list.append(precision)
        recall=float(len(df[(df['eluted']==1) & (df['nM']<=thresh)]))/float(len(df[df['eluted']==1]))
        recall_list.append(recall)
        f1.append(2*precision*recall/(recall+precision))

    df_report.loc[:,'threshold'] = thresholds_ic50
    df_report.loc[:,'precision'] = precision_list
    df_report.loc[:,'recall'] = recall_list
    df_report.loc[:,'f1'] = f1
    
    allele_result_netmhc[allele] = df_report.loc[df_report['f1'] == df_report['f1'].max()].values.tolist()[0]
    cutoff_netmhc.append(df_report.f1.tolist())

    nM_norm = 1 - np.log10(df.nM.tolist())/np.log10([50000.0]*len(df.nM))
    allele_result_netmhc[allele].append(metrics.roc_auc_score(df.eluted.tolist(), nM_norm))

    # Get label lists for NetMHC prediction
    netmhc_best_cutoff = allele_result_netmhc[allele][0]
    netmhc_labels = df.loc[df['nM'] <= netmhc_best_cutoff].Peptide.values

    # NP results
    df = pd.read_table(allele+"_elute_monotrained_MSscorei_5xNeg", header=0, sep="\t")
    df_report = pd.DataFrame()
    precision_list = []
    recall_list = []
    f1 = []
    auc = []
    for thresh in thresholds_np:

        try:
        	precision=float(len(df[(df['eluted']==1) & (df['eluted_prob']>=thresh)]))/float(len(df[df['eluted_prob']>= thresh]))
        	precision_list.append(precision)
        except:
        	precision = 0
        	precision_list.append(precision)
        recall=float(len(df[(df['eluted']==1) & (df['eluted_prob']>=thresh)]))/float(len(df[df['eluted']==1]))
        recall_list.append(recall)
        try:
        	f1.append(2*precision*recall/(recall+precision))
        except:
        	f1.append(0.0)

    df_report.loc[:,'threshold'] = thresholds_np
    df_report.loc[:,'precision'] = precision_list
    df_report.loc[:,'recall'] = recall_list
    df_report.loc[:,'f1'] = f1
    
    allele_result_np[allele] = df_report.loc[df_report['f1'] == df_report['f1'].max()].values.tolist()[0]
    cutoff_np.append(df_report.f1.tolist())
    
    allele_result_np[allele].append(metrics.roc_auc_score(df.eluted.tolist(), df.eluted_prob.tolist()))

    # Get label lists for NP prediction
    np_best_cutoff = allele_result_np[allele][0]
    np_labels = df.loc[df['eluted_prob'] >= np_best_cutoff].peptide.values

    # Assemble label lists: [NP, NetMHCpan, Exp.True]
    true_labels = df.loc[df['eluted'] == 1].peptide.values
    pred_labels[allele] = [np_labels, netmhc_labels, true_labels]


# Make plots 
sns.set_context("paper")
sns.set(style="whitegrid", rc={"grid.linewidth": 0.5})
assert len(allele_result_netmhc.keys()) == len(allele_result_np.keys()), "%d\t\t%d" \
	% (len(allele_result_netmhc.keys()),len(allele_result_np.keys()))
hla_list = []
threshs_map = []
ppv_map = []
f1_map = []
auc_map = []
for k1 in allele_result_netmhc:
	if k1 in allele_result_np:
		hla_list.append(k1)
		threshs_map.append([allele_result_netmhc[k1][0], allele_result_np[k1][0]])
		ppv_map.append([allele_result_netmhc[k1][1], allele_result_np[k1][1]])
		f1_map.append([allele_result_netmhc[k1][3], allele_result_np[k1][3]])
		auc_map.append([allele_result_netmhc[k1][4], allele_result_np[k1][4]])
stats_df = pd.DataFrame({"HLA.Type": hla_list+hla_list,
						"Specificity": np.append(np.array(ppv_map)[:,0],np.array(ppv_map)[:,1]),
						"F1": np.append(np.array(f1_map)[:,0],np.array(f1_map)[:,1]),
						"AUC": np.append(np.array(auc_map)[:,0],np.array(auc_map)[:,1])})
stats_df.loc[:,"Predictor"] = ["NetMHCpan"] * len(hla_list) + ["NP"] * len(hla_list)
stats_melt = pd.melt(stats_df, id_vars="Predictor", value_vars=["Specificity","F1","AUC"])

# data for plotting cutoff heatmap
cutoff_df = pd.DataFrame(cutoff_np, columns=thresholds_np, index=alleles)
# cutoff_df = pd.DataFrame(cutoff_netmhc, columns=thresholds_ic50, index=alleles)

# Initialize the figure
f, axes = plt.subplots(1, 2, figsize=(14,6))
# f, ax1 = plt.subplots(figsize=(7.5, 6))

# Show change of F1 as varying threshold
sns.heatmap(cutoff_df, cmap="coolwarm", ax=axes[0], cbar_kws={'label': 'F1 Score'})

# Set labels
axes[0].set_yticklabels(alleles, rotation=0)
axes[0].set_xticklabels(thresholds_np)
axes[0].set_ylabel("HLA Allele", fontsize=14)
axes[0].set_xlabel("Prediction Cutoff", fontsize=14)
# ax1.set_title("Heatmap of F1 Score at Different NetMHCpan Cutoff")
axes[0].tick_params(labelsize=10)

# plt.savefig("f1_s.png", dpi=600)
# sys.exit()

# Second figure
# sns.despine(left=True)
# Show each observation with a scatterplot
sns.stripplot(y="value", x="variable", hue="Predictor",
              data=stats_melt, dodge=True, jitter=True,
              alpha=.5, zorder=1, s=10, ax=axes[1])

# Show the conditional means
sns.pointplot(y="value", x="variable", hue="Predictor",
              data=stats_melt, dodge=.532, join=False, palette="dark",
              markers="d", scale=1.4, ci=None, ax=axes[1])

# Improve the legend 
handles, labels = axes[1].get_legend_handles_labels()
axes[1].legend(handles[2:], labels[2:], title="Predictor",
          handletextpad=0, columnspacing=1,
          loc="upper left", ncol=2, frameon=True)

# Set labels
axes[1].set_xlabel("Metrics", fontsize=18)
axes[1].set_ylabel("Value", fontsize=18)
axes[1].set_ylim(0.4,1.1)
axes[1].tick_params(labelsize=14)

plt.tight_layout(w_pad = 2.5)
plt.savefig("f1ab_new.png", dpi=600)
plt.close()

labels_A0101 = venn.get_labels(pred_labels['A0101'],fill=['number', 'logic'])
labels_A0201 = venn.get_labels(pred_labels['A0201'],fill=['number', 'logic'])
labels_B3501 = venn.get_labels(pred_labels['B3501'],fill=['number', 'logic'])

# Format label a bit
for key in labels_A0101:
	labels_A0101[key] = labels_A0101[key][4:]
for key in labels_A0201:
	labels_A0201[key] = labels_A0201[key][4:]	
for key in labels_B3501:
	labels_B3501[key] = labels_B3501[key][4:]

# Plot venn
#matplotlib.rcParams.update({"font.size": 40})
sns.set_context("talk")
ax3 = venn.venn3(labels_A0101, names=['NP', 'NetMHCpan', 'MS.True'])
plt.savefig("A0101.venn_new.png", dpi=300, bbox_inches='tight')
plt.close()
ax4 = venn.venn3(labels_A0201, names=['NP', 'NetMHCpan', 'MS.True'])
plt.savefig("A0201.venn_new.png", dpi=300, bbox_inches='tight')
plt.close()
ax5 = venn.venn3(labels_B3501, names=['NP', 'NetMHCpan', 'MS.True'])
plt.savefig("B3501.venn_new.png", dpi=300, bbox_inches='tight')
plt.close()



		


