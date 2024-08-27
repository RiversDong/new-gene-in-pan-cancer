import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib as mpl
mpl.rcParams["pdf.fonttype"]=42
mpl.rcParams["font.family"] = 'Arial'


def id_g():
	d = {}
	file = open(r'D:\essential_gene_evo\data\Homo_sapiens.GRCh38.108.gtf')
	for line in file.readlines()[5:]:
		if line.strip().split('\t')[2] == 'CDS':
			protein_id = ''
			gene_id = ''
			for i in line.strip().split('\t')[8].split(';'):
				if 'ENSG' in i:
					gene_id = i.split('"')[1]
				if 'gene_name' in i:
					name = i.split('"')[1]
					d[gene_id] = name
	return d

file2 = r'D:\essential_gene_evo\new_bin\13059_2022_2821_MOESM2_ESM_Ma_gb_2022.xlsx'
df_AGE = pd.read_excel(file2,sheet_name='S1',header = 3,index_col=1)
UC_gene = df_AGE[(df_AGE.age == 'UC')].index.tolist()
EM_gene = df_AGE[(df_AGE.age == 'EM')].index.tolist()
MM_gene = df_AGE[(df_AGE.age == 'MM')].index.tolist()
PSG_gene = df_AGE[(df_AGE.age == 'PSG')].index.tolist()


#gnomAD
df = pd.read_csv(r'D:\essential_gene_evo\data\release_2.1.1_constraint_gnomad.v2.1.1.lof_metrics.by_gene.txt\release_2.1.1_constraint_gnomad.v2.1.1.lof_metrics.by_gene.txt',header = 0,index_col =0,sep = '\t')
df = df.oe_lof_upper
df1 = df[df<0.35] 
gnomeAD_geme = df1.index.tolist()

file = r'D:\essential_gene_evo\new_bin\pan_score.csv'
df_panscore = pd.read_csv(file,index_col=0)
df_panscore.columns = ['panscore']

index_panscore=list(df_panscore.index)
dict_pan = df_panscore.panscore.to_dict()
sorted_index = df_panscore.sort_values('panscore', ascending=False).head(len(gnomeAD_geme)).index.tolist()
print(sorted_index)

UC1 = set(UC_gene) & set(gnomeAD_geme)
UC2 = set(UC_gene) & set(sorted_index)
EM1 = set(EM_gene) & set(gnomeAD_geme)
EM2 = set(EM_gene) & set(sorted_index)
MM1 = set(MM_gene) & set(gnomeAD_geme)
MM2 = set(MM_gene) & set(sorted_index)
PSG1 = set(PSG_gene) & set(gnomeAD_geme)
PSG2 = set(PSG_gene) & set(sorted_index)
print(len(UC1))
print(len(UC2))
print(len(EM1))
print(len(EM2))
print(len(MM1))
print(len(MM2))
print(len(PSG1))
print(len(PSG2))

fig, axes = plt.subplots(1, 4, figsize=(13, 8))
# 绘制韦恩图
venn2([UC1, UC2], set_labels=('gnomeAD', 'Sanger'),ax = axes[0], set_colors=("#1aafd0", "#3be8b0"))
axes[0].set_title('Comparision of UC',fontsize=18)
venn2([EM1, EM2], set_labels=('gnomeAD', 'Sanger'),ax = axes[1],set_colors=("#1aafd0", "#3be8b0"))
axes[1].set_title('Comparision of EM',fontsize=18)
venn2([MM1, MM2], set_labels=('gnomeAD', 'Sanger'),ax = axes[2],set_colors=("#1aafd0", "#3be8b0"))
axes[2].set_title('Comparision of MM',fontsize=18)
venn2([PSG1,PSG2], set_labels=('gnomeAD', 'Sanger'),ax = axes[3],set_colors=("#1aafd0", "#3be8b0"))
axes[3].set_title('Comparision of PSG',fontsize=18)
# 显示图形
plt.tight_layout()
#plt.savefig(r"D:\essential_gene_evo\new_figure_1\ess_divergent_cutoff.pdf")
plt.show()
