import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib as mlp

mlp.rcParams["pdf.fonttype"]=42
plt.rcParams.update({'font.family':'arial'})

#plt.rc('font', size=12, family='Times New Roman')
#plt.rc('font', weight='bold')
#plt.rc('axes', labelsize=12, titlesize=12, titleweight='bold')

data = pd.read_csv(r"D:\essential_gene_evo\data\release_2.1.1_constraint_gnomad.v2.1.1.lof_metrics.by_gene.txt\release_2.1.1_constraint_gnomad.v2.1.1.lof_metrics.by_gene.txt", sep="\t")
gnomad_genes = data[data["oe_lof_upper"]<0.35]["gene"]

sanger_data = pd.read_excel(r"D:\essential_gene_evo\result\Sup.xlsx", sheet_name="pan-score", )
sanger_genes = sanger_data["Gene"]
print("panscore", len(sanger_genes))

# age data
def get_age(age_file, age_type):
    data = pd.read_csv(age_file, sep="\t", low_memory="False")
    gene_age = data[["Symbol", age_type]]
    gene_age = gene_age.dropna(axis=0, how="any")
    return gene_age
gene2age_4 = get_age(r"D:\essential_gene_evo\data\13059_2022_2821_MOESM2_ESM_Ma_gb_2022.txt", "4 merged age groups")

agegroup=["UC", "EM", "MM", "PSG"]

# 读取pan-score文件
panscore = pd.read_excel(r"D:\essential_gene_evo\result\Sup.xlsx", sheet_name="pan-score")
fig, axes = plt.subplots(1, 4, figsize=(13, 8))
fig_index = 0
for i in agegroup:
    i_age_genes = gene2age_4[gene2age_4["4 merged age groups"] == i]["Symbol"]
    i_gnomad_genes = set(i_age_genes).intersection(set(gnomad_genes))
    print(i , "gnomad >0.35", len(gnomad_genes))
    print("nomad+age", i, len(i_gnomad_genes))

    i_sanger_genes = set(i_age_genes).intersection(set(sanger_genes))
    print("nomad+age", i, len(i_sanger_genes))
    print(i_sanger_genes)
    venn2([i_gnomad_genes, i_sanger_genes], set_labels=('Gnomad', 'Sanger'), ax=axes[fig_index], set_colors=("#1aafd0", "#3be8b0"))
    axes[fig_index].set_title('Comparision of '+i)
    fig_index+=1
plt.tight_layout()
#plt.savefig(r"D:\essential_gene_evo\new_figure_1\ess_divergent.pdf")
plt.show()






