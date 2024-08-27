import pandas as pd
from random import sample
from scipy.stats import pearsonr
import numpy as np
from scipy.stats import sem
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family':'times new roman', 'font.weight':"bold", "font.size":15})

fessential = r"D:\essential_gene_evo\data\41586_2019_1103_MOESM4_ESM\Supplementary Table 2.xlsx"
ori_data = pd.read_excel(fessential, sheet_name="a. fitness gene summary", skiprows=2)
tissue = ori_data.columns[2:]; tissue_index = list(range(2, len(tissue)))

# 选择100进行计算
mean_r_list = []
sem_r_list = []
x = []
x_index = 1
xticklabels = []
for dada in range(10, 300, 10):
    rs = []
    x.append(x_index)
    xticklabels.append(str(x_index)+"-"+str(dada))
    for lala in range(10):
        selected = [0,1]
        selected.extend(sample(tissue_index, dada))
        data = ori_data.iloc[:,selected]
        new_tissue = data.columns[2:]
        tissue_number = len(new_tissue)
        cancer_type = data.iloc[0][2:]
        tissue2cancer_type = {}
        for i in range(0, tissue_number):
            i_tissue = tissue[i].split(".")[0]; i_cancer = cancer_type[i]
            if i_tissue not in tissue2cancer_type:
                tissue2cancer_type[i_tissue] = [i_cancer]
            else:
                tissue2cancer_type[i_tissue].append(i_cancer)

        #< 计算findex
        tmp_pd = pd.DataFrame()
        tmp_pd_col = 1
        findex_pd = pd.DataFrame()
        findex_pd["Gene"] = data["Tissue"][3:]
        for i in tissue2cancer_type:
            cancer_type = set(tissue2cancer_type[i])
            for j in cancer_type:
                #print(j)
                # data.loc[:,data.iloc[0] == "Prostate Carcinoma"]
                #data.loc[:,data.iloc[0] == j]
                j_data = data.loc[:,data.iloc[0] == j]
                j_data = j_data[3:]
                gene_num, cell_line_num = j_data.shape
                tmp_pd[str(tmp_pd_col)] = j_data.sum(axis=1)/cell_line_num # ni/Ni
                tmp_pd_col = tmp_pd_col + 1
        gene_num, cancer_type_num =  tmp_pd.shape
        findex = tmp_pd.sum(axis=1)/cancer_type_num
        findex_pd["findex"] = findex

        all_panscore = pd.read_excel(r'D:\essential_gene_evo\result\Sup.xlsx', sheet_name="pan-score")
        ori_findex = all_panscore["findex"]
        r, pvalue = pearsonr(findex, ori_findex)
        rs.append(r)
    min_rs = min(rs)
    max_rs = max(rs)
    mean_r = np.mean(rs)
    sem_r = sem(rs)
    print(dada, min_rs, mean_r, max_rs, sem_r)
    mean_r_list.append(mean_r)
    sem_r_list.append(sem_r)
    x_index = x_index + 1
fig,axes = plt.subplots(figsize=(6,4.5))
plt.errorbar(x, mean_r_list, sem_r_list,  fmt="o",linewidth=1.5, markersize = 3, elinewidth=2, capsize=4, ecolor="red")
plt.plot(x, mean_r_list, 'k-', linewidth = 1.2, color="black")
#sns.regplot(x, mean_r_list, ci=95, scatter_kws={"color": "black", "s": 22,"facecolor":"black"}, line_kws={"color": "red", 'alpha':0.5}, marker=".", ax=axes)
plt.xticks(ticks=x, labels=xticklabels, rotation=90, fontsize=13)
plt.ylabel("Correlations between pan-scores")
plt.xlabel("The number of samples")
plt.tight_layout()
#plt.savefig(r"D:\essential_gene_evo\figures\Monte_Carlo.pdf")
plt.show()
