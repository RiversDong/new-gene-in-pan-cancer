import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy.stats as stats
import numpy as np

plt.rcParams.update({'font.family': 'times new roman', 'font.weight': "bold"})
plt.rcParams.update({'font.size': 13})

def tao(expression_profile_pd, tao_type="Yanai"):
    row_num, col_num = expression_profile_pd.shape
    tao_list = []
    gene_list = []
    for i in range(0, row_num):
        expression_profile = list(expression_profile_pd.iloc[i, 1:-1].astype(float))
        gene = expression_profile_pd.iloc[i, 0].split("(")[0].strip()
        if tao_type == "Yanai":
            tpm_max = max(expression_profile)
            gene_list.append(gene)
            if tpm_max != 0:
                tmp_lis = []
                for j in expression_profile:
                    tmp_lis.append(1 - (j / tpm_max))
                sum_tmp_lis = sum(tmp_lis)
                tao = sum_tmp_lis / (len(expression_profile) - 1)
                tao_list.append(tao)
            else:
                tao_list.append("NA")
        else:
            if tao_type == "Liao":
                new_expression_profile = []
                for j in expression_profile:
                    if j == 0:
                        new_expression_profile.append(2)
                    else:
                        new_expression_profile.append(j)
                tpm_max = max(new_expression_profile)
                tpm_max = np.log(tpm_max)
                tmp_lis = []
                for j in new_expression_profile:
                    tmp_lis.append(1 - (np.log(j)) / tpm_max)
                sum_tmp_lis = sum(tmp_lis)
                tao = sum_tmp_lis / (len(expression_profile) - 1)
                tao_list.append(tao)
    pd_tao = pd.DataFrame()
    pd_tao["Gene"] = gene_list
    pd_tao["tao"] = tao_list
    return pd_tao

pan_score_data = pd.read_excel(r"D:\essential_gene_evo\result\Sup.xlsx")
data = pd.read_csv(r"D:\essential_gene_evo\data\OmicsExpressionProteinCodingGenesTPMLogp1.csv", sep=",", header=None)
data = data.transpose()
data.columns = data.iloc[0]
new_data = data.iloc[1:]
pd_tao = tao(new_data, tao_type="Yanai")

merged_df = pd.merge(pan_score_data, pd_tao, how='outer', on="Gene")
merged_df = merged_df.dropna(how='any')

merged_df.to_csv(r"D:\essential_gene_evo\result\panscore_tao.csv", sep=",", index=False)

merged_df['panscore_level'] = pd.qcut(merged_df['findex'], 8, labels=False)

# 计算每个 panscore_level 的 tao 均值
level_mean_tao = merged_df.groupby('panscore_level')['tao'].mean().reset_index()


slope, intercept, r_value, p_value, std_err = stats.linregress(merged_df["findex"], merged_df['tao'])
print(f"相关系数: {r_value}, p值: {p_value}")