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
    for index, row in expression_profile_pd.iterrows():
        gene = row["Description"]
        expression_profile = list(row.iloc[2:])
        
        if tao_type == "Yanai":
            tpm_max = max(expression_profile)
            gene_list.append(gene)
            if tpm_max != 0:
                tmp_lis = []
                for i in expression_profile:
                    tmp_lis.append(1 - (i / tpm_max))
                sum_tmp_lis = sum(tmp_lis)
                tao = sum_tmp_lis / (len(expression_profile) - 1)
                tao_list.append(tao)
            else:
                tao_list.append("NA")
        else:
            if tao_type == "Liao":
                new_expression_profile = []
                for i in expression_profile:
                    if i == 0:
                        new_expression_profile.append(2)
                    else:
                        new_expression_profile.append(i)
                tpm_max = max(new_expression_profile)
                tpm_max = np.log(tpm_max)
                tmp_lis = []
                for i in new_expression_profile:
                    tmp_lis.append(1 - (np.log(i)) / tpm_max)
                sum_tmp_lis = sum(tmp_lis)
                tao = sum_tmp_lis / (len(expression_profile) - 1)
                tao_list.append(tao)
    pd_tao = pd.DataFrame()
    pd_tao["Gene"] = gene_list
    pd_tao["tao"] = tao_list
    return pd_tao

pan_score_data = pd.read_excel(r"D:\essential_gene_evo\result\Sup.xlsx")

# 替换成GTEx的表达量
data = pd.read_csv(r"D:\essential_gene_evo\data\GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", sep="\t", skiprows=2)
pd_tao = tao(data, tao_type="Yanai")
merged_df = pd.merge(pan_score_data, pd_tao, how='outer', on="Gene")

merged_df.to_csv(r"D:\essential_gene_evo\result\gtex_panscore_tao.csv", sep=",", index=False)

merged_df = pd.read_csv(r"D:\essential_gene_evo\result\gtex_panscore_tao.csv", sep=",")
merged_df = merged_df.dropna()

slope, intercept, r_value, p_value, std_err = stats.linregress(merged_df["findex"], merged_df["tao"])
print(r_value, p_value)
