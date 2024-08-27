import pandas as pd
from gtfparse import read_gtf
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 13

def get_age(data, age_type):
    gene_age = data[["Gene", age_type]]
    gene_age = gene_age.dropna(axis=0, how="any")
    return gene_age

def get_interaction(interaction, young, old, outfile):
    young = list(young); old=list(old)
    out = open(outfile, "w")
    data_interaction = pd.read_csv(interaction, sep="\t")
    out.write("gene1\tgene2\trank\ttype\n")
    for index, row in data_interaction.iterrows():
        gene1 = row["gene1"]; gene2 = row["gene2"]; rank = row["rank"]
        if (gene1 in young) and (gene2 in young):
            gene1_gene2_type = "YY"
            info = f"{gene1}\t{gene2}\t{rank}\t{gene1_gene2_type}\n"
            out.write(info)
        if (gene1 in young and gene2 in old) or (gene1 in old and gene2 in young):
            gene1_gene2_type = "YO"
            info = f"{gene1}\t{gene2}\t{rank}\t{gene1_gene2_type}\n"
            out.write(info)
        if (gene1 in old) and (gene2 in old):
            gene1_gene2_type = "OO"
            info = f"{gene1}\t{gene2}\t{rank}\t{gene1_gene2_type}\n"
            out.write(info)
    out.close()

def get_type_num(in_type):
    data = pd.read_csv(in_type, sep="\t")
    YY_row, col = data[data["type"] == "YY"].shape
    YO_row, col = data[data["type"] == "YO"].shape
    OO_row, col = data[data["type"] == "OO"].shape

    return YY_row, YO_row, OO_row

def read_file(infile):
    data = pd.read_csv(infile, sep="\t")
    gene1 = data["gene1"]; gene2 = data["gene2"]
    gene1 = set(gene1);  gene2=set(gene2)
    genes = gene1.union(gene2)
    gene_num = len(genes)
    link_num, column = data.shape
    
    return gene_num, link_num



if __name__ == "__main__":

    files = [("top_tcga_mutual_15.1", "top_gtex_mutual_15.1"), ("top_tcga_mutual_20", "top_gtex_mutual_20")]
    
    res_out = open("/data/chuand/essential_gene_evo/result/tcga_net_statistics", "w")
    
    res_out.write("\tNYY\tNYO\tNOO\n")
    
    labels = ['PSG-PSG', 'PSG-Other', 'Other-Other']
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
    width = 0.35
    axes_flat = axes.flatten()
    axes_index = 0
    for i, j in files:

        interaction_TCGA = f"/data/chuand/essential_gene_evo/data/{i}"
        interaction_GTEx = f"/data/chuand/essential_gene_evo/data/{j}"
        
        gene_num, link_num = read_file(interaction_TCGA)
        print(i, gene_num, link_num, link_num/gene_num)
        gene_num, link_num = read_file(interaction_GTEx)
        print(j, gene_num, link_num, link_num/gene_num)

        TCGA_OUT = f"/data/chuand/essential_gene_evo/data/TCGA_interaction_type_{i}"
        GTEx_OUT = f"/data/chuand/essential_gene_evo/data/GTEx_interaction_type_{j}" 

        age = "/data/chuand/essential_gene_evo/data/13059_2022_2821_MOESM2_ESM_Ma_gb_2022.txt"
        age_data = pd.read_csv(age, sep="\t", low_memory="False")
        gene2age_14 = get_age(age_data, "14 age groups")
        gene2age_4 = get_age(age_data, "4 merged age groups")

        all_list = []
        young = gene2age_4[(gene2age_4["4 merged age groups"]=="PSG")]["Gene"]
        all_list.extend(young)
        old = gene2age_4[(gene2age_4["4 merged age groups"]!="PSG")]["Gene"]
        all_list.extend(old)
        
        get_interaction(interaction_TCGA, young, old, TCGA_OUT)
        get_interaction(interaction_GTEx, young, old, GTEx_OUT)

        
        expected_young = len(list(combinations(young, 2)))/len(list(combinations(all_list, 2)))
        expected_young = round(expected_young, 6)

        expected_old = len(list(combinations(old, 2)))/len(list(combinations(all_list, 2)))
        expected_old = round(expected_old, 6)
        
        expected_young_old = (len(young)*len(old))/len(list(combinations(all_list, 2)))
        expected_young_old = round(expected_young_old, 6)

        tcga_YY_num, tcga_YO_num, tcga_OO_num = get_type_num(TCGA_OUT)
        gtex_YY_num, gtex_YO_num, gtex_OO_num = get_type_num(GTEx_OUT)

        #print(f"YY: {expected_young}")
        #print(f"YO: {expected_young_old}")
        #print(f"OO: {expected_old}")

        tcga_num = tcga_YY_num + tcga_YO_num + tcga_OO_num
        gtex_num = gtex_YY_num + gtex_YO_num + gtex_OO_num

        tcga_e_yy, tcga_e_yo, tcga_e_oo = tcga_num*expected_young, tcga_num*expected_young_old, tcga_num*expected_old
        gtex_e_yy, gtex_e_yo, gtex_e_oo = gtex_num*expected_young, gtex_num*expected_young_old, gtex_num*expected_old
        
        # tcga中观察和期望
        tcga_yy_ei = (tcga_YY_num - tcga_e_yy)/tcga_e_yy; tcga_yy_ei = round(tcga_yy_ei, 4)*100
        tcga_yo_ei = (tcga_YO_num - tcga_e_yo)/tcga_e_yo; tcga_yo_ei = round(tcga_yo_ei, 4)*100
        tcga_oo_ei = (tcga_OO_num - tcga_e_oo)/tcga_e_oo; tcga_oo_ei = round(tcga_oo_ei, 4)*100
        
        # GTEx中的观察和期望
        gtex_yy_ei = (gtex_YY_num-gtex_e_yy)/gtex_e_yy; gtex_yy_ei = round(gtex_yy_ei, 4)*100
        gtex_yo_ei = (gtex_YO_num-gtex_e_yo)/gtex_e_yo; gtex_yo_ei = round(gtex_yo_ei, 4)*100
        gtex_oo_ei = (gtex_OO_num-gtex_e_oo)/gtex_e_yo; gtex_yo_ei = round(gtex_yo_ei, 4)*100
        
        print(i, tcga_yy_ei, tcga_yo_ei, tcga_oo_ei)
        print(j, gtex_yy_ei, gtex_yo_ei, gtex_oo_ei)
        
        # 绘制图片 tcga_YY_num + tcga_YO_num + tcga_OO_num， gtex_YY_num + gtex_YO_num + gtex_OO_num
        tcga_obs = [tcga_YY_num, tcga_YO_num, tcga_OO_num]
        gtex_obs = [gtex_YY_num, gtex_YO_num, gtex_OO_num]
        x = np.arange(len(labels))
        # Aetna
        rects1 = axes_flat[axes_index].bar(x - width/2, tcga_obs, width, label='Cancer', color="#d20962")
        rects2 = axes_flat[axes_index].bar(x + width/2, gtex_obs, width, label='GTEx', color="#f47721")

        for rect in rects1 + rects2:
            height = rect.get_height()
            axes_flat[axes_index].annotate(f'{height:.0f}',
                        xy=(rect.get_x() + rect.get_width()/2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

        # 对数坐标
        axes_flat[axes_index].set_yscale('log')

        # 添加标签和标题
        axes_flat[axes_index].set_ylabel('Number of coexpression pairs')
        #axes[axes_index].set_title('Comparison of Counts between TCGA and GTEx')
        axes_flat[axes_index].set_xticks(x)
        axes_flat[axes_index].set_xticklabels(labels, rotation=45)
        axes_flat[axes_index].legend()

        axes_flat[axes_index].spines['right'].set_visible(False)
        axes_flat[axes_index].spines['top'].set_visible(False)
        
        # 绘制ei图
        tcga_ei = [tcga_yy_ei, tcga_yo_ei, tcga_oo_ei]
        gtex_ei = [gtex_yy_ei, gtex_yo_ei, gtex_oo_ei]
        rects3 = axes_flat[axes_index+1].bar(x - width/2, tcga_ei, width, label='Cancer', color="#d20962")
        rects4 = axes_flat[axes_index+1].bar(x + width/2, gtex_ei, width, label='GTEx', color="#f47721")
        for rect in rects3 + rects4:
            height = rect.get_height()
            axes_flat[axes_index + 1].annotate(f'{height:.0f}',
                        xy=(rect.get_x() + rect.get_width()/2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')
        axes_flat[axes_index + 1].set_ylabel('The excessed percentage')
        axes_flat[axes_index + 1].set_xticks(x)
        axes_flat[axes_index + 1].set_xticklabels(labels, rotation=45)
        axes_flat[axes_index + 1].legend()
        axes_flat[axes_index + 1].spines['right'].set_visible(False)
        axes_flat[axes_index + 1].spines['top'].set_visible(False)
        
        
        axes_index += 2
        res_out.write(f"E_TCGA\t{tcga_e_yy}\t{tcga_e_yo}\t{tcga_e_oo}\t\n")
        res_out.write(f"O_TCGA\t{tcga_YY_num}\t{tcga_YO_num}\t{tcga_OO_num}\t\n")
        
        res_out.write(f"E_GTEx\t{gtex_e_yy}\t{gtex_e_yo}\t{gtex_e_oo}\t\n")
        res_out.write(f"O_GTEx\t{gtex_YY_num}\t{gtex_YO_num}\t{gtex_OO_num}\t\n")

    res_out.close()

    plt.tight_layout()
    plt.savefig("/data/chuand/essential_gene_evo/result/get_interaction_type.png", bbox_inches='tight', transparent=True, dpi=300)



