import pandas as pd
from gtfparse import read_gtf
import igraph as ig

top_rank = 20

def get_rank14(infile, outfile):
    data = pd.read_csv(infile, sep="\t", chunksize=5000)
    out = open(outfile, "w")
    for i in data:
        data1 = i
        for index, row in data1.iterrows(): 
                gene=index
                corr = row[row<top_rank]
                corr = list(corr.index)
                if gene in corr:
                    corr.remove(gene)
                out.write("\t".join(corr) + "\t" + gene + "\n")
    data = pd.read_csv(infile, sep="\t")
    return data
    out.close()

def getIntersection(infile, outfile, rank):
    rank14 = infile
    data = open(rank14).read().split("\n")
    genes = set()
    gtf = read_gtf("/data/chuand/essential_gene_evo/data/Homo_sapiens.GRCh38.102.gtf")
    proteins = set(gtf[gtf["gene_biotype"]=="protein_coding"]["gene_id"])
    print(len(proteins))
    genes1 = set()
    OUT = open(outfile, "w")
    OUT.write("gene1" + "\t" + "gene2" +  "rank" + "\n")
    tmp = set()
    for i in data:
        row = i.split("\t")
        tgene = row[-1]
        pgene = row[0:-1]
        genes1 = genes1.union(set(row))
        if pgene:
            pgene_protein = set(pgene).intersection(proteins)
            if len(pgene_protein) != 0 and tgene in proteins:
                genes = genes.union(pgene_protein)
                genes.add(tgene)
                for j in pgene_protein:
                    interaction = j + tgene
                    if interaction not in tmp:
                        interaction_rank = rank.loc[tgene, j]
                        OUT.write(tgene + "\t" + j + "\t" + str(interaction_rank) + "\n")
                        tmp.add(tgene+j)
    OUT.close()
    print(len(genes))
    print(len(genes1.intersection(proteins)))

# 针对TCGA
rank = get_rank14("/data/chuand/essential_gene_evo/data/fromChen/mutual_rank_matrix_tcga.tsv", "/data/chuand/essential_gene_evo/data/tcga_mutual")
getIntersection("/data/chuand/essential_gene_evo/data/tcga_mutual", f"/data/chuand/essential_gene_evo/data/top_tcga_mutual_{top_rank}", rank)

# 针对GTEx
rank = get_rank14("/data/chuand/essential_gene_evo/data/fromChen/mutual_rank_matrix_gtex.tsv", "/data/chuand/essential_gene_evo/data/gtex_mutual")
getIntersection("/data/chuand/essential_gene_evo/data/gtex_mutual", f"/data/chuand/essential_gene_evo/data/top_gtex_mutual_{top_rank}", rank)
