import networkx as nx
import igraph as ig
import leidenalg as la
import matplotlib.pyplot as plt

def read_interaction_file(file_path):
    interactions = []
    with open(file_path, 'r') as file:
        next(file)  # 跳过标题行
        for line in file:
            proteins = line.strip().split()
            #if len(proteins) == 2:
            interactions.append((proteins[0], proteins[1]))
    return interactions

def create_interaction_network(interactions):
    G = nx.Graph()
    G.add_edges_from(interactions)
    return G

def detect_communities_leiden(G):
    # Convert NetworkX graph to iGraph graph
    g = ig.Graph.TupleList(G.edges(), directed=False)
    # Perform Leiden community detection
    partition = la.find_partition(g, la.ModularityVertexPartition)
    return g, partition

def get_submodules(g, partition):
    submodules = {}
    for community_number, nodes in enumerate(partition):
        submodules[community_number] = [g.vs[node]['name'] for node in nodes]
    return submodules

def visualize_communities(G, g, partition, figure):
    pos = nx.spring_layout(G)
    cmap = plt.get_cmap('viridis')
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    
    for community_number, nodes in enumerate(partition):
        node_names = [g.vs[node]['name'] for node in nodes]
        nx.draw_networkx_nodes(G, pos, nodelist=node_names, node_size=50, node_color=[cmap(community_number / len(partition))])
    #plt.show()
    plt.tight_layout()
    plt.savefig(figure)

if __name__ == "__main__":
    
    files = ["GTEx_interaction_type_top_gtex_mutual_20", "TCGA_interaction_type_top_tcga_mutual_20"]
    for file_path in files:
        interactions = read_interaction_file(file_path)
        G = create_interaction_network(interactions)
        g, partition = detect_communities_leiden(G)
        submodules = get_submodules(g, partition)
        OUT = open(file_path+".clu", "w")
        for community, nodes in submodules.items():
            print(f"Community {community}: {nodes}")
            OUT.write(f"Community {community}: {nodes}\n")
        visualize_communities(G, g, partition, file_path+".pdf")
        OUT.close()