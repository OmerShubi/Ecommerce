from typing import List, Set, Dict
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import networkx
import copy

INFECTED = 'infected'
SUSPECTIBLE = 'suspectible'

def LTM(graph: networkx.Graph, patients_0: List, iterations: int) -> Set:

    for node in graph.nodes:
        if node in patients_0:
            graph.nodes[node]['status'] = INFECTED
        else:
            graph.nodes[node]['status'] = SUSPECTIBLE

        neighbors = G.neighbors(node)
        num_neighbors = 0
        for _ in neighbors:
            num_neighbors += 1
        G.nodes[node]['num_neighbors'] = num_neighbors
        G.nodes[node]['concern'] = 0

    for _ in range(1, iterations+1):
        to_be_infected = []
        for node in G.nodes():
            neighbors = G.neighbors(node)
            node_threshold = 0
            num_infected_neighbors = 0
            for neighbor_id in neighbors:
                if G.nodes[neighbor_id]['status'] == INFECTED:
                    num_infected_neighbors += 1
                    node_threshold += G.get_edge_data(node, neighbor_id)['weight']
            if CONTAGION*node_threshold >= 1 + G.nodes[node]['concern']:
                to_be_infected.append(node)
            G.nodes[node]['num_infected_neighbors'] = num_infected_neighbors

        for node in G.nodes():
            if node in to_be_infected:
                G.nodes[node]['status'] = INFECTED

        for node in G.nodes():
            if G.nodes[node]['num_neighbors'] > 0:
                concern = G.nodes[node]['num_infected_neighbors'] / G.nodes[node]['num_neighbors']
            else:
                concern = 0
            G.nodes[node]['concern'] = concern

    total_infected = [n for n in G.nodes() if G.nodes[n]['status']==INFECTED]
    return set(total_infected)


def ICM(graph: networkx.Graph, patients_0: List, iterations: int) -> [Set, Set]:
    total_infected = set(patients_0)
    total_deceased = set()
    # TODO implement your code here
    return total_infected, total_deceased


def plot_degree_histogram(histogram: Dict):
    plt.bar(histogram.keys(), histogram.values(), width=0.80, color="b")

    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    plt.show()


def calc_degree_histogram(graph: networkx.Graph) -> Dict:
    """
    Example:
    if histogram[1] = 10 -> 10 nodes have only 1 friend
    """
    res = nx.degree_histogram(graph)
    histogram = {deg: num_nodes for deg, num_nodes in enumerate(res)}
    return histogram


def build_graph(filename: str) -> networkx.Graph:
    G = networkx.Graph()
    df = pd.read_csv(filename)
    G.add_nodes_from(np.unique(df.loc[:, ['from', 'to']].values))
    if 'w' in df.columns:
        G.add_weighted_edges_from(df.loc[:, ['from', 'to', 'w']].values)
    else:
        G.add_edges_from(df.loc[:, ['from', 'to']].values)
    return G


def clustering_coefficient(graph: networkx.Graph) -> float:
    triangles = sum(nx.triangles(graph).values()) / 3
    if triangles == 0:
        return 0
    # TODO did we miss any edge cases? because of devision by zero!
    transitivity = nx.transitivity(graph)
    triads = (3 * triangles) / transitivity
    open_triangles = triads - 3 * triangles
    cc = 3 * triangles / (3 * triangles + open_triangles)

    return cc


def compute_lethality_effect(graph: networkx.Graph, t: int) -> [Dict, Dict]:
    global LETHALITY
    mean_deaths = {}
    mean_infected = {}
    for l in (.05, .15, .3, .5, .7):
        LETHALITY = l
        for iteration in range(30):
            G = copy.deepcopy(graph)
            patients_0 = np.random.choice(list(G.nodes), size=50, replace=False, p=None)
            # TODO implement your code here

    return mean_deaths, mean_infected


def plot_lethality_effect(mean_deaths: Dict, mean_infected: Dict):
    # TODO implement your code here
    ...


def choose_who_to_vaccinate(graph: networkx.Graph) -> List:
    people_to_vaccinate = []
    # TODO implement your code here
    return people_to_vaccinate


def choose_who_to_vaccinate_example(graph: networkx.Graph) -> List:
    """
    The following heuristic for Part C is simply taking the top 50 friendly people;
     that is, it returns the top 50 nodes in the graph with the highest degree.
    """
    node2degree = dict(graph.degree)
    sorted_nodes = sorted(node2degree.items(), key=lambda item: item[1], reverse=True)[:50]
    people_to_vaccinate = [node[0] for node in sorted_nodes]
    return people_to_vaccinate


"Global Hyper-parameters"
CONTAGION = 1.05
LETHALITY = .15

if __name__ == "__main__":
    # sample = "sample.csv"
    # plt.subplot(121)
    # G = build_graph(filename=filename)
    # nx.draw(G, with_labels=True, font_weight='bold')
    # plt.show()
    PART_A = False
    if PART_A:
        filenames = ["PartA1.csv", "PartA2.csv"]
        for filename in filenames:
            G = build_graph(filename=filename)
            G_hist = calc_degree_histogram(G)
            plot_degree_histogram(G_hist)
            cc = clustering_coefficient(G)
            print(f"{filename} - {cc}")

    patients0 = pd.read_csv("patients0.csv", header=None)
    G = build_graph(filename="PartB-C.csv")
    infected = LTM(graph=G, patients_0=patients0[:20].values, iterations=6)
    print(len(infected))