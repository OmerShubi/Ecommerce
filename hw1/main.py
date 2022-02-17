import math
import random
from typing import List, Set, Dict
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import networkx
import copy
from collections import Counter
from math import factorial
import time

INFECTED = 'infected'
SUSPECTIBLE = 'suspectible'
REMOVED = 'removed'


def LTM(graph: networkx.Graph, patients_0: List, iterations: int) -> Set:
    global CONTAGION
    print(CONTAGION)
    for node in graph.nodes:
        if node in patients_0:
            graph.nodes[node]['status'] = INFECTED
        else:
            graph.nodes[node]['status'] = SUSPECTIBLE

        neighbors = graph.neighbors(node)
        num_neighbors = 0
        for _ in neighbors:
            num_neighbors += 1
        graph.nodes[node]['num_neighbors'] = num_neighbors
        graph.nodes[node]['concern'] = 0

    for _ in range(1, iterations + 1):
        to_be_infected = []
        for node in graph.nodes():
            neighbors = graph.neighbors(node)
            node_threshold = 0
            num_infected_neighbors = 0
            for neighbor_id in neighbors:
                if graph.nodes[neighbor_id]['status'] == INFECTED:
                    num_infected_neighbors += 1
                    node_threshold += graph.get_edge_data(node, neighbor_id)['weight']
            if graph.nodes[node]['status'] == SUSPECTIBLE and (
                    CONTAGION * node_threshold >= 1 + graph.nodes[node]['concern']):
                to_be_infected.append(node)
            graph.nodes[node]['num_infected_neighbors'] = num_infected_neighbors

        for node in graph.nodes():
            if node in to_be_infected:
                graph.nodes[node]['status'] = INFECTED

        for node in graph.nodes():
            if graph.nodes[node]['num_neighbors'] > 0:
                concern = graph.nodes[node]['num_infected_neighbors'] / graph.nodes[node]['num_neighbors']
            else:
                concern = 0
            graph.nodes[node]['concern'] = concern

    total_infected = [n for n in graph.nodes() if graph.nodes[n]['status'] == INFECTED]
    return set(total_infected)


def ICM(graph: networkx.Graph, patients_0: List, iterations: int) -> [Set, Set]:
    startTime = time.time()

    # Initialization
    for node in graph.nodes():
        # Update infected status for patients0
        if node in patients_0:
            # Remove patients based on lethality
            if np.random.binomial(1, LETHALITY):
                graph.nodes[node]['status'] = REMOVED
                graph.nodes[node]['NI'] = False

            else:
                graph.nodes[node]['status'] = INFECTED
                graph.nodes[node]['NI'] = True
        else:
            graph.nodes[node]['NI'] = False
            graph.nodes[node]['status'] = SUSPECTIBLE

        # Count number of neighbors for each node
        neighbors = graph.neighbors(node)
        num_neighbors = 0
        for _ in neighbors:
            num_neighbors += 1
        graph.nodes[node]['num_neighbors'] = num_neighbors
        # Initialize concern with 0
        graph.nodes[node]['concern'] = 0

    for _ in range(1, iterations + 1):
        to_be_infected = []
        to_be_removed = []
        for node in graph.nodes():
            if graph.nodes[node]['status'] == SUSPECTIBLE:
                num_infected_neighbors = 0
                num_removed_neighbors = 0
                for neighbor in graph.neighbors(node):
                    if graph.nodes[neighbor]['status'] == INFECTED:
                        num_infected_neighbors += 1
                    elif graph.nodes[neighbor]['status'] == REMOVED:
                        num_removed_neighbors += 1

                    if graph.nodes[neighbor]['NI'] and graph.nodes[neighbor]['status'] == INFECTED:
                        infect_prob = CONTAGION * graph.get_edge_data(neighbor, node)['weight'] * (
                                1 - graph.nodes[node]['concern'])
                        infect_prob = min(1, infect_prob)
                        if np.random.binomial(1, infect_prob):
                            if np.random.binomial(1, LETHALITY):
                                to_be_removed.append(node)
                            else:
                                to_be_infected.append(node)
                            break

                graph.nodes[node]['num_infected_neighbors'] = num_infected_neighbors
                graph.nodes[node]['num_removed_neighbors'] = num_removed_neighbors

        for node in graph.nodes():
            if node in set(to_be_infected):
                graph.nodes[node]['status'] = INFECTED
                graph.nodes[node]['NI'] = True
            elif node in set(to_be_removed):
                graph.nodes[node]['status'] = REMOVED
                graph.nodes[node]['NI'] = False
            else:
                graph.nodes[node]['NI'] = False

        for node in graph.nodes():
            if graph.nodes[node]['status'] == SUSPECTIBLE:
                if graph.nodes[node]['num_neighbors'] > 0:
                    concern = (graph.nodes[node]['num_infected_neighbors'] + 3 * graph.nodes[node]['num_removed_neighbors']) \
                              / graph.nodes[node]['num_neighbors']
                    concern = min(concern, 1)
                else:
                    concern = 0
                graph.nodes[node]['concern'] = concern

    total_infected = [n for n in graph.nodes() if graph.nodes[n]['status'] == INFECTED]
    total_deceased = [n for n in graph.nodes() if graph.nodes[n]['status'] == REMOVED]

    executionTime = (time.time() - startTime)
    print(f"ICM took {executionTime}")

    return set(total_infected), set(total_deceased)


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
    graph = networkx.Graph()
    df = pd.read_csv(filename)
    graph.add_nodes_from(np.unique(df.loc[:, ['from', 'to']].values))
    if 'w' in df.columns:
        graph.add_weighted_edges_from(df.loc[:, ['from', 'to', 'w']].values)
    else:
        graph.add_edges_from(df.loc[:, ['from', 'to']].values)
    return graph


def combinations(n, r):
    return factorial(n) // factorial(r) // factorial(n - r)


def clustering_coefficient(graph: networkx.Graph) -> float:
    triangles = 0
    for node in graph.nodes():
        node_neighbors = set(graph[node])
        common_neighbors_hist = Counter(len(node_neighbors & set(graph[w])) for w in node_neighbors)
        triangles += sum(num_common_neighbors * num_neighbors for num_common_neighbors, num_neighbors in common_neighbors_hist.items()) // 2

    triangles /= 3

    if triangles == 0:
        return 0

    triads = 0
    for node in graph.nodes():
        n = len(graph[node])
        if n >= 2:
            triads += combinations(n, 2)
    cc = 3 * triangles / triads

    return cc


def compute_lethality_effect(graph: networkx.Graph, t: int) -> [Dict, Dict]:
    global LETHALITY
    mean_deaths = {}
    mean_infected = {}
    for l in (.05, .15, .3, .5, .7):
        LETHALITY = l
        mean_infected[l] = 0
        mean_deaths[l] = 0
        for iteration in range(30):

            graph_copy = copy.deepcopy(graph)
            patients_0 = np.random.choice(list(graph_copy.nodes), size=50, replace=False, p=None)
            infected, removed = ICM(graph=graph_copy, patients_0=patients_0, iterations=t)
            mean_infected[l] += len(infected)
            mean_deaths[l] += len(removed)
            print(f"l-{l}, iteration-{iteration}")
        mean_infected[l] /= 30
        mean_deaths[l] /= 30

    return mean_infected, mean_deaths


def plot_lethality_effect(mean_deaths: Dict, mean_infected: Dict):
    plt.plot(mean_infected.keys(), mean_infected.values(), label='mean infected')
    plt.plot(mean_deaths.keys(), mean_deaths.values(), label='mean deaths')
    plt.title("Infections and death vs lethality")
    plt.ylabel("Mean Counts")
    plt.xlabel("Lethality")
    plt.legend()
    plt.show()
    plt.savefig("Infections and death vs lethality.png")


def vac_pagerank(graph: networkx.Graph):
    startTime = time.time()

    a = nx.pagerank(graph, weight='weight')
    executionTime = (time.time() - startTime)
    print(f"choose who to vac_pagerank took {executionTime}")
    return np.array(list(a.values()))


def vac_bridges(graph: networkx.Graph):
    startTime = time.time()
    x = np.array(list(nx.bridges(graph)))
    u, counts = np.unique(x.reshape((-1, 1)), return_counts=True)
    executionTime = (time.time() - startTime)
    print(f"choose who to vac_bridges took {executionTime}")
    return u.astype(int), counts


def vac_degree(graph: networkx.Graph):
    """
    The following heuristic for Part C is simply taking the top 50 friendly people;
     that is, it returns the top 50 nodes in the graph with the highest degree.
    """
    startTime = time.time()

    node2degree = dict(graph.degree)
    d = np.array(list(node2degree.values()))
    executionTime = (time.time() - startTime)
    print(f"choose who to vac_degree took {executionTime}")
    return d


def choose_who_to_vaccinate(graph: networkx.Graph) -> List:

    to_vac_degree = vac_degree(graph)
    to_vac_bridges_nodes, to_vac_bridges_vals = vac_bridges(graph)
    to_vac_pagerank = vac_pagerank(graph)
    to_vac_bridges = np.zeros_like(to_vac_degree)
    to_vac_bridges[to_vac_bridges_nodes] = to_vac_bridges_nodes
    joint = np.stack((to_vac_degree, to_vac_pagerank, to_vac_bridges))
    joint /= joint.max(axis=1).reshape((-1, 1))
    joint_sing_val = joint.sum(axis=0)
    num_nodes = 50
    people_to_vaccinate = np.argpartition(joint_sing_val, -num_nodes)[-num_nodes:]

    return list(people_to_vaccinate)


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
CONTAGION = 0.8
LETHALITY = .15

def run_all_three_parts():
    # sample = "sample.csv"
    # plt.subplot(121)
    # G = build_graph(filename=filename)
    # nx.draw(G, with_labels=True, font_weight='bold')
    # plt.show()
    PART_A = False
    if PART_A:
        filenames = ["PartA1.csv", "PartA2.csv", "PartB-C.csv"]
        for filename in filenames:
            graph = build_graph(filename=filename)
            G_hist = calc_degree_histogram(graph)
            plot_degree_histogram(G_hist)
            cc = clustering_coefficient(graph)
            print(f"{filename} - {cc}")

    PART_B = False
    if PART_B:
        patients0 = pd.read_csv("patients0.csv", header=None)
        graph = build_graph(filename="PartB-C.csv")

        infected = LTM(graph=copy.deepcopy(graph), patients_0=patients0[:50].values, iterations=6)
        print(f"LTM num infected - {len(infected)}")
        # infected = LTM(graph=copy.deepcopy(graph), patients_0=patients0[:48].values, iterations=6)
        # print(f"LTM num infected - {len(infected)}")
        # infected = LTM(graph=copy.deepcopy(graph), patients_0=patients0[:30].values, iterations=6)
        # print(f"LTM num infected - {len(infected)}")
        # infected = LTM(graph=copy.deepcopy(graph), patients_0=patients0[:20].values, iterations=6)
        # print(f"LTM num infected - {len(infected)}")

        # run_ICM(graph)
        mean_infected, mean_deaths = compute_lethality_effect(graph=copy.deepcopy(graph), t=6)
        plot_lethality_effect(mean_deaths=mean_deaths, mean_infected=mean_infected)

    COMPETITION = False
    if COMPETITION:
        graph = build_graph(filename="PartB-C.csv")
        graph_vac = copy.deepcopy(graph)
        startTime = time.time()

        to_vaccinate = choose_who_to_vaccinate(graph_vac)
        executionTime = (time.time() - startTime)
        print(f"choose who to vac took {executionTime}")
        graph_vac.remove_nodes_from(to_vaccinate)
        # run_ICM(graph_copy)
        np.random.seed(534)
        patients_0 = np.random.choice(list(graph_vac.nodes), size=50, replace=False, p=None)
        np.random.seed(534)
        infected, removed = ICM(graph=copy.deepcopy(graph_vac), patients_0=patients_0, iterations=6)
        print(f"ICM vac infected - {len(infected)},  removed - {len(removed)}")

        np.random.seed(534)
        infected, removed = ICM(graph=copy.deepcopy(graph), patients_0=patients_0, iterations=6)
        print(f"ICM  infected - {len(infected)},  removed - {len(removed)}")
