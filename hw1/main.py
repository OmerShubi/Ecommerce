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
            num_infected_neighbors = 0
            num_removed_neighbors = 0
            for neighbor in graph.neighbors(node):
                if graph.nodes[neighbor]['status'] == INFECTED:
                    num_infected_neighbors += 1
                elif graph.nodes[neighbor]['status'] == REMOVED:
                    num_removed_neighbors += 1

                if graph.nodes[neighbor]['NI'] and (graph.nodes[node]['status'] == SUSPECTIBLE):
                    infect_prob = CONTAGION * graph.get_edge_data(neighbor, node)['weight'] * (
                            1 - graph.nodes[node]['concern'])
                    infect_prob = min(1, infect_prob)
                    if np.random.binomial(1, infect_prob):
                        if np.random.binomial(1, LETHALITY):
                            to_be_removed.append(node)
                        else:
                            to_be_infected.append(node)

            graph.nodes[node]['num_infected_neighbors'] = num_infected_neighbors
            graph.nodes[node]['num_removed_neighbors'] = num_removed_neighbors

        for node in graph.nodes():
            if node in to_be_infected:
                graph.nodes[node]['status'] = INFECTED
                graph.nodes[node]['NI'] = True
            elif node in to_be_removed:
                graph.nodes[node]['status'] = REMOVED
                graph.nodes[node]['NI'] = False
            else:
                graph.nodes[node]['NI'] = False

        for node in graph.nodes():
            if graph.nodes[node]['num_neighbors'] > 0:
                concern = (graph.nodes[node]['num_infected_neighbors'] + 3 * graph.nodes[node]['num_removed_neighbors']) \
                          / graph.nodes[node]['num_neighbors']
                concern = min(concern, 1)
            else:
                concern = 0
            graph.nodes[node]['concern'] = concern

    total_infected = [n for n in graph.nodes() if graph.nodes[n]['status'] == INFECTED]
    total_deceased = [n for n in graph.nodes() if graph.nodes[n]['status'] == REMOVED]

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
    G = networkx.Graph()
    df = pd.read_csv(filename)
    G.add_nodes_from(np.unique(df.loc[:, ['from', 'to']].values))
    if 'w' in df.columns:
        G.add_weighted_edges_from(df.loc[:, ['from', 'to', 'w']].values)
    else:
        G.add_edges_from(df.loc[:, ['from', 'to']].values)
    return G


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
        n = len(G[node])
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
            startTime = time.time()

            G = copy.deepcopy(graph)
            patients_0 = np.random.choice(list(G.nodes), size=50, replace=False, p=None)
            infected, removed = ICM(graph=G, patients_0=patients_0, iterations=t)
            mean_infected[l] += len(infected)
            mean_deaths[l] += len(removed)
            executionTime = (time.time() - startTime)
            print(f"l-{l}, iteration-{iteration} took {executionTime}")
        mean_infected[l] /= 30
        mean_deaths[l] /= 30

    return mean_deaths, mean_infected


def plot_lethality_effect(mean_deaths: Dict, mean_infected: Dict):
    plt.plot(mean_infected.keys(), mean_infected.values(), label='mean infected')
    plt.plot(mean_deaths.keys(), mean_deaths.values(), label='mean deaths')
    plt.title("Infections and death vs lethality")
    plt.ylabel("Mean Counts")
    plt.xlabel("Lethality")
    plt.legend()
    plt.show()
    plt.savefig("Infections and death vs lethality.png")


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
CONTAGION = 0.8
LETHALITY = .2

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

    PART_B = True
    if PART_B:
        patients0 = pd.read_csv("patients0.csv", header=None)
        G = build_graph(filename="PartB-C.csv")
        infected = LTM(graph=copy.deepcopy(G), patients_0=patients0[:20].values, iterations=6)
        print(f"LTM num infected - {len(infected)}")
        i = 0
        r = 0
        trials = 4
        num_patients0 = 50
        iterations = 6
        for _ in range(trials):
            infected, removed = ICM(graph=copy.deepcopy(G), patients_0=patients0[:num_patients0].values, iterations=iterations)
            i += len(infected)
            r += len(removed)
        print(f"ICM patients0-{num_patients0} iterations-{iterations}, mean infected - {i / trials}, mean removed - {r / trials}")

        # mean_deaths, mean_infected = compute_lethality_effect(graph=copy.deepcopy(G), t=6)
        # plot_lethality_effect(mean_deaths=mean_deaths, mean_infected=mean_infected)
