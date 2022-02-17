import copy
import time

import numpy as np
import pandas as pd
#
# from main import LTM, ICM, compute_lethality_effect, choose_who_to_vaccinate, build_graph, calc_degree_histogram, \
#     plot_degree_histogram, clustering_coefficient, plot_lethality_effect, CONTAGION
import main
if __name__ == '__main__':

    filenames = ["PartA1.csv", "PartA2.csv", "PartB-C.csv"]
    for filename in filenames:
        graph = main.build_graph(filename=filename)
        G_hist = main.calc_degree_histogram(graph)
        main.plot_degree_histogram(G_hist)
        cc = main.clustering_coefficient(graph)
        print(f"{filename} - {cc}")


    patients0 = pd.read_csv("patients0.csv", header=None)
    graph = main.build_graph(filename="PartB-C.csv")
    main.CONTAGION = 1.0
    infected = main.LTM(graph=copy.deepcopy(graph), patients_0=patients0[:50].values, iterations=6)
    print(f"LTM 50 num infected - {len(infected)}")
    infected = main.LTM(graph=copy.deepcopy(graph), patients_0=patients0[:48].values, iterations=6)
    print(f"LTM 48 num infected - {len(infected)}")
    infected = main.LTM(graph=copy.deepcopy(graph), patients_0=patients0[:30].values, iterations=6)
    print(f"LTM 30 num infected - {len(infected)}")
    infected = main.LTM(graph=copy.deepcopy(graph), patients_0=patients0[:20].values, iterations=6)
    print(f"LTM 20 num infected - {len(infected)}")
    main.CONTAGION = 1.05
    infected = main.LTM(graph=copy.deepcopy(graph), patients_0=patients0[:30].values, iterations=6)
    print(f"LTM 1.05 30 num infected - {len(infected)}")
    infected = main.LTM(graph=copy.deepcopy(graph), patients_0=patients0[:20].values, iterations=6)
    print(f"LTM 1.05 20 num infected - {len(infected)}")


    main.CONTAGION = 0.8
    main.LETHALITY = 0.2
    infected, removed = main.ICM(graph=copy.deepcopy(graph), patients_0=patients0[:50], iterations=6)
    print(f"ICM 6 infected - {len(infected)},  removed - {len(removed)}")
    infected, removed = main.ICM(graph=copy.deepcopy(graph), patients_0=patients0[:20], iterations=4)
    print(f"ICM 4 infected - {len(infected)},  removed - {len(removed)}")

    main.CONTAGION = 0.8
    mean_infected, mean_deaths = main.compute_lethality_effect(graph=copy.deepcopy(graph), t=6)
    main.plot_lethality_effect(mean_deaths=mean_deaths, mean_infected=mean_infected)

    main.CONTAGION = 0.8
    main.LETHALITY = 0.15

    graph = main.build_graph(filename="PartB-C.csv")
    graph_vac = copy.deepcopy(graph)
    startTime = time.time()

    to_vaccinate = main.choose_who_to_vaccinate(graph_vac)
    executionTime = (time.time() - startTime)
    print(f"choose who to vac took {executionTime}")
    graph_vac.remove_nodes_from(to_vaccinate)
    patients_0 = np.random.choice(list(graph_vac.nodes), size=50, replace=False, p=None)

    np.random.seed(534)
    infected, removed = main.ICM(graph=copy.deepcopy(graph_vac), patients_0=patients_0, iterations=6)
    print(f"ICM vac infected - {len(infected)},  removed - {len(removed)}")

    np.random.seed(534)
    infected, removed = main.ICM(graph=copy.deepcopy(graph), patients_0=patients_0, iterations=6)
    print(f"ICM  infected - {len(infected)},  removed - {len(removed)}")

