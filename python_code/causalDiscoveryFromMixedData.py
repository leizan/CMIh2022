import sys

sys.path.append('..')
from src.lib.causal_discovery_from_mixed_data.graphResult import GraphResult
from src.lib.PCGCE.baselines.scripts_python.pcgce import pcgce_mixed

import networkx as nx
# import matplotlib.pyplot as plt
import numpy as np


class CausalDiscoveryFromMixedData:
    def __init__(self):
        self.tau_max = None
        self.pc_alpha = None
        self.sig_samples = None
        self.not_continuous_columns = None

        self.position_of_nodes = {'CpuG': (3, 0.2),
                                  'CpuH': (2.5 * np.cos(36 / 180 * np.pi), 2.5 * np.sin(36 / 180 * np.pi)),
                                  'NPH': (2.5 * np.cos(76 / 180 * np.pi), 2.5 * np.sin(76 / 180 * np.pi)),
                                  'RamH': (3 * np.cos(116 / 180 * np.pi), 2.5 * np.sin(116 / 180 * np.pi)),
                                  'NetOut': (3 * np.cos(208 / 180 * np.pi), -0.5),
                                  'NCM': (2.5 * np.cos(248 / 180 * np.pi), -0.8),
                                  'NPP': (-1, 0.8),
                                  'CpuP': (1, 0.5),
                                  'DiskW': (1.2, -0.8),
                                  'NetIn': (-3, 0.8),
                                  'code200': (2.5 * np.cos(76 / 180 * np.pi), 2.5 * np.sin(76 / 180 * np.pi) + 1.5),
                                  'code206': (2.5 * np.cos(76 / 180 * np.pi) + 1.5, 2.5 * np.sin(76 / 180 * np.pi) + 1)}

        self.position_of_nodes_ordi = {'CpuG': (3, 0.2),
                                       'CpuH': (2.5 * np.cos(36 / 180 * np.pi), 2.5 * np.sin(36 / 180 * np.pi)),
                                       'NPH': (2.5 * np.cos(76 / 180 * np.pi), 2.5 * np.sin(76 / 180 * np.pi)),
                                       'RamH': (3 * np.cos(116 / 180 * np.pi), 2.5 * np.sin(116 / 180 * np.pi)),
                                       'NetOut': (3 * np.cos(208 / 180 * np.pi), -0.5),
                                       'NCM': (2.5 * np.cos(248 / 180 * np.pi), -0.8),
                                       'NPP': (-1, 0.8),
                                       'CpuP': (1, 0.5),
                                       'DiskW': (1.2, -0.8),
                                       'NetIn': (-3, 0.8),
                                       'code_200_206': (
                                           2.5 * np.cos(76 / 180 * np.pi) + 1.5, 2.5 * np.sin(76 / 180 * np.pi) + 1)}

    def _draw_summary_causal_graph(self, graph, data_columns, graph_path, position_of_nodes):
        num_nodes = graph.shape[0]
        lag_max = graph.shape[-1]

        # return edge list
        oriented_edge_lists = []
        unoriented_edge_lists = []
        for lag in range(lag_max):
            if lag == 0:
                for n in range(num_nodes):
                    for m in range(n, num_nodes):
                        if graph[n, m, lag] == '-->':
                            oriented_edge_lists.append((data_columns[n], data_columns[m]))
                        if graph[n, m, lag] == '<--':
                            oriented_edge_lists.append((data_columns[m], data_columns[n]))
                        if graph[n, m, lag] == 'o-o' or graph[n, m, lag] == 'x-x':
                            unoriented_edge_lists.append((data_columns[m], data_columns[n]))
            else:
                for n in range(num_nodes):
                    for m in range(num_nodes):
                        if graph[n, m, lag] == '-->':
                            oriented_edge_lists.append((data_columns[n], data_columns[m]))
                        if graph[n, m, lag] == '<--':
                            oriented_edge_lists.append((data_columns[m], data_columns[n]))
                        if graph[n, m, lag] == 'o-o' or graph[n, m, lag] == 'x-x':
                            unoriented_edge_lists.append((data_columns[m], data_columns[n]))

        # remove duplicates in the unoriented edge list
        new_unoriented_edge_list = []
        for edge in unoriented_edge_lists:
            if (edge[0], edge[1]) not in new_unoriented_edge_list and (
                    edge[1], edge[0]) not in new_unoriented_edge_list and (
            edge[0], edge[1]) not in oriented_edge_lists and (
                    edge[1], edge[0]) not in oriented_edge_lists:
                new_unoriented_edge_list.append(edge)
        unoriented_edge_lists = new_unoriented_edge_list

        # true_pred_directed_edges = [edge for edge in oriented_edge_lists if edge in ground_truth]
        # faut_pred_directed_edges = [edge for edge in oriented_edge_lists if edge not in true_pred_directed_edges]
        #
        #
        # true_pred_undirected_edges = [edge for edge in unoriented_edge_lists if (edge[0], edge[1]) in ground_truth or (edge[1], edge[0]) in ground_truth]
        # faut_pred_undirected_edges = [edge for edge in unoriented_edge_lists if edge not in true_pred_undirected_edges]

        G = nx.DiGraph()
        G.add_nodes_from(data_columns)
        nx.draw_networkx_nodes(G, pos=position_of_nodes, node_size=500)
        nx.draw_networkx_labels(G, pos=position_of_nodes)
        # nx.draw_networkx_edges(G, pos=position_of_nodes, edgelist=true_pred_directed_edges, edge_color='r', arrows=True)
        # nx.draw_networkx_edges(G, pos=position_of_nodes, edgelist=faut_pred_directed_edges, arrows=True)
        # nx.draw_networkx_edges(G, pos=position_of_nodes, edgelist=true_pred_undirected_edges, edge_color='g', arrows=False)
        # nx.draw_networkx_edges(G, pos=position_of_nodes, edgelist=faut_pred_undirected_edges, arrows=False)
        nx.draw_networkx_edges(G, pos=position_of_nodes, edgelist=oriented_edge_lists, arrows=True)
        nx.draw_networkx_edges(G, pos=position_of_nodes, edgelist=unoriented_edge_lists, arrows=False)
        # plt.title(graph_path.split('/')[-1].split('.')[0])
        # plt.savefig(graph_path)
        # plt.show()

    def generate_causal_graph(self, data, tau_max=3, pc_alpha=0.05,
                              not_continuous_columns=[], length_of_data=50, sig_samples=20):

        self.tau_max = tau_max
        self.pc_alpha = pc_alpha
        self.not_continuous_columns = not_continuous_columns
        self.sig_samples = sig_samples

        # No NaN in data
        if np.isnan(data.values).sum() != 0:
            raise ValueError("NaNs in the data.")
        # No constant column
        if len(np.where(np.std(data.values, axis=1)==0)[0]) != 0:
            raise ValueError("Constant in the data.")


        if length_of_data < 0 or length_of_data > data.shape[0]:
            print('All data are used to establish the graph.')
        else:
            data = data.iloc[:length_of_data]

        graphs = pcgce_mixed(data=data, pc_alpha=self.pc_alpha, tau_max=self.tau_max,
                             verbose=False, not_continuous_columns=self.not_continuous_columns,
                             sig_samples=self.sig_samples)

        graph_result = GraphResult(causal_links=list(graphs.ghat.edges))

        return graph_result

    def plot_save_causal_graph(self, res, data):

        data_columns = data.columns
        graph = res['graph']

        if 'code200' in data.columns and 'code206' in data.columns:
            self._draw_summary_causal_graph(graph=graph, data_columns=data_columns, graph_path=self.save_graph_path,
                                            position_of_nodes=self.position_of_nodes)
        if 'code_200_206' in data.columns:
            self._draw_summary_causal_graph(graph=graph, data_columns=data_columns, graph_path=self.save_graph_path,
                                            position_of_nodes=self.position_of_nodes_ordi)


if __name__ == '__main__':
    CDM = CausalDiscoveryFromMixedData()
    # CDM.plot_save_causal_graph(res_file_path='results_09_06/res_strat_1_ordi_after.npy', save_graph_path='graph/strat_1_ordi_after.png', data_file_path='source_data/strat_1_ordi_after.csv')
    # CDM.plot_save_causal_graph(res_file_path='results_09_06/res_strat_1_bina_middle.npy', save_graph_path='graph/strat_1_bina_middle.png', data_file_path='source_data/strat_1_bina_middle.csv')

    # CDM.generate_causal_graph(pc_alpha=0.2, data_file_path='source_data/strat_1_ordi_after.csv', save_res_path='test_strat_1_ordi_after.npy')
    # CDM.plot_save_causal_graph(res_file_path='test_strat_1_ordi_after.npy', save_graph_path='test_strat_1_ordi_after.png', data_file_path='source_data/strat_1_ordi_after.csv')

