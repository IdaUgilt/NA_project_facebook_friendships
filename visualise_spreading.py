# Imports
import random
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from sir_spreading import SIR_classical
from sir_spreading import SIR_threshold
from sir_spreading import SIR_cascade

import warnings
warnings.filterwarnings('ignore')

# Get starting nodes
def get_starting_nodes(G, metric, n):
    '''
    Return starting nodes, which should be used as input to spreading proces.

    Args:
        G (nx.graph): Graph object to choose starting nodes from
        metric (str): Metrics used to pick starting nodes
        n (int): Number of starting nodes to return

    Returns:
        starting_nodes
    '''
    if metric == 'random':
        starting_nodes_idx = random.sample(sorted(G.nodes()), n)

    else:

        if metric == 'degree':
            node_dict = nx.degree_centrality(G)
        
        elif metric  == 'closeness':
            node_dict = nx.closeness_centrality(G)

        elif metric == 'betweenness':
            node_dict = nx.betweenness_centrality(G)
        
        elif metric  == 'eigenvector':
            node_dict = nx.eigenvector_centrality(G)
        
        elif metric == 'katz':
            node_dict = nx.katz_centrality(G)
        
        elif metric == 'harmonic':
            node_dict = nx.harmonic_centrality(G)
        
        else:
            print('metric is not supported!')


        # Sort the dictionary items in descending order based on centrality values
        sorted_nodes = sorted(node_dict.items(), key=lambda x: x[1], reverse=True)
            
        top_n_nodes = sorted_nodes[:n]
            
        # Extract only the node IDs from the sorted list
        starting_nodes_idx = [node_id for node_id, _ in top_n_nodes]

    return starting_nodes_idx

def get_spreading_data(G, beta, gamma, starting_nodes, trigger, metrics, kappa=None):
    '''
    Get data from the spreading proces.

    Args:
        trigger (str): Trigger mechanism.
        metrics (list): List of metrics.

    Returns:
        df (DataFram): Data frame with results
    '''

    # Dictionary connecting trigger mechanism and spreading model
    d_metrics = {'classical': SIR_classical, 'threshold': SIR_threshold, 'cascade': SIR_cascade}

    # Get data for each specified metric
    if trigger == 'classical' or trigger == 'cascade':
        d = {f'data_{metric}': d_metrics[trigger](G, beta, gamma, starting_nodes[metric]) for metric in metrics}
    elif trigger == 'threshold':
        d = {f'data_{metric}': d_metrics[trigger](G, kappa, beta, gamma, starting_nodes[metric]) for metric in metrics}
    else:
        print('Trigger type is not supported!')
        return

    return d

def visualise_spreading(data, custom_palette, ax):

    metrics = list(data.keys())

    # Create graph
    for i, metric in enumerate(metrics):
        sns.lineplot(data=data[metric], x='t', y='IR', label=f'{metric}', color=custom_palette[i], ax=ax)


