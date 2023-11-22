# Imports

import pandas as pd
import networkx as nx
from random import random

def sir_classical_spreading(G, beta = 1, gamma = None, starting_nodes = None):
    '''
    Run a SIR spreading process through a given network using classical triggering determined by parameter beta.
    Nodes stop attempting to infect neighbors after gamma time steps of being infected.
    
    Args:
        G (nx.Graph): Graph object to run spreading process on
        beta (float): probability of infected node infection neighbors
        gamma (int): time steps taken before node in I transitions to R (stops attempting to infect further)
        starting_nodes (int or list of ints): starting node(s) to use for spreading process

    Returns:
        results (pd.DataFrame): Data Frama containing time steps and corresponding fraction of nodes infected.
    '''

    # If no specific starting nodes given, use all nodes in G
    if starting_nodes == None:
        starting_nodes = list(G.nodes)

    # If only one starting node and given as int or str, contain in list
    elif type(starting_nodes) == int or type(starting_nodes) == str:
        starting_nodes = [starting_nodes]

    # Save amount of nodes in variable
    V = len(G.nodes)

    # Keep dictionary of fraction of infected and revomed nodes for each time step
    time_step_infection_rates = {0: [1/V]} # Time step 0 will always only have one infected node (the starting node)

    # Run complete spreading process for each starting node
    for s in starting_nodes:
        
        # Keep sets of nodes in each state
        S = set(G.nodes) - {s} # Initialize S as all nodes except starting node s
        I = {s} # Initialize I as only containing starting node
        R = set() # Initialize R as empty

        # Keep dict of nodes infected at each time step determine transition from I to R
        infected_at = {0: {s}} # Only starting node s is infected at time step 0

        # Time step iterator
        t = 1

        # Continue spreading process as long as I is not empty
        # Each iteration of this loop is the actions taken in time step t
        while I:

            # Keep set of newly infected nodes
            new_infected = set()

            # Iterate over all nodes in S to find new infections
            for u in S:

                # Check if u gets infected using classical trigger logic
                n = len(set(G.neighbors(u)) & I) # Amount of infected neighbors of u
                # print(u, 'infected neighbors:', set(G.neighbors(u)) & I)
                if random() < 1 - (1 - beta) ** n: 
                    new_infected.add(u)

            # Add new infected nodes to I and remove from S
            I = I.union(new_infected)
            S -= new_infected

            # Record infection time of newly infected nodes
            infected_at[t] = new_infected

            # Transition infected nodes from I to R after gamma time steps
            if t - gamma >= 0:
                to_transition = infected_at[t - gamma]
                R = R.union(to_transition)
                I -= to_transition

            # Add fraction of infected or removed nodes to time step record
            time_step_infection_rates.setdefault(t, []).append((len(I) + len(R)) / V)

            # Debug printing
            print(f't: {t}, |S|: {len(S)}, |I|: {len(I)}, |R|: {len(R)}')

            # Increment t
            t += 1

    # Average out values in dict of infection rates for each time step t
    time_step_infection_rates = {t: sum(i_r) / len(i_r) for t, i_r in time_step_infection_rates.items()}

    # Create Data Frame of infection rates for each time step t
    results = pd.DataFrame(list(time_step_infection_rates.items()), columns = ['t', 'IR'])

    return results