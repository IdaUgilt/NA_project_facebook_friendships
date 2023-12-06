# Imports

import pandas as pd
import networkx as nx
from random import random

def SIR_classical(G: nx.Graph, beta: float = .5, gamma: int = 1, starting_nodes = None):
    '''
    Run a SIR spreading process through a given network using classical triggering 
    determined by constant parameter beta.
    Nodes stop attempting to infect neighbors after gamma time steps of being infected.
    
    Args:
        G (nx.Graph): Graph object to run spreading process on
        beta (float): probability of infected node infecting neighbors
        gamma (int): time steps taken before node in I transitions to R (stops attempting to infect further)
        starting_nodes (int or list of ints): starting node(s) to use for spreading process

    Returns:
        results (pd.DataFrame): Data Frama containing time steps and corresponding fraction of nodes infected.
    '''

    # Save amount of nodes in variable
    V = len(G.nodes)

    # If no specific starting nodes given, use all nodes in G
    if starting_nodes == None:
        starting_nodes = list(G.nodes)

    # If only one starting node and given as int or str, contain in list
    elif type(starting_nodes) == int or type(starting_nodes) == str:
        starting_nodes = [starting_nodes]

    # Keep dictionary of shares of nodes in each state for each time step t
    time_step_rates = {0: { # Initialize with only t = 0
        'S': [(V-1)/V], # List of shares of S nodes, one entry for each starting node
        'I': [1/V], # List of shares of I nodes, one entry for each starting node
        'R': [0.0] # List of shares of R nodes, one entry for each starting node
    }}

    # Run complete spreading process for each starting node
    for s in starting_nodes:
        
        # Keep sets of nodes in each state
        S = set(G.nodes) - {s} # Initialize S as all nodes except starting node s
        I = {s} # Initialize I as only containing starting node
        R = set() # Initialize R as empty

        # Keep dict of nodes infected at each time step determine transition from I to R
        infected_at = {0: {s}} # Only starting node s is infected at time step 0

        # Time step iterator
        t = 0

        # Continue spreading process as long as I is not empty
        # Each iteration of this loop is the actions taken in time step t
        while I:

            # Increment t
            t += 1

            # Keep set of newly infected nodes
            new_infected = set()

            # Iterate over all nodes in S to find new infections
            for u in S:

                # Check if u gets infected using classical trigger logic
                n = len(set(G.neighbors(u)) & I) # Amount of infected neighbors of u
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

            # Add fraction of S, I, and R nodes to time step record
            # Set defualt if time step t not reached before
            time_step_rates.setdefault(t, {'S': [], 'I': [], 'R': []})['S'].append(len(S) / V) # Add share of nodes in S
            time_step_rates[t]['I'].append(len(I) / V) # Add share of nodes in I
            time_step_rates[t]['R'].append(len(R) / V) # Add share of nodes in R

    # Average out values in dict of S, I, and R rates for each time step t
    for t, rates in time_step_rates.items():
        time_step_rates[t]['S'] = sum(rates['S']) / len(rates['S'])
        time_step_rates[t]['I'] = sum(rates['I']) / len(rates['I'])
        time_step_rates[t]['R'] = sum(rates['R']) / len(rates['R'])

    # Create Data Frame of S, I, and R rates for each time step t
    results = pd.DataFrame([{'t': t, 
                             'SR': rates['S'], 
                             'IR': rates['I'], 
                             'RR': rates['R'], 
                             'IR+RR': rates['I'] + rates['R']} for t, rates in time_step_rates.items()])

    return results

def SIR_threshold(G: nx.Graph, kappa: int = 1, beta: float = .5, gamma: int = 1, starting_nodes = None):
    '''
    Run a SIR spreading process through a given network using threshold triggering 
    determined by constant parameter kappa and beta.
    Nodes stop attempting to infect neighbors after gamma time steps of being infected.
    
    Args:
        G (nx.Graph): Graph object to run spreading process on
        kappa (int): amount of infected neighbors necessary for eventual spreading
        beta (float): probability of infection once kappa threshold has been passed
        gamma (int): time steps taken before node in I transitions to R (stops attempting to infect further)
        starting_nodes (int or list of ints): starting node(s) to use for spreading process

    Returns:
        results (pd.DataFrame): Data Frama containing time steps and corresponding fraction of nodes infected.
    '''

    # Save amount of nodes in variable
    V = len(G.nodes)

    # If no specific starting nodes given, use all nodes in G
    if starting_nodes == None:
        starting_nodes = list(G.nodes)

    # If only one starting node and given as int or str, contain in list
    elif type(starting_nodes) == int or type(starting_nodes) == str:
        starting_nodes = [starting_nodes]

    # Keep dictionary of shares of nodes in each state for each time step t
    time_step_rates = {0: { # Initialize with only t = 0
        'S': [(V-1)/V], # List of shares of S nodes, one entry for each starting node
        'I': [1/V], # List of shares of I nodes, one entry for each starting node
        'R': [0.0] # List of shares of R nodes, one entry for each starting node
    }}

    # Run complete spreading process for each starting node
    for s in starting_nodes:
        
        # Keep sets of nodes in each state
        S = set(G.nodes) - {s} # Initialize S as all nodes except starting node s
        I = {s} # Initialize I as only containing starting node
        R = set() # Initialize R as empty

        # Keep dict of nodes infected at each time step determine transition from I to R
        infected_at = {0: {s}} # Only starting node s is infected at time step 0

        # Time step iterator
        t = 0

        # Continue spreading process as long as I is not empty
        # Each iteration of this loop is the actions taken in time step t
        while I:

            # Increment t
            t += 1

            # Keep set of newly infected nodes
            new_infected = set()

            # Iterate over all nodes in S to find new infections
            for u in S:
                
                if t == 1:
                    # Check if u gets infected using threshold trigger logic
                    n = len(set(G.neighbors(u)) & I) # Amount of infected neighbors of u
                    if n >= 1 and random() <= beta: # Beta probability of infection if kappa threshold passed
                        new_infected.add(u)
                else:
                    # Check if u gets infected using threshold trigger logic
                    n = len(set(G.neighbors(u)) & I) # Amount of infected neighbors of u
                    if n >= kappa and random() <= beta: # Beta probability of infection if kappa threshold passed
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

            # Add fraction of S, I, and R nodes to time step record
            # Set defualt if time step t not reached before
            time_step_rates.setdefault(t, {'S': [], 'I': [], 'R': []})['S'].append(len(S) / V) # Add share of nodes in S
            time_step_rates[t]['I'].append(len(I) / V) # Add share of nodes in I
            time_step_rates[t]['R'].append(len(R) / V) # Add share of nodes in R

    # Average out values in dict of S, I, and R rates for each time step t
    for t, rates in time_step_rates.items():
        time_step_rates[t]['S'] = sum(rates['S']) / len(rates['S'])
        time_step_rates[t]['I'] = sum(rates['I']) / len(rates['I'])
        time_step_rates[t]['R'] = sum(rates['R']) / len(rates['R'])

    results = pd.DataFrame([{'t': t, 
                             'SR': rates['S'], 
                             'IR': rates['I'], 
                             'RR': rates['R'], 
                             'IR+RR': rates['I'] + rates['R']} for t, rates in time_step_rates.items()])

    return results

def SIR_cascade(G: nx.Graph, beta: float = .5, gamma: int = 1, starting_nodes = None):
    '''
    Run a SIR spreading process through a given network using cascade triggering 
    determined by constant parameter beta. Note beta differs in meaning from classical and
    threshold models, as here it is the fraction of infected neighbors needed to be infected.
    Nodes stop attempting to infect neighbors after gamma time steps of being infected.
    
    Args:
        G (nx.Graph): Graph object to run spreading process on
        beta (float): fraction of neighbors being infected needed to infect a node
        gamma (int): time steps taken before node in I transitions to R (stops attempting to infect further)
        starting_nodes (int or list of ints): starting node(s) to use for spreading process

    Returns:
        results (pd.DataFrame): Data Frama containing time steps and corresponding fraction of nodes infected.
    '''

    # Save amount of nodes in variable
    V = len(G.nodes)

    # If no specific starting nodes given, use all nodes in G
    if starting_nodes == None:
        starting_nodes = list(G.nodes)

    # If only one starting node and given as int or str, contain in list
    elif type(starting_nodes) == int or type(starting_nodes) == str:
        starting_nodes = [starting_nodes]

    # Keep dictionary of shares of nodes in each state for each time step t
    time_step_rates = {0: { # Initialize with only t = 0
        'S': [(V-1)/V], # List of shares of S nodes, one entry for each starting node
        'I': [1/V], # List of shares of I nodes, one entry for each starting node
        'R': [0.0] # List of shares of R nodes, one entry for each starting node
    }}

    # Run complete spreading process for each starting node
    for s in starting_nodes:
        
        # Keep sets of nodes in each state
        S = set(G.nodes) - {s} # Initialize S as all nodes except starting node s
        I = {s} # Initialize I as only containing starting node
        R = set() # Initialize R as empty

        # Keep dict of nodes infected at each time step determine transition from I to R
        infected_at = {0: {s}} # Only starting node s is infected at time step 0

        # Time step iterator
        t = 0

        # Continue spreading process as long as I is not empty
        # Each iteration of this loop is the actions taken in time step t
        while I:

            # Increment t
            t += 1

            # Keep set of newly infected nodes
            new_infected = set()

            # Iterate over all nodes in S to find new infections
            for u in S:

                # Check if u gets infected using cascade trigger logic
                n_total = len(list(G.neighbors(u))) # Amount of total neighbors of u
                n_infected = len(set(G.neighbors(u)) & I) # Amount of infected neighbors of u
                if n_infected / n_total >= beta:
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

            # Add fraction of S, I, and R nodes to time step record
            # Set defualt if time step t not reached before
            time_step_rates.setdefault(t, {'S': [], 'I': [], 'R': []})['S'].append(len(S) / V) # Add share of nodes in S
            time_step_rates[t]['I'].append(len(I) / V) # Add share of nodes in I
            time_step_rates[t]['R'].append(len(R) / V) # Add share of nodes in R

    # Average out values in dict of S, I, and R rates for each time step t
    for t, rates in time_step_rates.items():
        time_step_rates[t]['S'] = sum(rates['S']) / len(rates['S'])
        time_step_rates[t]['I'] = sum(rates['I']) / len(rates['I'])
        time_step_rates[t]['R'] = sum(rates['R']) / len(rates['R'])

    results = pd.DataFrame([{'t': t, 
                             'SR': rates['S'], 
                             'IR': rates['I'], 
                             'RR': rates['R'], 
                             'IR+RR': rates['I'] + rates['R']} for t, rates in time_step_rates.items()])

    return results