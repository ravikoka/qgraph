import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from AndersonGraph import AndersonGraph


def psi_to_ipr(psi):
    '''
    Convert the wavefunction in the site basis into the Inverse Participation Ratio (IPR).
    '''

    return np.sum(np.abs(psi)**4)


def plot_ipr_evolution(anderson_graph, t_max, nt):
    '''
    Plot IPR of the wavefunction as function of time. The time domain runs from 0 to t_max in nt steps. 

    Args
        anderson_graph (AndersonGraph): anderson graph to time evolve the wavefunction on. 
        t_max (float): the final time. 
        nt (int): number of points in the time domain.  
    '''

    history = anderson_graph.solve(t_max=t_max, nt=nt)
    times = np.linspace(0, t_max, nt)

    ipr_history = [psi_to_ipr(psi) for psi in history]
    
    plt.title('Inverse Participation Ratio of the Wave Function over Time')
    plt.xlabel('time')
    plt.ylabel('Inverse Participation Ratio (IPR)')
    plt.plot(times, ipr_history)


def plot_ipr_vs_W(graph, time, site_num, t_hop, W_max, num_trials, graph_name):
    '''
    Plot IPR, a measure of wavefuntion localization, for various values of W, the disorder parameter. 
    Each IPR value is averaged over multiple trials, and the time is held fixed.
    
    Args
        graph (networkx graph): the underlying graph to time evolve on. 
        time (float): time at which IPR is calculated. 
        site_num (int): a value between 0 and the number of nodes in graph. Sets the intial site of the wavefunction. 
            The intial wavefunction is set to have probability density of 1 at the site indexed by site_num, and 0 elsewhere. 
        t_hop (float): hopping parameter. 
        W_max (float): the maximum value of W. The W domain runs between 0 and W_max in steps of 1. 
        num_trials (int): number of trials to average over. 
        graph_name (str): name of the graph. Used in the plot title.
    '''

    W_vals = np.linspace(0, W_max, W_max)
    n = graph.number_of_nodes()
    psi_0 = np.zeros(n)
    #psi_0[np.random.randint(0, n)] = 1
    psi_0[site_num] = 1

    ipr_averages = []
    for W in W_vals:

        ipr_vals = []
        for _ in range(num_trials):

            anderson_graph = AndersonGraph(graph=graph, psi_0=psi_0, eps_range=[-W, W], t_hop=t_hop) 
            psi_t = anderson_graph.psi_at_t(time)
            ipr_vals.append(psi_to_ipr(psi_t))

        ipr_averages.append(np.mean(ipr_vals))
    
    plt.plot(W_vals, ipr_averages)
    plt.xlabel('W')
    plt.ylabel('Average IPR')
    plt.title(f'IPR of Wave Function vs. W for the {graph_name}')


def plot_ipr_vs_p(num_pts, num_sites, time, W_vals, psi_0, t_hop, num_trials):
    '''
    Plot multiple IPR vs. p curves, each for a value of W.  
    Here the underlying graph is assumed to be a random_graph = nx.erdos_renyi_graph(n, p). 
    Each IPR value is averaged over multiple trials, and the time is held fixed.    

    Args
        num_pts (int): number of points in the p domain. The p domain runs from 0 to 1 in num_pts steps. 
        num_sites (int): number of sites in random graph
        time (float): time at which IPR is calculated. 
        W_vals (array-like[float]): values of W for which an IPR vs. p curve is calculated.
        psi_0 (1D nd array): initial wavefunction in the lattice site basis. 
        t_hop (float): hopping parameter.
        num_trials (int): number of trials to average over. 
    '''

    p_vals = np.linspace(0, 1, num_pts)
    ipr_curves_for_each_W = np.zeros(shape=(len(W_vals), len(p_vals)))
    
    for i, W in enumerate(W_vals):
        
        ipr_curve = []
        for p in p_vals:
            
            ipr_vals = []
            for _ in range(num_trials):
                graph = nx.erdos_renyi_graph(num_sites, p)
                anderson_graph = AndersonGraph(graph=graph, psi_0=psi_0, eps_range=[-W, W], t_hop=t_hop)
                psi_t = anderson_graph.psi_at_t(time)
                ipr_vals.append(psi_to_ipr(psi_t))
            
            ipr_curve.append(np.mean(ipr_vals))
        ipr_curves_for_each_W[i] = ipr_curve
    
    for ipr_curve in ipr_curves_for_each_W:
        plt.plot(p_vals, ipr_curve)
    
    plt.legend(W_vals)
    
    plt.xlabel('$p$')
    plt.ylabel('IPR average')
    plt.title('Inverse Participation Ratio vs. $p$')