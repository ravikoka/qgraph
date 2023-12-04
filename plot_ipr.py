import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from classes.AndersonGraph import AndersonGraph


def psi_to_ipr(psi):

    #newpsi = []
    #for x in psi:
    #    newpsi.append((abs(x))**4)
    #return sum(newpsi)
    return np.sum(np.abs(psi)**4)


def plot_ipr_evolution(anderson_graph, t_max, nt):
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
    œ
    Args
        graph (networkx graph):
        time ():
        site_num ():
        t_hop (float):
        num_trials (int):
        graph_name (str):
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
    Plot an IPR vs. p curve for multiple W values. 

    Args
        num_pts (int):
        num_sites (int): number of sites in random graph
        time (float):
        W_vals (array-like):
        psi_0 (1D nd array):
        t_hop (float):
        num_trials (int):
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