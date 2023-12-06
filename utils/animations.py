import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from matplotlib.animation import FuncAnimation
from AndersonGraph import AndersonGraph

plt.style.use('dark_background')


def animate_lattice_pdf(anderson_lattice, t_max):
    '''
    Animate the wavefunction PDF on a periodic square lattice vs. time. 
    To do: make this work with any 'rectangular' periodic lattice. 

    Args
        anderson_lattice (AndersonGraph): the anderson graph to time evolve. anderson_lattice.graph must be a periodic square graph.
            To make a periodic square lattice, use anderson_lattice.graph = nx.grid_graph(dim=[n, n], periodic=True).
        t_max (int): the final time. The time domain proceeds from 0 to t_max in steps of 1.
    '''

    n = int(np.sqrt(anderson_lattice.num_sites))

    fig = plt.figure(figsize = (12,10))
    
    x = np.arange(0, n, 1)
    y = np.arange(0, n, 1)
    X, Y = np.meshgrid(x, y)

    ax = plt.axes(projection='3d')
    ax.set_xlabel('x', labelpad=20)
    ax.set_ylabel('y', labelpad=20)
    ax.set_zlabel('P', labelpad=20)

    def update(frame):
        ax.clear()
        psi_t = anderson_lattice.psi_at_t(frame)
        density = np.real(np.multiply(psi_t.conj(), psi_t))
        density = np.reshape(density, (n, n))
        
        ax.plot_surface(X, Y, density, cmap = plt.cm.cividis)
        ax.set_zlim(0, 1)
        
        ax.set_title(f'Wave function probability density at time {frame}')
    
    ani = FuncAnimation(fig, update, frames=range(0, t_max), interval=100)
    path = f'plots/lattice_animation_size_{n}_tmax_{t_max}_W_{anderson_lattice.eps_range[1]}.gif'
    print(f'saving animation to {path}')
    ani.save(path, writer='ffmpeg', fps=5)
    
    plt.show()
    

def animate_random_graph_pdf(anderson_random_graph, t_max, p):
    '''
    Animate the wavefunction PDF on a random graph vs. time. 
    To do: determine how to get p out of a networkx graph. 

    Args
        anderson_random_graph (AndersonGraph): the anderson graph to time evolve. anderson_random_graph.graph must be a random graph.
            To create a random graph, use random_graph = nx.erdos_renyi_graph(n, p).
        t_max (int): the final time. The time domain proceeds from 0 to t_max in steps of 1.
        p (float, between 0.0 and 1.0): this is p associated with anderson_random_graph.graph, the probability an edge exists.
            Disclaimer: this is a workaround to save the animation with a title that contains p. 
    '''
    
    fig = plt.figure(figsize = (12,10))
    ax = plt.gca()
    
    def update(frame):
        ax.clear()
        psi_t = anderson_random_graph.psi_at_t(frame)
        density = np.real(np.multiply(psi_t.conj(), psi_t))
        ax.set_title("Wave function probability density at time " + str(frame))
        nx.draw(anderson_random_graph.graph, anderson_random_graph.pos, node_color = density, 
                cmap = plt.cm.cividis, ax = ax, edge_color = 'white')

    ani = FuncAnimation(fig, update, frames=range(0, t_max), interval=100)

    path = f'plots/random_graph_animation_n_{anderson_random_graph.num_sites}_tmax_{t_max}_W_{anderson_random_graph.eps_range[1]}_p_{p}.gif'
    print(f'saving animation to {path}')
    ani.save(path, writer='ffmpeg', fps=5)
    plt.show()