import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from matplotlib.animation import FuncAnimation
from classes.AndersonGraph import AndersonGraph

plt.style.use('dark_background')


def animate_lattice_pdf(anderson_lattice, t_max, W):
    
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
    path = f'plots/lattice_animation_size_{n}_tmax_{t_max}_W_{W}.gif'
    ani.save(path, writer='ffmpeg', fps=5)
    
    plt.show()
    

def animate_random_graph_pdf(anderson_random_graph, t_max, W):
    '''
    Animate
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

    path = f'plots/random_graph_animation_n_{anderson_random_graph.num_sites}_tmax_{t_max}_W_{W}_p_{anderson_random_graph.p}.gif'
    ani.save(path, writer='ffmpeg', fps=5)
    plt.show()