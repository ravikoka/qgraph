import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import expm


class AndersonGraph:
    '''
    Construct a graph to time-evolve the Anderson hamiltonian on. 
    We represent all operators and state vectors in the lattice site basis. 

    Attributes
        graph (): underlying graph
        psi_0 (1D nd array):
        eps_range (array-like):
        t_hop (float):
        num_sites (int):
        binding (2D nd array):
        pos ():
    '''

    def __init__(self, graph, psi_0, eps_range, t_hop):
        self.graph = graph
        self.psi_0 = psi_0
        self.eps_range = eps_range
        self.t_hop = t_hop

        self.num_sites = self.graph.number_of_nodes()
        self.binding = np.diagflat(np.random.uniform(*self.eps_range, size=self.num_sites))
        self.pos = nx.spring_layout(self.graph)


    def _hamiltonian(self):
        '''
        Construct the hamiltonian for the Anderson tight-binding model, a matrix representation in the occupancy site basis. 
        but with a random graph instead with probability p
        '''
        adjacency = nx.to_numpy_array(self.graph)
        hopping = -self.t_hop * adjacency

        return self.binding + hopping
    
    
    def _time_evolution_operator(self, time):
        '''
        Calculate the unitary time evolution operator for the given hamiltonian.

        Args
            time (float): time when time evolution operator is calculated

        Returns
            U(t) (ndarray of size num_sites x num_sites) 
        '''

        return expm(-1j * self._hamiltonian() * time)
    
    
    def psi_at_t(self, time):

        return self._time_evolution_operator(time) @ self.psi_0
    

    def solve(self, t, nt): #t_steps):
        '''
        Calculate psi(t). 

        Args
            t (array-like): time range of form (t_initial, t_final)
            nt (int): number of time steps
        '''        

        times = np.linspace(0, t, nt)
        
        history = []

        for time in times:
            psi_t = self._time_evolution_operator(time) @ self.psi_0
            history.append(psi_t) 
        
        return history
    

    def plot_graph(self):
        
        return nx.draw(self.graph)


    def plot_density(self, t):# axisstabilized = False):
        '''
        Plot the probability density, aka |psi(t)|^2

        Args:
            t (float): time to plot it
        '''
        
        #fig = plt.figure(figsize = (12,10))
        #psi_t = self._time_evolution(t) @ self.psi_0
        psi_t = self.psi_at_t(t)
        density = np.real(np.multiply(psi_t.conj(), psi_t))
        
        #plt.title("Wave function probability density at time " + str(t) + "\n p = " + str(self.p)) 
      
        nx.draw(self.graph, self.pos, node_color = density, cmap = plt.cm.cividis)
        
        plt.show()