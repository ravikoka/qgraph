import jax.numpy as jnp
from jax import random
import networkx as nx
import matplotlib.pyplot as plt
from jax.scipy.linalg import expm

class JAXAndersonGraph:
    '''
    Construct a graph to time-evolve the Anderson hamiltonian on. 
    We represent all operators and state vectors in the lattice site basis. 

    Attributes
        graph (periodic networkx graph): underlying graph, representing lattice sites we simulate on.
        psi_0 (1D nd array): initial wavefunction in the site basis.
        eps_range (array-like): range to draw random values of epsilon from to create a random diagonal on the hamiltonian.
        t_hop (float): hopping parameter.
        num_sites (int): number of lattice sites.
        binding (2D nd array): binding term of the hamiltonian (ie. not the hoppoing term).
        pos (dict): dict with nodes of self.graph as keys and positions as values.
    '''

    def __init__(self, graph, psi_0, eps_range, t_hop=1, seed=0):
        self.graph = graph
        self.psi_0 = psi_0

        self.num_sites = self.graph.number_of_nodes()
        
        # Checks if number of nodes on graph matches number of nodes for wave function
        if self.num_sites != len(self.psi_0):
            raise ValueError("Number of nodes on the graph does not match the length of the wave function")

        self.eps_range = eps_range
        self.t_hop = t_hop

        self.binding = jnp.diagflat(random.uniform(key=random.PRNGKey(seed), minval=self.eps_range[0], maxval=self.eps_range[1], shape=(self.num_sites,)))
        self.pos = nx.spring_layout(self.graph)

        self.hamiltonian = self.binding - self.t_hop * self._get_adjacency()


    def _get_adjacency(self): 
        '''
        Construct the adjacency matrix for the graph, a matrix representation in the occupancy site basis. 

        Returns
            adjacency (2d ndarray): matrix representation of the graph in the occupancy basis.
        '''
        adjacency = nx.to_numpy_array(self.graph)
        adjacency[-1, 0] = 1
        adjacency[0, -1] = 1

        return jnp.asarray(adjacency)
        
    def _time_evolution_operator(self, time):
        '''
        Calculate the unitary time evolution operator for the given hamiltonian.

        Args
            time (float): time when time evolution operator is calculated

        Returns
            U(t) (ndarray of size num_sites x num_sites): unitary time evolution operator. 
        '''

        return expm(-1j * self.hamiltonian * time)
    
    
    def psi_at_t(self, time):
        '''
        Returns the wavefucntion, in the lattice site basis, at a given time t. 

        Args
            time (float): time at which the wavefunction is calculated.
        
        Returns
            psi (1D ndarray of size num_sites): wavefunction at given time.
        '''
        return jnp.dot(self._time_evolution_operator(time), self.psi_0).block_until_ready()
        
    

    def simulate(self, t_max, nt): #t_steps):
        '''
        Calculate psi(t) for a sequence of times. 

        Args
            t_max (float): final time.
            nt (int): number of time steps.

        Returns
            history (ndarray of size nt x num_sites): wavefunctions as a function of time from 0 to t_max.
        '''        

        times = jnp.linspace(0, t_max, nt)
        history = []

        for time in times:
            #self._time_evolution_operator(time) @ self.psi_0
            history.append(self.psi_at_t(time))
        return history
    

    def plot_graph(self):
        '''
        Draw graph attribute, using the networkx draw method.
        '''
        
        return nx.draw(self.graph)


    def plot_density(self, t, node_size=10, line_width=1, layout=None):# axisstabilized = False):
        '''
        Plot the probability density, aka |psi(t)|^2

        Args
            t (float): time when the probability density is plotted
        '''
        
        #fig = plt.figure(figsize = (12,10))
        #psi_t = self._time_evolution(t) @ self.psi_0
        plt.style.use('dark_background')
        psi_t = self.psi_at_t(t)
        density = jnp.real(jnp.multiply(psi_t.conj(), psi_t))
        
        #plt.title("Wave function probability density at time " + str(t) + "\n p = " + str(self.p)) 
        if layout == None:
            self.pos=nx.spring_layout(self.graph)
        else:
            self.pos = layout
        nx.draw(self.graph, self.pos, node_color = density, cmap = plt.cm.cividis, node_size = node_size, width = line_width)
        
        plt.show()