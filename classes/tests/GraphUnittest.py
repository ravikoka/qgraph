# Created by Anish, 12/5/2023
import unittest
import numpy as np
import networkx as nx
import scipy.linalg as la
import matplotlib.pyplot as plt
from classes.AndersonGraph import AndersonGraph

class GraphUnittestClass(unittest.TestCase):
    def setUp(self) -> None:
        self.G = nx.StarGraph(5)
        self.A = nx.to_numpy_array(self.G)
        self.psi_0 = np.zeros(5)
        self.psi_0[5//2] = 1
        self.eps = [10, 20, 30, 40, 50, 60]

        self.t_hop = 1

        self.AndersonGraph = AndersonGraph(self.G, self.psi_0, self.t_hop, eps=self.eps, store_history=True)
    
    def test_time_evolution(self):
        self.assertEqual(self.AndersonGraph._time_evolution(0), self.psi_0)

    def test_propogate(self):
        pass
    def test_simulate(self):
        pass

if __name__ == '__main__':
    unittest.main()