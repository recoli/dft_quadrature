import unittest
from partition import BeckeVoronoi
import numpy as np

class BeckeVoronoiTest(unittest.TestCase):

    def test_becke_voronoi(self):

        atoms = [np.array((0,0,0)),np.array((0,0,.74))]
        bv_0 = BeckeVoronoi(atoms,0)
        self.assertEqual(bv_0.evaluate(np.array((0,0,-1))),1)
        self.assertEqual(bv_0.evaluate(np.array((0,0,1))),0)

        bv_1 = BeckeVoronoi(atoms,1)
        self.assertEqual(bv_1.evaluate(np.array((0,0,-1))),0)
        self.assertEqual(bv_1.evaluate(np.array((0,0,1))),1)