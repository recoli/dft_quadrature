import unittest
from molecular_quadrature import MolecularQuadrature
from slater_test_f import slater_1s
import numpy as np
import pdb
import math

class TestIntegrateRhoH2(unittest.TestCase):

    def test_integrate_rho_h1(self):
    	density =lambda r: 0.5*pow(slater_1s(1,np.array((0,0,0)))(r),2)
        atoms = map(np.array,[(0,0,0)])
        quad = MolecularQuadrature(atoms)
        n = quad.integrate(density)
        self.assertEqual(n,1)

    def test_integrate_rho_h2(self):
    	density =lambda r: pow(slater_1s(1,np.array((0,0,0)))(r) + slater_1s(1,np.array((0,0,.74)))(r),2)
        atoms = map(np.array,[(0,0,0),(0,0,.74)])
        quad = MolecularQuadrature(atoms)
        n = quad.integrate(density)
        # self.assertEqual(n,2)

if __name__ == '__main__':
    unittest.main()
