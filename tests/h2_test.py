import unittest
from molecular_quadrature import MolecularQuadrature
from slater_test_f import slater_1s
import numpy as np
import pdb
import math
# import cProfile
# import pstats

class TestIntegrateRhoH2(unittest.TestCase):

    def test_integrate_rho_h1(self):
        pass
    	# density = lambda r: pow(slater_1s(1,np.array((0,0,0)))(r),2)
        # h_becke_rad = 0.35
        # r_quad = rq.euler_maclaurin_radial_quad(189,h_becke_rad)
        # atoms = map(np.array,[(0,0,0)])
        # quad = MolecularQuadrature(atoms)
        # n = quad.integrate(density)
        # print (1-n)
        # self.assertLess(abs(1-n),1e-14)

    def test_integrate_rho_h2(self):
        pass
        # Moar moar moar digits
        # S12 = 0.917267406386372645191084997515436658741779275938969500240635
        # N = 1/math.sqrt(2*(1+S12))
    	# density =lambda r: 2*pow(N*(slater_1s(1,np.array((0,0,0)))(r) + slater_1s(1,np.array((0,0,.74)))(r)),2)
        # atoms = map(np.array,[(0,0,0),(0,0,.74)])
        # quad = MolecularQuadrature(atoms)
        # n = quad.integrate(density)
        # self.assertLess(abs(2-n),1e-14)

if __name__ == '__main__':
    unittest.main()
