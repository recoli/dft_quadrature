import unittest
from molecular_quadrature import MolecularQuadrature
from slater_test_f import slater_1s
import numpy as np
import math
import radial_quadrature as rq
import angular_quadrature as aq

class TestIntegrateRhoH2(unittest.TestCase):

    def test_integrate_rho_h1(self):
        density = lambda r: pow(slater_1s(1,np.array((0,0,0)))(r),2)
        atoms = map(np.array,[(0,0,0)])
        h_becke_rad = 0.35
        r_weights,r_roots = rq.euler_maclaurin_radial_quad(125,h_becke_rad)
        theta_weights,theta_roots,phi_weights,phi_roots = aq.theta_phi_product_quad(10,2*10)
        quad = MolecularQuadrature(atoms,r_weights,r_roots,theta_weights,theta_roots,phi_weights,phi_roots)
        n = quad.integrate(density)
        self.assertLess(abs(1-n),1e-12)

    def test_integrate_rho_h2(self):
        # Moar moar moar digits
        S12 = 0.917267406386372645191084997515436658741779275938969500240635
        N = 1/math.sqrt(2*(1+S12))
    	density =lambda r: 2*pow(N*(slater_1s(1,np.array((0,0,0)))(r) + slater_1s(1,np.array((0,0,.74)))(r)),2)
        atoms = map(np.array,[(0,0,0),(0,0,.74)])
        h_becke_rad = 0.35
        r_weights,r_roots = rq.euler_maclaurin_radial_quad(100,h_becke_rad)
        theta_weights,theta_roots,phi_weights,phi_roots = aq.theta_phi_product_quad(13,2*13)
        quad = MolecularQuadrature(atoms,r_weights,r_roots,theta_weights,theta_roots,phi_weights,phi_roots)
        n = quad.integrate(density)
        print 2-n
        self.assertLess(abs(2-n),1e-10)

if __name__ == '__main__':
    unittest.main()
