import unittest
import math
from molecular_quadrature import MolecularQuadrature

class TestIntegrateRhoH2(unittest.TestCase):

    def slater_1s(zeta,R):
        def slater_functor(r):
            norm = math.sqrt(sum(pow(r[i] - R[i],2) for i in range(0,len(r))))
            return (pow(zeta,3)/3)*math.exp(-zeta*norm)

    def density(r):
        return 2*pow(slater_functor(1,0) + slater_functor(1,.74),2)


    def test_integrate_rho_h2(self):
        atoms = [(0,0,0),(0,0,.74)]
        quad = MolecularQuadrature(atoms)
        n = quad.integrate(self.density)
        self.assertEqual(n,2)

if __name__ == '__main__':
    unittest.main()
