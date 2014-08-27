import unittest
from molecular_quadrature import MolecularQuadrature
from slater_test_f import slater_1s

class TestIntegrateRhoH2(unittest.TestCase):

    def density(r):
        return 2*pow(slater_1s(1,0)(r) + slater_1s(1,.74)(r),2)


    def test_integrate_rho_h2(self):
        atoms = [(0,0,0),(0,0,.74)]
        #quad = MolecularQuadrature(atoms)
        #n = quad.integrate(self.density)
        #self.assertEqual(n,2)

if __name__ == '__main__':
    unittest.main()
