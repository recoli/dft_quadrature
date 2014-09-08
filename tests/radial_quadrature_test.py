import unittest
import radial_quadrature as rq
from slater_test_f import slater_1s
import numpy as np
import math

class RadialQuadratureTest(unittest.TestCase):

    def test_euler_maclaurin_radial_quad(self):
        # Becke used this instead of Bragg-Slater value of 0.25A
        h_becke_rad = 0.35
        weights,roots = rq.euler_maclaurin_radial_quad(1000,h_becke_rad)
        sf = slater_1s(1,np.array((0,0,0)))
        quad_point_fn = lambda i: pow(sf(np.array((0,0,roots[i]))),2)
        integral = 4*math.pi*sum(weights[i] * quad_point_fn(i) for i in range(0,len(weights)))
        self.assertLess(abs(1-integral),1e-14)