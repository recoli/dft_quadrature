import angular_quadrature as aq
import unittest
import math
import itertools

class AngularQuadratureTest(unittest.TestCase):

    def test_product_grid(self): 
        # Angular points - gauss-legendre theta with equally spaced phi
        n_theta_points = 9
        n_phi_points = 2*n_theta_points
        theta_weights,theta_roots,phi_weights,phi_roots  = aq.theta_phi_product_quad(n_theta_points,n_phi_points)
        
        # We transform the theta coordinate with a change of variables to gauss legendre qudrature
        integral = sum(p[0]*p[1] for p in itertools.product(theta_weights,phi_weights))
        self.assertLess(abs(integral - 4*math.pi),1e-14)
