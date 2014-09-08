import quadrature as quad
import math

def theta_phi_product_quad(n_theta_points,n_phi_points):        
        theta_weights,theta_roots_trans = quad.gauss_legendre_quad(n_theta_points)
        theta_roots = [math.acos(theta_trans) for theta_trans in theta_roots_trans]
        phi_roots = ((2*math.pi*i)/(n_phi_points-1) for i in range(0,n_phi_points))
        phi_weights =  [2*math.pi/n_phi_points]*n_phi_points
        return theta_weights,theta_roots,phi_weights,phi_roots
