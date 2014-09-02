from partition import BeckeVoronoi
from quadrature import gauss_legendre_quad,euler_maclaurin_quad
import math
import numpy as np
import pdb
# import matplotplib.pyplot as plt

class MolecularQuadrature:

    def __init__(self,atoms):
        self.atoms = atoms
        self.weights = []
        self.xyz = []
        self.weight_fns = [BeckeVoronoi(atoms,i) for i in range(0,len(atoms))]

        # Set up radial points
        h_becke_rad = 0.35
        r_weights,r_roots = euler_maclaurin_quad(0,1,75)

        # Our integrands are assumed to vanish at r = 0 and r = inf
        r_weights = r_weights[1:-2]
        r_roots_transf = r_roots[1:-2]

        # Change of variables to improve quadrature accuracy from Murray/Handy/Laming
        m = 2
        r_of_q = lambda q: h_becke_rad*pow(q,m)*pow(1 - q,-m)
        jacob_of_q = lambda q: h_becke_rad*m*pow(q,m-1)*pow(1 - q,(-m)-1)
        jacobian_vals = map(jacob_of_q,r_roots_transf)
        r_roots = map(r_of_q,r_roots_transf)

        # Angular points - gauss-legendre theta equally spaced phi
        n_theta_points = 11
        theta_weights,theta_roots_trans = gauss_legendre_quad(n_theta_points)
        theta_roots = [math.acos(val) for val in theta_roots_trans]

        phi_weights,phi_roots = euler_maclaurin_quad(0,2*math.pi,2*n_theta_points) 
        phi_weights[0] *= 2
        phi_weights[-1] *= 2

        # Now translate to absolute cartesians
        for atom in atoms:
            for r_weight,r_root,jacobian_val in zip(r_weights,r_roots,jacobian_vals):
                for theta_weight,theta_root in zip(theta_weights,theta_roots):
                    for phi_weight,phi_root in zip(phi_weights,phi_roots):
                        weight = r_weight*theta_weight*phi_weight*pow(r_root,2)*jacobian_val
                        self.weights.append(weight)
                        x = r_root*math.sin(theta_root)*math.cos(phi_root)
                        y = r_root*math.sin(theta_root)*math.sin(phi_root)
                        z = r_root*math.cos(theta_root)
                        self.xyz.append(np.array((x,y,z)))

    def integrate(self,integrand):
        integral = 0
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.scatter(*zip(*self.xyz))
        # plt.show()
        for atom_idx in range(0,len(self.atoms)):
            for weight,point in zip(self.weights,self.xyz):
                abs_xyz = point + self.atoms[atom_idx]
                weight_fns_norm_fac = sum(weight_fn.evaluate(abs_xyz) for weight_fn in self.weight_fns)
                #(self.weight_fns[atom_idx].evaluate(abs_xyz)/float(weight_fns_norm_fac))
                integral += weight*integrand(abs_xyz)
        return integral