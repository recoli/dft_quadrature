from partition import BeckeVoronoi
from quadrature import gauss_legendre_quad,euler_maclaurin_quad
import math
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D
import radial_quadrature as rq

class MolecularQuadrature:

    def __init__(self,atoms,r_weights,r_roots,theta_weights,theta_roots,phi_weights,phi_roots):
        self.atoms = atoms
        self.weights = []
        self.xyz = []
        self.weight_fns = [BeckeVoronoi(atoms,i) for i in range(0,len(atoms))]

        # Convert into atom-relative cartesian coordinates
        for r_weight,r_root in zip(r_weights,r_roots):
            for theta_weight,theta_root in zip(theta_weights,theta_roots):
                for phi_weight,phi_root in zip(phi_weights,phi_roots):
                    weight = r_weight*theta_weight*phi_weight
                    self.weights.append(weight)
                    x = r_root*math.sin(theta_root)*math.cos(phi_root)
                    y = r_root*math.sin(theta_root)*math.sin(phi_root)
                    z = r_root*math.cos(theta_root)
                    self.xyz.append(np.array((x,y,z)))

    def integrate(self,integrand):
        integral = 0
        abs_xyzs = []
        for atom_idx in range(0,len(self.atoms)):
            for quad_weight,point in zip(self.weights,self.xyz):
                abs_xyz = point + self.atoms[atom_idx]
                abs_xyzs.append(abs_xyz)
                becke_weight_norm_fac = sum(weight_fn.evaluate(abs_xyz) for weight_fn in self.weight_fns)
                becke_weight = self.weight_fns[atom_idx].evaluate(abs_xyz)/float(becke_weight_norm_fac)
                integral += quad_weight*becke_weight*integrand(abs_xyz)
        return integral