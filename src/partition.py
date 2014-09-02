import numpy as np
import operator
import math

def p(mu):
    return (float(3)/2)*mu - (float(1)/2)*pow(mu,3)

def make_cutoff_fn(RA,RB):
    def new_cutoff_fn(r):
        r0 = np.linalg.norm(r - RA)
        r1 = np.linalg.norm(r - RB)
        R01 = np.linalg.norm(RA - RB)
        return 0.5*(1 - p(p(p((r0 - r1)/R01))))
    return new_cutoff_fn

class BeckeVoronoi:

    def __init__(self,atoms,atom_idx):

        self.cutoff_fns = []
        this_atom = atoms[atom_idx]
        for other_atom_idx in range(0,len(atoms)):
            if(other_atom_idx == atom_idx):
                continue
            self.cutoff_fns.append(make_cutoff_fn(this_atom,atoms[other_atom_idx]))

    def evaluate(self,r):
        return reduce(operator.mul,(fn(r) for fn in self.cutoff_fns),1.0)