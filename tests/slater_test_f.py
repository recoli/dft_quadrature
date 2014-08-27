import math
import numpy as np

def slater_1s(zeta,RA):
    def slater_functor(r):
        norm = np.linalg.norm(r - RA)
        return math.sqrt(pow(float(zeta),3)/math.pi)*math.exp(-zeta*norm)
    return slater_functor
