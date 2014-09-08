import quadrature as quad

def euler_maclaurin_radial_quad(n_points,alpha):
    q_weights,q_roots = quad.euler_maclaurin_quad(0,1,n_points)

    # Our integrands are assumed to vanish at r = 0 and r = inf
    q_weights = q_weights[1:-2]
    q_roots = q_roots[1:-2]

    # Change of variables to improve quadrature accuracy from Murray/Handy/Laming
    m = 2
    r_of_q = lambda q: alpha*pow(q,m)*pow(1 - q,-m)
    jacob_of_q = lambda q: alpha*m*pow(q,m-1)*pow(1 - q,(-m)-1)
    jacobian_vals = map(jacob_of_q,q_roots)
    r_roots = map(r_of_q,q_roots)
    scaled_weights = [q_weights[i]*pow(r_roots[i],2)*jacobian_vals[i] for i in range(0,len(q_roots))]
    return scaled_weights,r_roots