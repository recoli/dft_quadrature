def euler_maclaurin_quad(start,end,n_points):
    roots = [start + i*float(end - start)/n_points for i in range(0,n_points+1)]
    weights = [(end-start)*float(1)/n_points]*n_points
    weights[0] = weights[0]/2
    weights[-1] = weights[-1]/2
    return weights,roots
