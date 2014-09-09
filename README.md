This repository implements numerical integration of arbitrary functions over an
atom-centered molecular grid using Becke partitioning. Becke’s crucial discovery
allows for the subdivision of space into individual atom-centered grids so that
nuclear cusps can be accurately integrated, while simultaneously respecting
normalization over all space. This is accomplished by using a weighting function
for each atom that decays rapidly in the interatomic regions, thus dividing
space into Voronoi cells in a smooth fashion. You can see the cool way that the
weight functions for two hydrogen molecules partition the grid points between
them in the image below, where the red points have non-negligible weight values
with respect to the lower hydrogen and the blue points have non-negligible
weight values with respect to the upper hydrogen (this is with an Euler-Maclaurin radial grid and a product angular grid):

[!H2 atom-centered grid with Becke Partitioning](images/becke_partitioning.png)
