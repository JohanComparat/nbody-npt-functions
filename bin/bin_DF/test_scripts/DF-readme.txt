

They are written in  unformatted fortran  format. 

I enclose a  demo code on how to access the data. 

The density fields are given in units of mean density = Omega_M * rhoc_crit(z). SO  the mean value of the cells should be 1. 

The velocity fields are in km/s.  To degrade the mesh to coarser  grids, you need to weight by the density in each cell. It is basically  the same as computing the center of mass velocity of  each  block of cells  you collapse to degrade the resolution. 