# MNPDynamics

This MATLAB package provides a framework to simulate nanoparticle dynamics obeying either Brownian or Néel rotation.

Two discretization methods are contained in this tool box: A spherical harmonic decomposition as well as a Finite Volume method relying on a triagonal grid on the unit sphere. This grid should be as uniform as possible and must have the property that each triangle's circumcenter is contained in its interior. A different toolbox that produces such grids is cited in the source code of the file FV_matrix.m.

If Néel rotation with a time-dependent easy axis is of interest, the files contained in the folder matlab_dynamic_n should be used. For stationary easy axes, the methods in the directory matlab are superior.

A simple example is presented in the corresponding directory. The documentation of each function can be found in the form of comments directly in its source code.

A research article analyzing the performance and presenting applications of this method is in development.
