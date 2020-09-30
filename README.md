# MNPDynamics

This MATLAB package provides a framework to simulate nanoparticle dynamics obeying either Brownian or Néel rotation.

Two discretization methods are contained in this tool box: A spherical harmonic decomposition as well as a Finite Volume method relying on a triagonal grid on the unit sphere. This grid should be as uniform as possible and must have the property that each triangle's circumcenter is contained in its interior. A different toolbox that produces such grids is cited in the source code of the file FV_matrix.m.

A research article analyzing the performance and presenting applications of this method is in development.
