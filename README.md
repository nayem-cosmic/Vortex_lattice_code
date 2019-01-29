# Analysis of Potential Flow around a 3D Hydro-wing by Vortex Lattice Method
The thesis deals with the potential flow problem around the three-dimensional hydrofoil which is constructed by specifying the Laplace equation with the body boundary condition, the boundary condition for the flow disturbance far from the body and the Kutta condition at the trailing edge of the lifting hydrofoil.

In order to solve the problem, the Laplace equation is transformed into an integral equation in terms of a distribution of singular solutions on the boundaries. After satisfying the boundary conditions the integral equation can be written into a matrix form and this matrix is solved by Gaussian Elimination procedure.

The Vortex Lattice Method represents the wing as a surface on which a grid of vortex rings is superimposed. The velocities induced by each ring vortex at a specified control point are calculated using the law of Biot-Savart. A summation is performed for all control points on the wing to produce a set of linear algebraic equations for the ring vortex strengths that satisfy the boundary condition of no flow through the wing. The vortex strengths are related to the wing circulation and the pressure differential between the upper and lower wing surfaces. The pressure differentials are integrated to yield the total forces and moments.

Based on the mathematical theory a computer program is developed in Fortran 95 to analyze the hydrodynamic characteristics of the 3-D hydrofoil and the planar wing. Python’s matplotlib module is used to plot the graphs.
