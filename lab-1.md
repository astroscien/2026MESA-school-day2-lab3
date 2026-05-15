## A Short Introduction to How GYRE Finds Oscillation Modes

GYRE solves the stellar pulsation equations as a boundary-value eigenvalue problem.  
The unknowns are the oscillation variables inside the star, and the oscillation frequency is the eigenvalue. Only special frequencies satisfy the physical boundary conditions at both the center and the surface. Those special frequencies are the stellar modes. 

GYRE uses a method called **Magnus Multiple Shooting**. Instead of solving the whole star in one step, it divides the star into many small radial intervals. In each interval, it computes how the oscillation solution moves from one grid point to the next, and then it requires all neighboring pieces to match smoothly.

This matching problem can be written in matrix form as

$$
S u = 0,
$$

where $u$ contains the oscillation variables on the shooting grid, and $S$ is the global matrix that enforces the boundary conditions and the matching conditions between adjacent intervals. A non-trivial solution exists only when

$$
\det(S) = 0.
$$

GYRE searches for frequencies at which the matrix $S$ becomes singular. Those frequencies are the eigenfrequencies of the stellar model.

Once a root is found, GYRE reconstructs the corresponding eigenfunction, that is, the radial shape of the mode inside the star. This useful tool lets us compare the mode frequencies predicted by a stellar model with the frequencies observed in real stars.

## Reference
- [Townsend, R. H. D., & Teitler, S. A. 2013, *MNRAS*, 435, 3406](https://ui.adsabs.harvard.edu/abs/2013MNRAS.435.3406T/abstract)
