## A Short Introduction to How GYRE Finds Oscillation Modes (Proofread by Rich)

GYRE solves the stellar oscillation equations, a set of differential equations and boundary conditions that describe small, periodic perturbations about a star's equilibrium state. Consistent solutions to these equations (known as "modes") can be found only for certain specific choices of the perturbation frequency, and so the frequency takes on the mathematical role of an eigenvalue.

To calculate the frequency eigenvalues (or "eigenfrequencies") of a stellar model (obtained, for instance, from MESA), GYRE sets up a large system of algebraic equations. These equations are derived from finite-difference approximations to the oscillation differential equations, taken between pairs of adjacent spatial grid points. In symbolic form, the algebraic equations can be written as

$$
S u = 0,
$$

where $S$ is a matrix of coefficients and $u$ is a vector of unknowns representing the perturbations at each grid point. Solutions to this equation only exist when

$$
\det(S) = 0.
$$

and so GYRE's task is to search for the frequencies at which the determinant of $S$ vanishes. Once these eigenfrequencies are found, GYRE reconstructs the associated eigenfunctions (describing the spatial dependence of perturbations) from the vector $u$.

The eigenfrequencies of a star depend on the detailed internal structure of the star. Therefore, by comparing a set of eigenfrequencies for a given stellar model against those observed in a real star, we can test how well the model represents the real star --- a technique known as asteroseismology.

## Reference
- [Townsend, R. H. D., & Teitler, S. A. 2013, *MNRAS*, 435, 3406](https://ui.adsabs.harvard.edu/abs/2013MNRAS.435.3406T/abstract)
