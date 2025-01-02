# BEM module

## Interior and exterior studies
The boundary element method (BEM), as explained in {cite:ps}`Mechel2004FormulasOA`, aims to solve the following equation for exterior field studies:

$$

\iint_S \left[ p(\textbf{y}) \frac{\partial G(\textbf{x}, \textbf{y})}{\partial n(\textbf{y})} - \frac{\partial p(\textbf{y})}{\partial n(\textbf{y})}G(\textbf{x}, \textbf{y}) \right]dS = \begin{cases} p(\textbf{x}), & \textbf{x} \in B_e \\ \frac{1}{2}p(\textbf{x}), & \textbf{x} \in S \\ 0,& \textbf{x} \in B_i \end{cases}

$$

where **x** is a point in space, **y** is a point on the system boundary $S$; and where $B_e$ and $B_i$ the exterior and interior domains. In the case of bempp-cl (and these equations) the propagation convention is $e^{+jkr}$[^conventionStuff]. The green function $G(\textbf{x}, \textbf{y})$ used is:

[^conventionStuff]: ElectroacPy uses $-k$ in the boundary and potential equations to change the convention to $e^{-jkr}$. Better results were obtained with it. Wouldn't be able to explain why.

$$
G(\textbf{x}, \textbf{y}) = \frac{1}{4\pi ||\textbf{x}-\textbf{y}||} e^{+jk||\textbf{x}-\textbf{y}||}.
$$

The BEM is done in two steps:
- calculation of the pressure on the system’s boundaries,
- calculation of the pressure at external points (or grids). Because it is closely related to the **evaluation** class, this point is discussed in the next section.

By writing the integral equation under the matrix form $Ax=b$, we can solve for $x$ --- the acoustic pressure:

$$
    \left[ D - \frac{1}{2} I \right] P_s = j\omega \rho S U,
$$

Where $P_s$​ represents the acoustic pressure on the boundaries, $I$ is the identity matrix, $D$ and $S$ denote the double and single layer terms, and $U$ is the source velocity.

The single and double layer terms correspond to the configuration of direct (acoustic radiator) and indirect sources (system boundaries). Simply put, the single layer represents monopole sources, while the double layer corresponds to dipole sources.

```{figure} ./BEM_backend_images/single_double_layer.svg
    :name: SD-layers

Single layer (blue) and double layer (black). The external pressure is the sum of the two projected fields. This is a simplified representation of BEM — see {cite:ps}`Mechel2004FormulasOA` for more complete explanation.
```

In the case of interior studies, the previous equation becomes
 
$$
    \left[ D + \frac{1}{2} I \right] P_s = j\omega \rho S U.
$$

The following code shows how to compute the acoustic pressure over the boundary of a system. To make it easier to read, a single frequency step is computed. The full frequency response can be computed by looping over a frequency array.

In the case of interior studies, the previous equation becomes
 
$$
    \left[ D + \frac{1}{2} I \right] P_s = j\omega \rho S U.
$$

The following code shows how to compute the acoustic pressure over the boundary of a system. To make it easier to read, a single frequency step is computed. The full frequency response can be computed by looping over a frequency array.

```python
import bempp.api
from bempp.api.operators.boundary import helmholtz, sparse
from bempp.api.linalg import gmres

f = 5000          # frequency (Hz)
w = 2 * np.pi * f # angular frequency (rad/s)
k = w / 343       # wavenumber

# importing system "grid"
mesh = 'path/to/meshFile.msh'      # we can take the studio monitor
grid = bempp.api.import_grid(mesh)

# domain where scattered pressure is evaluated
spaceP       = bempp.api.function_space(grid, "P", 1)
identity     = sparse.identity(spaceP, spaceP, spaceP)
double_layer = helmholtz.double_layer(spaceP, spaceP, spaceP)

# domain where source velocity is defined (known quantity)
spaceU       = bempp.api.function_space(grid, "DP", 0, segments=[1]) # segment=[1]: physical group corresponding to the driver
dof_count    = spaceU.grid_dof_count  # to apply the same velocity on all elements of the radiator
u_in         = bempp.api.GridFunction(spaceU, 
                                      coefficients=np.ones(dof_count))
single_layer = helmholtz.single_layer(spaceU, spaceP, spaceP, k)

# get total scattered field
p_total, _ = gmres(double_layer - 0.5 * identity, 
                1j * w * rho * single_layer * u_in, tol=1e-5)
# if Gmsh is installed and linked to Python:
p_total.plot(transformation="real")
```

If Gmsh is installed --- and accessible through Python --- `p_total.plot(transformation="real")` will display the pressure over the system's mesh:

```{figure} ./BEM_backend_images/pressure_sphere.png
    :name: pressure-sphere
    :scale: 75 %

Pressure distribution over a spherical loudspeaker.
```

