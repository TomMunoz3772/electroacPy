# BEM module

## Interior and exterior studies
The boundary element method (BEM), as explained in {cite:ps}`Mechel2004FormulasOA`, aims to solve the following equation:

$$

\iint_S \left[ p(\textbf{y}) \frac{\partial G(\textbf{x}, \textbf{y})}{\partial n(\textbf{y})} - \frac{\partial p(\textbf{y})}{\partial n(\textbf{y})}G(\textbf{x}, \textbf{y}) \right]dS = \begin{cases} p(\textbf{x}), & \textbf{x} \in B_e \\ \frac{1}{2}p(\textbf{x}), & \textbf{x} \in S \\ 0,& \textbf{x} \in B_i \end{cases}

$$

where **x** is a point in space, **y** is a point on the system boundary $S$; and where $B_e$ and $B_i$ the exterior and interior domains. In the case of bempp-cl (and these equations) the propagation convention is $e^{+jkr}$[^conventionStuff]. The green function $G(\textbf{x}, \textbf{y})$ used is:

[^conventionStuff]: ElectroacPy uses $-k$ in the boundary and potential equations to change the convention to $e^{-jkr}$. Better results were obtained with it. Wouldn't be able to explain why.

$$
G(\textbf{x}, \textbf{y}) = \frac{1}{4\pi ||\textbf{x}-\textbf{y}||} e^{+jk||\textbf{x}-\textbf{y}||}.
$$

Generally, a such system is solved in two steps:
- calculation of the pressure on the system’s boundaries (scattered field),
- calculation of the pressure at external points or grids. This second step will be discussed in the next section, because it is closely related to the **evaluation** class.

By writing the integral equations in a matrix form we can rewrite the integral equation under the form $Ax=b$, where we solve for $x$. For exterior-field calculations, the equation is the following:

$$
    \left[ D - \frac{1}{2} I \right] P_s = j\omega \rho S U,
$$

while it is 

$$
    \left[ D + \frac{1}{2} I \right] P_s = j\omega \rho S U,
$$

for interior-field calculation. With $P_s$ the scattered field, $D$ the double layer term, $I$ the identity matrix, $S$ the single layer term and $U$ the source velocity. Single and Double layer means that the scattered pressure is projected as dipoles, while the source $U$ is described as monopoles --- see ???. 


