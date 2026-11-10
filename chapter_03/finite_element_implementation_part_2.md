# Assembly of the stiffness matrix

% CLAUDE: Briefly discuss background material that the students need to know for solving @project/project_part_2.md. 
% In particular:
% 1. Explain that the assembly of the stiffness matrix is very similar to the assembly of the mass matrix except for the following
% - need to compute local stiffness matrix $A_T$ instead 
% - matrix is assembled without taking Dirichlet b.c. into account; % Thus we need to postprocess the matrix and the right-hand side vector to account for Dirichlet b.c.
% For Neumann and Robin b.c., we need to compute the corresponding contributions to the right-hand side vector $b$.
% Focus the discussion on piecewise linear functions, including the implicated simplications.
% Explain how we arrive at the formula for the coefficients of
% the linear shape function $\lambda_i^{(T)} = a_i + b_i x + c_i y$. Derive the corresponding analytical expression similar to the one mention in @project/project_part_2.md when the students are asked to implement the "compute_hat_gradients".
%
% Afterwards explain how to incorporate Dirichlet b.c. into the linear system in accordance with the description in the apply_bcs function in @project/project_part_2_solutions.md

This section collects the background needed for
[part 2 of the project](../project/project_part_2.ipynb), where we build a finite element
solver for the Poisson problem. We continue to work with continuous, piecewise linear
functions and use the notation of the previous section, except that we write $T$ for a
generic triangle of $\mathcal{T}_h$ and drop the index $k$.

(ssec:from-mass-to-stiffness)=
## From the mass matrix to the stiffness matrix

We consider the Poisson problem with mixed boundary conditions,

$$
\begin{aligned}
-\Delta u &= f && \text{in } \Omega, \\
u &= g_D && \text{on } \partial\Omega_D, \\
\nabla u \cdot n &= g_N && \text{on } \partial\Omega_N,
\end{aligned}
$$ (eq:poisson-mixed-p1)

where $n$ is the outward pointing unit normal. Multiplying with a test function $v$ that
vanishes on $\partial\Omega_D$ and integrating by parts gives

$$
\int_\Omega \nabla u \cdot \nabla v \,\mathrm{d}x
= \int_\Omega f v \,\mathrm{d}x + \int_{\partial\Omega_N} g_N v \,\mathrm{d}S .
$$

Inserting $u_h = \sum_{j=0}^{N_p-1} U_j \phi_j$ and testing with $v = \phi_i$ leads to a
linear system $A \mathbf{U} = \mathbf{b}$ with the **stiffness matrix** $A$ and the load
vector $\mathbf{b}$,

$$
A_{ij} = \int_\Omega \nabla \phi_j \cdot \nabla \phi_i \,\mathrm{d}x,
\qquad
b_i = \int_\Omega f \phi_i \,\mathrm{d}x + \int_{\partial\Omega_N} g_N \phi_i \,\mathrm{d}S,
\qquad i, j = 0, \ldots, N_p - 1.
$$ (eq:stiffness-load-p1)

The assembly of $A$ is very similar to the assembly of the mass matrix in
{prf:ref}`alg:assemble-mass-2d`: we loop over the elements, compute a local matrix and add
it to the global matrix via the local-to-global map. There are only three differences.

1. On each element, we compute the **element stiffness matrix** $A_T$ instead of the
   element mass matrix.
2. We assemble $A$ for **all** vertices, without taking the Dirichlet boundary conditions
   into account. Afterwards, we postprocess the matrix $A$ and the vector $\mathbf{b}$ to
   account for them.
3. Neumann and Robin boundary conditions lead to additional integrals over boundary edges.

(ssec:element-stiffness-matrix-p1)=
## The element stiffness matrix

Let $T$ be a triangle with the vertices $N_0, N_1, N_2$ in local numbering, and let
$\lambda^{(T)}_\alpha$, $\alpha = 0, 1, 2$, be its barycentric coordinates, so that
$\phi_i|_T = \lambda^{(T)}_\alpha$ if $i = \iota(T, \alpha)$. All non-zero contributions of
$T$ to $A$ are collected in the matrix $A_T \in \mathbb{R}^{3 \times 3}$ with the entries

$$
(A_T)_{\alpha\beta}
= \int_T \nabla \lambda^{(T)}_\beta \cdot \nabla \lambda^{(T)}_\alpha \,\mathrm{d}x,
\qquad \alpha, \beta = 0, 1, 2.
$$

For $\mathbb{P}_1$-elements, the shape functions are affine, that is,

$$
\lambda^{(T)}_\alpha(x, y) = a_\alpha + b_\alpha x + c_\alpha y,
$$ (eq:p1-shape-function)

and their gradients $\nabla \lambda^{(T)}_\alpha = (b_\alpha, c_\alpha)^\top$ are constant
on $T$. The integrand is therefore constant as well, and the integral is simply the value
of the integrand times the area. With $b = (b_0, b_1, b_2)^\top$ and
$c = (c_0, c_1, c_2)^\top$, we obtain

$$
(A_T)_{\alpha\beta} = |T| \, (b_\alpha b_\beta + c_\alpha c_\beta),
\qquad \text{that is,} \qquad
A_T = |T| \, \big( b \, b^\top + c \, c^\top \big).
$$ (eq:element-stiffness-matrix-p1)

No quadrature rule is needed. In contrast to the mass matrix, however, $A_T$ is not a
multiple of one fixed reference matrix, since the coefficients $b$ and $c$ depend on the
shape of $T$. It remains to compute them.

(ssec:p1-shape-function-coefficients)=
## The coefficients of the linear shape functions

Let $(x_\alpha, y_\alpha)$ be the coordinates of $N_\alpha$. The shape functions
$\lambda^{(T)}_\alpha$ are the barycentric coordinates of $T$
({prf:ref}`def:barycentric-simplex`), so we can start from the representation
@eq:barycentric-cramer, which we obtained with Cramer's rule. For $d = 2$, the barycentric coordinates are the
solution of the linear system

$$
\underbrace{
\begin{pmatrix}
x_0 & x_1 & x_2 \\
y_0 & y_1 & y_2 \\
1 & 1 & 1
\end{pmatrix}
}_{=: B}
\begin{pmatrix}
\lambda^{(T)}_0(x, y) \\ \lambda^{(T)}_1(x, y) \\ \lambda^{(T)}_2(x, y)
\end{pmatrix}
=
\begin{pmatrix} x \\ y \\ 1 \end{pmatrix},
$$

and Cramer's rule gives

$$
\lambda^{(T)}_\alpha(x, y) = \frac{\det B_\alpha(x, y)}{\det B},
\qquad \alpha = 0, 1, 2,
$$ (eq:p1-barycentric-cramer)

where $B_\alpha(x, y)$ is obtained from $B$ by replacing its column $\alpha$ with
$(x, y, 1)^\top$. For the denominator, we subtract the first column of $B$ from the other
two and expand along the last row, which yields

$$
D := \det B = (x_1 - x_0)(y_2 - y_0) - (x_2 - x_0)(y_1 - y_0) ,
$$ (eq:p1-signed-determinant)

and $|D| = 2|T|$ by {prf:ref}`cor:simplex-volume`. For the numerator, we take local
indices modulo $3$. A cyclic permutation of the columns does not change the determinant,
so we can move the column $(x, y, 1)^\top$ to the front and expand along it:

$$
\begin{aligned}
\det B_\alpha(x, y)
&= \det \begin{pmatrix}
x & x_{\alpha+1} & x_{\alpha+2} \\
y & y_{\alpha+1} & y_{\alpha+2} \\
1 & 1 & 1
\end{pmatrix} \\
&= \big( x_{\alpha+1} y_{\alpha+2} - x_{\alpha+2} y_{\alpha+1} \big)
+ \big( y_{\alpha+1} - y_{\alpha+2} \big) \, x
+ \big( x_{\alpha+2} - x_{\alpha+1} \big) \, y .
\end{aligned}
$$

Dividing by $D$ and comparing with @eq:p1-shape-function gives

$$
a_\alpha = \frac{x_{\alpha+1} y_{\alpha+2} - x_{\alpha+2} y_{\alpha+1}}{D},
\qquad
b_\alpha = \frac{y_{\alpha+1} - y_{\alpha+2}}{D},
\qquad
c_\alpha = \frac{x_{\alpha+2} - x_{\alpha+1}}{D}.
$$ (eq:p1-shape-function-coefficients)

Only $b_\alpha$ and $c_\alpha$ enter the stiffness matrix. Collecting them for
$\alpha = 0, 1, 2$, we recognize the components of two cross products in $\mathbb{R}^3$:

$$
b = \frac{1}{D} \, (y_0, y_1, y_2)^\top \times (1, 1, 1)^\top,
\qquad
c = \frac{1}{D} \, (1, 1, 1)^\top \times (x_0, x_1, x_2)^\top .
$$ (eq:p1-gradients-cross-product)

The sign of $D$ matters here. We have $D = 2|T| > 0$ if the vertices are numbered
counter-clockwise and $D = -2|T| < 0$ if they are numbered clockwise. Dividing by the
*signed* determinant $D$ gives the correct gradients for both orderings, whereas the area
in @eq:element-stiffness-matrix-p1 is $|T| = |D|/2$. This is worth keeping in mind, since
not every mesh generator returns consistently oriented triangles.

(ssec:neumann-robin-assembly-p1)=
## Neumann and Robin boundary conditions

The Neumann term in @eq:stiffness-load-p1 only changes the load vector. As for the volume
term, we split the integral into contributions from the boundary edges $E$ of the mesh
that make up $\partial\Omega_N$,

$$
\int_{\partial\Omega_N} g_N \phi_i \,\mathrm{d}S
= \sum_{E \subset \partial\Omega_N} \int_E g_N \phi_i \,\mathrm{d}S .
$$

On an edge $E$, only the two basis functions associated with the end points of $E$ are
non-zero, so every edge contributes to exactly two entries of $\mathbf{b}$. What we need
is a list of these edges together with their end points, which plays the role of the
connectivity matrix $C$ for the boundary.

A Robin condition is treated in the same way and contributes a similar boundary integral of
its data to $\mathbf{b}$. In addition, since a Robin condition involves the unknown $u$
itself, it also contributes a boundary integral to the matrix $A$, which couples the two
end points of every Robin edge.

(ssec:dirichlet-postprocessing-p1)=
## Incorporating Dirichlet boundary conditions

After the assembly, the system $A \mathbf{U} = \mathbf{b}$ contains one equation for every
vertex. Let $\mathcal{I}_D$ be the set of the indices $i$ of all vertices $N_i$ on
$\partial\Omega_D$. For $i \in \mathcal{I}_D$, the basis function $\phi_i$ does not vanish
on $\partial\Omega_D$, so it is not an admissible test function and the $i$-th equation
has no justification. We use these rows instead to prescribe the boundary values, and we
replace the equation by

$$
U_i = g_D(N_i), \qquad i \in \mathcal{I}_D,
$$

which means that $g_D$ is approximated by its nodal interpolant. In terms of $A$ and
$\mathbf{b}$, we perform the following steps for all $i \in \mathcal{I}_D$:

1. set $b_i \gets g_D(N_i)$,
2. set the entire row $i$ of $A$ to zero, $A_{ij} \gets 0$ for all $j$,
3. set the diagonal entry to one, $A_{ii} \gets 1$.

To see the effect, let us number the free vertices first, which we do only for this
illustration. The subscripts $F$ and $D$ denote the free and the Dirichlet vertices. The
modified system then reads

$$
\begin{pmatrix} A_{FF} & A_{FD} \\ 0 & I \end{pmatrix}
\begin{pmatrix} \mathbf{U}_F \\ \mathbf{U}_D \end{pmatrix}
=
\begin{pmatrix} \mathbf{b}_F \\ \mathbf{g}_D \end{pmatrix}.
$$ (eq:dirichlet-modified-system)

The second block row gives $\mathbf{U}_D = \mathbf{g}_D$, and the first block row then
becomes $A_{FF} \mathbf{U}_F = \mathbf{b}_F - A_{FD} \mathbf{g}_D$. These are exactly the
Galerkin equations for the free vertices, with the known boundary values moved to the
right-hand side. Since the columns belonging to $\mathcal{I}_D$ are left untouched, the
modified matrix is no longer symmetric, which poses no problem for a sparse direct solver.

Two practical remarks are in order. The Dirichlet step has to come **after** all other
contributions have been added, since a vertex at the end of a Neumann edge may also lie on
$\partial\Omega_D$, and its load vector entry must be overwritten. Moreover, the step
changes individual matrix entries, so it should be carried out in a format that supports
this, such as DOK, before the matrix is converted to CSR or CSC for the solve, see
[](#ssec:sparse-storage-formats).