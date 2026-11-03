---
jupytext:
  cell_metadata_filter: -jupyter
  formats: ipynb,md:myst
  main_language: python
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.5
kernelspec:
  display_name: tma4220-finite-elements (3.14.10)
  language: python
  name: python3
---

# Some aspects of the implementation of finite element methods

In this section, we discuss the implementation of a finite element solver for
a Poisson problem 
based on $\mathbb{P}_1$-Lagrange elements. Here, we rather focus on some prevalent implementation aspects and principles without aiming for the most efficient realization,
acknowledging that 
high-performance computing capabilities of FEM software are indispensable
for simulating larger and realistic problems.

+++

## Data structures for meshes

+++

The simplest way to store a triangulation $\mathcal{T}_h$ consisting of $N_p$
vertices and $N_t$ triangles uses two arrays.

* The **point matrix** $P \in \mathbb{R}^{N_p \times 2}$ stores the vertex
  coordinates. Row $i$ contains the coordinates $(x_i, y_i)$ of the vertex
  with global index $i$.
* The **connectivity matrix** $C \in \mathbb{N}^{N_t \times 3}$ stores the
  triangles. Row $k$ contains the global indices of the three vertices of
  triangle $T_k$, which we list in counter-clockwise order.

All geometric information is therefore kept in $P$, while $C$ only describes
how the vertices are connected to form triangles. To get the coordinates of the
vertices of $T_k$, we look up the rows $C_{k0}, C_{k1}, C_{k2}$ of $P$. Following
Python's convention, we start counting vertices and triangles at $0$.

```{figure} figures/mesh.png
:label: fig-mesh-point-connectivity
:alt: The unit square divided into eight triangles. Nine vertices, numbered 0 to 8 in red, lie on a uniform 3 by 3 grid with spacing 0.5, numbered row by row from the bottom left. The triangles are labelled #0 to #7. Four triangles, #0, #2, #4 and #6, form a diamond around the centre vertex 4, and the remaining four, #1, #3, #5 and #7, fill the corners of the square.
:width: 70%

A triangulation of the unit square with $N_p = 9$ vertices (red numbers) and
$N_t = 8$ triangles (labelled #0–#7).
```

For the triangulation in [](#fig-mesh-point-connectivity), the point and
connectivity matrices are

$$
P =
\begin{pmatrix}
0   & 0   \\
0.5 & 0   \\
1   & 0   \\
0   & 0.5 \\
0.5 & 0.5 \\
1   & 0.5 \\
0   & 1   \\
0.5 & 1   \\
1   & 1
\end{pmatrix},
\qquad
C =
\begin{pmatrix}
1 & 4 & 3 \\
0 & 1 & 3 \\
1 & 5 & 4 \\
1 & 2 & 5 \\
3 & 4 & 7 \\
3 & 7 & 6 \\
4 & 5 & 7 \\
5 & 8 & 7
\end{pmatrix}.
$$

For example, the second row of $C$ says that triangle #1 has the vertices
$0$, $1$ and $3$, which, according to rows $0$, $1$ and $3$ of $P$, lie at
$(0,0)$, $(0.5,0)$ and $(0,0.5)$.

+++

Next, let us define and plot some simple functions on a given triangulation.
For that, let's get a bit acquainted with the local **meshtools** and the **plottools** package, which provides a convenient interface to the mesh data structures and allows for easy plotting of functions defined on the mesh.

First, we create a mesh for the unit square $\Omega = (0,1)^2$ with a uniform triangulation based on $N$ subdivision in each 
axis.

```{code-cell} ipython3
%matplotlib widget
import numpy as np
import meshtools as mt
import plottools as pt
```

```{code-cell} ipython3
N = 10
P, T = mt.unitsquare_mesh(N)
pt.plot_mesh_2d(P, T)
```

Next, we plot a function on top of the mesh. For that, we define a function $f(x,y)$ and evaluate it at the vertices of the mesh.
Note that this corresponds to a piecewise linear interpolation of the function on the triangulation.

```{code-cell} ipython3
def f(x,y):
    return np.sin(2*np.pi*x)*np.cos(np.pi*y)

X = P[:,0]
Y = P[:,1]
F = f(X,Y)
pt.plot2D(X, Y, F, T)
```

## Outline of the assembly algorithm

As a first guiding example in 2D, we consider the $L^2$-projection onto continuous,
piecewise linear functions. It requires the assembly of a mass matrix and a load vector,
which already contains all the key steps of the general assembly process, while avoiding
the derivatives and boundary terms of a full boundary value problem.

### The $L^2$-projection onto $\mathbb{P}_1$-elements

Let $\mathcal{T}_h = \{T_k\}_{k=0}^{N_t-1}$ be a triangulation of $\Omega \subset \mathbb{R}^2$
with vertices $\{a_i\}_{i=0}^{N_p-1}$, stored in the point matrix $P$ and the
connectivity matrix $C$ as described above, and let
$V_h = \mathbb{P}_1^c(\mathcal{T}_h)$ be the space of continuous, piecewise linear
functions. For a given $f \in L^2(\Omega)$, the $L^2$-projection $P_h f \in V_h$ is
defined by

$$
\int_\Omega P_h f \, v_h \,\mathrm{d}x
= \int_\Omega f \, v_h \,\mathrm{d}x
\quad \forall\, v_h \in V_h.
$$ (eq:l2-projection-p1)

Let $\{\phi_i\}_{i=0}^{N_p-1}$ be the nodal basis of $V_h$, that is,
$\phi_i \in V_h$ with $\phi_i(a_j) = \delta_{ij}$. Writing
$P_h f = \sum_{j=0}^{N_p-1} U_j \phi_j$ and testing with $v_h = \phi_i$, we see that
@eq:l2-projection-p1 is equivalent to the linear system

$$
M \mathbf{U} = \mathbf{b}
$$

for the coefficient vector $\mathbf{U} = (U_j)_{j=0}^{N_p-1} \in \mathbb{R}^{N_p}$, where the
**mass matrix** $M \in \mathbb{R}^{N_p \times N_p}$ and the **load vector**
$\mathbf{b} \in \mathbb{R}^{N_p}$ are given by

$$
M_{ij} = \int_\Omega \phi_j \phi_i \,\mathrm{d}x,
\qquad
b_i = \int_\Omega f \phi_i \,\mathrm{d}x,
\qquad i, j = 0, \ldots, N_p - 1.
$$

As in 1D, we compute these quantities element by element and exploit that each basis
function $\phi_i$ is non-zero only on the triangles that share the vertex $a_i$.

+++

### Assembly of the mass matrix

Decomposing the integral into a sum over the elements gives

$$
M_{ij} = \sum_{k=0}^{N_t-1} \int_{T_k} \phi_j \phi_i \,\mathrm{d}x.
$$

The integral over $T_k$ **can only be non-zero if both $a_i$ and $a_j$ are vertices
of $T_k$**. Let $a^{(k)}_0, a^{(k)}_1, a^{(k)}_2$ denote the vertices of $T_k$ in the
order given by the $k$-th row of $C$. If $a_i = a^{(k)}_\alpha$, the restriction
$\phi_i|_{T_k}$ is the linear function that equals $1$ at $a^{(k)}_\alpha$ and $0$ at the
two other vertices, so it is precisely the barycentric coordinate $\lambda^{(k)}_\alpha$
of $T_k$:

$$
\phi_i|_{T_k} = \lambda^{(k)}_\alpha .
$$

Here $\alpha \in \{0, 1, 2\}$ is a **local index** on the element $T_k$, while $i$ is a
**global index**. The **local-to-global map** $i = \iota(k, \alpha)$, which in 1D had to be
constructed separately, is for $\mathbb{P}_1$-elements simply read off from the connectivity
matrix:

$$
\iota(k, \alpha) = C_{k\alpha}, \qquad \alpha = 0, 1, 2, \quad k = 0, \ldots, N_t - 1.
$$

All non-zero contributions of $T_k$ to $M$ are thus collected in the
**element mass matrix** $M^{(k)} \in \mathbb{R}^{3 \times 3}$,

$$
M^{(k)}_{\alpha\beta} = \int_{T_k} \lambda^{(k)}_\beta \lambda^{(k)}_\alpha \,\mathrm{d}x,
\qquad \alpha, \beta = 0, 1, 2.
$$ (eq:element-mass-matrix-p1)

To compute $M^{(k)}$, we again pass to the reference triangle
$\widehat{T} = \operatorname{conv}\{(0,0), (1,0), (0,1)\}$ with the barycentric coordinates

$$
\widehat{\lambda}_0(\hat{x}, \hat{y}) = 1 - \hat{x} - \hat{y}, \qquad
\widehat{\lambda}_1(\hat{x}, \hat{y}) = \hat{x}, \qquad
\widehat{\lambda}_2(\hat{x}, \hat{y}) = \hat{y},
$$

see @eq:barycentric-unit-simplex. The affine map

$$
F_k(\hat{x}) = B_k \hat{x} + a^{(k)}_0,
\qquad
B_k = \big( a^{(k)}_1 - a^{(k)}_0 \;\big|\; a^{(k)}_2 - a^{(k)}_0 \big) \in \mathbb{R}^{2 \times 2},
$$ (eq:affine-map-p1)

maps the vertices of $\widehat{T}$ to the vertices of $T_k$, and since affine maps preserve
barycentric coordinates ({prf:ref}`lem:barycentric-affine-invariance`), we have
$\lambda^{(k)}_\alpha = \widehat{\lambda}_\alpha \circ F_k^{-1}$. A change of variables
with $\mathrm{d}x = |\det B_k| \,\mathrm{d}\hat{x}$ and
$|\det B_k| = |T_k| / |\widehat{T}|$ then yields

$$
M^{(k)}_{\alpha\beta}
= \frac{|T_k|}{|\widehat{T}|} \int_{\widehat{T}} \widehat{\lambda}_\beta \widehat{\lambda}_\alpha \,\mathrm{d}\hat{x}
= |T_k| \, \widehat{M}_{\alpha\beta},
\qquad
\widehat{M}_{\alpha\beta} := \frac{1}{|\widehat{T}|} \int_{\widehat{T}} \widehat{\lambda}_\beta \widehat{\lambda}_\alpha \,\mathrm{d}\hat{x}.
$$ (eq:element-mass-matrix-ref-p1)

As in 1D, the **reference mass matrix** $\widehat{M}$ is independent of the element and
is computed once and for all. We include the factor $1/|\widehat{T}|$ in its definition,
so that $\widehat{M}_{\alpha\beta}$ is the *mean value* of
$\widehat{\lambda}_\beta \widehat{\lambda}_\alpha$ over $\widehat{T}$, in the same way as
the area is factored out of the quadrature rules in @eq:quadrature-reference-triangle.
For $\widehat{T} = [0,1]$ in 1D, where $|\widehat{T}| = 1$, this agrees with the earlier
definition. The only element-dependent quantity is then the area

$$
|T_k| = \tfrac{1}{2} |\det B_k|
= \tfrac{1}{2} \big| (x_1 - x_0)(y_2 - y_0) - (x_2 - x_0)(y_1 - y_0) \big|,
$$

where $(x_\alpha, y_\alpha)$ are the coordinates of $a^{(k)}_\alpha$. Since $\widehat{M}$ is
the integral of a quadratic polynomial, it can be computed exactly, either symbolically or
with any quadrature rule of degree $2$ from @tab:quadrature-rules-triangle.

:::{prf:algorithm} Assembly of the $\mathbb{P}_1$ mass matrix
:label: alg:assemble-mass-2d

**Input:**
- point matrix $P \in \mathbb{R}^{N_p \times 2}$
- connectivity matrix $C \in \mathbb{N}^{N_t \times 3}$

**Output:** Mass matrix $M \in \mathbb{R}^{N_p \times N_p}$ with
$M_{ij} = \int_\Omega \phi_j \phi_i \,\mathrm{d}x$ for $i, j = 0, \ldots, N_p - 1$.

1. $M \gets 0$
2. **for** $\alpha, \beta = 0, 1, 2$ **do**
   1. $\widehat{M}_{\alpha\beta} \gets \frac{1}{|\widehat{T}|} \int_{\widehat{T}} \widehat{\lambda}_\beta \, \widehat{\lambda}_\alpha \,\mathrm{d}\hat{x}$
3. **for** $k = 0, \ldots, N_t - 1$ **do**
   1. $a^{(k)}_\alpha \gets$ row $C_{k\alpha}$ of $P$ for $\alpha = 0, 1, 2$
   2. $|T_k| \gets \frac{1}{2} \big| \det \big( a^{(k)}_1 - a^{(k)}_0 \;\big|\; a^{(k)}_2 - a^{(k)}_0 \big) \big|$
   3. $M^{(k)} \gets |T_k| \, \widehat{M}$
   4. **for** $\alpha, \beta = 0, 1, 2$ **do**
      1. $i \gets C_{k\alpha}$, $\; j \gets C_{k\beta}$
      2. $M_{ij} \gets M_{ij} + M^{(k)}_{\alpha\beta}$
4. **return** $M$
:::

+++

### Assembly of the load vector

The load vector is decomposed in the same way,

$$
b_i = \sum_{k=0}^{N_t-1} \int_{T_k} f \phi_i \,\mathrm{d}x,
$$

and only the elements with $a_i = a^{(k)}_\alpha$ contribute, through the
**element load vector**

$$
b^{(k)}_\alpha = \int_{T_k} f \lambda^{(k)}_\alpha \,\mathrm{d}x, \qquad \alpha = 0, 1, 2.
$$

For a general $f$, this integral cannot be computed exactly, and we approximate it with a
quadrature rule $\{(\widehat{\zeta}_q, \omega_q)\}_{q=1}^{n_q}$ on $\widehat{T}$, see
[](#ssec:quadrature-rules-triangle). By
{prf:ref}`rem:quadrature-physical-triangle`, the change of variables gives

$$
b^{(k)}_\alpha
= \frac{|T_k|}{|\widehat{T}|} \int_{\widehat{T}} f\big(F_k(\hat{x})\big)\, \widehat{\lambda}_\alpha(\hat{x}) \,\mathrm{d}\hat{x}
\approx |T_k| \sum_{q=1}^{n_q} \omega_q \, \widehat{\lambda}_\alpha(\widehat{\zeta}_q) \, f\big(F_k(\widehat{\zeta}_q)\big).
$$ (eq:element-load-vector-p1)

The values of the reference basis functions at the quadrature points,

$$
\widehat{\Lambda}_{\alpha q} = \widehat{\lambda}_\alpha(\widehat{\zeta}_q),
\qquad \alpha = 0, 1, 2, \quad q = 1, \ldots, n_q,
$$

are independent of the element, and so are the weighted values
$\widehat{L}_{\alpha q} = \omega_q \widehat{\Lambda}_{\alpha q}$. Both are precomputed once
for the chosen quadrature rule. For $\mathbb{P}_1$-elements, the column
$(\widehat{\Lambda}_{0q}, \widehat{\Lambda}_{1q}, \widehat{\Lambda}_{2q})$ is simply the
vector of barycentric coordinates of $\widehat{\zeta}_q$, as listed in
@tab:quadrature-rules-triangle. Element by element, it only remains to evaluate $f$ at the
mapped quadrature points and to scale with $|T_k|$. The mapped points are again computed
from the barycentric coordinates,

$$
F_k(\widehat{\zeta}_q) = \sum_{\beta=0}^{2} \widehat{\Lambda}_{\beta q} \, a^{(k)}_\beta ,
$$

which follows from @eq:affine-map-p1 and @eq:barycentric-unit-simplex.

:::{prf:algorithm} Assembly of the $\mathbb{P}_1$ load vector
:label: alg:assemble-load-2d

**Input:**
- point matrix $P \in \mathbb{R}^{N_p \times 2}$
- connectivity matrix $C \in \mathbb{N}^{N_t \times 3}$
- right-hand side $f$
- quadrature rule $\{(\widehat{\zeta}_q, \omega_q)\}_{q=1}^{n_q}$ on $\widehat{T}$ with $\sum_q \omega_q = 1$

**Output:** Load vector $\mathbf{b} \in \mathbb{R}^{N_p}$ with
$b_i \approx \int_\Omega f \phi_i \,\mathrm{d}x$ for $i = 0, \ldots, N_p - 1$.

1. $\mathbf{b} \gets 0$
2. **for** $\alpha = 0, 1, 2$, $\; q = 1, \ldots, n_q$ **do**
   1. $\widehat{\Lambda}_{\alpha q} \gets \widehat{\lambda}_\alpha(\widehat{\zeta}_q)$
   2. $\widehat{L}_{\alpha q} \gets \omega_q \, \widehat{\Lambda}_{\alpha q}$
3. **for** $k = 0, \ldots, N_t - 1$ **do**
   1. $a^{(k)}_\alpha \gets$ row $C_{k\alpha}$ of $P$ for $\alpha = 0, 1, 2$
   2. $|T_k| \gets \frac{1}{2} \big| \det \big( a^{(k)}_1 - a^{(k)}_0 \;\big|\; a^{(k)}_2 - a^{(k)}_0 \big) \big|$
   3. **for** $q = 1, \ldots, n_q$ **do**
      1. $\widehat{f}_q \gets f\big(F_k(\widehat{\zeta}_q)\big) = f\big(\sum_{\beta=0}^{2} \widehat{\Lambda}_{\beta q} \, a^{(k)}_\beta\big)$
   4. **for** $\alpha = 0, 1, 2$ **do**
      1. $b^{(k)}_\alpha \gets |T_k| \sum_{q=1}^{n_q} \widehat{L}_{\alpha q} \, \widehat{f}_q$
      2. $i \gets C_{k\alpha}$
      3. $b_i \gets b_i + b^{(k)}_\alpha$
4. **return** $\mathbf{b}$
:::

```{prf:remark} Matrix form of the element computations
:label: rem:load-vector-matrix-form

The loop body of {prf:ref}`alg:assemble-load-2d` can be written compactly with matrices.
Let $X^{(k)} \in \mathbb{R}^{3 \times 2}$ contain the coordinates of the vertices
$a^{(k)}_0, a^{(k)}_1, a^{(k)}_2$ as rows, that is, the rows $C_{k0}, C_{k1}, C_{k2}$ of $P$.
Then the $q$-th row of

$$
\widehat{\Lambda}^\top X^{(k)} \in \mathbb{R}^{n_q \times 2}
$$

is the mapped quadrature point $F_k(\widehat{\zeta}_q)$. Applying $f$ row-wise gives the
vector $\widehat{\mathbf{f}}^{(k)} = (\widehat{f}_q)_{q=1}^{n_q}$, and the element load
vector is the matrix-vector product

$$
\mathbf{b}^{(k)} = |T_k| \, \widehat{L} \, \widehat{\mathbf{f}}^{(k)} \in \mathbb{R}^3.
$$

Similarly, the reference mass matrix of @eq:element-mass-matrix-ref-p1 is obtained from any
quadrature rule of degree at least $2$ as $\widehat{M} = \widehat{L} \, \widehat{\Lambda}^\top$.
In NumPy, with `Lam` $= \widehat{\Lambda}$ and `L` $= \widehat{L}$, these read
`X_q = Lam.T @ P[C[k]]`, `b_k = area * L @ f(X_q[:, 0], X_q[:, 1])` and `M_ref = L @ Lam.T`.
```

Once $M$ and $\mathbf{b}$ are assembled, the coefficients of $P_h f$ are obtained by solving
$M \mathbf{U} = \mathbf{b}$. Note that both algorithms share the same structure: a
precomputation on the reference element, followed by a single loop over the elements that
computes a geometric factor, forms the local contribution, and adds it to the global
matrix or vector via the local-to-global map. The same structure carries over to the
stiffness matrix of the Poisson problem, where the gradients of the basis functions enter
as an additional geometric factor.

+++

(ssec:quadrature-rules-triangle)=
## Quadrature rules on the reference triangle

+++

When we assemble the load vector, or a mass matrix with variable coefficients, the
element integrals are approximated by a quadrature rule. As in the one-dimensional
case, we state the rules once on the reference triangle

$$
\widehat{T} = \operatorname{conv}\{(0,0), (1,0), (0,1)\}
= \{ (\hat{x}, \hat{y}) : \hat{x}, \hat{y} \geqslant 0,\ \hat{x} + \hat{y} \leqslant 1 \},
$$

and approximate

$$
\int_{\widehat{T}} \hat{f}(\hat{x}) \,\mathrm{d}\hat{x}
\approx |\widehat{T}| \sum_{q=1}^{n_q} \omega_q \, \hat{f}(\widehat{\zeta}_q)
$$ (eq:quadrature-reference-triangle)

with points $\widehat{\zeta}_q \in \widehat{T}$ and weights $\omega_q$, where
$|\widehat{T}| = \frac{1}{2}$. With the area factored out, the weights of every rule
below sum to $1$. A rule has
*degree of exactness* $p$ if @eq:quadrature-reference-triangle holds with equality for
all $\hat{f} \in \mathbb{P}_p(\widehat{T})$ but not for all $\hat{f} \in \mathbb{P}_{p+1}(\widehat{T})$.
To check this, it suffices to test the monomials, for which

$$
\int_{\widehat{T}} \hat{x}^i \hat{y}^j \,\mathrm{d}\hat{x} = \frac{i!\, j!}{(i+j+2)!},
$$ (eq:monomial-integral-triangle)

see @exer-monomial-integral-triangle.

Quadrature points on triangles are most conveniently given in barycentric
coordinates $(\widehat{\lambda}_0, \widehat{\lambda}_1, \widehat{\lambda}_2)$. By
@eq:barycentric-unit-simplex, the point with barycentric coordinates
$(\lambda_0, \lambda_1, \lambda_2)$ has the Cartesian coordinates
$(\hat{x}, \hat{y}) = (\lambda_1, \lambda_2)$ in $\widehat{T}$.
The rules in common use are *symmetric*: whenever a point is used, so are all points
obtained by permuting its barycentric coordinates, and all of them carry the same
weight. The *multiplicity* of a point is the number of distinct points obtained this way.
For instance, $(a, a, 1-2a)$ with $a \neq \frac{1}{3}$ has multiplicity $3$ and stands for
the points $(a, a, 1-2a)$, $(a, 1-2a, a)$ and $(1-2a, a, a)$, that is, for
$(\hat{x}, \hat{y}) = (a, 1-2a)$, $(1-2a, a)$ and $(a, a)$.
A symmetric rule treats the three vertices alike, so it does not depend on how the
vertices of a triangle are numbered.

```{table} Quadrature rules on the reference triangle $\widehat{T}$. Except for the vertex rule, they follow Table 8.2 in Ern and Guermond (2004). Points are given in barycentric coordinates, each together with its multiplicity, and the weights $\omega_q$ are relative to the area as in @eq:quadrature-reference-triangle.
:label: tab:quadrature-rules-triangle

| Rule | Degree | $n_q$ | Barycentric coordinates | Multiplicity | Weights $\omega_q$ |
|:--|:-:|:-:|:--|:-:|:--|
| Centroid rule | $1$ | $1$ | $(\frac{1}{3}, \frac{1}{3}, \frac{1}{3})$ | $1$ | $1$ |
| Vertex rule | $1$ | $3$ | $(1, 0, 0)$ | $3$ | $\frac{1}{3}$ |
| Interior three-point rule | $2$ | $3$ | $(\frac{1}{6}, \frac{1}{6}, \frac{2}{3})$ | $3$ | $\frac{1}{3}$ |
| Edge midpoint rule | $2$ | $3$ | $(\frac{1}{2}, \frac{1}{2}, 0)$ | $3$ | $\frac{1}{3}$ |
| Four-point rule | $3$ | $4$ | $(\frac{1}{3}, \frac{1}{3}, \frac{1}{3})$ <br> $(\frac{1}{5}, \frac{1}{5}, \frac{3}{5})$ | $1$ <br> $3$ | $-\frac{9}{16}$ <br> $\frac{25}{48}$ |
| Seven-point rule | $3$ | $7$ | $(\frac{1}{3}, \frac{1}{3}, \frac{1}{3})$ <br> $(\frac{1}{2}, \frac{1}{2}, 0)$ <br> $(1, 0, 0)$ | $1$ <br> $3$ <br> $3$ | $\frac{9}{20}$ <br> $\frac{2}{15}$ <br> $\frac{1}{20}$ |
| Six-point rule | $4$ | $6$ | $(a_i, a_i, 1-2a_i)$, $i = 1, 2$ | $3$ | $\omega_i$, $i = 1, 2$ |
| Radon's seven-point rule | $5$ | $7$ | $(\frac{1}{3}, \frac{1}{3}, \frac{1}{3})$ <br> $(a_i, a_i, 1-2a_i)$, $i = 1, 2$ | $1$ <br> $3$ | $\frac{9}{40}$ <br> $\omega_i$, $i = 1, 2$ |
| Twelve-point rule | $6$ | $12$ | $(a_i, a_i, 1-2a_i)$, $i = 1, 2$ <br> $(a, b, 1-a-b)$ | $3$ <br> $6$ | $\omega_i$, $i = 1, 2$ <br> $\omega_3$ |
```

The six-point rule, Radon's seven-point rule and the twelve-point rule use the
following parameters. Only Radon's rule has a closed form; for the others, we give the
values to 15 digits.

$$
\begin{aligned}
&\text{Six-point}: &
a_1 &= 0.445948490915965, & \omega_1 &= 0.223381589678011, \\
&& a_2 &= 0.091576213509771, & \omega_2 &= 0.109951743655322, \\[1ex]
&\text{Radon}: &
a_1 &= \tfrac{6 - \sqrt{15}}{21}, & \omega_1 &= \tfrac{155 - \sqrt{15}}{1200}, \\
&& a_2 &= \tfrac{6 + \sqrt{15}}{21}, & \omega_2 &= \tfrac{155 + \sqrt{15}}{1200}, \\[1ex]
&\text{Twelve-point}: &
a_1 &= 0.063089014491502, & \omega_1 &= 0.050844906370206, \\
&& a_2 &= 0.249286745170910, & \omega_2 &= 0.116786275726378, \\
&& a &= 0.310352451033785, \quad b = 0.053145049844816, & \omega_3 &= 0.082851075618374.
\end{aligned}
$$

Apart from the vertex rule, the table is taken from {cite}`ErnGuermond2004`. The
interior three-point rule and the four-point rule are due to
{cite}`HammerMarloweStroud1956`, Radon's seven-point rule to {cite}`Radon1948`, and the
six-point and twelve-point rules to {cite}`Dunavant1985`, who tabulates symmetric rules
up to degree $20$. The remaining rules are classical and are collected, together with
many others, in {cite}`Stroud1971`.

```{prf:remark} Quadrature on a physical triangle
:label: rem:quadrature-physical-triangle

Let $F_T(\hat{x}) = B_T \hat{x} + b_T$ be the affine map from $\widehat{T}$ onto a triangle
$T$. Since $|\det B_T| = |T| / |\widehat{T}|$, the change of variables formula turns
@eq:quadrature-reference-triangle into

$$
\int_T f(x) \,\mathrm{d}x
= \frac{|T|}{|\widehat{T}|} \int_{\widehat{T}} f(F_T(\hat{x})) \,\mathrm{d}\hat{x}
\approx |T| \sum_{q=1}^{n_q} \omega_q \, f\big(F_T(\widehat{\zeta}_q)\big).
$$

Since $F_T$ maps vertices to vertices and preserves barycentric coordinates
({prf:ref}`lem:barycentric-affine-invariance`), the mapped points
$F_T(\widehat{\zeta}_q)$ have the same barycentric coordinates in $T$ as listed in
@tab:quadrature-rules-triangle, and the weights are unchanged.
As $F_T$ is affine, $f \circ F_T$ is a polynomial of the same degree as $f$, so the
mapped rule has the same degree of exactness as the reference rule. For
$\mathbb{P}_1$-elements, the gradients of the basis functions are constant and the
element stiffness matrix needs no quadrature beyond the rule of degree $1$. The entries of
the element mass matrix, $\int_T \phi_i \phi_j \,\mathrm{d}x$, are quadratic and are
integrated exactly by both rules of degree $2$. For the load vector,
$\int_T f \phi_i \,\mathrm{d}x$ is integrated exactly by a rule of degree $p+1$ whenever
$f$ is a polynomial of degree $p$.
```

```{prf:remark} Negative weights
:label: rem:quadrature-negative-weights

The four-point rule has the smallest number of points among the listed rules of degree
$3$, but its centroid weight $-\frac{9}{16}$ is negative. As a consequence, the
quadrature approximation of $\int_T f \,\mathrm{d}x$ can be negative for a function
$f \geqslant 0$, and a mass matrix assembled with this rule need not be positive
definite for higher-order elements. Moreover, the sum $\sum_q |\omega_q| = \frac{17}{8}$
exceeds $\sum_q \omega_q = 1$, so rounding errors in the values $f(\widehat{\zeta}_q)$
are amplified. Rules with positive weights and interior points are therefore preferred,
and the six-point rule, of degree $4$ with only two points more, is a common choice even
when degree $3$ would suffice.
```

```{exercise} Integrals of monomials over the reference triangle
:label: exer-monomial-integral-triangle

Show that for all integers $i, j \geqslant 0$, the monomial $\hat{x}^i \hat{y}^j$
satisfies @eq:monomial-integral-triangle.
```

```{admonition} Hint 1: Reduce to one dimension
:class: hint dropdown

For fixed $\hat{x} \in [0,1]$, the variable $\hat{y}$ runs from $0$ to $1 - \hat{x}$.
Carry out the inner integration over $\hat{y}$ to obtain

$$
\int_{\widehat{T}} \hat{x}^i \hat{y}^j \,\mathrm{d}\hat{x}
= \frac{1}{j+1} \int_0^1 \hat{x}^i (1-\hat{x})^{j+1} \,\mathrm{d}\hat{x}.
$$
```

```{admonition} Hint 2: Integrate by parts
:class: hint dropdown

Set $I(i,k) := \int_0^1 x^i (1-x)^k \,\mathrm{d}x$. For $k \geqslant 1$, integrate by
parts, differentiating $(1-x)^k$ and integrating $x^i$, to show that

$$
I(i,k) = \frac{k}{i+1}\, I(i+1, k-1).
$$

Each step moves one power from $(1-x)$ to $x$. Apply the recursion $k$ times, evaluate
$I(i+k, 0)$, and conclude that $I(i,k) = \frac{i!\, k!}{(i+k+1)!}$. Then combine this
with Hint 1 for $k = j+1$.
```
