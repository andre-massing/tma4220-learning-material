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

(ssec:sparse-storage-formats)=
## Storage formats for sparse matrices

So far, we have written the global matrix $M$ as if it were an ordinary
$N_p \times N_p$ array. For realistic meshes this is not an option, and how the matrix
is stored matters as much as how its entries are computed.

### Sparsity of finite element matrices

By construction, $M_{ij} = \int_\Omega \phi_j \phi_i \,\mathrm{d}x$ can only be non-zero
if the supports of $\phi_i$ and $\phi_j$ overlap in a set of positive area, that is, if
$a_i$ and $a_j$ are vertices of a common triangle. Row $i$ therefore contains one diagonal
entry and one entry for each vertex connected to $a_i$ by an edge. Counting the edges gives
the number of non-zero entries,

$$
\operatorname{nnz}(M) = N_p + 2 N_e,
$$ (eq:nnz-p1)

where $N_e$ is the number of edges of $\mathcal{T}_h$. For a triangulation of a simply
connected domain, Euler's formula $N_p - N_e + N_t = 1$ together with $3 N_t \approx 2 N_e$
(almost every edge is shared by two triangles) gives $N_e \approx 3 N_p$ and
$N_t \approx 2 N_p$ for large meshes. Hence

$$
\operatorname{nnz}(M) \approx 7 N_p,
\qquad
\frac{\operatorname{nnz}(M)}{N_p^2} \approx \frac{7}{N_p} \to 0
\quad \text{as } h \to 0 .
$$

Each row has a bounded number of non-zeros, independent of the mesh size, and the fraction
of non-zero entries, the **density**, tends to zero. A matrix with this property is called
**sparse**. The stiffness matrix of the Poisson problem has exactly the same sparsity
pattern, since it also couples only vertices that share a triangle.

Let us look at the mass matrix on the uniform mesh of the unit square with $N = 10$,
generated by `mt.unitsquare_mesh(10)` and assembled with {prf:ref}`alg:assemble-mass-2d`.
The mesh has $N_p = 121$ vertices and $N_t = 200$ triangles.

```{figure} figures/mass-matrix-sparsity-n10.png
:label: fig-mass-matrix-sparsity
:width: 60%

Sparsity pattern of the $\mathbb{P}_1$ mass matrix on the uniform mesh of the unit square
with $N = 10$. Each dot marks a non-zero entry.
```

The figure confirms @eq:nnz-p1: with $N_e = N_p + N_t - 1 = 320$ we get
$\operatorname{nnz}(M) = 121 + 2 \cdot 320 = 761$. The rows have between $3$ and $9$
non-zeros, $6.29$ on average. The mean is still below $7$ because a coarse mesh has many boundary vertices,
which have fewer neighbors. The pattern also reflects the **numbering** of the vertices. Since
`unitsquare_mesh` numbers them row by row, $N + 1$ vertices per row, the neighbors of
vertex $i$ have the indices $i \pm 1$ (horizontal edges), $i \pm (N+1)$ (vertical edges)
and $i \pm N$ or $i \pm (N+2)$ (diagonal edges). All non-zeros thus lie within a band of
width $N + 2$ around the diagonal. A different numbering of the same mesh gives a
different pattern, which is exploited by reordering algorithms in direct solvers.

The table below shows what is at stake when $M$ is stored as a dense array of
$8$-byte floating point numbers instead of storing only its non-zeros (in the CSR format
introduced below).

```{table} Memory required for the P1 mass matrix on the uniform unit square mesh.
:label: tab:sparse-vs-dense-memory

| $N$ | $N_p$ | $\operatorname{nnz}(M)$ | density | dense storage | sparse storage (CSR) |
|---:|---:|---:|---:|---:|---:|
| 10 | 121 | 761 | $5.2 \cdot 10^{-2}$ | 0.12 MB | 0.01 MB |
| 100 | 10 201 | 70 601 | $6.8 \cdot 10^{-4}$ | 832 MB | 0.89 MB |
| 300 | 90 601 | 631 801 | $7.7 \cdot 10^{-5}$ | 66 GB | 7.9 MB |
| 1000 | 1 002 001 | 7 006 001 | $7.0 \cdot 10^{-6}$ | 8.0 TB | 88 MB |
```

Dense storage costs $\mathcal{O}(N_p^2)$ memory, sparse storage $\mathcal{O}(N_p)$. The
same holds for the work: a matrix-vector product costs $\mathcal{O}(N_p^2)$ operations in
dense and $\mathcal{O}(\operatorname{nnz}) = \mathcal{O}(N_p)$ in sparse storage, and
Gaussian elimination on a dense matrix costs $\mathcal{O}(N_p^3)$, about $10^{18}$
operations for the last row of the table. Solvers that exploit the sparsity do much
better, see {prf:ref}`rem:sparse-solver-complexity`, but only if the matrix is stored in a
format that lets them skip the zeros.

:::{prf:remark} Outlook: how fast can sparse linear systems be solved?
:label: rem:sparse-solver-complexity

This remark is only meant to give an idea of what is possible. The analysis of sparse
solvers is beyond the scope of this course.

**Sparse direct solvers** compute an LU or Cholesky factorization of the matrix. Their cost
is governed by the **fill-in**: entries that are zero in the matrix but become non-zero in
the factors. The fill-in depends strongly on the order in which the unknowns are
eliminated. With the row-by-row numbering of `unitsquare_mesh`, for instance, the matrix
is a band matrix of bandwidth $\sim \sqrt{N_p}$, the band fills up completely, and the
factorization costs $\mathcal{O}(N_p^2)$ operations. Practical solvers such as SuperLU
(used by SciPy's `spsolve` and `splu`), UMFPACK, MUMPS or PARDISO therefore reorder the
unknowns first. For meshes like ours, the **nested dissection** ordering
{cite}`George1973` achieves

| | factorization | memory for the factors |
|---|---|---|
| 2D | $\mathcal{O}(N_p^{3/2})$ | $\mathcal{O}(N_p \log N_p)$ |
| 3D | $\mathcal{O}(N_p^{2})$ | $\mathcal{O}(N_p^{4/3})$ |

and no elimination order can do asymptotically better {cite}`HoffmanMartinRose1973`. Once
computed, the factors can be reused: every further right-hand side costs only two
triangular solves, that is, $\mathcal{O}(N_p \log N_p)$ operations in 2D.

**Iterative solvers** such as the conjugate gradient (CG) method only need matrix-vector
products, at $\mathcal{O}(N_p)$ operations each, but the number of iterations grows with
the condition number $\kappa$ of the matrix, roughly like $\sqrt{\kappa}$ for CG. For the
mass matrix on a quasi-uniform mesh, $\kappa$ is bounded independently of $h$, so CG
solves $M \mathbf{U} = \mathbf{b}$ in $\mathcal{O}(N_p)$ operations. For the stiffness
matrix of the Poisson problem, $\kappa \sim h^{-2}$, and plain CG needs
$\mathcal{O}(N_p^{3/2})$ operations in 2D. Multigrid methods, used as solvers or as
preconditioners for CG, reach the optimal complexity $\mathcal{O}(N_p)$.

As a rule of thumb, sparse direct solvers are robust and often the method of choice for
2D problems with up to about $10^6$ unknowns, while large 3D problems require
preconditioned iterative methods.
:::

### Common storage formats

There are several ways to store only the non-zero entries, each favoring different
operations. We illustrate the four most common ones with the small matrix

$$
A =
\begin{pmatrix}
1 & 0 & 2 & 0 \\
0 & 3 & 0 & 0 \\
4 & 5 & 6 & 0 \\
0 & 0 & 7 & 8
\end{pmatrix},
\qquad \operatorname{nnz}(A) = 8,
$$

and give the corresponding class in SciPy's
[`scipy.sparse`](https://docs.scipy.org/doc/scipy/reference/sparse.html) module. Indices
start at $0$ as in Python.

#### Coordinate format (COO)

The coordinate format stores one triplet $(i, j, a_{ij})$ per entry in three arrays `row`,
`col` and `data` of equal length. The triplets may come in **any order**, and the same
position $(i, j)$ may appear **several times**; the values are then summed when the
matrix is used or converted.

```{figure} figures/sparse-format-coo.svg
:label: fig-sparse-format-coo
:width: 90%

COO storage of $A$. The triplets are stored in no particular order, and the entry
$a_{22} = 6$ is stored as two triplets with the values $4$ and $2$ (orange).
```

This is exactly the output of the assembly loop: element $T_k$ produces the nine triplets
$\big(\iota(k,\alpha), \iota(k,\beta), M^{(k)}_{\alpha\beta}\big)$, and contributions of
neighboring elements to the same entry are summed automatically. The update
$M_{ij} \gets M_{ij} + M^{(k)}_{\alpha\beta}$ in {prf:ref}`alg:assemble-mass-2d` thus
becomes "append a triplet", and the triplet arrays of length $9 N_t$ can be allocated in
advance.

```python
A = sp.coo_array((data, (row, col)), shape=(n, n))
```

COO is the format of choice for **constructing** a matrix, but it supports neither
efficient access to a row or column nor a direct solve. Before computing with the matrix,
we convert it to CSR or CSC. With $32$-bit indices, COO needs $16$ bytes per triplet, and
since there are $9 N_t \approx 18 N_p$ triplets but only $\operatorname{nnz} \approx 7 N_p$
distinct entries, the triplet arrays are about $2.5$ times larger than the final matrix.

#### Compressed sparse row format (CSR)

CSR sorts the entries row by row and stores the row index only implicitly. The arrays
`data` and `indices` hold the values and column indices of the non-zeros of row $0$, then
row $1$, and so on, and the **row pointer** array `indptr` of length $n + 1$ records where
each row starts:

$$
\text{row } i \text{ occupies the positions } k = \texttt{indptr}[i], \ldots, \texttt{indptr}[i+1] - 1 .
$$

```{figure} figures/sparse-format-csr.svg
:label: fig-sparse-format-csr
:width: 90%

CSR storage of $A$. The colors mark the rows. `indptr[i]` points to the first entry of
row $i$, and `indptr[n] = nnz` marks the end.
```

Compared with COO, the array `row` of length $\operatorname{nnz}$ is replaced by `indptr`
of length $n + 1$, so CSR needs about $12$ bytes per non-zero. More importantly, all
entries of a row are stored contiguously, which makes **extracting rows** and computing
**matrix-vector products** fast:

:::{prf:algorithm} Matrix-vector product in CSR format
:label: alg:csr-matvec

**Input:** $A \in \mathbb{R}^{n \times m}$ in CSR format (`data`, `indices`, `indptr`),
vector $\mathbf{x} \in \mathbb{R}^m$

**Output:** $\mathbf{y} = A \mathbf{x}$

1. **for** $i = 0, \ldots, n - 1$ **do**
   1. $y_i \gets 0$
   2. **for** $k = \texttt{indptr}[i], \ldots, \texttt{indptr}[i+1] - 1$ **do**
      1. $y_i \gets y_i + \texttt{data}[k] \, x_{\texttt{indices}[k]}$
2. **return** $\mathbf{y}$
:::

The cost is $\mathcal{O}(\operatorname{nnz} + n)$. Iterative solvers such as the
conjugate gradient method access the matrix only through such products, so CSR is the
standard format for them. The price of the compression is that **changing the sparsity
pattern** is expensive: inserting a new non-zero shifts all subsequent entries of `data`
and `indices`. SciPy issues a `SparseEfficiencyWarning` when this happens.

```python
A = sp.csr_array(A_coo)   # or A_coo.tocsr(); duplicates are summed
```

#### Compressed sparse column format (CSC)

CSC is the same idea with the roles of rows and columns exchanged: the entries are sorted
column by column, `indices` holds row indices, and `indptr` points to the start of each
column.

```{figure} figures/sparse-format-csc.svg
:label: fig-sparse-format-csc
:width: 90%

CSC storage of $A$. The colors mark the columns.
```

The three CSC arrays of $A$ are the three CSR arrays of $A^\top$. For a symmetric matrix
such as $M$, both formats therefore store identical arrays. CSC is efficient for
**extracting columns**, and it is the format used internally by the sparse direct solver
SuperLU behind `scipy.sparse.linalg.spsolve` and `splu`. `spsolve` accepts CSR and CSC,
while `splu` converts any other input to CSC with a warning.

#### Dictionary of keys format (DOK)

DOK stores the matrix as a Python dictionary (a hash table) that maps each key $(i, j)$ to
the value $a_{ij}$.

```{figure} figures/sparse-format-dok.svg
:label: fig-sparse-format-dok
:width: 75%

DOK storage of $A$: a hash table with keys $(i, j)$, here listed in the order in which the
entries were inserted.
```

Reading, inserting or updating a single entry costs $\mathcal{O}(1)$ on average, so DOK
allows the assembly loop of {prf:ref}`alg:assemble-mass-2d` to be written literally,
`M[i, j] += M_k[alpha, beta]`. The convenience comes at a high price. Every entry is stored
as separate Python objects (a tuple of two integers as key, a float as value, plus the
hash table slot), roughly $200$ bytes per non-zero compared with $12$ to $16$ bytes in
CSR or COO, and every update runs through the Python interpreter. DOK is therefore useful
for small matrices and for modifying a few entries, but not for assembling large ones.
For computations it has to be converted to CSR or CSC as well.

The `scipy.sparse` module provides further formats, notably LIL (list of lists, similar in
spirit to DOK), BSR (block sparse rows, useful for vector-valued problems) and DIA
(diagonal storage, for matrices with few non-zero diagonals).

### Choosing a format

```{table} Comparison of sparse storage formats from the point of view of finite element computations ($n$ is the number of rows, $32$-bit indices).
:label: tab:sparse-formats

| Format | SciPy class | Memory | Assembly | Matrix-vector product | Solvers |
|---|---|---|---|---|---|
| COO | `coo_array` | $16 \operatorname{nnz}$ bytes | very good: append triplets, duplicates summed | possible, but convert first | convert first |
| CSR | `csr_array` | $12 \operatorname{nnz} + 4(n+1)$ bytes | poor: changing the pattern is costly | very good | iterative solvers; `spsolve` |
| CSC | `csc_array` | $12 \operatorname{nnz} + 4(n+1)$ bytes | poor: changing the pattern is costly | good | direct solvers (SuperLU: `spsolve`, `splu`) |
| DOK | `dok_array` | $\approx 200 \operatorname{nnz}$ bytes | simple but slow: entry-wise in Python | convert first | convert first |
```

(ssec:sparse-assembly-workflow)=
### A workflow for finite element matrices

In practice, the formats are combined: each one is used for what it does well.

1. **Assemble** in COO. Allocate `row`, `col` and `data` of length $9 N_t$ and fill them
   element by element, or compute all triplets at once with NumPy: `np.repeat(T, 3, axis=1)` and `np.tile(T, (1, 3))` list
   $\iota(k, \alpha)$ and $\iota(k, \beta)$ for all $k, \alpha, \beta$, and
   `area[:, None, None] * M_ref` is the stack of all element matrices $M^{(k)}$.
2. **Convert** once with `.tocsr()` (or `.tocsc()`), which sorts the entries and sums the
   duplicates.
3. **Solve** with a sparse solver from
   [`scipy.sparse.linalg`](https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html),
   for instance the direct solver `spsolve(M, b)` or the conjugate gradient method
   `cg(M, b)` for symmetric positive definite matrices.

@tab:assembly-timings compares this workflow with the entry-wise assembly into a
DOK matrix and into a dense array, both written as in {prf:ref}`alg:assemble-mass-2d`.

```{table} Assembly of the $\mathbb{P}_1$ mass matrix on the uniform unit square mesh with $N = 40$ ($N_p = 1681$, $\operatorname{nnz} = 11\,441$), measured on a laptop.
:label: tab:assembly-timings

| Assembly | Time | Memory of the result |
|---|---:|---:|
| element loop, adding into a DOK matrix | 0.21 s | 2.3 MB |
| element loop, adding into a dense array | 0.024 s | 22.6 MB |
| vectorized COO triplets, converted to CSR | 0.0007 s | 0.14 MB |
```

For this small mesh, the dense loop is faster than the DOK loop, as indexing a NumPy array
is much cheaper than updating a dictionary through SciPy. The dense array, however, already
takes up more than $100$ times as much memory as the CSR matrix, and this ratio grows linearly with
$N_p$. The vectorized COO assembly is both the fastest and the most memory-efficient
option.

```{prf:remark} Sparse arrays and sparse matrices in SciPy
:label: rem:scipy-sparse-array-vs-matrix

SciPy offers every format in two variants, for instance `csr_array` and the older
`csr_matrix`. The array variants behave like NumPy arrays and are recommended for new code.
The difference that matters most in practice is the operator `*`: it is the **matrix
product** for `csr_matrix` but the **entry-wise product** for `csr_array`. Writing `A @ x`
for matrix products avoids this pitfall with both, see the SciPy
[migration guide](https://docs.scipy.org/doc/scipy/reference/sparse.migration_to_sparray.html).
```

For more details on the formats, their conversion and the available operations, see the
SciPy [user guide on sparse arrays](https://docs.scipy.org/doc/scipy/tutorial/sparse.html)
and the reference pages of
[`coo_array`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_array.html),
[`csr_array`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_array.html),
[`csc_array`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_array.html) and
[`dok_array`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.dok_array.html).

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
