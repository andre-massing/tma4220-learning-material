---
jupytext:
  cell_metadata_filter: -jupyter
  formats: md:myst,ipynb
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.5
---

# Implementation of the finite element method in 1D

In this section, we discuss the implementation of the finite element
method in 1D.  We will strive to abstract some of the key ideas in
such a way that they can be reused and adapted when we will turn the
implementation of the finite element method in 2D and 3D.  That being
said, the code we will develop here is by no means optimized, but
rather presents a simple and straightforward implementation that is
easy to understand.  Possible adaption to improve efficiency and
performance will considered when we discuss the implementation of the
finite element method in 2D and 3D.

## Solving the Poisson problem with Neumann boundary conditions

As first guiding example, we will consider the
the Poisson problem with Neumann boundary conditions,

$$
\begin{aligned}
- u'' + u &=  f \quad \text{on } I 
\\
    -u'(a) &= g_N(a)
    \\
     u'(b) &= g_N(b)
\end{aligned}
$$ (eq:poisson-neumann-1d)

where we have added a low order term to the differential operator to ensure that the problem is well-posed.

As before, we introduce a mesh $\mathcal{T}_h = \{T_l\}_{l=0}^{M}$
where $T_l = [x_l, x_{l+1}]$ and $h = \max_{T\in \mathcal{T}_h} h_T$ with $h_{T_l} = x_{l+1}-x_l$.
For simplicity, we will often assume that the mesh is uniform, i.e. $h_T = h$ for all $T \in \mathcal{T}_h$.
We discretize @eq:poisson-neumann-1d using
piecewise polynomials of order $k$; that is,
we seek to find $u_h \in V_{h} = \mathbb{P}_k^c(\mathcal{T}_h)$
such that

$$
a(u_h, v_h)
= l(v_h) \quad \forall v_h \in V_{h,0}
$$

where 

$$
\begin{aligned}
a(u_h, v_h) 
&=\int_a^b (u_h' v_h' + u_h v_h),
\\
l(v_h)
&=\int_a^b f v_h + g_N(a)v_h(a) + g_N(b)v_h(b)
\end{aligned}
$$

Now, let $\{\phi_i\}_{i=0}^{N+1}$ be the set of global basis functions we constructed
from the element-wise Lagrange cardinal basis function in 
[previous section](#ssec:finite-element-1d-method), where $\phi_0$ and $\phi_{N+1}$ are the basis functions associated with the boundary nodes $x_0 = a$ and $x_{N+1} = b$, respectively. 
Then clearly $V_h = \mathbb{P}_k^c(\mathcal{T}_h) = \text{span}\{\phi_i\}_{i=0}^{N+1}$.

Finding the discrete solution $u_h = \sum_{i=0}^{N+1} U_i \phi_i$ amounts
to solve a linear algebra system

$$
\mathcal{A} \mathbf{U} = \mathbf{b}
$$

where $\mathbf{U} = (U_i)_{i=0}^{N+1}\in \mathbb{R}^{N+1}$ is the 
**coefficient vector** associated with $u_h$, and
the elements of the matrix $\mathcal{A}\in \mathbb{R}^{(N+2)\times(N+2)}$ 
(often also referred to as the **stiffness matrix**)
are given by $a_{ij} = a(\phi_j, \phi_i), i,j= 0,\ldots,N+1$.
Similarly, the components of the 
vector $\mathbf{b} \in \mathbb{R}^{N+2}$ (the **load vector**)
are given by $b_i = l(\phi_i), i = 0,\ldots, N+1$.

The process of computing the matrix $\mathcal{A}$ and the vector $\mathbf{b}$ is called **assembly**. 
An efficient implementation of the assembly process exploits the fact that the basis functions $\phi_i$ have local support, and on a given element $T$, there are only a few basis functions that do not vanish.
This will allow us to compute the contributions to $A$ and $\mathbf{b}$ element-wise and then assemble them into the global matrix and vector.  

### Assembly of the stiffness matrix
The first step is to break-up the computation of the element entry $A_{ij} = a(\phi_j, \phi_i)$ by decomposing the integral into a sum of integrals over the elements $T$:

$$
a(\phi_j, \phi_i) = \sum_{l=0}^{N} a^{T}(\phi_j, \phi_i) = \sum_{l=0}^{N} \int_{T} ( \phi_j' \phi_i' + \phi_j \phi_i)
$$

Now, the integral $a^{T}(\phi_j, \phi_i)= \int_{T} (\phi_j' \phi_i' + \phi_j \phi_i)$ **can only lead to a non-zero contribution if both $\phi_i$ and $\phi_j$ are non-zero on $T$**.
Further, we know that $\phi_i$ is non-zero on $T$ only if its restriction to $T$
corresponds to a Lagrange cardinal basis function $\ell_\alpha^{(T)}$
spanning the space of polynomials of order $k$
on $T$, i.e. if

$$
\phi_i|_{T} = \ell_\alpha^{(T)} \text{ for some }\alpha \in \{0, \ldots, k\}.
$$

Note that the index $\alpha$ is a **local index** on the element $T$,
which only depends on the chosen polynomial order $k$ and therefore
only runs from $0$ to $k$ (independent of the number of elements $M$),
while the index $i$ is a **global index** which runs from $1$ to $N$
which in turn depends on both the chosen polynomial order $k$
and the number of elements $M$.
Thus to identify which local index $\alpha$ corresponds to a given global index $i$, we need to introduce a so-called **local-to-global** mapping $i = \iota(\alpha, l)$ which maps the local index $\alpha$
on the element $T$ to the global index $i$. 

Then all non-zero contributions to the matrix $A$ from the element $T$ can be encoded in the local 
(stiffness) matrix $A^{(T)}$ associated with the element $T$ which is defined as follows:

$$
A^{(T)}_{\alpha \beta} = a^{T}(\ell_\beta^{(T)}, \ell_\alpha^{(T)}) 
= 
\int_{T} 
\Big(
  (\ell_\beta^{(T)})' (\ell_\alpha^{(T)})'
+ \ell_\beta^{(T)} \ell_\alpha^{(T)} 
\Big)
\quad \forall \alpha, \beta = 0, \ldots, k
$$ (eq:local-stiffness-matrix)

As we realize, the computation of the local stiffness matrix $A^{(T)}$
leads for each element $T$ to very similar integrals, and it might be
worthwhile to try to precompute parts of these integrals a priori.
More precisely, we introduce a reference element $\widehat{T} = [0,1]$
and define the reference Lagrange cardinal basis functions
$\widehat{\ell}_\alpha(\widehat{x})$ on $\widehat{T}$ for a given
polynomial order $k$. For each "physical" element $T = [x_l,
x_{l+1}]$, we can then define an affine mapping $\Phi_l: \widehat{T}
\to T$ which maps the reference element to the physical element $T$,
with corresponding inverse mapping $\Phi_T^{-1}: T \to \widehat{T}$:

$$
\begin{aligned}
\Phi_T(\widehat{x}) &= x(\widehat{x}) = x_l + h_l \widehat{x} \\
\Phi_T^{-1}(x) &=  \widehat{x}(x) = \frac{x - x_l}{h_l}
\end{aligned}
$$

Next, we observe that the Lagrange cardinal basis functions on the physical element $T$ can be expressed in terms of the reference basis functions on $\widehat{T}$, i.e. $\ell_\alpha^{(T)}(x) = \widehat{\ell}_\alpha(\Phi_T^{-1}(x))$ for all $x \in T$. So
$\ell_\alpha^{(T)}(x) = \widehat{\ell}_\alpha(\widehat{x})$ with $\widehat{x} = \Phi_T^{-1}(x)$.
Using the chain rule, we can compute the derivative of $\ell_\alpha^{(T)}$ with respect to $x$ as follows:

$$
\begin{aligned}
\frac{d}{dx} \ell_\alpha^{(T)}(x)
&= \frac{d}{dx} \widehat{\ell}_\alpha(\widehat{x}(x)) \\
&= \frac{d}{d\widehat{x}} \widehat{\ell}_\alpha(\widehat{x}) \cdot \frac{d}{dx} \widehat{x}(x) \\
&= \frac{1}{h_l} \frac{d}{d\widehat{x}} \widehat{\ell}_\alpha(\widehat{x})
\end{aligned} 
$$

Similarly, a change of variables in the integral, we can also express the integral over $T$ in terms of an integral over $\widehat{T}$ via the change of the integral measure

$$
\mathrm{d}x = |(\Phi_T^{-1})'(\widehat{x})| \,\mathrm{d}\widehat{x} = h_T \,\mathrm{d}\widehat{x}.
$$

Combining these two observations, we can express the integral in @eq:local-stiffness-matrix in terms of an integral over the reference element $\widehat{T}$ as follows:

$$
\begin{aligned}
a^T(\ell_{\beta}^{(T)}, \ell_{\alpha}^{(T)})
= &\int_{T} 
\Big(
  (\ell_\beta^{(T)})'(x) (\ell_\alpha^{(T)})'(x)
+ \ell_\beta^{(T)}(x) \ell_\alpha^{(T)}(x)
\Big)
 \, \mathrm{d}x
 \\
& =
\int_{\widehat{T}}
\Big( \frac{1}{h_T} (\widehat{\ell}_\beta)'(\widehat{x}) \cdot \frac{1}{h_T} (\widehat{\ell}_\alpha)'(\widehat{x}) 
+
\widehat{\ell}_\beta(\widehat{x}) \widehat{\ell}_\alpha(\widehat{x}) 
\Big)
\cdot h_T \,\mathrm{d}\widehat{x}
\\
&=
\frac{1}{h_T} \int_{\widehat{T}} (\widehat{\ell}_\beta)'(\widehat{x}) \cdot (\widehat{\ell}_\alpha)'(\widehat{x})  \,\mathrm{d}\widehat{x}
+
h_T \int_{\widehat{T}} \widehat{\ell}_\beta(\widehat{x}) \cdot \widehat{\ell}_\alpha(\widehat{x})  \,\mathrm{d}\widehat{x}
% \quad \forall \alpha, \beta = 0, \ldots, k
\\
&\eqqcolon
\dfrac{1}{h_T} \widehat{A}_{\alpha\beta} + h_T \widehat{M}_{\alpha\beta}
\end{aligned}
$$ (eq:local-stiffness-matrix-ref)

Except for the factors $1/h_T$ and $h_T$,
the integrals on the right-hand side of @eq:local-stiffness-matrix-ref is independent of the element $T$ and can be precomputed once and for all for a given polynomial order $k$.
Since the two integrals contributions scale differently with $h_T$, we need to compute both integrals separately, and we collect their contributions in the matrices $\widehat{A}$ and $\widehat{M}$. 

```{prf:remark}
This will be a recurring theme in the implementation of finite element methods, where we will
use a change of variables to represent integrals on physical elements in terms of *precomputable integrals* on a reference element and geometric factors, that can computed rather cheaply
on a per-element basis.
```

```{prf:remark}
In the considerations above, we also used the fact that the Lagrange cardinal basis functions on the physical element $T$ can be expressed in terms of Lagrange cardinal basis functions on the reference element $\widehat{T}$, i.e. $\ell_\alpha^{(T)}(x) = \widehat{\ell}_\alpha(\Phi_T^{-1}(x))$ for all $x \in T$. This is a special case of so-called **affine-equivalence** of finite elements,
where the basis functions on a physical element can be generated by composing the basis functions on a reference element with an affine mapping. 
```

:::{prf:algorithm} Assembly of the stiffness matrix
:label: alg:assemble-stiffness-1d

**Input:** 
- mesh vertices $x_0 < x_1 < \dots < x_{M+1}$
- mesh $\mathcal{T}_h = \{T_l\}_{l=0}^{M}$ with $T_l = [x_l, x_{l+1}]$
- polynomial order $k$
- local-to-global map $\iota(\alpha, l)$

**Output:** Stiffness matrix $\mathcal{A} \in \mathbb{R}^{(N+2) \times (N+2)}$ with $A_{ij} = a(\phi_j, \phi_i)$ \
for $i, j = 0, \ldots, N+1$.

1. $\mathcal{A} \gets 0$
2. **for** $\alpha, \beta = 0, \ldots, k$ **do**
   1. $\widehat{A}_{\alpha\beta} \gets \int_{\widehat{T}} (\widehat{\ell}_\beta)' \, (\widehat{\ell}_\alpha)' \,\mathrm{d}\widehat{x}$
   2. $\widehat{M}_{\alpha\beta} \gets \int_{\widehat{T}} \widehat{\ell}_\beta \, \widehat{\ell}_\alpha \,\mathrm{d}\widehat{x}$
3. **for** $l = 0, \ldots, M$ **do**
   1. $h_l \gets x_{l+1} - x_l$
   2. $A^{(T_l)} \gets \frac{1}{h_l} \widehat{A}$ + $h_l \widehat{M}$
   3. **for** $\alpha, \beta = 0, \ldots, k$ **do**
      1. $i \gets \iota(\alpha, l)$, $\; j \gets \iota(\beta, l)$
      2. $A_{ij} \gets A_{ij} + A^{(T_l)}_{\alpha\beta}$
4. **return** $A$
:::


### Assembly of the load vector
Next, we consider the assembly of the load vector $\mathbf{b}$.
Again, we can decompose the integral into a sum of integrals over the elements $T$:

$$
b_i = \sum_{T \in \mathcal{T}_h} \int_T f \phi_i \,\mathrm{d}x + 
g_N(a) \phi_i(a) + g_N(b) \phi_i(b).
$$

In contrast to stiffness matrix assembly, where we could precompute the integrals of the basis functions **exactly** on the reference element,
the integral contribution $\int_T f \phi_i \,\mathrm{d}x$ cannot be computed exactly in general
for an arbitrary function $f$. Instead, the integral is approximated using suitable quadrature rules $\{(\zeta_q^T, w_q^T)\}_{q=1}^{d}$ on each element $T$:

$$
\int_T f \phi_i \,\mathrm{d}x \approx \sum_{q=1}^{d} w_q^T f(\zeta_q^T) \phi_i(\zeta_q^T)
$$
As before, we only need to consider the elements $T$ for which $\phi_i$ is non-zero, i.e. the elements $T$ for which $\phi_i|_T = \ell_\alpha^{(T)}$ for some $\alpha = 0, \ldots, k$.
Then again, we can use a change of variables to express the integral on the physical element $T$ in terms of an integral on the reference element $\widehat{T}$  and the corresponding quadrature rule $\{(\widehat{\zeta}_q, \widehat{w}_q)\}_{q=1}^{d}$ on $\widehat{T}$:

$$
\begin{aligned}
\int_T f(x) \ell_{\alpha}^{(T)}(x) \,\mathrm{d}x
&=
\int_{\widehat{T}} f (\Phi_T(\widehat{x})) {\ell}_{\alpha}^T(\Phi_T(\widehat{x})) |(\Phi_T^{-1})'(\widehat{x})| \,\mathrm{d}\widehat{x}
\\
&=\int_{\widehat{T}} \widehat{f} (\widehat{x}) {\ell}_{\alpha}(\widehat{x}) h_T \,\mathrm{d}\widehat{x}
\\
&\approx h_T \sum_{q=1}^{d} \widehat{w}_q \widehat{f}(\widehat{\zeta}_q) {\ell}_{\alpha}(\widehat{\zeta}_q)
\end{aligned}
$$
where $(\{\widehat{\zeta}_q, \widehat{w}_q\})_{q=1}^{d}$ are the quadrature points and weights on the reference element $\widehat{T}$, and $\widehat{f}(\widehat{x}) = f(\Phi_T(\widehat{x}))$ is the function $f$ evaluated at the quadrature points mapped to the physical element $T$.

Note that the products $\widehat{w}_q \widehat{\ell}_\alpha(\widehat{\zeta}_q)$ are again independent of
the element $T$, so they can be precomputed once for a given polynomial order $k$ and quadrature rule,
exactly as the reference matrices $\widehat{A}$ and $\widehat{M}$ were.
Only the evaluation of $f$ at the mapped quadrature points and the scaling by $h_T$ remain
to be done element by element.
Collecting these observations, we arrive at the following algorithm.

:::{prf:algorithm} Assembly of the load vector
:label: alg:assemble-load-1d

**Input:** 
- mesh vertices $x_0 < x_1 < \dots < x_{M+1}$
- mesh $\mathcal{T}_h = \{T_l\}_{l=0}^{M}$ with $T_l = [x_l, x_{l+1}]$
- polynomial order $k$
- local-to-global map $\iota(\alpha, l)$
- right-hand side $f$ and Neumann data $g_N(a)$, $g_N(b)$
- quadrature rule $\{(\widehat{\zeta}_q, \widehat{w}_q)\}_{q=1}^{d}$ on $\widehat{T} = [0,1]$

**Output:** Load vector $\mathbf{b} \in \mathbb{R}^{N+2}$ with $b_i = l(\phi_i)$ \
for $i = 0, \ldots, N+1$.

1. $\mathbf{b} \gets 0$
2. **for** $\alpha = 0, \ldots, k$, $\; q = 1, \ldots, d$ **do**
   1. $\widehat{L}_{\alpha q} \gets \widehat{w}_q \, \widehat{\ell}_\alpha(\widehat{\zeta}_q)$
3. **for** $l = 0, \ldots, M$ **do**
   1. $h_l \gets x_{l+1} - x_l$
   2. **for** $q = 1, \ldots, d$ **do**
      1. $\widehat{f}_q \gets f(\Phi_{T_l}(\widehat{\zeta}_q)) = f(x_l + h_l \widehat{\zeta}_q)$
   3. **for** $\alpha = 0, \ldots, k$ **do**
      1. $b^{(T_l)}_\alpha \gets h_l \sum_{q=1}^{d} \widehat{L}_{\alpha q} \, \widehat{f}_q$
      2. $i \gets \iota(\alpha, l)$
      3. $b_i \gets b_i + b^{(T_l)}_\alpha$
4. *Neumann boundary contributions,* using $\phi_{\iota(0,0)}(a) = 1$ and $\phi_{\iota(k,M)}(b) = 1$,
   while all other basis functions vanish at $a$ and $b$:
   1. $b_{\iota(0,0)} \gets b_{\iota(0,0)} + g_N(a)$
   2. $b_{\iota(k,M)} \gets b_{\iota(k,M)} + g_N(b)$
5. **return** $\mathbf{b}$
:::

For convenience, we recall some popular quadrature rules on the reference element $\widehat{T} = [0,1]$ in the following table.

```{table} Quadrature rules $\{(\widehat{\zeta}_q, \widehat{w}_q)\}_{q=1}^{d}$ on the reference element $\widehat{T} = [0,1]$, together with the polynomial degree they integrate exactly.
:label: tab:quadrature-rules-ref

| Rule | $d$ | Points $\widehat{\zeta}_q$ | Weights $\widehat{w}_q$ | Exact for degree |
|:--|:-:|:--|:--|:-:|
| Midpoint rule | $1$ | $\frac{1}{2}$ | $1$ | $1$ |
| Trapezoidal rule | $2$ | $0,\; 1$ | $\frac{1}{2},\; \frac{1}{2}$ | $1$ |
| Simpson's rule | $3$ | $0,\; \frac{1}{2},\; 1$ | $\frac{1}{6},\; \frac{2}{3},\; \frac{1}{6}$ | $3$ |
| Gauss–Legendre, $2$ points | $2$ | $\frac{1}{2} - \frac{1}{2\sqrt{3}},\; \frac{1}{2} + \frac{1}{2\sqrt{3}}$ | $\frac{1}{2},\; \frac{1}{2}$ | $3$ |
| Gauss–Legendre, $3$ points | $3$ | $\frac{1}{2} - \frac{1}{2}\sqrt{\frac{3}{5}},\; \frac{1}{2},\; \frac{1}{2} + \frac{1}{2}\sqrt{\frac{3}{5}}$ | $\frac{5}{18},\; \frac{4}{9},\; \frac{5}{18}$ | $5$ |
```

```{prf:remark}
Note that the weights are normalized such that $\sum_{q=1}^{d} \widehat{w}_q = 1 = |\widehat{T}|$;
the factor $h_T$ accounting for the size of the physical element $T$ appears separately
in the assembly, see step 3.3.1 of @alg:assemble-load-1d.
When choosing a rule, recall that the integrand $\widehat{f} \widehat{\ell}_\alpha$ contains the
basis function of degree $k$, so a rule which is exact for degree $p + k$ integrates the
element contribution exactly whenever $\widehat{f}$ is a polynomial of degree $p$.
```

### Example: Realization with linear elements

We now carry out the assembly explicitly for linear elements, $k = 1$, on a uniform mesh
with $h_{T_l} = h$ for all $l = 0, \ldots, M$. The reference basis functions on $\widehat{T} = [0,1]$ and their derivatives are

\begin{equation}
\widehat{\ell}_0(\widehat{x}) = 1 - \widehat{x}, \quad
\widehat{\ell}_1(\widehat{x}) = \widehat{x},
\qquad
(\widehat{\ell}_0)' = -1, \quad (\widehat{\ell}_1)' = 1 ,
\end{equation}

and the local-to-global map is simply $\iota(\alpha, l) = l + \alpha$. Since each element
carries $k+1 = 2$ basis functions, all local matrices are of size $2 \times 2$, and the
global system has $N + 2 = M + 2$ unknowns, one per mesh vertex.

#### The reference stiffness matrix

The derivatives are constant, so the first integral in @eq:local-stiffness-matrix-ref is
evaluated without effort:

$$
\widehat{A}_{\alpha\beta}
= \int_0^1 (\widehat{\ell}_\beta)'(\widehat{\ell}_\alpha)' \,\mathrm{d}\widehat{x}
\quad \Longrightarrow \quad
\widehat{A} =
\begin{pmatrix}
1 & -1 \\
-1 & 1
\end{pmatrix} .
$$ (eq:ref-stiffness-p1)

Each element $T_l$ contributes $\frac{1}{h}\widehat{A}$ to the rows and columns $l$ and $l+1$.
Every interior vertex $x_i$ is shared by the two elements $T_{i-1}$ and $T_i$ and therefore
receives the diagonal entry twice, whereas the two boundary vertices belong to one element only.
Collecting the contributions of all elements, the stiffness part of $\mathcal{A}$ is the
tridiagonal matrix

$$
\frac{1}{h}
\begin{pmatrix}
 1 & -1 &        &    &    \\
-1 &  2 & -1     &    &    \\
   & \ddots & \ddots & \ddots &  \\
   &    & -1     &  2 & -1 \\
   &    &        & -1 &  1
\end{pmatrix}
\in \mathbb{R}^{(N+2) \times (N+2)} .
$$ (eq:global-stiffness-p1)

#### The reference mass matrix

The second integral in @eq:local-stiffness-matrix-ref gives, with
$\int_0^1 (1-\widehat{x})^2 = \int_0^1 \widehat{x}^2 = \frac{1}{3}$ and
$\int_0^1 \widehat{x}(1-\widehat{x}) = \frac{1}{6}$,

$$
\widehat{M}_{\alpha\beta}
= \int_0^1 \widehat{\ell}_\beta \widehat{\ell}_\alpha \,\mathrm{d}\widehat{x}
\quad \Longrightarrow \quad
\widehat{M} =
\begin{pmatrix}
\frac{1}{3} & \frac{1}{6} \\[2pt]
\frac{1}{6} & \frac{1}{3}
\end{pmatrix}
= \frac{1}{6}
\begin{pmatrix}
2 & 1 \\
1 & 2
\end{pmatrix} .
$$ (eq:ref-mass-p1)

Assembled in the same way, the mass part $h \widehat{M}$ yields

$$
\frac{h}{6}
\begin{pmatrix}
 2 & 1 &        &    &    \\
 1 & 4 & 1      &    &    \\
   & \ddots & \ddots & \ddots &  \\
   &    & 1      &  4 & 1 \\
   &    &        & 1 &  2
\end{pmatrix} ,
$$

so that the local matrix on each element is
$A^{(T_l)} = \frac{1}{h}\widehat{A} + h \widehat{M}$ and the full stiffness matrix
$\mathcal{A}$ is the sum of the two tridiagonal matrices above.

#### The load vector with the trapezoidal rule

The trapezoidal rule from @tab:quadrature-rules-ref uses the endpoints of $\widehat{T}$, where
the reference basis functions satisfy $\widehat{\ell}_\alpha(\widehat{\zeta}_q) = \delta_{\alpha q}$.
Only one term survives in each local contribution:

$$
b^{(T_l)}_0 = \frac{h}{2} f(x_l),
\qquad
b^{(T_l)}_1 = \frac{h}{2} f(x_{l+1}) .
$$

Adding the two contributions at every interior vertex and the Neumann data at the two boundary
vertices, the load vector becomes

$$
\mathbf{b}
= h
\begin{pmatrix}
\frac{1}{2} f(x_0) \\
f(x_1) \\
\vdots \\
f(x_{N}) \\
\frac{1}{2} f(x_{N+1})
\end{pmatrix}
+
\begin{pmatrix}
g_N(a) \\ 0 \\ \vdots \\ 0 \\ g_N(b)
\end{pmatrix} .
$$

This is the load vector one would also obtain from a finite difference discretization.
Note that the integrand $f \ell_\alpha^{(T)}$ is of degree $p + 1$ for a polynomial $f$ of
degree $p$, so this choice integrates the element contributions exactly only for constant $f$.

#### The load vector with Simpson's rule

Simpson's rule adds the midpoint $x_{l+1/2} = \frac{1}{2}(x_l + x_{l+1})$, where both reference
basis functions take the value $\frac{1}{2}$. With the weights
$\frac{1}{6}, \frac{2}{3}, \frac{1}{6}$ we obtain

$$
b^{(T_l)}_0 = \frac{h}{3}\Big( \frac{1}{2} f(x_l) + f(x_{l+1/2}) \Big),
\qquad
b^{(T_l)}_1 = \frac{h}{3}\Big( f(x_{l+1/2}) + \frac{1}{2} f(x_{l+1}) \Big) ,
$$

and therefore

$$
\mathbf{b}
= \frac{h}{3}
\begin{pmatrix}
\frac{1}{2} f(x_0) + f(x_{1/2}) \\
f(x_{1/2}) + f(x_1) + f(x_{3/2}) \\
\vdots \\
f(x_{N-1/2}) + f(x_{N}) + f(x_{N+1/2}) \\
f(x_{N+1/2}) + \frac{1}{2} f(x_{N+1})
\end{pmatrix}
+
\begin{pmatrix}
g_N(a) \\ 0 \\ \vdots \\ 0 \\ g_N(b)
\end{pmatrix} .
$$

Now the midpoint values of $f$ enter as well, and since Simpson's rule is exact for degree $3$,
the element contributions are computed exactly for every $f$ of degree $p \leqslant 2$.
