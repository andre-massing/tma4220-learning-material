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
    -u(b) &= g_N(b)
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
This will allow us to compute the contributions to $A$ and $\widetilde{b}$ element-wise and then assemble them into the global matrix and vector.  

More precisely, we compute the element entry $A_{ij} = a(\phi_j, \phi_i)$ by decomposing the integral into a sum of integrals over the elements $T$:

$$
a(\phi_j, \phi_i) = \sum_{l=0}^{N} a^{T}(\phi_j, \phi_i) = \sum_{l=0}^{N} \int_{T} ( \phi_j' \phi_i' + \phi_j \phi_i)
$$

Now, the integral $a^{T}(\phi_j, \phi_i)= \int_{T} (\phi_j' \phi_i' + \phi_j \phi_i)$ can only lead to a non-zero contribution if both $\phi_i$ and $\phi_j$ are non-zero on $T$.
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
- mesh $\mcT_h = \{T_l\}_{l=0}^{M}$ with $T_l = [x_l, x_{l+1}]$
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
