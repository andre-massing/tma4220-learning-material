(sec:fem-general)=
# The finite element method in multiple dimensions

So far, we have constructed finite element spaces on an interval by patching
together local Lagrange interpolation polynomials. We now take a step back and
formulate what a finite element *is*, in a way which is independent of the space
dimension and of the particular polynomial space used. This abstraction is what
allows the same code, and the same theory, to serve the plethora of finite element
families in use today.

In general, the construction of a global finite element space on a domain $\Omega$ proceeds in three steps:

1. Provide a subdivision of the domain $\Omega$ into a finite number of non-overlapping subdomains, called **elements**. The collection of elements is called a **mesh** or **triangulation**.

2. On each element, define a finite-dimensional space of functions, called the **local finite element space**. 

3. Define a global finite element space by patching together the local spaces, in a way that ensures that the resulting functions are globally continuous (or satisfy other desired global compatibility conditions).

```{figure} figures/fem-construction-steps.svg
:label: fig-fem-construction-steps
:alt: Four panels: an L-shaped domain subdivided into triangles with one element highlighted, that element enlarged with its three degrees of freedom, the mesh with all nodes marked and the support of one basis function shaded, and a three-dimensional plot of that basis function as a pyramid over the mesh.
:width: 90%

The three steps, and what they produce. Left: the domain $\Omega$ is subdivided into elements $T$, here
triangles, one of which is highlighted. Middle: on that element a local space with its
degrees of freedom is chosen, here $\mathbb{P}_1(T)$ with the values at the three
vertices. Third: the local spaces are patched together into a global space $V_h$, whose
degrees of freedom are the values at all mesh nodes; the shaded elements form the support
of the global basis function belonging to the ringed node. Right: that basis function
itself, the piecewise linear "hat" which equals $1$ at the ringed node, vanishes at every
other node, and is supported on the shaded patch.
```

## Construction of local finite element spaces
First, we introduce a powerful generalization of the finite elements
according to Ciarlet, which extracts the main components we  have
already encounter in the case of Lagrange element in one dimension.
Recall that in 1D, we constructed finite element basis functions 
by glueing together the Lagrange interpolation basis/cardinal functions
for polynomials of order $k$ on each subinterval $T$ of a given interval
$\Omega = (a,b)$. More precisely, starting
from a mesh element For $T = (x_l, x_{l+1})$, we introduced $k+1$
interpolation nodes $x_l = \xi_0 < \xi_1 < \ldots < \xi_k = x_{l+1}$,
together with their associated basis functions/cardinal functions 
$\lambda_i(x) \in \mathbb{P}_k(T)$ satisfying $\lambda_i(\xi_j) = \delta_{ij}$. Thanks to the last property, we can write any polynomial $p \in \mathbb{P}_k(T)$ as 

$$
p(x) = \sum_{i=0}^k p(\xi_i) \lambda_i(x)
$$

and correspondingly, the interpolation polynomial $I_T f \in \mathbb{P}_k(T)$ for $f \in
C(\overline{T})$ can be calculated via

$$
I_T f (x) = \sum_{i=0}^k f(\xi_i) \lambda_i(x)
$$

Note that the evaluations and $p(\xi)$ and $f(\xi_i)$ can be abstractly 
written as linear functionals acting on $p$ and $f$. More precisely,
we define the linear functionals $\{\sigma_i\}_{i=0}^k$ associated
with the nodes $\{\xi_i\}_{i=0}^k$ by

$$
\mathbb{P}_{k}(T) \ni p \mapsto \sigma(p) := p(\xi_i)
$$

and similar for $f \in C(\overline{T})$.
Then $p$ and $I_T f$ can be written as

$$
p(x) = \sum_{i=0}^k \sigma_i(p) \lambda_i(x),
\quad
I_T f (x) = \sum_{i=0}^k \sigma_i(f) \lambda_i(x)
$$

This illustrates already all the important ideas of general finite
element construction: We will need a **mesh element** $T$, and
**finite dimensional function space** $\mathcal{P}$ of dimension
$n_{sh}$, and a set of linear functionals $\{\sigma\}_{j=1}^{n_{sh}}$
(**the degree of freedoms**) which determines functions $p \in
\mathcal{P}$ in a sense that will be made more precise further down
below.  From this the set of linear functionals, we will construct a
set of basis functions (the **nodal basis** or **shape functions**
$\{\lambda_{i}\}_{i=1}^{n_{sh}}$ for $\mathcal{P}$, which satisfy
$\sigma_{j}(\lambda_i) = \delta_{ij}$.
With this picture in mind, it is time to get started with the precise mathematical definitions.


```{prf:definition} Ciarlet's construction
:label: def:ciarlet-finite-element

A **finite element** consists of a triple $(T, \mathcal{P}, \Sigma)$ where

i) $T \subseteq \mathbb{R}^d$ is a bounded, closed domain, the **element**, such that
   the interior $\mathring{T}$ is non-empty and $T$ is a $C^1$ polyhedron;

ii) $\mathcal{P}$ is a non-trivial, finite-dimensional vector space of functions
    $p : T \to \mathbb{R}^q$, the space of **shape functions**, for some integer
    $q \geqslant 1$, typically $q \in \{1, d\}$;

iii) $\Sigma = \{\sigma_1, \ldots, \sigma_{n_{\mathrm{sh}}}\}$ is a set of linear
     functionals on $\mathcal{P}$ which is a basis for the dual space $\mathcal{P}'$.
     The linear forms $\sigma_i$ are called **degrees of freedom**, or **nodal
     variables**.
```

```{prf:remark}
:label: rem:dim-shape-functions

Since $\mathcal{P}$ is a finite-dimensional vector space,
$\dim \mathcal{P} = \dim \mathcal{P}' = n_{\mathrm{sh}}$.
```

Often, we will need to check whether a given set of linear functionals actually
constitutes a basis of $\mathcal{P}'$. The following lemma provides two criteria
which are easier to verify than the definition.

```{prf:lemma} Unisolvence
:label: lem:unisolvence

Let $\mathcal{P}$ be a finite-dimensional space of dimension $n_{\mathrm{sh}}$ and let
$\Sigma = \{\sigma_1, \ldots, \sigma_{n_{\mathrm{sh}}}\} \subseteq \mathcal{P}'$ be a
set of linear functionals. Then the following statements are equivalent.

a) $\Sigma$ is a basis of $\mathcal{P}'$.

b) For all $p \in \mathcal{P}$: if $\sigma_i(p) = 0$ for $i = 1, \ldots, n_{\mathrm{sh}}$,
   then $p \equiv 0$.

c) The mapping
   $\Phi_\Sigma : \mathcal{P} \to \mathbb{R}^{n_{\mathrm{sh}}}$ given by
   $\Phi_\Sigma(p) = (\sigma_i(p))_{i=1}^{n_{\mathrm{sh}}}$ is an isomorphism.
```

```{prf:proof}
Let $\{\varphi_i\}_{i=1}^{n_{\mathrm{sh}}}$ be a basis of $\mathcal{P}$.

**Reformulation of a).** If $\Sigma$ is a basis of $\mathcal{P}'$, then for every
$l \in \mathcal{P}'$ there is a unique
$(\alpha_1, \ldots, \alpha_{n_{\mathrm{sh}}}) \in \mathbb{R}^{n_{\mathrm{sh}}}$ such that

$$
l = \sum_{j=1}^{n_{\mathrm{sh}}} \alpha_j \sigma_j .
$$ (eq:dual-basis-expansion)

Since $l$ is a linear functional, {eq}`eq:dual-basis-expansion` is equivalent to

$$
b_i := l(\varphi_i)
= \sum_{j=1}^{n_{\mathrm{sh}}} \alpha_j \underbrace{\sigma_j(\varphi_i)}_{=: A_{ij}} ,
\qquad i = 1, \ldots, n_{\mathrm{sh}} ,
$$

or equivalently, that for any given $\mathbf{b} \in \mathbb{R}^{n_{\mathrm{sh}}}$ there
is a unique $\boldsymbol{\alpha} \in \mathbb{R}^{n_{\mathrm{sh}}}$ with
$A \boldsymbol{\alpha} = \mathbf{b}$. This in turn is equivalent to
$A \in \mathbb{R}^{n_{\mathrm{sh}} \times n_{\mathrm{sh}}}$ being invertible.

**Reformulation of b).** For a given $p \in \mathcal{P}$, write
$p = \sum_{i=1}^{n_{\mathrm{sh}}} \beta_i \varphi_i$. Then the assumption in b) reads

$$
0 = \sigma_j(p)
= \sum_{i=1}^{n_{\mathrm{sh}}} \beta_i \underbrace{\sigma_j(\varphi_i)}_{=: B_{ji}} ,
\qquad j = 1, \ldots, n_{\mathrm{sh}} ,
$$

or in other words $\mathbf{0} = B \boldsymbol{\beta}$. Now $p \equiv 0$ if and only if
$\beta_1 = \ldots = \beta_{n_{\mathrm{sh}}} = 0$, so statement b) is equivalent to
$B \boldsymbol{\beta} = \mathbf{0} \Rightarrow \boldsymbol{\beta} = \mathbf{0}$, that is,
$\ker B = \{0\}$. Since $B$ is a square matrix, injectivity implies bijectivity, so b) is
equivalent to $B$ being invertible. But $B = A^{\top}$, and hence a) and b) are
equivalent.

**a) $\Leftrightarrow$ c).** The matrix representation of the linear mapping
$\Phi_\Sigma$ with respect to the basis
$\{\varphi_1, \ldots, \varphi_{n_{\mathrm{sh}}}\}$ of $\mathcal{P}$ and the standard
basis $\{\mathbf{e}_1, \ldots, \mathbf{e}_{n_{\mathrm{sh}}}\}$ of
$\mathbb{R}^{n_{\mathrm{sh}}}$ is exactly the matrix $B = A^{\top}$ from the previous
step. Thus $\Phi_\Sigma$ is an isomorphism if and only if $B$, and therefore $A$, is
invertible.
```

```{prf:remark} Unisolvence and the generalized Vandermonde matrix
:label: rem:generalized-vandermonde

Property b) or c) is often referred to as **unisolvence**.

The matrix $\big( \sigma_j(\varphi_i) \big)_{i,j=1}^{n_{\mathrm{sh}}}$ is often called
the **generalized Vandermonde matrix**. In 1D, if $\varphi_i(x) = x^i$ and
$\xi_0 < \ldots < \xi_k$ are $k+1$ distinct interpolation points, so that
$\sigma_j(v) = v(\xi_j)$, then
$\big( \sigma_j(\varphi_i) \big)_{ij} = \big( \xi_j^{\,i} \big)_{ij}$
is the usual Vandermonde matrix.
```

```{prf:lemma} Existence of the nodal basis/shape functions
:label: lem:nodal-basis

i) There is a basis $\{\theta_i\}_{i=1}^{n_{\mathrm{sh}}}$ of $\mathcal{P}$ such that

$$
\sigma_j(\theta_i) = \delta_{ij} =
\begin{cases}
1 & i = j, \\
0 & \text{else.}
\end{cases}
$$ (eq:nodal-basis-property)

This basis is called the **nodal basis**, and the $\theta_i$ are called **shape
functions**.

ii) Let $\{\varphi_i\}_{i=1}^{n_{\mathrm{sh}}}$ be a basis of $\mathcal{P}$ and define
$V \in \mathbb{R}^{n_{\mathrm{sh}} \times n_{\mathrm{sh}}}$ by
$V_{ij} = \sigma_j(\varphi_i)$. Then the shape functions can be computed as

$$
\theta_i = \sum_{j=1}^{n_{\mathrm{sh}}} (V^{-1})_{ij} \, \varphi_j ,
\qquad i = 1, \ldots, n_{\mathrm{sh}} .
$$ (eq:shape-functions-from-vandermonde)
```

```{exercise} Existence and computation of the nodal basis
:label: exer-nodal-basis

Prove {prf:ref}`lem:nodal-basis`. For part ii), insert the ansatz
$\theta_i = \sum_j c_{ij} \varphi_j$ into {eq}`eq:nodal-basis-property` and identify the
resulting linear system for the coefficients $c_{ij}$.
```

```{prf:proof}
We make the ansatz $\theta_k = \sum_{j=1}^{n_{\mathrm{sh}}} c_{ik} \varphi_k$ and evaluate $\sigma_k$ on both sides 
find a system of equations for the coefficients $c_{ij}$,

$$
\delta_{ij} =
\sigma_j(\theta_i) = \sum_{k=1}^{n_{\mathrm{sh}}} c_{ik} \sigma_j(\varphi_k) \quad \text{for } i,j = 1, \ldots, n_{\mathrm{sh}},
$$

or in matrix form $I = C V$, where $C = (c_{ij})_{ij}$ and $V = (\sigma_j(\varphi_k))_{kj}$.
Since the invertibility of $V$ is guaranteed by {prf:ref}`lem:unisolvence`, we find that
$C = V^{-1}$ and {eq}`eq:shape-functions-from-vandermonde` follows. 
```

We conclude this section by introducing the concept of interpolation for a general finite element, which generalizes the computation
of the Lagrange interpolation of a continuous function, we discussed
earlier,

$$
I_Tf(x) = \sum_{i=1}^k f(\xi_i) \lambda_i(x)
\sum_{i=1}^k \sigma_i(f) \lambda_i(x)
$$

In general, we wish to interpolate (in a broad sense) function
$f$ from some Banach function space $V(T)$ which contains $P$ and
which is contained in $L^1(T, \mathbb{R}^q)$.
Depending on $V(T)$ the evaluation $\sigma_i(f)$ might not be well-defined
for $f\in V(T)$, e.g. think of point evaluation of Sobolev functions, which
is not well-defined if the regularity of $f$ is too low,

Consequently, to make the evaluation $\sigma_i(f)$ well-defined,
we for a given set $\sigma = \{\sigma_j\}_{j=1}^{n_{sh}} \subset \mathcal{P}'$ associated
with a finite element triple $(T,\mathcal{P}, \Sigma)$, there
is a corresponding sets of extended bounded functionals 

$\widetilde{\sigma} = \{\widetilde{\sigma}_j\}_{j=1}^{n_{sh}} 
\subset \mathcal{V(T)}'$. 
More precisely, we require for the linear functionals that

$$
\begin{aligned}
|\widetilde{\sigma}_i(f)| &\leqslant C_i \| f \|_{V(T)} && \forall f \in V(T) ,
\\
\widetilde{\sigma}_i(p) &= \sigma_i(p) && \forall p \in \mathcal{P},
\end{aligned}
$$

for some constants $C_i$. With such an extension in mind,
we usually skip the $\widetilde{\cdot}$ and just write $\sigma_i(f)$
instead.
The **local interpolation operator**
can now be defined as before, with the nodal basis $\{\theta_i\}_{i=1}^{n_{\mathrm{sh}}}$ of {prf:ref}`lem:nodal-basis`, by setting

$$
I_Tf(x) = \sum_{i=1}^{n_{\mathrm{sh}}} \sigma_i(f) \, \theta_i(x)
$$ (eq:def-local-inter-op)

which gives a well-defined and bounded operator
$I_T : V(T) \to \mathcal{P}$.

```{exercise} Boundedness of the local interpolation operator
:label: exer-local-interpolation-bounded

Prove that the local interpolation operator $I_T$ defined in
{eq}`eq:def-local-inter-op` is a well-defined, bounded linear operator
$I_T : V(T) \to \mathcal{P}$.
```

## Polynomial spaces in several variables

```{prf:lemma} Dimension of $\mathbb{P}_k(\mathbb{R}^d)$
:label: lem:dim-pk

Let $\mathbb{P}_k(\mathbb{R}^d)$ denote the space of polynomials of degree $k$ in $d$
variables, that is,

$$
\mathbb{P}_k(\mathbb{R}^d)
= \Big\{ p : \mathbb{R}^d \to \mathbb{R} \;\Big|\;
p(x) = \sum_{|\alpha| \leqslant k} c_\alpha x^\alpha \Big\} .
$$

Then

$$
\dim \mathbb{P}_k(\mathbb{R}^d) = \binom{d+k}{k} .
$$ (eq:dim-pk)
```

```{exercise} Dimension of the polynomial space
:label: exer-dim-pk

Prove {prf:ref}`lem:dim-pk` by determining the number 
of distinct monomials $x^\alpha$ with $|\alpha| \leqslant k$, or equivalently, by counting
of distinct multi-indices in the set 
$\mathcal{A}_{k,d} = \{ \alpha \in \mathbb{N}^d : |\alpha| \leqslant k \}$.
```

```{admonition} Hint
:class: hint dropdown

Sort the monomials by their total degree. Writing 

$$
\mathbb{H}_l(\mathbb{R}^d) = \Big\{ p : \mathbb{R}^d \to \mathbb{R} \;\Big|\;
p(x) = \sum_{|\alpha| = l} c_\alpha x^\alpha, \ c_\alpha \in \mathbb{R} \Big\}
$$
 for the
space of **homogeneous** polynomials of degree exactly $l$, spanned by the monomials
$x^\alpha$ with $|\alpha| = l$, every polynomial splits uniquely into its homogeneous
parts, so that

$$
\mathbb{P}_k(\mathbb{R}^d) = \bigoplus_{l=0}^{k} \mathbb{H}_l(\mathbb{R}^d),
\qquad \text{and hence} \qquad
\dim \mathbb{P}_k(\mathbb{R}^d) = \sum_{l=0}^{k} \dim \mathbb{H}_l(\mathbb{R}^d) .
$$

It therefore suffices to count the multi-indices with $|\alpha| = l$, that is, the ways
of distributing $l$ units among $d$ variables. Show that

$$
\dim \mathbb{H}_l(\mathbb{R}^d) = \binom{d + l - 1}{l} ,
$$

for instance by induction on $d$, or by the following combinatorial argument: a
multi-index $\alpha$ with $|\alpha| = l$ is the same as an arrangement of $l$ stars and
$d-1$ bars in a row, the bars splitting the stars into the $d$ groups
$\alpha_1, \ldots, \alpha_d$. Such an arrangement consists of $l + d - 1$ symbols and is
fixed as soon as one decides which $l$ of them are the stars, which leaves
$\binom{d+l-1}{l}$ possibilities.

Summing over $l$ then gives

$$
\dim \mathbb{P}_k(\mathbb{R}^d)
= \sum_{l=0}^{k} \binom{d+l-1}{l}
= \binom{d+k}{k} ,
$$

where the last identity follows from repeated application of Pascal's rule
$\binom{m}{j} = \binom{m-1}{j-1} + \binom{m-1}{j}$.
```

Next, we derive a decomposition lemma for multivariate polynomials that vanish on a
hyperplane. Recall that an **affine function** $\lambda : \mathbb{R}^d \to \mathbb{R}$ is
a linear function plus a constant,

$$
\lambda(x) = w \cdot x + c ,
\qquad w \in \mathbb{R}^d, \quad c \in \mathbb{R} ,
$$ (eq:affine-function)

so that $\lambda$ is linear precisely when $c = 0$.
For $w \neq 0$, the **zero set** of $\lambda$

$$
L = \{ x \in \mathbb{R}^d \mid \lambda(x) = 0 \}
$$

is a hyperplane, and conversely every hyperplane is the zero set of such a $\lambda$; it
passes through the origin if and only if $\lambda$ is linear. Often, we will denote the
hyperplane simply by $\lambda$. The following lemma will play an important role when
deriving unisolvence for certain finite elements discussed afterwards.

```{prf:lemma} Decomposition of polynomials vanishing on a hyperplane
:label: lem:hyperplane-decomposition

Let $\lambda$ be an affine function as in {eq}`eq:affine-function` and let
$p \in \mathbb{P}_k(\mathbb{R}^d)$ vanish on the hyperplane $L = \{ \lambda = 0 \}$. Then we can write $p = \lambda q$ where
$q \in \mathbb{P}_{k-1}(\mathbb{R}^d)$.
```

```{prf:proof}
Since $\mathbb{P}_k(\mathbb{R}^d)$ is invariant under composition with affine mappings,
we may, after a suitable affine change of coordinates, assume that

$$
L = \{ x = (\widetilde{x}, x_d) \in \mathbb{R}^d \mid x_d = 0 \} ,
$$

in other words $\lambda(x) = x_d$; rescaling $\lambda$ only rescales $q$, so this
is no restriction. Rewrite a given multi-index
$\alpha = (\alpha_1, \ldots, \alpha_{d-1}, \alpha_d) \in \mathbb{N}^d$
as $\alpha = (\widetilde{\alpha}, \alpha_d)$ 
then we can write any polynomial $p \in \mathbb{P}_k(\mathbb{R}^d)$ as

$$
\begin{aligned}
p(x)
&= \sum_{|\alpha| \leqslant k} c_\alpha x^\alpha
\\
&= \sum_{\alpha_d = 0}^{k} \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_d}
  c_{\widetilde{\alpha}, \alpha_d} \, \widetilde{x}^{\widetilde{\alpha}} x_d^{\alpha_d}
  \\
&= 
\sum_{|\widetilde{\alpha}| \leqslant k}
  c_{\widetilde{\alpha}, 0} \, \widetilde{x}^{\widetilde{\alpha}}
+
\sum_{\alpha_d = 1}^{k} \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_d}
  c_{\widetilde{\alpha}, \alpha_d} \, \widetilde{x}^{\widetilde{\alpha}} x_d^{\alpha_d} ,
\end{aligned}
$$

Evaluating on the hyperplane then gives

$$
0 = p(\widetilde{x}, 0)
= \underbrace{\sum_{|\widetilde{\alpha}| \leqslant k}
  c_{\widetilde{\alpha}, 0} \, \widetilde{x}^{\widetilde{\alpha}}}_{=: p_L(\widetilde{x})}
+ \underbrace{\sum_{\alpha_d = 1}^k \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_d}
  c_{\widetilde{\alpha}, \alpha_d} \, \widetilde{x}^{\widetilde{\alpha}} \, 0^{\alpha_d}}_{= 0} .
$$

Since $p_L(\widetilde{x}) = 0$ for all $\widetilde{x} \in \mathbb{R}^{d-1}$, we conclude
that $c_{\widetilde{\alpha}, 0} = 0$ for all $|\widetilde{\alpha}| \leqslant k$.
Therefore

$$
p(x)
= \sum_{\alpha_d = 1}^{k} \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_d}
  c_{\widetilde{\alpha}, \alpha_d} \, \widetilde{x}^{\widetilde{\alpha}} x_d^{\alpha_d}
= \underbrace{x_d}_{= \lambda(x)} \cdot
  \underbrace{\sum_{\alpha_d = 1}^{k} \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_d}
  c_{\widetilde{\alpha}, \alpha_d} \, \widetilde{x}^{\widetilde{\alpha}} x_d^{\alpha_d - 1}}_{=: q(x)} ,
$$

and $q \in \mathbb{P}_{k-1}(\mathbb{R}^d)$, since every remaining monomial has degree
$|\widetilde{\alpha}| + \alpha_d - 1 \leqslant k - 1$.
```

## Examples of finite element on triangles

### Lagrange elements on triangles

Throughout this section, $T \subset \mathbb{R}^2$ is a non-degenerate triangle with
vertices $a_1, a_2, a_3$, and $\lambda_1, \lambda_2, \lambda_3$ denote its **barycentric
coordinates**: $\lambda_i$ is the affine function, in the sense of
{eq}`eq:affine-function`, with $\lambda_i(a_j) = \delta_{ij}$, so
that $\lambda_i$ vanishes precisely on the line $L_i$ through the two vertices other
than $a_i$.

```{prf:definition} The $\mathbb{P}_1$ Lagrange element on a triangle
:label: def:p1-triangle

The **$\mathbb{P}_1$ Lagrange element** on $T$ is the triple $(T, \mathcal{P}, \Sigma)$ with
$\mathcal{P} = \mathbb{P}_1(T)$ and

$$
\Sigma = \{ \sigma_1, \sigma_2, \sigma_3 \},
\qquad \sigma_i(p) = p(a_i), \quad i = 1, 2, 3 .
$$
```

```{figure} figures/p1-triangle.svg
:label: fig-p1-triangle
:alt: A triangle with its three vertices marked as the degrees of freedom of the P1 element.
:width: 45%

The $\mathbb{P}_1$ element: one degree of freedom per vertex, matching
$\dim \mathbb{P}_1(\mathbb{R}^2) = 3$.
```

By {prf:ref}`lem:dim-pk`, $\dim \mathbb{P}_1(\mathbb{R}^2) = \binom{3}{1} = 3$, which
equals the number of degrees of freedom. It remains to verify unisolvence.

```{prf:lemma} Unisolvence of the $\mathbb{P}_1$ element
:label: lem:p1-unisolvence

The triple of {prf:ref}`def:p1-triangle` is a finite element in the sense of
{prf:ref}`def:ciarlet-finite-element`, that is, $\Sigma$ is a basis of $\mathcal{P}'$.
```

```{prf:proof}
Since $\dim \mathbb{P}_1(T) = 3 = \# \Sigma$, {prf:ref}`lem:unisolvence` reduces the claim
to statement b): every $p \in \mathbb{P}_1(T)$ with $p(a_1) = p(a_2) = p(a_3) = 0$
vanishes identically.

So let $p$ be such a polynomial and consider the line $L_3$ through $a_1$ and $a_2$, the
zero set of $\lambda_3$. Parametrising $L_3$ affinely turns $p|_{L_3}$ into a univariate
polynomial of degree at most $1$ which vanishes at the two distinct points $a_1$ and
$a_2$. A non-trivial univariate polynomial of degree at most $1$ has at most one root,
hence $p$ vanishes identically on $L_3$.

Since $\lambda_3$ is an affine function vanishing exactly on $L_3$,
{prf:ref}`lem:hyperplane-decomposition` applies and we may write

$$
p = \lambda_3 \, q, \qquad q \in \mathbb{P}_0(\mathbb{R}^2),
$$

so $q$ is a constant. Evaluating at the remaining vertex gives

$$
0 = p(a_3) = \lambda_3(a_3) \, q = q ,
$$

because $\lambda_3(a_3) = 1 \neq 0$: the triangle is non-degenerate, so $a_3 \notin L_3$.
Hence $q = 0$ and $p \equiv 0$.
```

```{prf:remark} The nodal basis of the $\mathbb{P}_1$ element
:label: rem:p1-nodal-basis

The barycentric coordinates are exactly the shape functions of this element: they satisfy
$\lambda_i \in \mathbb{P}_1(T)$ and $\sigma_j(\lambda_i) = \lambda_i(a_j) = \delta_{ij}$,
which is the defining property {eq}`eq:nodal-basis-property` of the nodal basis. In the
language of {prf:ref}`lem:nodal-basis`, no Vandermonde system has to be solved here — the
nodal basis is available in closed form.
```

```{prf:definition} The $\mathbb{P}_2$ Lagrange finite element on a triangle
:label: def:p2-triangle

Let $a_{ij} = \frac{1}{2}(a_i + a_j)$ denote the midpoint of the edge joining $a_i$ and
$a_j$. The **$\mathbb{P}_2$ Lagrange element** on $T$ is the triple $(T, \mathcal{P}, \Sigma)$ with
$\mathcal{P} = \mathbb{P}_2(T)$ and

$$
\Sigma = \{ p \mapsto p(a_i) \}_{i=1}^{3}
\; \cup \;
\{ p \mapsto p(a_{ij}) \}_{1 \leqslant i < j \leqslant 3} .
$$
```

```{figure} figures/p2-triangle.svg
:label: fig-p2-triangle
:alt: A triangle with its three vertices and three edge midpoints marked as the degrees of freedom of the P2 element.
:width: 45%

The $\mathbb{P}_2$ element: one degree of freedom per vertex and per edge midpoint,
matching $\dim \mathbb{P}_2(\mathbb{R}^2) = 6$.
```

Again the counts agree: $\dim \mathbb{P}_2(\mathbb{R}^2) = \binom{4}{2} = 6 = \# \Sigma$.

```{exercise} Unisolvence of the $\mathbb{P}_2$ element
:label: exer-p2-unisolvence

Show that the triple of {prf:ref}`def:p2-triangle` is a finite element.
```

```{admonition} Hint
:class: hint dropdown

Follow the proof of {prf:ref}`lem:p1-unisolvence`. Each edge of $T$ now carries three
nodes, two vertices and one midpoint, so the restriction of $p \in \mathbb{P}_2(T)$ to
the line $L_i$ is a univariate quadratic with three distinct roots and vanishes
identically. Apply {prf:ref}`lem:hyperplane-decomposition` twice on two distinct
edges $L_i$ and $L_j$ to show that $p = c \lambda_i \lambda_j
$ with $c$ being a constant and finally deduce that $c = 0$.
```


```{prf:definition} The $\mathbb{P}_3$ Lagrange finite element on a triangle
:label: def:p3-triangle

Let $a_{iij} = \frac{1}{3}(2 a_i + a_j)$ for $i \neq j$ denote the two points dividing the
edge from $a_i$ to $a_j$ into three equal parts, and let
$a_{123} = \frac{1}{3}(a_1 + a_2 + a_3)$ be the barycenter of $T$. The **$\mathbb{P}_3$
element** on $T$ is the triple $(T, \mathcal{P}, \Sigma)$ with
$\mathcal{P} = \mathbb{P}_3(T)$ and

$$
\Sigma = \{ p \mapsto p(a_i) \}_{i=1}^{3}
\; \cup \;
\{ p \mapsto p(a_{iij}) \}_{i \neq j}
\; \cup \;
\{ p \mapsto p(a_{123}) \} .
$$
```

```{figure} figures/p3-triangle.svg
:label: fig-p3-triangle
:alt: A triangle with three vertices, two points on each edge and the barycenter marked as the degrees of freedom of the P3 element.
:width: 45%

The $\mathbb{P}_3$ element: one degree of freedom per vertex, two per edge and one in the
interior, matching $\dim \mathbb{P}_3(\mathbb{R}^2) = 10$.
```

Here $\dim \mathbb{P}_3(\mathbb{R}^2) = \binom{5}{3} = 10 = 3 + 6 + 1 = \# \Sigma$. Note
that a single point per edge would not suffice: the interior degree of freedom and the
second point on each edge are exactly what makes the counts match.

```{exercise} Unisolvence of the $\mathbb{P}_3$ element
:label: exer-p3-unisolvence

Show that the triple of {prf:ref}`def:p3-triangle` is a finite element.
```

```{admonition} Hint
:class: hint dropdown

Each edge now carries four nodes, so the restriction of $p \in \mathbb{P}_3(T)$ to each
line $L_i$ is a univariate cubic with four distinct roots. Peeling off all three
barycentric coordinates with {prf:ref}`lem:hyperplane-decomposition` leaves
$p = c \, \lambda_1 \lambda_2 \lambda_3$ with a constant $c$, since the degrees already
match. What does the remaining degree of freedom at the barycenter give you?
```

### The cubic Hermite element on a triangle

In all elements so far, the degrees of freedom were point evaluations, and such elements
are called **Lagrange elements**. {prf:ref}`def:ciarlet-finite-element` allows any linear
functionals on $\mathcal{P}$, however, and prescribing derivatives is the natural next
choice. Elements whose degrees of freedom involve derivatives are called **Hermite
elements**.

```{prf:definition} The cubic Hermite triangle
:label: def:p3-hermite-triangle

Let $a_{123} = \frac{1}{3}(a_1 + a_2 + a_3)$ be the barycenter of $T$. The **cubic
Hermite element** on $T$ is the triple $(T, \mathcal{P}, \Sigma)$ with
$\mathcal{P} = \mathbb{P}_3(T)$ and

$$
\Sigma =
\{ p \mapsto p(a_i) \}_{i=1}^{3}
\; \cup \;
\{ p \mapsto \partial_1 p(a_i), \; p \mapsto \partial_2 p(a_i) \}_{i=1}^{3}
\; \cup \;
\{ p \mapsto p(a_{123}) \} ,
$$

where $\partial_1, \partial_2$ denote the partial derivatives with respect to the
coordinates of $\mathbb{R}^2$.
```

```{figure} figures/p3-hermite-triangle.svg
:label: fig-p3-hermite-triangle
:alt: A triangle whose three vertices carry a function value and a gradient, drawn as a dot inside a ring, together with the barycenter carrying a function value.
:width: 48%

The cubic Hermite element: at each vertex the value and the two first derivatives, and
the value at the barycenter. The ring around a vertex stands for the two gradient
degrees of freedom.
```

The counts agree again: $3 + 6 + 1 = 10 = \dim \mathbb{P}_3(\mathbb{R}^2)$. The element
uses the same polynomial space as the cubic Lagrange triangle of
{prf:ref}`def:p3-triangle`, but distributes its degrees of freedom differently.

```{prf:lemma} Unisolvence of the cubic Hermite element
:label: lem:hermite-unisolvence

The triple of {prf:ref}`def:p3-hermite-triangle` is a finite element in the sense of
{prf:ref}`def:ciarlet-finite-element`.
```

```{prf:proof}
Since $\# \Sigma = 10 = \dim \mathbb{P}_3(T)$, {prf:ref}`lem:unisolvence` again reduces
the claim to statement b). So let $p \in \mathbb{P}_3(T)$ satisfy

$$
p(a_i) = 0, \quad \nabla p(a_i) = 0 \quad (i = 1, 2, 3),
\qquad p(a_{123}) = 0 ,
$$

and let us show that $p$ vanishes identically.

**Step 1: $p$ vanishes on each edge line.** Consider the edge joining $a_1$ and $a_2$,
which lies on the line $L_3 = \{ \lambda_3 = 0 \}$, and parametrise that line by
$\gamma(t) = a_1 + t (a_2 - a_1)$. Then $\varphi := p \circ \gamma$ is a univariate
polynomial of degree at most $3$ with

$$
\varphi(0) = p(a_1) = 0, \qquad \varphi(1) = p(a_2) = 0 ,
$$

and, by the chain rule,

$$
\varphi'(0) = \nabla p(a_1) \cdot (a_2 - a_1) = 0 ,
\qquad
\varphi'(1) = \nabla p(a_2) \cdot (a_2 - a_1) = 0 .
$$

Thus $\varphi$ has roots of multiplicity at least two at $t = 0$ and at $t = 1$, that is,
at least four roots counted with multiplicity. A non-trivial univariate polynomial of
degree at most $3$ has at most three, so $\varphi \equiv 0$ and $p$ vanishes on all of
$L_3$. The same argument applied to the other two edges shows that $p$ vanishes on
$L_1$, $L_2$ and $L_3$.

**Step 2: peeling off the barycentric coordinates.** By
{prf:ref}`lem:hyperplane-decomposition`, applied to the affine function $\lambda_3$,

$$
p = \lambda_3 \, q_1, \qquad q_1 \in \mathbb{P}_2(\mathbb{R}^2) .
$$

Now $p$ also vanishes on $L_1$, while $\lambda_3$ vanishes at exactly one point of $L_1$,
namely the vertex $a_2 = L_1 \cap L_3$. Hence $q_1$ vanishes on $L_1$ with the possible
exception of that single point, and therefore on all of $L_1$ by continuity. Applying the
lemma twice more gives

$$
q_1 = \lambda_1 \, q_2, \quad q_2 \in \mathbb{P}_1(\mathbb{R}^2),
\qquad
q_2 = \lambda_2 \, c, \quad c \in \mathbb{P}_0(\mathbb{R}^2) ,
$$

so that $p = c \, \lambda_1 \lambda_2 \lambda_3$ with a constant $c$.

**Step 3: the interior degree of freedom.** Since
$\lambda_i(a_{123}) = \frac{1}{3}$ for $i = 1, 2, 3$,

$$
0 = p(a_{123}) = c \, \Big( \frac{1}{3} \Big)^3 = \frac{c}{27} ,
$$

hence $c = 0$ and $p \equiv 0$.
```

```{prf:remark} Why the derivative degrees of freedom need care
:label: rem:hermite-affine

The gradient degrees of freedom do not transform as simply as point values. If
$\Phi_T(\widehat{x}) = B \widehat{x} + b$ maps a reference triangle onto $T$ and
$\widehat{p} = p \circ \Phi_T$, then
$\nabla \widehat{p}(\widehat{x}) = B^{\top} \nabla p(\Phi_T(\widehat{x}))$, so the
Jacobian enters whenever degrees of freedom are transferred between the reference element
and a physical element. For the Lagrange elements above, whose degrees of freedom are
point values, the transfer is immediate.
```

## Examples of finite element on tetrahedra

```{admonition} TODO
:class: warning
Add examples of finite elements on tetrahedra, for instance the $\mathbb{P}_1$
elements.
```

## $\mathbb{P}_k$ Lagrange elements on d-simplices

Everything so far was carried out on a triangle. The construction does not depend on the
space dimension, and we now formulate it for simplices in $\mathbb{R}^d$: intervals for
$d = 1$, triangles for $d = 2$, tetrahedra for $d = 3$.

### Simplices and barycentric coordinates

```{prf:definition} Regular $d$-simplex
:label: def:simplex

Let $a_0, a_1, \ldots, a_d \in \mathbb{R}^d$ be **affinely independent**, that is, let the
$d$ vectors $a_1 - a_0, \ldots, a_d - a_0$ be linearly independent. The **$d$-simplex**
with vertices $a_0, \ldots, a_d$ is their convex hull

$$
\begin{aligned}
T = \operatorname{conv}\{a_0, \ldots, a_d\}
&= \Big\{ \sum_{i=0}^{d} \lambda_i a_i \;\Big|\;
\lambda_i \geqslant 0, \; \sum_{i=0}^{d} \lambda_i = 1 \Big\}
\\
&= \Big\{ a_0 + \sum_{i=1}^{d} \lambda_i (a_i-a_0) \;\Big|\;
\lambda_i \geqslant 0, \; \sum_{i=1}^{d} \lambda_i \leqslant 1 \Big\},
\end{aligned}
$$ (eq:simplex)
that is, the set of all **convex combinations** of its vertices. 

Such a simplex is called **regular**, or non-degenerate. Affine independence is equivalent
to the invertibility of the matrix
$B = [\, a_1 - a_0, \; \ldots, \; a_d - a_0 \,] \in \mathbb{R}^{d \times d}$, and to $T$
having non-empty interior.

More generally, a **$k$-simplex in $\mathbb{R}^d$** is the convex hull of $k+1$ 
affinely independent points in $\mathbb{R}^d$.

For each $i$, the **face opposite $a_i$** is the $(d-1)$-simplex
$F_i = \operatorname{conv}\{ a_j : j \neq i \}$.
of the $d$-simplex $T = \operatorname{conv}\{a_0, \ldots, a_d\}$.

A $k$-face of $T$ is the convex hull of any $k+1$ of its vertices (and thus a $k$-simplex). 
A $d-1$-face is called a **facet**, $1$-faces are **edges** and $0$-faces are **vertices**.
```

A regular $1$-simplex is an interval, a regular $2$-simplex a triangle and a regular
$3$-simplex a tetrahedron. The whole construction below rests on one set of functions
attached to $T$, which we already met for triangles.

Thanks to the linear independence of the $d$ vectors $a_1 - a_0, \ldots, a_d - a_0$, we see
that every point $x \in T$ can be written as **unique** convex combination
$\sum_{i=0}^{d} \lambda_i a_i$
 of the vertices $a_0, \ldots, a_d$, 
 thus defining a mapping
 $\boldsymbol{x} \mapsto \boldsymbol{\lambda}(\boldsymbol{x}) = (\lambda_0(\boldsymbol{x}), \ldots, \lambda_d(\boldsymbol{x})) \in \mathbb{R}^{d+1}$.
If we allow $\lambda_i$ to take negative values with the total sum
still being 1, then the same argument shows that every point $x \in
\mathbb{R}^d$ can be written as a **unique affine combination** of the
vertices.

```{prf:definition} Barycentric coordinates
:label: def:barycentric-simplex

Let $T \operatorname{conv}\{a_0, \ldots, a_d\}$ be a regular $d$-simplex. 
For $x \in \mathbb{R}^d$ the system

$$
\sum_{i=0}^{d} \lambda_i(x) = 1,
\qquad
\sum_{i=0}^{d} \lambda_i(x) \, a_i = x
$$ (eq:barycentric-simplex-system)

has exactly one solution $(\lambda_0(x), \ldots, \lambda_d(x)) \in \mathbb{R}^{d+1}$. The
functions $\lambda_0, \ldots, \lambda_d$ are called the **barycentric coordinates** of $T$.
```

For the unit simplex $\widehat{T} = \operatorname{conv}\{0, e_1, \ldots, e_d\}$, the barycentric coordinates are particularly simple to write down:

$$
\widehat{x} = \sum_{i=1}^{d} \widehat{x}_i e_i 
  = \Big( 1 - \sum_{i=1}^{d} \widehat{x}_i \Big) 0
  + \sum_{i=1}^{d} \widehat{x}_i e_i,
$$

and therefore,

$$
\widehat{\lambda}_0(\widehat{x}) = 1 - \sum_{i=1}^{d} \widehat{x}_i, \quad \text{and } \quad \widehat{\lambda}_i(\widehat{x}) = \widehat{x}_i \quad \text{ for } i = 1, \ldots, d.
$$ (eq:barycentric-unit-simplex)

In particular, we see that the barycentric coordinates 
are affine functions of $x$,
that constitutes the basis of the $\mathbb{P}_1(\mathbb{R}^d)$
$$
\mathbb{P}_1(\mathbb{R}^d)
= \operatorname{span}\{ \widehat{\lambda}_i\}_{i=0}^d.
$$

Let us investigate some properties of the barycentric coordinates for a general regular $d$-simplex $T = \operatorname{conv}\{a_0, \ldots, a_d\}$.
We can rewrite {eq}`eq:barycentric-simplex-system` in matrix form as

$$
\underbrace{
\begin{pmatrix} 
  a_0 & \cdots & a_d 
  \\ 
  1 & \cdots & 1 
  \end{pmatrix}
}_{=: B}
\begin{pmatrix} \lambda_0(x) \\ \vdots \\ \lambda_d(x) \end{pmatrix}
= \begin{pmatrix} x \\ 1 \end{pmatrix} ,
$$

The invertibility of the $(d+1) \times (d+1)$ matrix $B$ is equivalent to the affine
independence of the vertices, since a column-wise subtraction column 1 containing 
$(a_0^{\top}, 1)^{\top}$ from the subsequent columns does not change the determinant, and yields the matrix

$$B \sim \widetilde{B} = 
\begin{pmatrix} a_0 & a_1 - a_0 & \cdots & a_d - a_0 \\ 
1 & 0 & 0 & \cdots & 0 
\end{pmatrix}
$$

Writing $C = B^{-1}$, we see that the barycentric coordinates are affine functions of $x$:

$$
\begin{pmatrix} \lambda_0(x) \\ \vdots \\ \lambda_d(x) \end{pmatrix}
= C \begin{pmatrix} x \\ 1 \end{pmatrix} .
$$ (eq:barycentric-coord-matrix-system)

or equivalently,

$$
\lambda_i(x) =  \sum_{j=1}^{d} c_{i,j} x_j + c_{i, d+1} ,
$$

for $i = 0, \ldots, d$,

Since the affine functions 
$p_i(x) = x_i$ for $i = 1, \ldots, d$ and $p_0(x) = 1$ form a basis of the space of affine functions on $\mathbb{R}^d$, 
so do the barycentric coordinates $\lambda_i(x)$ for $i = 0, \ldots, d$,
thanks to @eq:barycentric-coord-matrix-system
and invertibility of the matrix $C$.
Further, since $a_i = 1 a_i = \sum_{j=0}^{d} \lambda_j(a_i) a_j$, we have that $\lambda_i(a_i) = 1$ and $\lambda_j(a_i) = 0$ for $j \neq i.$
In other words, the barycentric coordinates are the (uniquely defined) affine functions on $\mathbb{R}^d$ that take the value $1$ at one vertex and $0$ at all other vertices.
This leads us to define our first finite element on a $d$-simplex.

```{prf:proposition} $\mathbb{P}_1$ Lagrange element on a $d$-simplex
:label: prop:p1-simplex
We let
- $T = \operatorname{conv}\{a_0, \ldots, a_d\}$ be a regular $d$-simplex
- $\mathcal{P} = \mathbb{P}_1(T)$ be the space of affine functions on $T$
- $\Sigma = \{ \sigma_0, \ldots, \sigma_d \}$ with $\sigma_i(p) = p(a_i)$ for $i = 0, \ldots, d$.


Then the triple $(T, \mathcal{P}, \Sigma)$ is a finite element.
Moreover, the barycentric coordinates $\lambda_i$ for $i = 0, \ldots, d$ are the shape functions of this element and thus for
$p \in \mathcal{P}$, we have that the representation

$$
p(x) = \sum_{i=0}^{d} p(a_i) \, \lambda_i(x) .
$$ (eq:p1-lagrange-representation)

```

```{prf:proof}
We start from the back and show @eq:p1-lagrange-representation first. Since we already concluded that the barycentric coordinates form a basis of $\mathbb{P}_1(\mathbb{R}^d)$, we can write any $p \in \mathbb{P}_1(T)$ as

$$
p(x) = \sum_{i=0}^{d} c_i \lambda_i(x)
$$

Since $\lambda_i(a_j) = \delta_{ij}$, evaluating at $x = a_j$ gives $p(a_j) = c_j$, so that $c_j = p(a_j)$.
From here, it easy to conclude that $\Sigma$ 
satisfies the unisolvence property of {prf:ref}`def:ciarlet-finite-element`, 
since 
$0 = p(a_i)$ for $i = 0, \ldots, d$ implies that $p \equiv 0$.
```

Since $\lambda_i$ are the affine functions characterized by the property $\lambda_i(a_j) = \delta_{ij}$, 
we can derive a few alternative formulas for the barycentric coordinates.

```{prf:proposition} Alternative formulas for the barycentric coordinates
:label: prop:barycentric-alternative-formulas
For a regular $d$-simplex $T = \operatorname{conv}\{a_0, \ldots, a_d\}$, the barycentric coordinates $\lambda_i$ for $i = 0, \ldots, d$ can be expressed as

$$
\begin{aligned}
\lambda_i(x) &= 1 - \dfrac{|F_i|_{d-1}}{d |T|_{d} } (x-a_i)\cdot n_i  \\
             &= 1 - \dfrac{(x - a_i)\cdot n_i}{(a_l - a_i) \cdot n_i} \quad \text{for any } l \neq i.
\end{aligned}
$$

```

```{prf:proof}
Both formulas follow from the fact that $\lambda_i$ is an affine function, and thus we need only to verify that satifies
$\lambda_i(a_j) = \delta_{ij}$.
For the second formula, we can verify as follow. Obviously, $\lambda_i(a_i) = 1$.
For $j \neq i$, we have that 

$$
(a_j - a_i)\cdot n_i = \underbrace{(a_j - a_l)\cdot n_i}_{= 0} + 
(a_l - a_i)\cdot n_i =
(a_l - a_i)\cdot n_i 
$$
where the first term is zero since $a_l, a_j \in F_i$ and thus (a_l - a_j) is parallel to $F_i$ while $n_i$ is orthogonal to the face $F_i$.
and
$(a_l - a_i)\cdot n_i \neq 0$ for any $l \neq i$ since the $a_i$ is not included in the hyperplane containing the face $F_i$.

Thus

$$
\lambda_i(a_j) = 1 - \dfrac{(a_j - a_i)\cdot n_i}{(a_l - a_i) \cdot n_i} = 1 - \dfrac{(a_l - a_i)\cdot n_i}{(a_l - a_i) \cdot n_i} = 0
$$

The proof of the first formula is left as {prf:ref}`exer-barycentric-normal-formula`.
```

```{exercise} The normal-vector formula for the barycentric coordinates
:label: exer-barycentric-normal-formula

Prove the first formula of {prf:ref}`prop:barycentric-alternative-formulas`,

$$
\lambda_i(x) = 1 - \frac{|F_i|_{d-1}}{d \, |T|_d} \, (x - a_i) \cdot n_i ,
$$

where $n_i$ denotes the outward unit normal of $T$ on the face $F_i$ opposite $a_i$.
```

```{admonition} Hint
:class: hint dropdown

Let $h_i$ be the distance from the vertex $a_i$ to the hyperplane containing $F_i$.

1. Show that $(a_l - a_i) \cdot n_i = h_i$ for every $l \neq i$.
2. Show that $|T|_d = \frac{1}{d} |F_i|_{d-1} \, h_i$. Slice $T$ by the hyperplanes
   parallel to $F_i$, as in the proof of {prf:ref}`lem:unit-simplex-volume`. At distance
   $t \in [0, h_i]$ from $a_i$, the slice is a copy of $F_i$ scaled by $t / h_i$.
3. Conclude that the first formula coincides with the second one.
```








Another important property of the barycentric coordinates is that they are 
invariant under affine transformations and that they allows to 
define bijective affine mappings between two given non-degenerate simplices
$\widehat{T}$ and $T$. More precisely, we have the following lemma:


```{prf:lemma} Affine invariance of barycentric coordinates
:label: lem:barycentric-affine-invariance

Given a regular $d$-simplex
$\widehat{T} = \operatorname{conv}\{\widehat{a}_0, \ldots, \widehat{a}_d\}$
and an regular affine mapping $F : \mathbb{R}^d \to \mathbb{R}^d$, 
$\widehat{x} \mapsto F(\widehat{x}) = B \widehat{x} + b$ with $B \in \mathbb{R}^{d \times d}$ invertible and $b \in \mathbb{R}^d$.
Then $F$ maps $\widehat{T}$ bijectively onto the regular $d$-simplex
$T = \operatorname{conv}\{a_0, \ldots, a_d\}$ with $a_i = F(\widehat{a}_i)$, and for $x = F(\widehat{x})$, we have that 
$\mathbf{\lambda}(x) = \widehat{\mathbf{\lambda}}(\widehat{x}) = \widehat{\lambda}(F^{-1}(x))$.

In particular, recalling @eq:barycentric-simplex-system, we have that

$$
\begin{aligned}
\lambda_i(x) &= 
\widehat{\lambda}_i(F^{-1}(x)) =
(F^{-1}(x))_i \quad i = 1, \ldots, d
\\
\lambda_0(x) &= 1 - \sum_{i=1}^{d} (F^{-1}(x))_i
\end{aligned}
$$

```


```{prf:proof}
Using that $\sum_{i=0}^{d} \widehat{\lambda}_i(\widehat{x}) = 1$ and thus
$b = \sum_{i=0}^{d} \widehat{\lambda}_i(\widehat{x}) b$, we deduce that
for $\widehat{x} = \sum_{i=0}^{d} \widehat{\lambda}_i(\widehat{x}) \, \widehat{a}_i$ we have that

$$
\begin{aligned}
x = F(\widehat{x}) 
  &= B \widehat{x} + b 
  \\
  &= B \sum_{i=0}^{d} \widehat{\lambda}_i(\widehat{x}) \, \widehat{a}_i + 
  \sum_{i=0}^{d} \widehat{\lambda}_i(\widehat{x}) b
  \\
  &= \sum_{i=0}^{d} \widehat{\lambda}_i(\widehat{x}) \, (B \widehat{a}_i + b)
  = \sum_{i=0}^{d} \widehat{\lambda}_i(\widehat{x}) \, a_i .
\end{aligned}
$$
so $\lambda_i(x) = \widehat{\lambda}_i(\widehat{x})$ for $i = 0, \ldots, d$.
```

The proof of the previous lemma gives us also a recipe for constructing a bijective affine mapping between two given non-degenerate simplices
$\widehat{T} = \operatorname{conv}\{\widehat{a}_0, \ldots, \widehat{a}_d\}$ and $T = \operatorname{conv}\{a_0, \ldots, a_d\}$.
Using the fact that $\widehat{\lambda}(\widehat{x})$ are affine functions
and satisfy $\widehat{\lambda}_i(\widehat{a}_j) = \delta_{ij}$, 
the affine mapping $F$  given by
$$
F(\widehat{x}) = \sum_{i=0}^{d} \widehat{\lambda}_i(\widehat{x}) \, a_i .
$$
satisfies $F(\widehat{a}_i) = a_i$ for $i = 0, \ldots, d$.

We proceed by investigating the geometric meaning of the barycentric
coordinates.  To do so, we will need a formula for the volume of a
simplex. As any regular $d$-simplex can be mapped onto the unit
simplex by an affine transformation, it is natural to compute the
volume of the unit simplex first and then use the change of variables
formula to compute the volume of any regular $d$-simplex. Let's start
with a formula for the volume of the unit simplex.
Here, the idea to use to combine induction over the dimension and Fubini's theorem
by slicing the unit simplex over the last coordinate.

```{figure} figures/unit-simplex-slicing.svg
:label: fig-unit-simplex-slicing
:alt: The unit tetrahedron with vertices at the origin and at 1 on each of the three coordinate axes. A horizontal triangle at height x3 = s is shaded and labelled (1 - s) times the unit triangle.
:width: 55%

Slicing the unit simplex $\widehat{T}_3$ at height $x_3 = s$. The slice is the unit
triangle $\widehat{T}_2$ scaled by $1 - s$, so its area is $(1 - s)^2 / 2$, and
integrating over $s \in [0, 1]$ gives $|\widehat{T}_3| = 1/6$.
```

```{prf:lemma} Volume of the unit simplex
:label: lem:unit-simplex-volume

The unit simplex
$\widehat{T}_d = \{ x \in \mathbb{R}^d : x_i \geqslant 0, \; \sum_{i=1}^{d} x_i \leqslant 1 \}$
has volume $|\widehat{T}_d|_d = 1/d!$.
```

```{prf:proof}
We use induction on $d$. For $d = 1$, $\widehat{T}_1 = [0, 1]$ has length $1 = 1/1!$.

Let $d \geqslant 2$. For fixed $x_d = s \in [0, 1]$, the slice of $\widehat{T}_d$ is

$$
\Big\{ x' \in \mathbb{R}^{d-1} : x'_i \geqslant 0, \; \sum_{i=1}^{d-1} x'_i \leqslant 1 - s \Big\}
= (1 - s) \, \widehat{T}_{d-1} ,
$$

and for $s \notin [0, 1]$ the slice is empty. Under the linear map $x' \mapsto t x'$ volumes
in $\mathbb{R}^{d-1}$ scale by $t^{d-1}$, so the slice has volume
$(1 - s)^{d-1} |\widehat{T}_{d-1}|_{d-1} = (1 - s)^{d-1} / (d-1)!$ by the induction
hypothesis. Fubini's theorem then gives

$$
|\widehat{T}_d|_d
= \int_0^1 (1 - s) |\widehat{T}_{d-1}|_{d-1} \, \mathrm{d}s
= \int_0^1 \frac{(1 - s)^{d-1}}{(d-1)!} \, \mathrm{d}s
= \frac{1}{(d-1)!} \cdot \frac{1}{d}
= \frac{1}{d!} .
$$
```

<!-- 
```{prf:remark} A proof by cutting the unit cube
:label: rem:unit-simplex-volume-cube

For a permutation $\sigma$ of $\{1, \ldots, d\}$ let
$C_\sigma = \{ x \in [0, 1]^d : x_{\sigma(1)} \leqslant x_{\sigma(2)} \leqslant \cdots \leqslant x_{\sigma(d)} \}$.
Every point of the cube lies in some $C_\sigma$, and two different $C_\sigma$ overlap only
in points with two equal coordinates, a set of measure zero. Permuting coordinates
preserves volume, so all $d!$ sets $C_\sigma$ have the same volume and

$$
|C_{\mathrm{id}}| = \frac{|[0,1]^d|}{d!} = \frac{1}{d!} .
$$

The partial-sum map $y_k = x_1 + \cdots + x_k$, $k = 1, \ldots, d$, is linear with a lower
triangular matrix of ones, so its determinant is $1$. It maps $\widehat{T}_d$ bijectively
onto $C_{\mathrm{id}} = \{ 0 \leqslant y_1 \leqslant \cdots \leqslant y_d \leqslant 1 \}$:
the conditions $x_i \geqslant 0$ become $y_{i-1} \leqslant y_i$ (with $y_0 = 0$), and
$\sum_i x_i \leqslant 1$ becomes $y_d \leqslant 1$. Hence
$|\widehat{T}_d| = |C_{\mathrm{id}}| = 1/d!$.
``` 
-->

```{prf:corollary} Volume of a simplex
:label: cor:simplex-volume

Let $T = \operatorname{conv}\{a_0, \ldots, a_d\}$ be a regular $d$-simplex and
$B_T = [\, a_1 - a_0, \; \ldots, \; a_d - a_0 \,] \in \mathbb{R}^{d \times d}$. Then

$$
|T|_d = \frac{|\det B_T|}{d!}
= \frac{1}{d!} \left| \det
\begin{pmatrix}
  a_0 & \cdots & a_d
  \\
  1 & \cdots & 1
\end{pmatrix}
\right| .
$$ (eq:simplex-volume)
```

```{prf:proof}
By the second representation in {eq}`eq:simplex`, the affine map
$\Phi_T(\widehat{x}) = B_T \widehat{x} + a_0$ maps $\widehat{T}_d$ bijectively onto $T$.
The change of variables formula and {prf:ref}`lem:unit-simplex-volume` give

$$
|T|_d = \int_{\widehat{T}_d} |\det B_T| \, \mathrm{d}\widehat{x}
= |\det B_T| \, |\widehat{T}_d|_d = \frac{|\det B_T|}{d!} .
$$

The second expression in {eq}`eq:simplex-volume` follows from the column operations above:
the determinant of $\widetilde{B}$ equals, up to sign, $\det B_T$, as one sees by expanding
along its last row.
```

Next, recall that Cramer's rule for computing the solution of a linear system $Ax = b$ states that the $i$-th component of the solution is given by

$$
x_i = \frac{\det A_i}{\det A} ,
$$
where $A_i$ is the matrix obtained from $A$ by replacing its $i$-th column with the right-hand side vector $b$. Applying this to our system,
we can express the barycentric coordinates in terms of determinants:

$$
\lambda_i(x) = \frac{\det B_i(x)}{\det B} , \qquad i = 0, \ldots, d ,
$$

where 

$$
B_i(x) =
\begin{pmatrix}
a_0 & \cdots & a_{i-1} & x & a_{i+1} & \cdots & a_d \\
1 & \cdots & 1 & 1 & 1 & \cdots & 1
\end{pmatrix}
$$

As for matrix $B$, the determinant of $B_i(x)$ 
is related to the volume of a simplex
$T_i(x) = \operatorname{conv}\{a_0, \ldots, a_{i-1}, x, a_{i+1}, \ldots, a_d\}$,

$$
|T_i(x)|_d = \frac{|\det B_i(x)|}{d!}.
$$
and consequently, we deduce that

$$
\lambda_i(x) = 
\frac{\det B_i(x)}{\det B} 
=
\frac{\det B_i(x)/d!}{\det B/d!} 
= \frac{|T_i(x)|_d}{|T|_d}
$$
where we used the fact that $\det B$
and $\det B_i(x)$ have the same sign (why?). 

In other words, the barycentric coordinate $\lambda_i(x)$ is the ratio of the volume of the simplex $T_i(x)$ which is obtained from $T$ by replacing the vertex $a_i$ with $x$,
and the volume of the original simplex $T$. 

```{figure} figures/barycentric-coordinates.svg
:label: fig-barycentric-coordinates
:alt: A triangle with vertices a0, a1, a2 and an interior point x joined to all three vertices, splitting it into three shaded sub-triangles labelled lambda_0, lambda_1 and lambda_2.
:width: 50%

Barycentric coordinates for $d = 2$. Joining $x$ to the vertices splits $T$ into three
sub-triangles; by part iv) of {prf:ref}`lem:barycentric-properties`, $\lambda_i(x)$ is the
area of the one opposite $a_i$, divided by the area of $T$. The coordinate $\lambda_i$
vanishes on the edge opposite $a_i$ and equals $1$ at $a_i$.
```


Let us summarize the properties of the barycentric coordinates
we have derived so far in the following proposition:

```{prf:proposition} Properties of the barycentric coordinates
:label: lem:barycentric-properties

Let $T = \operatorname{conv}\{a_0, \ldots, a_d\}$ be a regular $d$-simplex with barycentric coordinates $\lambda_0, \ldots, \lambda_d$.

i) Each $\lambda_i$ is the unique affine function satisfying $\lambda_i(a_j) = \delta_{ij}$.

ii) Every affine function $p$ satisfies
    $p(x) = \sum_{i=0}^{d} p(a_i) \, \lambda_i(x)$. In particular
    $\{\lambda_i\}_{i=0}^{d}$ is the nodal basis/shape functions for the $\mathbb{P}_1(T)$ Lagrange element.

iii) $T = \{ x \in \mathbb{R}^d : \lambda_i(x) \geqslant 0 \text{ for } i = 0, \ldots, d \}$,
     and the face opposite $a_i$ is $F_i = T \cap \{ \lambda_i = 0 \}$.

iv) $\lambda_i(x) = |T_i(x)|_d / |T|_d$, where $T_i(x)$ is the simplex obtained from $T$ by
    replacing the vertex $a_i$ with $x$, and $|\cdot|$ denotes $d$-dimensional volume.
```


### The $\mathbb{P}_k$ Lagrange element on a simplex

With the barycentric coordinates and their basic properties at hand, we can now define the $\mathbb{P}_k$ Lagrange on a regular $d$-simplex $T$.

The degrees of freedom are point values, taken on the points whose barycentric coordinates
are multiples of $1/k$.

```{prf:definition} Principal lattice and the $\mathbb{P}_k$ Lagrange element
:label: def:pk-simplex

Let $T = \operatorname{conv}\{a_0, \ldots, a_d\}$ be a regular $d$-simplex 

and $k \geqslant 1$. The **principal lattice** of order
$k$ is the set of points defined through their barycentric coordinates as follows:

$$
\begin{aligned}
\mathcal{L}_k(T) &=
\{ a_{\alpha}\}_{\alpha \in \mathcal{A}_{k,d}}
\\
&= \Big\{ a_\alpha \in T \mid  
\boldsymbol{\lambda}(a_\alpha) = 
\dfrac{1}{k} (\underbrace{k-|\alpha|}_{=:\alpha_0}, \underbrace{\alpha_1, \ldots, \alpha_d}_{\alpha}),
\alpha \in \mathcal{A}_{k,d} \Big\}
\\
&= 
\Big\{ a_\alpha 
 =  \sum_{i=0}^{d} \frac{\alpha_i}{k} \, a_i
\;:\; (\alpha_0, \alpha) \in \mathbb{N}_0^{1+d}, \; |(\alpha_0, \alpha)| = k \Big\},
\end{aligned}
$$ (eq:principal-lattice)
that is, the points of $T$ whose barycentric coordinates are integer multiples of $1/k$.
The **$\mathbb{P}_k$ Lagrange element** on $T$ is the triple $(T, \mathcal{P}, \Sigma)$
with $\mathcal{P} = \mathbb{P}_k(T)$ and

$$
\Sigma = \{ p \mapsto p(x_\alpha) \;:\; x_\alpha \in \mathcal{L}_k(T) \} .
$$
```

The number of degrees of freedom is right: counting the multi-indices
$\alpha \in \mathcal{A}_{k,d}$ is the same stars and bars count as in
{prf:ref}`exer-dim-pk`, with $k$ units distributed over $d+1$ slots, so

$$
\# \mathcal{A}_{k,d} = \binom{d+k}{k} = \dim \mathbb{P}_k(\mathbb{R}^d)
$$ (eq:lattice-count)

by {prf:ref}`lem:dim-pk`. For $d = 2$ the definition reproduces the three elements of the
previous section: the vertices for $k = 1$, the vertices and edge midpoints for $k = 2$,
and the vertices, two points per edge and the barycenter for $k = 3$.

Before we turn to the proof of unisolvence, let us record how Lagrange elements on
different simplices are related. The relation is an instance of the following notion.

```{prf:definition} Affine equivalent finite elements
:label: def:affine-equivalence

Two finite elements $(\widehat{T}, \widehat{\mathcal{P}}, \widehat{\Sigma})$ and
$(T, \mathcal{P}, \Sigma)$ with
$\widehat{\Sigma} = \{\widehat{\sigma}_1, \ldots, \widehat{\sigma}_{n_{\mathrm{sh}}}\}$ and
$\Sigma = \{\sigma_1, \ldots, \sigma_{n_{\mathrm{sh}}}\}$ are called **affine equivalent**
if there is an invertible affine mapping $F(\widehat{x}) = B \widehat{x} + b$, with
$B \in \mathbb{R}^{d \times d}$ invertible and $b \in \mathbb{R}^d$, such that

i) $T = F(\widehat{T})$;

ii) $\mathcal{P} = \{ \widehat{p} \circ F^{-1} : \widehat{p} \in \widehat{\mathcal{P}} \}$,
    that is, $p \in \mathcal{P}$ if and only if $p \circ F \in \widehat{\mathcal{P}}$;

iii) $\sigma_i(p) = \widehat{\sigma}_i(p \circ F)$ for all $p \in \mathcal{P}$ and
     $i = 1, \ldots, n_{\mathrm{sh}}$, after a suitable numbering of the degrees of
     freedom.
```

```{prf:lemma} Affine equivalence of the Lagrange elements
:label: lem:pk-affine-equivalence

Let $\widehat{T} = \operatorname{conv}\{\widehat{a}_0, \ldots, \widehat{a}_d\}$ be a regular
$d$-simplex, $F(\widehat{x}) = B \widehat{x} + b$ an invertible affine mapping, and
$T = F(\widehat{T})$ with vertices $a_i = F(\widehat{a}_i)$. Then $F$ maps the principal
lattice of $\widehat{T}$ onto that of $T$, with $F(\widehat{a}_\alpha) = a_\alpha$ for
all $\alpha \in \mathcal{A}_{k,d}$, and the $\mathbb{P}_k$ Lagrange elements on
$\widehat{T}$ and $T$ are affine equivalent through $F$.
```

```{prf:proof}
*Shape functions.* If $p$ is a polynomial of degree at most $k$, so is $p \circ F$: each
monomial $x^\beta$ becomes a product of $|\beta|$ affine functions of $\widehat{x}$. The
same holds for $F^{-1}$, which is affine as well. Hence $p \in \mathbb{P}_k(T)$ if and
only if $p \circ F \in \mathbb{P}_k(\widehat{T})$, which is ii) of
{prf:ref}`def:affine-equivalence`.

*Degrees of freedom.* By {prf:ref}`lem:barycentric-affine-invariance`, the barycentric
coordinates of $F(\widehat{a}_\alpha)$ with respect to $T$ coincide with those of
$\widehat{a}_\alpha$ with respect to $\widehat{T}$, namely $\frac{1}{k}(\alpha_0, \alpha)$.
By the uniqueness of barycentric coordinates, $F(\widehat{a}_\alpha) = a_\alpha$. Numbering
the degrees of freedom on both elements by $\alpha$, we obtain

$$
\sigma_\alpha(p) = p(a_\alpha) = (p \circ F)(\widehat{a}_\alpha)
= \widehat{\sigma}_\alpha(p \circ F)
\quad \text{for all } p \in \mathbb{P}_k(T) ,
$$

which is iii). Property i) holds by assumption.
```


```{figure} figures/principal-lattice-p4.svg
:label: fig-principal-lattice-p4
:alt: The unit triangle with vertices a0 at the origin, a1 at (1,0) and a2 at (0,1), subdivided by thin lines into 16 small triangles. Fifteen dots sit at the lattice points, each labelled by its multi-index; three vertex dots, nine edge dots and three interior dots are coloured differently.
:width: 55%

The principal lattice $\mathcal{L}_4(\widehat{T})$ on the unit triangle, with each point
$a_\alpha = (\alpha_1/4, \alpha_2/4)$ labelled by $(\alpha_1, \alpha_2)$. The thin lines are the
level sets $\lambda_i = j/4$, $j = 1, 2, 3$, of the barycentric coordinates, and the
lattice points are their intersections. There are $3$ vertices, $3 (k-1) = 9$ points
inside the edges and $\binom{k-1}{2} = 3$ interior points, in total
$\binom{k+2}{2} = 15 = \dim \mathbb{P}_4(\mathbb{R}^2)$.
```

```{prf:theorem} Unisolvence of the $\mathbb{P}_k$ Lagrange element
:label: thm:pk-simplex-unisolvence

The triple of {prf:ref}`def:pk-simplex` is a finite element in the sense of
{prf:ref}`def:ciarlet-finite-element`, for every regular $d$-simplex $T$ and every
$k \geqslant 1$.
```

```{prf:proof}
We follow here the presentation of @ErnGuermond2021 [Proof of Proposition 7.12]
According the previous observation,
$\# \mathcal{A}_{k,d} = \dim \mathbb{P}_k(\mathbb{R}^d)$, so it remains to show that the degrees of freedom are unisolvent.
This can by induction over the dimension $d$ and the polynomial degree $k$.
The statement that the unisolvence holds for the $\mathcal{P}_k$
Lagrange element on one (and thus by affine-invariance on any
non-degenerate) $d$-simplex will be simply abbreviated a
$[\mathcal{P}_{k,d}]$.

Induction start: $[\mathcal{P}_{k,1}]$ holds in 1 dimension for $k \geqslant 1$, see our discussion of the 1D case in 
[](#ssec:lagrange-interpolation).

Induction step: Assume that $[\mathcal{P}_{k,d-1}]$ holds for all $k \geqslant 1$ and $d \geqslant 2$,
 we want to show that $[\mathcal{P}_{k,d}]$ holds for all $k \geqslant 1$. Consider facet $F_0$ of $T$ opposite to vertex $a_0$. 
Now $a_{\alpha} \in \mathcal{L}_k(T) \cap F_0
\Leftrightarrow \lambda_0(a_\alpha) = 0 \Leftrightarrow \alpha_0 = k - |\alpha| = 0$, so the nodes that are located
on $F_0$ are precisely whose barycentric coordinates are made up from multi-indices
$\alpha = (\alpha_1, \ldots, \alpha_d)$ with $|\alpha| = k$.

Let $L_0 = \{ x \in \mathbb{R}^d \mid \lambda(x) = 0 \}$ be the
hyperplane that contains $F_0$.
Then there is a bijective affine mapping 
$F: \mathbb{R}^{d-1} \to L_0$
that maps the reference $d-1$ dimensional 
$\widehat{T} = \operatorname{conv}\{e_0, e_1, \ldots, e_{d-1}\}$ bijectively to $F = \operatorname{conv}\{a_1, \ldots, a_{d}\}$
with $F(e_i) = a_{i+1}$ for $i = 0, \ldots, d-1$.
For $\alpha \in \mathcal{A_{k,d}}$ with $|\alpha| = k$, we can 
rewrite $\alpha = (k - |\beta|, \beta_1, \ldots, \beta_{d-1})$ 
with $(\beta_1, \ldots, \beta_{d-1}) =: \beta \in \mathcal{A}_{k,d-1}$.

For $\beta \in \mathcal{A}_{k,d-1}$, we see that the points 
$b_{\beta} = F^{-1}(a_\alpha) \in \widehat{T}$ with $\alpha = (k - |\beta|, \beta_1, \ldots, \beta_{d-1})$
have barycentric coordinates
$\boldsymbol{\lambda}(b_\beta) = \frac{1}{k} (k - |\beta|, \beta_1,
\ldots, \beta_{d-1})$ as the barycentric coordinates are invariant under
affine transformations.
So $\{b_\beta\}_{\beta \in \mathcal{A}_{k,d-1}}$ is exactly the principal lattice of order $k$ of the $d-1$ dimensional reference simplex $\widehat{T}$.

Now $p_0 = p \circ F \in \mathbb{P}_k(\widehat{T})$ vanishes on the
principal lattice of order $k$ of $\widehat{T}$, so by the induction
hypothesis $[\mathcal{P}_{k,d-1}]$, we have that $p_0 \equiv 0$ and
thus $p|_{F_0} \equiv 0$. Thanks to lemma
{prf:ref}`lem:hyperplane-decomposition`, we can write $p = \lambda_0
q$ with $q \in \mathbb{P}_{k-1}(T)$. To conclude the proof, we need to
show that $q \equiv 0$, which we will do by induction on $k$.
For $k = 1$, we have that $q \in \mathbb{P}_0(T)$ is constant, and since $\lambda_0(a_0) = 1$ and $p(a_0) = 0$, 
we deduce that that $q \equiv 0$.

Assume not that $k \geqslant 2$ and that $[\mathcal{P}_{k-1,d}]$ holds.
We consider the remaining lattice points
$a_\alpha \in \mathcal{L}_k(T) \setminus F_0$ which are exactly those with
barycentric coordinates $\boldsymbol{\lambda}(a_\alpha) = \frac{1}{k} (\alpha_0, \ldots, \alpha_d)$ made up from multi-indices 
$(\alpha_0, \alpha) = (k - |\alpha|, \alpha)$ with 
$\alpha_0 \geqslant 1$ or equivalentlly $|\alpha| \leqslant k - 1$.
Our goal is to identify these points with the principal lattice of order $k-1$ of some regular $d$-simplex.
To see this, we rewrite the convex combination of $a_\alpha$ as follows,

$$
\begin{aligned}
a_\alpha 
= \sum_{i=0}^{d} \frac{\alpha_i}{k} a_i
&= a_0 + \sum_{i=1}^{d} \frac{\alpha_i}{k} (a_i - a_0)
\\
&= a_0 + \sum_{i=1}^{d} \frac{\alpha_i}{k-1} \dfrac{(k-1)}{k} (a_i - a_0)
\\
&= a_0 + \sum_{i=1}^{d} \frac{\alpha_i}{k-1} (\widetilde{a}_i - a_0)
\end{aligned}
$$

using

$$
\dfrac{(k-1)}{k} (a_i - a_0)  =
\underbrace{\dfrac{(k-1)}{k} a_i  + \dfrac{1}{k} a_0}_{=:\widetilde{a}_i} - a_0 = \widetilde{a}_i - a_0
$$


With  $\widetilde{T} = \operatorname{conv}\{a_0, \widetilde{a}_1, \ldots, \widetilde{a}_d\}$, 
we conclude 
we have that $\widetilde{T}$ is a regular $d$-simplex and the points $a_\alpha$ with $\alpha_0 \geqslant 1$ are exactly the points of the principal lattice of order $k-1$ of $\widetilde{T}.$
Since $p$ vanishes on these points, we have that $q$ vanishes on the principal lattice of order $k-1$ of $\widetilde{T}$, and by the induction hypothesis $[\mathcal{P}_{k-1,d}]$, we conclude that $q \equiv 0$ and thus $p \equiv 0$.
```

```{prf:remark} The shape functions for $k = 1$
:label: rem:pk-simplex-nodal-basis

For $k = 1$ part ii) of {prf:ref}`lem:barycentric-properties` already exhibits the nodal
basis: the shape functions are the barycentric coordinates themselves, in any dimension.
This generalises {prf:ref}`rem:p1-nodal-basis` from triangles to simplices.
```