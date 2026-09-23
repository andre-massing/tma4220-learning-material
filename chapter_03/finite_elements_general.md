(sec:fem-general)=
# The finite element method in multiple dimensions

So far, we have constructed finite element spaces on an interval by patching
together local Lagrange interpolation polynomials. We now take a step back and
formulate what a finite element *is*, in a way which is independent of the space
dimension and of the particular polynomial space used. This abstraction is what
allows the same code, and the same theory, to serve the plethora of finite element
families in use today.

```{admonition} TODO
:class: warning

Motivation from the 1D examples. The handwritten notes reserve a section here but
leave it empty.
```

## Construction of local finite element spaces

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

```{prf:lemma} Existence of the nodal basis
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

## Polynomial spaces in several variables

```{prf:lemma} Dimension of $\mathbb{P}_k(\mathbb{R}^n)$
:label: lem:dim-pk

Let $\mathbb{P}_k(\mathbb{R}^n)$ denote the space of polynomials of degree $k$ in $n$
variables, that is,

$$
\mathbb{P}_k(\mathbb{R}^n)
= \Big\{ p : \mathbb{R}^n \to \mathbb{R} \;\Big|\;
p(x) = \sum_{|\alpha| \leqslant k} c_\alpha x^\alpha \Big\} .
$$

Then

$$
\dim \mathbb{P}_k(\mathbb{R}^n) = \binom{n+k}{k} .
$$ (eq:dim-pk)
```

````{exercise} Dimension of the polynomial space
:label: exer-dim-pk

Prove {prf:ref}`lem:dim-pk`, that is, count the multi-indices
$\alpha \in \mathbb{N}^n$ with $|\alpha| \leqslant k$, one for each monomial
$x^\alpha$ spanning $\mathbb{P}_k(\mathbb{R}^n)$.

```{admonition} Hint
:class: hint dropdown

Sort the monomials by their total degree. Writing $\mathbb{H}_l(\mathbb{R}^n)$ for the
space of **homogeneous** polynomials of degree exactly $l$, spanned by the monomials
$x^\alpha$ with $|\alpha| = l$, every polynomial splits uniquely into its homogeneous
parts, so that

$$
\mathbb{P}_k(\mathbb{R}^n) = \bigoplus_{l=0}^{k} \mathbb{H}_l(\mathbb{R}^n),
\qquad \text{and hence} \qquad
\dim \mathbb{P}_k(\mathbb{R}^n) = \sum_{l=0}^{k} \dim \mathbb{H}_l(\mathbb{R}^n) .
$$

It therefore suffices to count the multi-indices with $|\alpha| = l$, that is, the ways
of distributing $l$ units among $n$ variables. Show that

$$
\dim \mathbb{H}_l(\mathbb{R}^n) = \binom{n + l - 1}{l} ,
$$

for instance by induction on $n$, or by the following combinatorial argument: a
multi-index $\alpha$ with $|\alpha| = l$ is the same as an arrangement of $l$ stars and
$n-1$ bars in a row, the bars splitting the stars into the $n$ groups
$\alpha_1, \ldots, \alpha_n$. Such an arrangement consists of $l + n - 1$ symbols and is
fixed as soon as one decides which $l$ of them are the stars, which leaves
$\binom{n+l-1}{l}$ possibilities.

Summing over $l$ then gives

$$
\dim \mathbb{P}_k(\mathbb{R}^n)
= \sum_{l=0}^{k} \binom{n+l-1}{l}
= \binom{n+k}{k} ,
$$

where the last identity follows from repeated application of Pascal's rule
$\binom{m}{j} = \binom{m-1}{j-1} + \binom{m-1}{j}$.
```
````

Next, we derive a decomposition lemma for multivariate polynomials that vanish on a
hyperplane. Recall that an **affine function** $\lambda : \mathbb{R}^n \to \mathbb{R}$ is
a linear function plus a constant,

$$
\lambda(x) = w \cdot x + c ,
\qquad w \in \mathbb{R}^n \setminus \{0\}, \quad c \in \mathbb{R} ,
$$ (eq:affine-function)

so that $\lambda$ is linear precisely when $c = 0$. Its zero set

$$
L = \{ x \in \mathbb{R}^n \mid \lambda(x) = 0 \}
$$

is a hyperplane, and conversely every hyperplane is the zero set of such a $\lambda$; it
passes through the origin if and only if $\lambda$ is linear. Often, we will denote the
hyperplane simply by $\lambda$. The following lemma will play an important role when
deriving unisolvence for certain finite elements discussed afterwards.

```{prf:lemma} Decomposition of polynomials vanishing on a hyperplane
:label: lem:hyperplane-decomposition

Let $\lambda$ be an affine function as in {eq}`eq:affine-function` and let
$p \in \mathbb{P}_k(\mathbb{R}^n)$ vanish on the hyperplane $L = \{ \lambda = 0 \}$. Then we can write $p = \lambda q$ where
$q \in \mathbb{P}_{k-1}(\mathbb{R}^n)$.
```

```{prf:proof}
Since $\mathbb{P}_k(\mathbb{R}^n)$ is invariant under composition with affine mappings,
we may, after a suitable affine change of coordinates, assume that

$$
L = \{ x = (\widetilde{x}, x_n) \in \mathbb{R}^n \mid x_n = 0 \} ,
$$

in other words $\lambda(x) = x_n$; rescaling $\lambda$ only rescales $q$, so this
is no restriction. Rewrite a given multi-index
$\alpha = (\alpha_1, \ldots, \alpha_{n-1}, \alpha_n) \in \mathbb{N}^n$
as $\alpha = (\widetilde{\alpha}, \alpha_n)$ 
then we can write any polynomial $p \in \mathbb{P}_k(\mathbb{R}^n)$ as

$$
\begin{aligned}
p(x)
&= \sum_{|\alpha| \leqslant k} c_\alpha x^\alpha
\\
&= \sum_{\alpha_n = 0}^{k} \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_n}
  c_{\widetilde{\alpha}, \alpha_n} \, \widetilde{x}^{\widetilde{\alpha}} x_n^{\alpha_n}
  \\
&= 
\sum_{|\widetilde{\alpha}| \leqslant k}
  c_{\widetilde{\alpha}, 0} \, \widetilde{x}^{\widetilde{\alpha}}
+
\sum_{\alpha_n = 1}^{k} \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_n}
  c_{\widetilde{\alpha}, \alpha_n} \, \widetilde{x}^{\widetilde{\alpha}} x_n^{\alpha_n} ,
\end{aligned}
$$

Evaluating on the hyperplane then gives

$$
0 = p(\widetilde{x}, 0)
= \underbrace{\sum_{|\widetilde{\alpha}| \leqslant k}
  c_{\widetilde{\alpha}, 0} \, \widetilde{x}^{\widetilde{\alpha}}}_{=: p_L(\widetilde{x})}
+ \underbrace{\sum_{\alpha_n = 1}^k \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_n}
  c_{\widetilde{\alpha}, \alpha_n} \, \widetilde{x}^{\widetilde{\alpha}} \, 0^{\alpha_n}}_{= 0} .
$$

Since $p_L(\widetilde{x}) = 0$ for all $\widetilde{x} \in \mathbb{R}^{n-1}$, we conclude
that $c_{\widetilde{\alpha}, 0} = 0$ for all $|\widetilde{\alpha}| \leqslant k$.
Therefore

$$
p(x)
= \sum_{\alpha_n = 1}^{k} \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_n}
  c_{\widetilde{\alpha}, \alpha_n} \, \widetilde{x}^{\widetilde{\alpha}} x_n^{\alpha_n}
= \underbrace{x_n}_{= \lambda(x)} \cdot
  \underbrace{\sum_{\alpha_n = 1}^{k} \; \sum_{|\widetilde{\alpha}| \leqslant k - \alpha_n}
  c_{\widetilde{\alpha}, \alpha_n} \, \widetilde{x}^{\widetilde{\alpha}} x_n^{\alpha_n - 1}}_{=: q(x)} ,
$$

and $q \in \mathbb{P}_{k-1}(\mathbb{R}^n)$, since every remaining monomial has degree
$|\widetilde{\alpha}| + \alpha_n - 1 \leqslant k - 1$.
```

## Examples of finite element on triangles

Throughout this section, $T \subset \mathbb{R}^2$ is a non-degenerate triangle with
vertices $a_1, a_2, a_3$, and $\lambda_1, \lambda_2, \lambda_3$ denote its **barycentric
coordinates**: $\lambda_i$ is the affine function, in the sense of
{eq}`eq:affine-function`, with $\lambda_i(a_j) = \delta_{ij}$, so
that $\lambda_i$ vanishes precisely on the line $L_i$ through the two vertices other
than $a_i$.

### The $\mathbb{P}_1$ element on a triangle

```{prf:definition} The linear Lagrange triangle
:label: def:p1-triangle

The **$\mathbb{P}_1$ element** on $T$ is the triple $(T, \mathcal{P}, \Sigma)$ with
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

### The $\mathbb{P}_2$ element on a triangle

```{prf:definition} The quadratic Lagrange triangle
:label: def:p2-triangle

Let $a_{ij} = \frac{1}{2}(a_i + a_j)$ denote the midpoint of the edge joining $a_i$ and
$a_j$. The **$\mathbb{P}_2$ element** on $T$ is the triple $(T, \mathcal{P}, \Sigma)$ with
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

````{exercise} Unisolvence of the $\mathbb{P}_2$ element
:label: exer-p2-unisolvence

Show that the triple of {prf:ref}`def:p2-triangle` is a finite element.

```{admonition} Hint
:class: hint dropdown

Follow the proof of {prf:ref}`lem:p1-unisolvence`. Each edge of $T$ now carries three
nodes, two vertices and one midpoint, so the restriction of $p \in \mathbb{P}_2(T)$ to
the line $L_i$ is a univariate quadratic with three distinct roots and vanishes
identically. Apply {prf:ref}`lem:hyperplane-decomposition` once for each edge to peel off
the factors $\lambda_1, \lambda_2, \lambda_3$ one at a time, and compare the resulting
degree with $2$.
```
````

### The $\mathbb{P}_3$ element on a triangle

```{prf:definition} The cubic Lagrange triangle
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

````{exercise} Unisolvence of the $\mathbb{P}_3$ element
:label: exer-p3-unisolvence

Show that the triple of {prf:ref}`def:p3-triangle` is a finite element.

```{admonition} Hint
:class: hint dropdown

Each edge now carries four nodes, so the restriction of $p \in \mathbb{P}_3(T)$ to each
line $L_i$ is a univariate cubic with four distinct roots. Peeling off all three
barycentric coordinates with {prf:ref}`lem:hyperplane-decomposition` leaves
$p = c \, \lambda_1 \lambda_2 \lambda_3$ with a constant $c$, since the degrees already
match. What does the remaining degree of freedom at the barycenter give you?
```
````

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
