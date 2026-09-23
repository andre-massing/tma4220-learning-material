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

```{admonition} TODO
:class: warning

Proof of {prf:ref}`lem:dim-pk`. The handwritten notes announce a proof but the page is
left blank.
```

Next, we derive a decomposition lemma for multivariate polynomials that vanish on a
hyperplane $L = \{ x \in \mathbb{R}^n \mid \lambda(x) = 0 \}$, where
$\lambda \in (\mathbb{R}^n)'$ is a linear functional. Often, we will denote the
hyperplane simply by $\lambda$. The following lemma will play an important role when
deriving unisolvence for certain finite elements discussed afterwards.

```{prf:lemma} Decomposition of polynomials vanishing on a hyperplane
:label: lem:hyperplane-decomposition

Let $p \in \mathbb{P}_k(\mathbb{R}^n)$ vanish on a hyperplane
$L = \{ \lambda = 0 \}$. Then we can write $p = \lambda q$ where
$q \in \mathbb{P}_{k-1}(\mathbb{R}^n)$.
```

```{prf:proof}
Since $\mathbb{P}_k(\mathbb{R}^n)$ is invariant under composition with affine mappings,
we may, after a suitable affine change of coordinates, assume that

$$
L = \{ x = (\widetilde{x}, x_n) \in \mathbb{R}^n \mid x_n = 0 \} ,
$$

in other words $\lambda(x) = x_n$. Rewrite a given multindex
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
\qquad
\alpha = (\alpha_1, \ldots, \alpha_{n-1}, \alpha_n) \in \mathbb{N}^n 
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
