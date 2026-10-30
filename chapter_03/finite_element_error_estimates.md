# Error estimates for the finite element approximation

We now combine the abstract error theory for Galerkin methods from
{prf:ref}`Céa's lemma<ceas-lemma>` with the
{prf:ref}`interpolation error estimate<thm:fe-interpolation-error>` from the previous
section. This gives convergence rates for the finite element approximation, first in
the natural norm of the weak formulation and then, by a duality argument, in the
$L^2$-norm.

Throughout this section, $\Omega \subset \mathbb{R}^d$, $d \in \{1, 2, 3\}$, is a bounded,
connected and polygonal ($d = 2$) or polyhedral ($d = 3$) domain, so that it can be
covered exactly by a simplicial mesh, $\overline{\Omega} = \bigcup_{T \in \mathcal{T}_h} T$.

## Error estimates in the energy norm

### An abstract error estimate

Let $V \subset H^1(\Omega)$ be a closed subspace, for instance $V = H^1_0(\Omega)$,
equipped with the norm $\|\cdot\|_{1,\Omega}$. We consider the weak formulation:
find $u \in V$ such that

$$
a(u, v) = l(v) \quad \forall v \in V,
$$ (eq:abstract-weak-form-err)

where $a(\cdot,\cdot)$ and $l(\cdot)$ satisfy the assumptions of the
{prf:ref}`Lax-Milgram theorem<thm-lax-milgram>` with boundedness constant $C_a$ and
coercivity constant $\alpha$. Given a finite element space $V_h$, the finite element
approximation $u_h \in V_h$ is defined by

$$
a(u_h, v_h) = l(v_h) \quad \forall v_h \in V_h.
$$ (eq:abstract-discrete-form-err)

:::{prf:theorem} Abstract error estimate in the energy norm
:label: thm:abstract-energy-error
Let $\{\mathcal{T}_h\}_{h>0}$ be a shape-regular family of simplicial meshes of $\Omega$,
and let $V_h$ be a finite element space on $\mathcal{T}_h$ such that

* $V_h \subset V$, i.e. the method is **conforming**;
* the elements of $V_h$ satisfy the assumptions of
  {prf:ref}`thm:fe-interpolation-error` for some $k \geqslant 1$;
* the global interpolation operator maps into the finite element space,
  $I_h v \in V_h$ for all $v \in V \cap H^{k+1}(\Omega)$.

If the solution $u$ of {eq}`eq:abstract-weak-form-err` satisfies
$u \in H^{k+1}(\Omega)$, then the finite element approximation $u_h$ satisfies

$$
\| u - u_h \|_{1,\Omega} \leqslant C \frac{C_a}{\alpha} h^{k} |u|_{k+1,\Omega},
$$ (eq:abstract-energy-error)

where $C$ depends only on $d$, $k$, the shape-regularity constant $\sigma$, the
reference element and $\operatorname{diam}(\Omega)$.
:::

:::{prf:proof}
Since $V_h \subset V$, {prf:ref}`Céa's lemma<ceas-lemma>` applies. Choosing
$v_h = I_h u \in V_h$ in {eq}`eq:ceas-lemma` gives

$$
\| u - u_h \|_{1,\Omega}
\leqslant \frac{C_a}{\alpha} \inf_{v_h \in V_h} \| u - v_h \|_{1,\Omega}
\leqslant \frac{C_a}{\alpha} \| u - I_h u \|_{1,\Omega}.
$$

The interpolation estimate {eq}`eq:global-interpolation-error` with $m = 0$ and $m = 1$
yields

$$
\| u - I_h u \|_{1,\Omega}^2
= \| u - I_h u \|_{0,\Omega}^2 + | u - I_h u |_{1,\Omega}^2
\leqslant C^2 \bigl( h^{2(k+1)} + h^{2k} \bigr) |u|_{k+1,\Omega}^2
\leqslant C^2 \bigl( 1 + \operatorname{diam}(\Omega)^2 \bigr) h^{2k} |u|_{k+1,\Omega}^2,
$$

where we used $h \leqslant \operatorname{diam}(\Omega)$. Combining both estimates proves
{eq}`eq:abstract-energy-error`.
:::

The proof shows the general pattern: Céa's lemma reduces the error of the
finite element approximation to a *best approximation error*, and the interpolation
operator provides a concrete element of $V_h$ for which this approximation error can be
estimated. The convergence rate $h^k$ is thus determined by the polynomial degree $k$
of the finite element space, provided the solution $u$ is smooth enough.

### The $\mathbb{P}_k$ Lagrange method for the Poisson problem

We now apply {prf:ref}`thm:abstract-energy-error` to the Poisson problem with
homogeneous Dirichlet boundary conditions,

$$
\left\{
\begin{alignedat}{2}
- \Delta u &= f & &\quad \text{in } \Omega, \\
         u &= 0 & &\quad \text{on } \Gamma = \partial \Omega,
\end{alignedat}
\right.
$$ (eq:poisson-dirichlet-err)

with $f \in L^2(\Omega)$. As discussed in {ref}`sec:weak-form-pdes`, the weak
formulation {eq}`eq:dirichlet-problem-II-weak` reads: find $u \in V := H^1_0(\Omega)$
such that

$$
a(u, v) := \int_{\Omega} \nabla u \cdot \nabla v = \int_{\Omega} f v =: l(v)
\quad \forall v \in V.
$$

For $k \geqslant 1$, the $\mathbb{P}_k$ Lagrange finite element space with homogeneous
Dirichlet boundary conditions is

$$
V_h := \bigl\{ v_h \in C^0(\overline{\Omega}) \mid
v_h|_T \in \mathbb{P}_k(T) \;\forall T \in \mathcal{T}_h, \;
v_h = 0 \text{ on } \Gamma \bigr\}.
$$ (eq:pk-lagrange-space-dirichlet)

Since $a(v, v) = |v|_{1,\Omega}^2$, the natural *energy norm* of this problem is
$\|v\|_a := a(v,v)^{1/2} = |v|_{1,\Omega}$, which by the
{prf:ref}`Poincaré inequality<cor:poincare>` is a norm on $H^1_0(\Omega)$, equivalent to
$\|\cdot\|_{1,\Omega}$.

:::{prf:corollary} Energy error estimate for the $\mathbb{P}_k$ Lagrange method
:label: cor:poisson-energy-error
Let $\{\mathcal{T}_h\}_{h>0}$ be a shape-regular family of simplicial meshes of $\Omega$,
let $V_h= \mathbb{P}_k^c(\mathcal{T}_h)$ be the space of continuous and elementwise polynomials of degree $k$ built from  $\mathbb{P}_k$ Lagrange space {eq}`eq:pk-lagrange-space-dirichlet`
with $k \geqslant 1$.  Let $u_h \in V_h$ be the finite element approximation of the
Poisson problem {eq}`eq:poisson-dirichlet-err`. Then

$$
| u - u_h |_{1,\Omega} = \min_{v_h \in V_h} | u - v_h |_{1,\Omega}.
$$ (eq:poisson-best-approximation)

If in addition $u \in H^{k+1}(\Omega)$, then

$$
| u - u_h |_{1,\Omega} \leqslant C h^{k} |u|_{k+1,\Omega}
\quad \text{and} \quad
\| u - u_h \|_{1,\Omega} \leqslant C (1 + C_P^2)^{1/2} h^{k} |u|_{k+1,\Omega},
$$ (eq:poisson-energy-error)

where $C$ depends only on $d$, $k$ and the shape-regularity constant $\sigma$.
:::

:::{prf:proof}
**Best approximation.** By {prf:ref}`cor:poincare`, $|\cdot|_{1,\Omega}$ is a norm on
$V = H^1_0(\Omega)$, and with respect to this norm $a(\cdot,\cdot)$ is bounded and
coercive with $C_a = \alpha = 1$, since $a(v, w) \leqslant |v|_{1,\Omega} |w|_{1,\Omega}$
by the Cauchy-Schwarz inequality and $a(v,v) = |v|_{1,\Omega}^2$. Moreover,
$V_h \subset H^1_0(\Omega)$, because continuous, piecewise polynomial functions belong
to $H^1(\Omega)$ and $v_h = 0$ on $\Gamma$. Thus {prf:ref}`Céa's lemma<ceas-lemma>`
applied with the norm $|\cdot|_{1,\Omega}$ gives
$|u - u_h|_{1,\Omega} \leqslant \inf_{v_h \in V_h} |u - v_h|_{1,\Omega}$. The reverse
inequality is trivial since $u_h \in V_h$, and the infimum is attained at $v_h = u_h$.

**Convergence rate.** We check the assumptions of {prf:ref}`thm:fe-interpolation-error`.
By {prf:ref}`lem:pk-affine-equivalence`, all $\mathbb{P}_k$ Lagrange elements are affine
equivalent to the one on a reference simplex, and $\widehat{\mathcal{P}} = \mathbb{P}_k(\widehat{T})$.
For $d \leqslant 3$ and $k \geqslant 1$, the Sobolev embedding theorem gives
$H^{k+1}(\widehat{T}) \subset H^2(\widehat{T}) \hookrightarrow C^0(\widehat{T})$
continuously, see e.g. {cite}`Evans2010` (Chapter 5). Hence the point evaluations
$\widehat{\sigma}_i(\widehat{v}) = \widehat{v}(\widehat{a}_i)$ are well-defined and bounded
on $H^{k+1}(\widehat{T})$, and so is $\widehat{I}$. The point evaluations on $T$ and
$\widehat{T}$ are related by $\sigma_i(v) = v(a_i) = \widehat{v}(\widehat{a}_i)$ for all
continuous $v$.

Next, $I_h u \in V_h$ for $u \in H^1_0(\Omega) \cap H^{k+1}(\Omega)$: on a common facet
of two neighbouring elements, both restrictions of $I_h u$ are polynomials of degree $k$
which coincide in the Lagrange nodes on that facet, and these nodes determine a polynomial
of degree $k$ on the facet uniquely, as shown in the proof of
{prf:ref}`thm:pk-simplex-unisolvence`. Hence $I_h u$ is continuous. Likewise, $I_h u$
vanishes in all Lagrange nodes on $\Gamma$, since $u = 0$ there, and thus on all of
$\Gamma$.

Now the best approximation property {eq}`eq:poisson-best-approximation` with
$v_h = I_h u$ and the interpolation estimate {eq}`eq:global-interpolation-error` with
$m = 1$ give

$$
| u - u_h |_{1,\Omega} \leqslant | u - I_h u |_{1,\Omega} \leqslant C h^{k} |u|_{k+1,\Omega}.
$$

The second estimate in {eq}`eq:poisson-energy-error` follows from the first one and the
norm equivalence {eq}`eq:poincare-norm-equivalence`.
:::

```{prf:remark}
In the energy norm, the finite element approximation is not only quasi-optimal but
*optimal*: by {eq}`eq:poisson-best-approximation`, $u_h$ is the best approximation of $u$
in $V_h$ with respect to $|\cdot|_{1,\Omega}$. In other words, $u_h$ is the orthogonal
projection of $u$ onto $V_h$ with respect to the inner product $a(\cdot,\cdot)$, which is
exactly what {prf:ref}`Galerkin orthogonality<lem-galerkin-ortho>` expresses.
```

## $L^2$ error estimates: the Aubin-Nitsche duality argument

The interpolation estimate {eq}`eq:global-interpolation-error` with $m = 0$ suggests
that the $L^2$-error of the finite element approximation should converge with the rate
$h^{k+1}$, one order higher than in the energy norm. Céa's lemma, however, only
controls the error in the $H^1$-norm, and estimating
$\|u - u_h\|_{0,\Omega} \leqslant \|u - u_h\|_{1,\Omega}$ gives away the extra power
of $h$. The **Aubin-Nitsche duality argument** recovers this extra order by testing the
error against the solution of an auxiliary problem. This requires that the auxiliary
problem has a smooth solution.

### The Poisson problem

```{prf:definition} Elliptic regularity
:label: def:elliptic-regularity-poisson

We say that the Poisson problem {eq}`eq:poisson-dirichlet-err` on $\Omega$ is
**$H^2$-regular** if there is a constant $C_{\mathrm{reg}}$ such that for every
$g \in L^2(\Omega)$, the solution $z \in H^1_0(\Omega)$ of

$$
a(v, z) = (g, v)_{\Omega} \quad \forall v \in H^1_0(\Omega)
$$ (eq:poisson-dual-problem)

satisfies $z \in H^2(\Omega)$ and

$$
|z|_{2,\Omega} \leqslant C_{\mathrm{reg}} \|g\|_{0,\Omega}.
$$ (eq:poisson-elliptic-regularity)
```

```{prf:remark}
The Poisson problem is $H^2$-regular if $\Omega$ is a **convex** polygon or polyhedron,
or if $\partial \Omega$ is smooth. It is not $H^2$-regular on non-convex polygons,
such as the L-shaped domain, where the solution has a singularity at the re-entrant
corner. Proving these statements is a topic of regularity theory for elliptic PDEs,
and we will simply assume $H^2$-regularity where needed.
```

Problem {eq}`eq:poisson-dual-problem` is called the **dual problem**. Since $a(\cdot,\cdot)$
is symmetric, it is again a Poisson problem, $-\Delta z = g$ in $\Omega$ and $z = 0$ on
$\Gamma$, now with right-hand side $g$.

:::{prf:theorem} $L^2$ error estimate for the Poisson problem
:label: thm:poisson-l2-error
Let the assumptions of {prf:ref}`cor:poisson-energy-error` hold, and assume that the
Poisson problem is $H^2$-regular. Then

$$
\| u - u_h \|_{0,\Omega} \leqslant C C_{\mathrm{reg}} \, h \, | u - u_h |_{1,\Omega}.
$$ (eq:poisson-aubin-nitsche)

In particular, if $u \in H^{k+1}(\Omega)$, then

$$
\| u - u_h \|_{0,\Omega} \leqslant C C_{\mathrm{reg}} \, h^{k+1} |u|_{k+1,\Omega}.
$$ (eq:poisson-l2-error)

Here $C$ depends only on $d$, $k$ and the shape-regularity constant $\sigma$.
:::

:::{prf:proof}
Let $e := u - u_h \in H^1_0(\Omega)$ and let $z \in H^1_0(\Omega) \cap H^2(\Omega)$ be the
solution of the dual problem {eq}`eq:poisson-dual-problem` with right-hand side $g = e$.
Testing the dual problem with $v = e$ gives

$$
\| e \|_{0,\Omega}^2 = (e, e)_{\Omega} = a(e, z).
$$

By {prf:ref}`Galerkin orthogonality<lem-galerkin-ortho>`, $a(e, v_h) = 0$ for all
$v_h \in V_h$, and in particular for $v_h = I_h z$. Together with the Cauchy-Schwarz
inequality, we obtain

$$
\| e \|_{0,\Omega}^2 = a(e, z - I_h z)
\leqslant | e |_{1,\Omega} \, | z - I_h z |_{1,\Omega}.
$$

Since $\mathbb{P}_1 \subset \mathbb{P}_k$ and $\widehat{I}$ is bounded on
$H^2(\widehat{T})$ for $d \leqslant 3$, we may apply
{prf:ref}`thm:fe-interpolation-error` with $k$ replaced by $1$ and $m = 1$, and then the
regularity estimate {eq}`eq:poisson-elliptic-regularity` with $g = e$:

$$
| z - I_h z |_{1,\Omega} \leqslant C h |z|_{2,\Omega}
\leqslant C C_{\mathrm{reg}} \, h \, \| e \|_{0,\Omega}.
$$

Hence $\| e \|_{0,\Omega}^2 \leqslant C C_{\mathrm{reg}} \, h \, |e|_{1,\Omega} \|e\|_{0,\Omega}$,
and dividing by $\|e\|_{0,\Omega}$ (the case $\|e\|_{0,\Omega} = 0$ being trivial) gives
{eq}`eq:poisson-aubin-nitsche`. Finally, {eq}`eq:poisson-l2-error` follows by inserting
the energy error estimate {eq}`eq:poisson-energy-error` into
{eq}`eq:poisson-aubin-nitsche`.
:::

The idea of the proof is that the $L^2$-error is represented as $a(e, z)$, and Galerkin
orthogonality allows us to replace $z$ by the interpolation error $z - I_h z$. The
error $e$ thus gets multiplied by a quantity which is itself of order $h$, as long as $z$
is smooth enough.

### General second order elliptic problems

The duality argument does not rely on the specific form of the Laplace operator. We now
consider the second order operator in divergence form
$\mathcal{A} u = - \nabla \cdot (A \nabla u)$ from {ref}`sec:weak-form-pdes`, with
coefficients $A \in (L^{\infty}(\Omega))^{d \times d}$ which are
{prf:ref}`uniformly elliptic<def:ellipticity>`, together with homogeneous Dirichlet
boundary conditions. The weak formulation reads: find $u \in V := H^1_0(\Omega)$ such that

$$
a(u, v) := \int_{\Omega} A \nabla u \cdot \nabla v = \int_{\Omega} f v
\quad \forall v \in V.
$$ (eq:general-elliptic-weak-err)

The bilinear form $a(\cdot,\cdot)$ is bounded and coercive on $H^1_0(\Omega)$, but it is
**not symmetric** unless $A$ is symmetric. For the duality argument, the roles of the two
arguments of $a(\cdot,\cdot)$ therefore matter.

```{prf:definition} Adjoint problem
:label: def:adjoint-problem

Given $g \in L^2(\Omega)$, the **dual** or **adjoint problem** of
{eq}`eq:general-elliptic-weak-err` is: find $z \in V$ such that

$$
a(v, z) = (g, v)_{\Omega} \quad \forall v \in V.
$$ (eq:general-dual-problem)

Since $A \nabla v \cdot \nabla z = \nabla v \cdot A^{T} \nabla z$, the dual problem
{eq}`eq:general-dual-problem` is the weak formulation of the **adjoint PDE**

$$
\left\{
\begin{alignedat}{2}
\mathcal{A}^{*} z := - \nabla \cdot (A^{T} \nabla z) &= g & &\quad \text{in } \Omega, \\
                                               z &= 0 & &\quad \text{on } \Gamma.
\end{alignedat}
\right.
$$ (eq:adjoint-pde)

We say that the adjoint problem is **$H^2$-regular** if there is a constant
$C_{\mathrm{reg}}$ such that for every $g \in L^2(\Omega)$, the solution $z$ satisfies
$z \in H^2(\Omega)$ and $|z|_{2,\Omega} \leqslant C_{\mathrm{reg}} \|g\|_{0,\Omega}$.
```

The adjoint bilinear form $a^{*}(z, v) := a(v, z)$ satisfies the assumptions of the
{prf:ref}`Lax-Milgram theorem<thm-lax-milgram>` with the same constants $C_a$ and
$\alpha$ as $a(\cdot,\cdot)$, so the adjoint problem {eq}`eq:general-dual-problem` is
well-posed. Since $A^T$ satisfies the same ellipticity condition as $A$, the adjoint
problem is again a second order elliptic problem. For a symmetric coefficient matrix
$A = A^{T}$, it coincides with the original problem, as in the Poisson case.

:::{prf:theorem} $L^2$ error estimate for general second order elliptic problems
:label: thm:general-l2-error
Let $V_h$ be the $\mathbb{P}_k$ Lagrange space {eq}`eq:pk-lagrange-space-dirichlet` on a
shape-regular family of simplicial meshes, $k \geqslant 1$, let $u_h \in V_h$ be the
finite element approximation of {eq}`eq:general-elliptic-weak-err`, and assume that the
adjoint problem is $H^2$-regular. Then

$$
\| u - u_h \|_{0,\Omega} \leqslant C \, C_a C_{\mathrm{reg}} \, h \, \| u - u_h \|_{1,\Omega}.
$$ (eq:general-aubin-nitsche)

In particular, if $u \in H^{k+1}(\Omega)$, then

$$
\| u - u_h \|_{0,\Omega} \leqslant C \frac{C_a^2}{\alpha} C_{\mathrm{reg}} \, h^{k+1} |u|_{k+1,\Omega}.
$$ (eq:general-l2-error)

Here $C$ depends only on $d$, $k$, the shape-regularity constant $\sigma$ and
$\operatorname{diam}(\Omega)$.
:::

:::{prf:proof}
Let $e := u - u_h$ and let $z \in V \cap H^2(\Omega)$ be the solution of the adjoint
problem {eq}`eq:general-dual-problem` with $g = e$. Note that the error $e$ enters the
*first* argument of $a(\cdot,\cdot)$ in the dual problem, which is exactly the argument
in which {prf:ref}`Galerkin orthogonality<lem-galerkin-ortho>` holds,
$a(e, v_h) = 0$ for all $v_h \in V_h$. Hence

$$
\| e \|_{0,\Omega}^2 = a(e, z) = a(e, z - I_h z)
\leqslant C_a \| e \|_{1,\Omega} \| z - I_h z \|_{1,\Omega}.
$$

As in the proof of {prf:ref}`thm:abstract-energy-error`, the interpolation estimate with
$k$ replaced by $1$ and $m \in \{0, 1\}$ gives
$\| z - I_h z \|_{1,\Omega} \leqslant C h |z|_{2,\Omega}$ with $C$ depending also on
$\operatorname{diam}(\Omega)$. Using the regularity of the adjoint problem, we obtain

$$
\| e \|_{0,\Omega}^2 \leqslant C \, C_a C_{\mathrm{reg}} \, h \, \| e \|_{1,\Omega} \| e \|_{0,\Omega},
$$

which proves {eq}`eq:general-aubin-nitsche`. Combining it with the energy error estimate
{eq}`eq:abstract-energy-error` from {prf:ref}`thm:abstract-energy-error`, whose
assumptions are verified exactly as in the proof of {prf:ref}`cor:poisson-energy-error`,
gives {eq}`eq:general-l2-error`.
:::

```{prf:remark}
The duality argument explains why the *adjoint* problem, and not the original one,
has to be regular: the error $e$ is placed in the first argument of $a(\cdot,\cdot)$ so
that Galerkin orthogonality can be used, which forces the dual solution $z$ into the
second argument. For the $L^2$-estimate, it is therefore the regularity of
$\mathcal{A}^{*}$ that matters. For instance, the adjoint problem is $H^2$-regular if
$\Omega$ is convex and the coefficients $a_{ij}$ are Lipschitz continuous.
```
