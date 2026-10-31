# Inverse estimates and applications 

The interpolation estimates of {prf:ref}`thm:fe-interpolation-error` bound how well
*smooth* functions can be approximated by finite element functions. In this section, we
consider the reverse question: how do different norms of a *finite element function*
compare with each other? The answer, the so-called inverse estimates, is a versatile
tool. We use it to analyse the $L^2$-projection onto finite element spaces and to estimate
the condition number of the stiffness matrix.

As in the previous sections, $\Omega \subset \mathbb{R}^d$, $d \in \{1,2,3\}$, is a bounded
polygonal or polyhedral domain, $\{\mathcal{T}_h\}_{h>0}$ is a shape-regular family of
simplicial meshes of $\Omega$ with shape-regularity constant $\sigma$, and
$h := \max_{T \in \mathcal{T}_h} h_T$.
For the $L^2$-norms on an element and on the domain, we use the shorthand notation
$\| v \|_{T}$ and $\| v \|_{\Omega}$ from {eq}`eq:l2-norm-shorthand`.

## Inverse estimates

### A local inverse estimate

:::{prf:lemma} Inverse estimate
:label: lem:inverse-estimate
For every $T \in \mathcal{T}_h$, let $F_T(\widehat{x}) = B_T \widehat{x} + b_T$ map a
fixed reference simplex $\widehat{T}$ onto $T$, and let $\mathcal{P}$ be a space of
functions on $T$ such that $v_h \circ F_T \in \widehat{\mathcal{P}}$ for all
$v_h \in \mathcal{P}$, for a fixed finite dimensional space
$\widehat{\mathcal{P}} \subset H^1(\widehat{T})$. This holds, for instance, for
$\mathcal{P} = \mathbb{P}_k(T)$ and $\widehat{\mathcal{P}} = \mathbb{P}_k(\widehat{T})$.
Then there is a constant $C_I$, depending only on $d$, $\sigma$, $\widehat{T}$ and
$\widehat{\mathcal{P}}$, such that

$$
\| \nabla v_h \|_{T} \leqslant C_I h_T^{-1} \| v_h \|_{T}
\quad \forall v_h \in \mathcal{P}, \; \forall T \in \mathcal{T}_h.
$$ (eq:inverse-estimate)

:::

:::{prf:proof}
**Step 1: Finite dimensional argument on $\widehat{T}$.**
Both $\|\cdot\|_{1,\widehat{T}}$ and $\|\cdot\|_{0,\widehat{T}}$ are norms on the
finite dimensional space $\widehat{\mathcal{P}}$. Since all norms on a finite dimensional
space are equivalent, there is a constant $\widehat{C}$, depending only on
$\widehat{T}$ and $\widehat{\mathcal{P}}$, such that

$$
| \widehat{v} |_{1,\widehat{T}}
\leqslant \| \widehat{v} \|_{1,\widehat{T}}
\leqslant \widehat{C} \| \widehat{v} \|_{0,\widehat{T}}
\quad \forall \widehat{v} \in \widehat{\mathcal{P}}.
$$

**Step 2: Transformation to $T$.**
Let $v_h \in \mathcal{P}$ and $\widehat{v} := v_h \circ F_T \in \widehat{\mathcal{P}}$.
By the {prf:ref}`transformation formula<lem:transformation-sobolev-seminorms>`
{eq}`eq:transformation-formula-sobolev-seminorms`, used with $l = 1$ for the
first inequality below and with $l = 0$ for the last one, we obtain

$$
\| \nabla v_h \|_{T} = | v_h |_{1,T}
\leqslant C \| B_T^{-1} \|_2 \, |\det B_T|^{1/2} \, | \widehat{v} |_{1,\widehat{T}}
\leqslant C \widehat{C} \| B_T^{-1} \|_2 \, |\det B_T|^{1/2} \, \| \widehat{v} \|_{0,\widehat{T}}
\leqslant C \widehat{C} \| B_T^{-1} \|_2 \, \| v_h \|_{T}.
$$

**Step 3: Shape regularity.**
By {prf:ref}`lem:affine-map-norm` and $h_T / \rho_T \leqslant \sigma$,

$$
\| B_T^{-1} \|_2 \leqslant \frac{h_{\widehat{T}}}{\rho_T}
= \frac{h_{\widehat{T}}}{h_T} \frac{h_T}{\rho_T}
\leqslant h_{\widehat{T}} \, \sigma \, h_T^{-1},
$$

which proves {eq}`eq:inverse-estimate` with $C_I := C \widehat{C} h_{\widehat{T}} \sigma$.
:::

```{prf:remark} Why "inverse" estimate?
:label: rem:inverse-estimate-name

Recall the {prf:ref}`Poincaré inequality<thm:poincare>`

$$
\| v \|_{\Omega} \leqslant C_P \| \nabla v \|_{\Omega}
\quad \forall v \in H^1_0(\Omega),
$$

where a scaling argument shows that $C_P \sim \operatorname{diam}(\Omega)$. The inverse
estimate {eq}`eq:inverse-estimate` looks like a local version of this inequality, read
in the *inverse* direction: the gradient is bounded by the function itself, and the
length scale $\operatorname{diam}(\Omega)$ is replaced by the local length scale $h_T$.

Such an estimate cannot hold for general $v \in H^1(T)$: the functions
$v_n(x) = \sin(n x_1)$ satisfy $\|v_n\|_T \leqslant |T|^{1/2}$, while
$\|\nabla v_n\|_T \to \infty$ as $n \to \infty$. The inverse estimate holds only
because $v_h$ is taken from a finite dimensional space, and its constant $C_I h_T^{-1}$
blows up as the mesh is refined.
```

### The notation $\lesssim$

Keeping track of all constants quickly becomes tedious, and their precise values are
rarely of interest. From now on, we therefore write

$$
a \lesssim b
\quad :\Longleftrightarrow \quad
a \leqslant C b
$$

for a constant $C > 0$ that is **independent of the mesh size $h$, the element $T$ and
the functions involved**, and we write $a \sim b$ if both $a \lesssim b$ and $b \lesssim a$
hold. The hidden constant $C$ may still depend on other quantities, and we will always
state on which ones. Unless said otherwise, the hidden constants in this section depend
only on the dimension $d$, the polynomial degree $k$, the shape-regularity constant
$\sigma$ and the reference element. For instance, the inverse estimate
{eq}`eq:inverse-estimate` now reads

$$
\| \nabla v_h \|_{T} \lesssim h_T^{-1} \| v_h \|_{T}
\quad \forall v_h \in \mathcal{P}.
$$

### General inverse estimates

The proof of {prf:ref}`lem:inverse-estimate` extends to higher derivatives and to
traces on facets. In the following exercises, we write
$\| D^j v \|_{T} := | v |_{j,T}$ for the $L^2$-norm of all derivatives of order $j$.

````{exercise} Inverse estimates for higher derivatives
:label: exer-inverse-estimate-higher

Let the assumptions of {prf:ref}`lem:inverse-estimate` hold with
$\widehat{\mathcal{P}} = \mathbb{P}_k(\widehat{T})$. Show that

$$
\| D^j v_h \|_{T} \lesssim h_T^{-i} \| D^{j-i} v_h \|_{T}
\quad \forall v_h \in \mathcal{P}, \; 0 \leqslant i \leqslant j.
$$

*Hint:* On $\widehat{T}$, the seminorm $|\cdot|_{j-i,\widehat{T}}$ is in general only a
seminorm on $\mathbb{P}_k(\widehat{T})$, since it vanishes on $\mathbb{P}_{j-i-1}$.
Show first that
$|\widehat{v}|_{j,\widehat{T}} \leqslant \widehat{C} |\widehat{v}|_{j-i,\widehat{T}}$ for
all $\widehat{v} \in \mathbb{P}_k(\widehat{T})$. To this end, note that both sides
vanish on $\mathbb{P}_{j-i-1}(\widehat{T})$, and that $|\cdot|_{j-i,\widehat{T}}$ is a
norm on any complement $W$ of this subspace in $\mathbb{P}_k(\widehat{T})$, so that the
finite dimensional argument applies on $W$. Then proceed as in the proof of {prf:ref}`lem:inverse-estimate`, using both
inequalities of {eq}`eq:transformation-formula-sobolev-seminorms`.
````

````{exercise} Inverse trace estimate
:label: exer-inverse-trace-estimate

Let the assumptions of {prf:ref}`lem:inverse-estimate` hold with
$\widehat{\mathcal{P}} = \mathbb{P}_k(\widehat{T})$, let $F$ be a facet of
$T \in \mathcal{T}_h$ with unit normal $n$, and write $\partial_n v_h := \nabla v_h \cdot n$.
Show that

$$
\| \partial_n v_h \|_{F} \lesssim h_T^{-1/2} \| \nabla v_h \|_{T}
\quad \forall v_h \in \mathcal{P}.
$$

*Hint:* Since $|\partial_n v_h| \leqslant |\nabla v_h|$, it suffices to show
$\| w \|_{F} \lesssim h_T^{-1/2} \| w \|_{T}$ for every component $w$ of $\nabla v_h$,
which is a polynomial of degree $k-1$. On the reference element, use that
$\|\cdot\|_{\widehat{F}}$ and $\|\cdot\|_{\widehat{T}}$ are a seminorm and a norm on
$\mathbb{P}_{k-1}(\widehat{T})$, respectively. Then transform the integrals over $F$ and
$T$ to $\widehat{F} = F_T^{-1}(F)$ and $\widehat{T}$, and use that
$|F| / |\widehat{F}| \lesssim h_T^{d-1}$ and $|T| / |\widehat{T}| \gtrsim h_T^{d}$
by shape regularity.
````

## The $L^2$-projection

### Definition and well-posedness

In this section, let $k \geqslant 1$ and let

$$
V_h := \bigl\{ v_h \in C^0(\overline{\Omega}) \mid
v_h|_T \in \mathbb{P}_k(T) \;\forall T \in \mathcal{T}_h \bigr\}
$$ (eq:pk-lagrange-space-no-bc)

be the $\mathbb{P}_k$ Lagrange space without boundary conditions.

:::{prf:definition} $L^2$-projection
:label: def:l2-projection
The **$L^2$-projection** $\Pi_h : L^2(\Omega) \to V_h$ is defined by: for
$u \in L^2(\Omega)$, find $\Pi_h u \in V_h$ such that

$$
(\Pi_h u, v_h)_{\Omega} = (u, v_h)_{\Omega} \quad \forall v_h \in V_h.
$$ (eq:l2-projection)

:::

:::{prf:lemma} Properties of the $L^2$-projection
:label: lem:l2-projection-properties
For every $u \in L^2(\Omega)$, problem {eq}`eq:l2-projection` has a unique solution
$\Pi_h u \in V_h$, and $\Pi_h : L^2(\Omega) \to V_h$ is linear. Moreover,

1. **Projection property:** $\Pi_h v_h = v_h$ for all $v_h \in V_h$.
2. **Orthogonality:** $(u - \Pi_h u, v_h)_{\Omega} = 0$ for all $v_h \in V_h$.
3. **Stability:** $\| \Pi_h u \|_{\Omega} \leqslant \| u \|_{\Omega}$.
4. **Best approximation:**
   $\| u - \Pi_h u \|_{\Omega} = \min_{v_h \in V_h} \| u - v_h \|_{\Omega}$.
:::

:::{prf:proof}
**Well-posedness.** We apply the {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>` on the
Hilbert space $V_h$ with the inner product and norm of $L^2(\Omega)$. Since $V_h$ is finite
dimensional, it is a closed subspace of $L^2(\Omega)$ and hence a Hilbert space. The
bilinear form $(w_h, v_h) \mapsto (w_h, v_h)_{\Omega}$ is bounded with $C_a = 1$ by the
Cauchy-Schwarz inequality and coercive with $\alpha = 1$. The linear form
$v_h \mapsto (u, v_h)_{\Omega}$ is bounded with constant $\|u\|_{\Omega}$. Hence
{eq}`eq:l2-projection` has a unique solution, and the stability estimate
{eq}`eq:lax-milgram-stab` gives 3. Linearity follows from uniqueness.

**Projection property and orthogonality.** If $u = v_h \in V_h$, then $v_h$ itself solves
{eq}`eq:l2-projection`, so $\Pi_h v_h = v_h$ by uniqueness. Property 2 is just a
reformulation of {eq}`eq:l2-projection`.

**Best approximation.** As in the proof of {prf:ref}`Céa's lemma<ceas-lemma>`, we add
and subtract an arbitrary $v_h \in V_h$ and use orthogonality with
$v_h - \Pi_h u \in V_h$:

$$
\| u - \Pi_h u \|_{\Omega}^2
= (u - \Pi_h u, u - v_h)_{\Omega} + \underbrace{(u - \Pi_h u, v_h - \Pi_h u)_{\Omega}}_{= 0}
\leqslant \| u - \Pi_h u \|_{\Omega} \| u - v_h \|_{\Omega}.
$$

Dividing by $\| u - \Pi_h u \|_{\Omega}$ (the case $\| u - \Pi_h u \|_{\Omega} = 0$
being trivial) gives $\| u - \Pi_h u \|_{\Omega} \leqslant \| u - v_h \|_{\Omega}$.

:::

In terms of the nodal basis $\{\phi_i\}_{i=1}^N$ of $V_h$, problem {eq}`eq:l2-projection`
is the linear system $M U = b$ with the **mass matrix** $M_{ij} = (\phi_j, \phi_i)_{\Omega}$
and $b_i = (u, \phi_i)_{\Omega}$, where $\Pi_h u = \sum_{j} U_j \phi_j$.

### Error estimates

For the $H^1$-estimate, we need a stronger assumption on the mesh than shape regularity,
namely that the family of meshes is
{prf:ref}`quasi-uniform<def:shape-regular-quasi-uniform>` with constant $c_q$.

On quasi-uniform meshes, all elements have comparable size, $h_T \sim h$, and the local
inverse estimate turns into a global one.

:::{prf:corollary} Global inverse estimate
:label: cor:global-inverse-estimate
Let $\{\mathcal{T}_h\}_{h>0}$ be quasi-uniform. Then

$$
| v_h |_{1,\Omega} \lesssim h^{-1} \| v_h \|_{\Omega} \quad \forall v_h \in V_h,
$$ (eq:global-inverse-estimate)

where the hidden constant depends only on $d$, $k$, $\sigma$ and $c_q$.
:::

:::{prf:proof}
Since $V_h \subset H^1(\Omega)$, we can compute $|v_h|_{1,\Omega}$ elementwise. Applying
{eq}`eq:inverse-estimate` on each element and using $h_T^{-1} \leqslant c_q^{-1} h^{-1}$
gives

$$
| v_h |_{1,\Omega}^2 = \sum_{T \in \mathcal{T}_h} \| \nabla v_h \|_{T}^2
\leqslant C_I^2 \sum_{T \in \mathcal{T}_h} h_T^{-2} \| v_h \|_{T}^2
\leqslant C_I^2 c_q^{-2} h^{-2} \| v_h \|_{\Omega}^2.
$$

:::

:::{prf:theorem} Error estimates for the $L^2$-projection
:label: thm:l2-projection-error
Let $d \leqslant 3$, $k \geqslant 1$ and $u \in H^{k+1}(\Omega)$. Then

$$
\| u - \Pi_h u \|_{\Omega} \lesssim h^{k+1} |u|_{k+1,\Omega}.
$$ (eq:l2-projection-error-l2)

If, in addition, the family of meshes is quasi-uniform, then

$$
| u - \Pi_h u |_{1,\Omega} \lesssim h^{k} |u|_{k+1,\Omega}.
$$ (eq:l2-projection-error-h1)

The hidden constants depend only on $d$, $k$, $\sigma$, and for
{eq}`eq:l2-projection-error-h1` also on $c_q$.
:::

:::{prf:proof}
Since $d \leqslant 3$ and $k \geqslant 1$, the Lagrange interpolant $I_h u \in V_h$ is
well-defined, and it satisfies the interpolation estimates
{eq}`eq:global-interpolation-error`, as shown in the proof of
{prf:ref}`cor:poisson-energy-error`.

**$L^2$-estimate.** We add and subtract the interpolant and use the projection property
$I_h u = \Pi_h (I_h u)$:

$$
u - \Pi_h u
= (u - I_h u) + (I_h u - \Pi_h u)
= (u - I_h u) + \Pi_h (I_h u - u).
$$ (eq:l2-projection-splitting)

Using the stability of $\Pi_h$ and the interpolation estimate with $m = 0$, we obtain

$$
\| u - \Pi_h u \|_{\Omega}
\leqslant \| u - I_h u \|_{\Omega} + \| \Pi_h (I_h u - u) \|_{\Omega}
\leqslant 2 \| u - I_h u \|_{\Omega}
\lesssim h^{k+1} |u|_{k+1,\Omega}.
$$

**$H^1$-estimate.** We use the same splitting {eq}`eq:l2-projection-splitting`. Since
$\Pi_h (I_h u - u) \in V_h$, the global inverse estimate
{eq}`eq:global-inverse-estimate`, the stability of $\Pi_h$ and the interpolation estimates
with $m = 1$ and $m = 0$ give

$$
\begin{aligned}
| u - \Pi_h u |_{1,\Omega}
&\leqslant | u - I_h u |_{1,\Omega} + | \Pi_h (I_h u - u) |_{1,\Omega}
\\
&\lesssim | u - I_h u |_{1,\Omega} + h^{-1} \| \Pi_h (I_h u - u) \|_{\Omega}
\\
&\leqslant | u - I_h u |_{1,\Omega} + h^{-1} \| I_h u - u \|_{\Omega}
\\
&\lesssim h^{k} |u|_{k+1,\Omega} + h^{-1} h^{k+1} |u|_{k+1,\Omega}
\lesssim h^{k} |u|_{k+1,\Omega}.
\end{aligned}
$$

:::

```{prf:remark} The $L^2$-projection is not local
:label: rem:l2-projection-nonlocal

The Lagrange interpolant is a *local* operator: $(I_h u)|_T$ depends only on the values of
$u$ at the nodes of $T$. The $L^2$-projection, in contrast, is in general *not* local.
Computing $\Pi_h u$ requires solving a linear system with the mass matrix $M$, and although
$M$ is sparse, its inverse $M^{-1}$ is in general a dense matrix. Changing $u$ on a single
element therefore changes $\Pi_h u$ on the whole domain, even though the effect decays
quickly away from that element.
```

```{prf:remark} On the regularity assumption
:label: rem:l2-projection-regularity

The $L^2$-projection is well-defined for every $u \in L^2(\Omega)$. The assumptions
$u \in H^{k+1}(\Omega)$ and $d \leqslant 3$ in {prf:ref}`thm:l2-projection-error` are only
needed because the proof uses the Lagrange interpolant, which requires point values of $u$.
Other proofs avoid this, for instance by replacing $I_h$ with a *quasi-interpolation
operator* such as the Clément or the Scott-Zhang operator, which only use local averages
of $u$. In this way one obtains, for instance,
$\| u - \Pi_h u \|_{\Omega} \lesssim h^{s} |u|_{s,\Omega}$ for all
$u \in H^{s}(\Omega)$, $1 \leqslant s \leqslant k+1$, see e.g. {cite}`ErnGuermond2021`.
```

## Condition number estimates for the stiffness matrix

### Stiffness matrix and condition number

We consider the $\mathbb{P}_k$ Lagrange method for the Poisson problem
{eq}`eq:poisson-dirichlet-err` with homogeneous Dirichlet boundary conditions, with the
finite element space $V_h \subset H^1_0(\Omega)$ from
{eq}`eq:pk-lagrange-space-dirichlet`. Let $\{\phi_i\}_{i=1}^N$ be the nodal basis of
$V_h$, which consists of the global Lagrange basis functions associated with the nodes
$x_i$ in the *interior* of $\Omega$. As in {eq}`discrete-system-derivation`, the
**stiffness matrix** $A \in \mathbb{R}^{N \times N}$ is given by

$$
A_{ij} = a(\phi_j, \phi_i) = \int_{\Omega} \nabla \phi_j \cdot \nabla \phi_i.
$$

Every $v_h \in V_h$ has a unique representation $v_h = \sum_{i=1}^N V_i \phi_i$ with
coefficient vector $V = (V_i)_{i=1}^N \in \mathbb{R}^N$, and

$$
V^{T} A V = a(v_h, v_h) = | v_h |_{1,\Omega}^2.
$$ (eq:stiffness-energy)

Hence $A$ is symmetric and positive definite.

The **condition number** of an invertible matrix $A$ with respect to the Euclidean norm
$|\cdot|$ on $\mathbb{R}^N$ is

$$
\kappa(A) := \| A \|_2 \, \| A^{-1} \|_2.
$$

It measures how much relative perturbations of the right-hand side $b$ can be amplified
in the solution of $A U = b$, and it governs the convergence speed of iterative solvers
such as the conjugate gradient method. For a symmetric positive definite matrix,
$\| A \|_2 = \lambda_{\max}(A)$ and $\| A^{-1} \|_2 = \lambda_{\min}(A)^{-1}$, and the
extreme eigenvalues are given by the Rayleigh quotients

$$
\lambda_{\max}(A) = \max_{V \neq 0} \frac{V^{T} A V}{|V|^2},
\qquad
\lambda_{\min}(A) = \min_{V \neq 0} \frac{V^{T} A V}{|V|^2},
\qquad
\kappa(A) = \frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}.
$$ (eq:rayleigh-quotients)

### The condition number estimate

The estimate is based on three ingredients. The first two are already available: the
Poincaré inequality, which holds in particular for all $v_h \in V_h \subset H^1_0(\Omega)$,
and the inverse estimate. The third one relates the $L^2$-norm of a finite element
function to the Euclidean norm of its coefficient vector.

:::{prf:lemma} Equivalence of $L^2$-norm and scaled Euclidean norm
:label: lem:l2-euclidean-equivalence
Let $\{\mathcal{T}_h\}_{h>0}$ be quasi-uniform. Then

$$
h^{d} |V|^2 \lesssim \| v_h \|_{\Omega}^2 \lesssim h^{d} |V|^2
\quad \forall v_h = \sum_{i=1}^N V_i \phi_i \in V_h,
$$ (eq:l2-euclidean-equivalence)

where the hidden constants depend only on $d$, $k$, $\sigma$ and $c_q$.
:::

````{exercise} Proof of the norm equivalence
:label: exer-l2-euclidean-equivalence

Prove {prf:ref}`lem:l2-euclidean-equivalence`. Proceed in the following steps.

1. **Reference element.** For $\widehat{v} = \sum_{j} \widehat{V}_j \widehat{\theta}_j$ in
   the nodal basis of $\mathbb{P}_k(\widehat{T})$, show that
   $\| \widehat{v} \|_{\widehat{T}}^2 \sim \sum_j \widehat{V}_j^2$.
   *Hint:* Both sides define a norm on the finite dimensional space $\mathbb{R}^{n_{\mathrm{sh}}}$
   of coefficient vectors.
2. **Single element.** For $T \in \mathcal{T}_h$, let $V_T$ denote the vector of the
   coefficients of $v_h$ associated with the nodes of $T$, where the coefficients of
   boundary nodes are $0$. Show that $\| v_h \|_{T}^2 \sim h_T^{d} |V_T|^2$.
   *Hint:* Use Step 1 together with
   $\| v_h \|_{T}^2 = |\det B_T| \, \| \widehat{v} \|_{\widehat{T}}^2$ and
   $|\det B_T| = |T| / |\widehat{T}| \sim h_T^{d}$, which follows from
   $\rho_T^{d} \lesssim |T| \leqslant h_T^{d}$ and shape regularity.
3. **Summation.** Sum over all elements. Use that every node belongs to at least one and,
   by shape regularity, at most a bounded number of elements, and that $h_T \sim h$ by
   quasi-uniformity.
````

:::{prf:theorem} Condition number of the stiffness matrix
:label: thm:stiffness-condition-number
Let $\{\mathcal{T}_h\}_{h>0}$ be quasi-uniform and let $A$ be the stiffness matrix of the
$\mathbb{P}_k$ Lagrange method for the Poisson problem with homogeneous Dirichlet boundary
conditions. Then

$$
\lambda_{\max}(A) \lesssim h^{d-2},
\qquad
\lambda_{\min}(A) \gtrsim h^{d},
\qquad
\kappa(A) \lesssim h^{-2},
$$ (eq:stiffness-condition-number)

where the hidden constants depend only on $d$, $k$, $\sigma$, $c_q$ and the Poincaré
constant $C_P$ of $\Omega$.
:::

:::{prf:proof}
Let $V \in \mathbb{R}^N \setminus \{0\}$ and $v_h = \sum_{i} V_i \phi_i \in V_h$. By
{eq}`eq:stiffness-energy`, $V^{T} A V = |v_h|_{1,\Omega}^2$.

**Upper bound.** The global inverse estimate {eq}`eq:global-inverse-estimate` and the
upper bound in {eq}`eq:l2-euclidean-equivalence` give

$$
V^{T} A V = | v_h |_{1,\Omega}^2 \lesssim h^{-2} \| v_h \|_{\Omega}^2
\lesssim h^{-2} h^{d} |V|^2 = h^{d-2} |V|^2.
$$

**Lower bound.** The {prf:ref}`Poincaré inequality<thm:poincare>`, applied to
$v_h \in H^1_0(\Omega)$, and the lower bound in {eq}`eq:l2-euclidean-equivalence` give

$$
V^{T} A V = | v_h |_{1,\Omega}^2 \geqslant C_P^{-2} \| v_h \|_{\Omega}^2
\gtrsim h^{d} |V|^2.
$$

By the Rayleigh quotients {eq}`eq:rayleigh-quotients`, $\lambda_{\max}(A) \lesssim h^{d-2}$
and $\lambda_{\min}(A) \gtrsim h^{d}$, and thus
$\kappa(A) = \lambda_{\max}(A) / \lambda_{\min}(A) \lesssim h^{-2}$.
:::

```{prf:remark}
:label: rem:condition-number-sharp

The estimate $\kappa(A) \lesssim h^{-2}$ cannot be improved in general. For $\mathbb{P}_1$
elements on a uniform mesh of $\Omega = (0,1)$ with $N+1$ elements of size
$h = 1/(N+1)$, the stiffness matrix is $A = h^{-1} \operatorname{tridiag}(-1, 2, -1)$, with
eigenvalues $\lambda_j = 4 h^{-1} \sin^2(j \pi h / 2)$, $j = 1, \ldots, N$. Hence

$$
\kappa(A) = \frac{\sin^2(N \pi h / 2)}{\sin^2(\pi h / 2)}
= \cot^2(\pi h / 2) \sim \frac{4}{\pi^2} h^{-2}.
$$

As a consequence, the number of conjugate gradient iterations needed to reach a fixed
accuracy grows like $\sqrt{\kappa(A)} \sim h^{-1}$ under mesh refinement. This is one
motivation for preconditioning.
```






