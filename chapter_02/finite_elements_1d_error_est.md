(sec:fem-1d-error-estimates)=
# Error estimates for the finite element method in 1D

So far, we have constructed finite element spaces and discussed how the
resulting linear systems are assembled. We now turn to the question of how
accurate the discrete solution actually is. By Céa's lemma, see
{prf:ref}`ceas-lemma`, the finite element solution is quasi-optimal, so it
suffices to estimate the error of one conveniently chosen approximation of $u$
in $V_h$. The natural candidate is the interpolant of $u$, and the main part of
this section is therefore devoted to interpolation error estimates.

Throughout this section, we write $\| w \|_{\omega} := \| w \|_{L^2(\omega)}$
for the $L^2$ norm over a subinterval $\omega \subseteq \Omega$.

## The interpolation operator

Let $\Omega = (a,b)$ and let
$a = x_0 < x_1 < \ldots < x_M = b$ be a uniform subdivision of $\Omega$ such that

$$
x_l = x_0 + l h, \quad l = 0, \ldots, M,
\qquad \text{with } h = \frac{b-a}{M} .
$$

The associated mesh is $\mathcal{T}_h = \{T_l\}_{l=0}^{M-1}$ with
$T_l = [x_l, x_{l+1}]$.

Recall from [the previous section](#ssec:finite-element-1d-method) that for
$\mathbb{P}^{\mathrm{c}}_k(\mathcal{T}_h)$ we constructed a set of basis functions
$\{\phi_i\}_{i=0}^{N+1}$ by constructing Lagrange interpolation polynomials on
each subinterval $T_l$ and by patching these polynomials together in such a way
that the resulting functions remain continuous on $(a,b)$. Using the properties
of the Lagrange basis functions, the very same idea defines an interpolation
operator.

```{prf:definition} Nodal interpolation operator
:label: def:nodal-interpolation-1d

The **interpolation operator**
$\mathcal{I}_h : C([a,b]) \to \mathbb{P}^{\mathrm{c}}_k(\mathcal{T}_h)$
is defined element by element. On each element $T_l$ we set

$$
(\mathcal{I}_h v)|_{T_l} := \mathcal{I}_{T_l} v ,
$$

where $\mathcal{I}_{T_l} v$ is the usual Lagrange interpolation polynomial
computed from $v$,

$$
\mathcal{I}_{T_l} v(x) = \sum_{\alpha = 0}^{k} v(\xi^l_\alpha) \, \ell^{(T_l)}_\alpha(x) .
$$ (eq:local-interpolant-1d)

Here $x_l = \xi^l_0 < \xi^l_1 < \ldots < \xi^l_k = x_{l+1}$ are the chosen
interpolation points on $T_l$, often uniformly distributed, and
$\ell^{(T_l)}_\alpha$ is the $\alpha$-th Lagrange cardinal function on $T_l$,
satisfying $\ell^{(T_l)}_\alpha(\xi^l_\beta) = \delta_{\alpha\beta}$.
```

It can be checked that the resulting global interpolation operator
$\mathcal{I}_h v$ is globally continuous, so that indeed
$\mathcal{I}_h v \in \mathbb{P}^{\mathrm{c}}_k(\mathcal{T}_h)$.

```{exercise} Continuity of the interpolant
:label: exer-interpolant-continuity

Prove the claim above, that is, show that $\mathcal{I}_h v$ is continuous
on $[a,b]$ for every $v \in C([a,b])$. What role does the distribution of the
interpolation points $\{\xi^l_\alpha\}_{\alpha=0}^{k}$ on each element play?
```

For the interpolation estimate below, we need one fact which is particular to
one space dimension.

```{prf:remark} A Sobolev embedding which holds only in 1D
:label: rem:sobolev-embedding-1d

For an interval $\Omega = (a,b) \subset \mathbb{R}$ we have

$$
H^{k+1}(\Omega) \subset C^{k}(\overline{\Omega}) .
$$

In particular, $v \in H^2(\Omega)$ implies $v \in C^1(\overline{\Omega})$, so
that the point values $v(x_l)$ used in the definition of
$\mathcal{I}_h v$ are well defined.
Remember that this embedding holds only in **one** dimension: in $\mathbb{R}^n$,
an embedding into $C^k(\overline{\Omega})$ requires $H^{s}(\Omega)$ with
$s > k + n/2$.
```

## The interpolation error estimate

```{prf:theorem} Interpolation estimate for $\mathbb{P}^{\mathrm{c}}_1(\mathcal{T}_h)$
:label: thm:interpolation-estimate-p1

Let $v \in H^2(\Omega)$ and let $\mathcal{I}_h v \in \mathbb{P}^{\mathrm{c}}_1(\mathcal{T}_h)$
be its interpolating function. Then

$$
\| (v - \mathcal{I}_h v)' \|_{\Omega} \leqslant C_1 \, h \, \| v'' \|_{\Omega} ,
$$ (eq:interpolation-estimate-h1)

$$
\| v - \mathcal{I}_h v \|_{\Omega} \leqslant C_0 \, h^2 \, \| v'' \|_{\Omega} .
$$ (eq:interpolation-estimate-l2)

Both statements can be condensed into the single estimate

$$
| v - \mathcal{I}_h v |_{H^r(\Omega)}
:= \Big\| \frac{\mathrm{d}^r}{\mathrm{d}x^r} (v - \mathcal{I}_h v) \Big\|_{\Omega}
\leqslant C_r \, h^{2-r} \, \| v'' \|_{\Omega} ,
\qquad r = 0, 1 .
$$ (eq:interpolation-estimate-hr)
```

Before we turn to the proof, we record an immediate but important corollary.

```{prf:corollary} Energy norm error estimate for linear elements
:label: cor:energy-error-p1

Let $u_h$ be the continuous, piecewise linear finite element solution to either

$$
\begin{aligned}
-u'' &= f \quad \text{in } I = (a,b), &
u(a) &= g_D(a), \quad u(b) = g_D(b),
\end{aligned}
$$

or

$$
\begin{aligned}
-u'' + u &= f \quad \text{in } I = (a,b), &
-u'(a) &= g_N(a), \quad u'(b) = g_N(b) .
\end{aligned}
$$

Assume furthermore that $u \in H^2(I)$. Then

$$
\| (u - u_h)' \|_{I} \leqslant C \, h \, \| u'' \|_{I} .
$$ (eq:energy-error-p1)
```

```{exercise} Proof of the energy norm error estimate
:label: exer-energy-error-p1

Prove {prf:ref}`cor:energy-error-p1`. Combine Céa's lemma {prf:ref}`ceas-lemma` with
{prf:ref}`thm:interpolation-estimate-p1`, using $v_h = \mathcal{I}_h u$ as the
comparison function. Convince yourself first that $\mathcal{I}_h u$ is an
admissible choice in each of the two cases.
```

Now we turn to the proof of the interpolation estimate.

```{prf:proof}
Since $v \in H^2(\Omega)$, we know from {prf:ref}`rem:sobolev-embedding-1d` that
$v \in C^1(\overline{\Omega})$. By the interpolation property, the interpolant
reproduces the nodal values of $v$,

$$
v(x_l) = \mathcal{I}_h v(x_l), \qquad l = 0, \ldots, M ,
$$

so the error function $e := v - \mathcal{I}_h v$ satisfies
$e(x_l) = 0$ for $l = 0, \ldots, M$.

Since $e$ vanishes at both endpoints of each element, Rolle's theorem provides
for every $l = 0, \ldots, M-1$ a point $\zeta_l \in (x_l, x_{l+1})$ such that

$$
e'(\zeta_l) = 0 .
$$

Using the fundamental theorem of calculus, we see that for $x \in (x_l, x_{l+1})$

$$
e'(x) = \underbrace{e'(\zeta_l)}_{= 0} + \int_{\zeta_l}^{x} e''(y) \,\mathrm{d}y
= \int_{\zeta_l}^{x} v''(y) \,\mathrm{d}y ,
$$

where we used that $(\mathcal{I}_h v)'' = 0$ on $T_l$, since
$(\mathcal{I}_h v)|_{T_l}$ is a linear polynomial. Applying the Cauchy-Schwarz
inequality to the product $1 \cdot v''(y)$, we obtain on each element the bound

$$
|e'(x)|
\leqslant \int_{x_l}^{x_{l+1}} 1 \cdot |v''(y)| \,\mathrm{d}y
\leqslant \Big( \int_{x_l}^{x_{l+1}} 1^2 \,\mathrm{d}y \Big)^{\frac{1}{2}}
          \Big( \int_{x_l}^{x_{l+1}} |v''(y)|^2 \,\mathrm{d}y \Big)^{\frac{1}{2}}
\leqslant h^{\frac{1}{2}} \| v'' \|_{T_l} .
$$ (eq:interpolation-proof-bound)

Squaring and integrating this inequality over $T_l$, we arrive at

$$
\int_{x_l}^{x_{l+1}} |e'(x)|^2 \,\mathrm{d}x
\leqslant \int_{x_l}^{x_{l+1}} h \, \| v'' \|^2_{T_l} \,\mathrm{d}x
= h^2 \, \| v'' \|^2_{T_l} ,
$$

and summing over all elements $T_l$ yields

$$
\| e' \|^2_{\Omega}
\leqslant \sum_{l=0}^{M-1} h^2 \, \| v'' \|^2_{T_l}
= h^2 \, \| v'' \|^2_{\Omega} ,
$$

which is precisely {eq}`eq:interpolation-estimate-h1` with $C_1 = 1$.

To derive {eq}`eq:interpolation-estimate-l2`, we combine the fundamental theorem
of calculus, the fact that $e(x_l) = 0$, and the bound
{eq}`eq:interpolation-proof-bound` to infer that

$$
|e(x)|
= \Big| \underbrace{e(x_l)}_{= 0} + \int_{x_l}^{x} e'(y) \,\mathrm{d}y \Big|
\leqslant \int_{x_l}^{x_{l+1}} |e'(y)| \,\mathrm{d}y
\leqslant h^{\frac{3}{2}} \| v'' \|_{T_l}
\qquad \forall x \in (x_l, x_{l+1}) .
$$

Squaring and integrating once more leads to

$$
\int_{x_l}^{x_{l+1}} |e(x)|^2 \,\mathrm{d}x
\leqslant \int_{x_l}^{x_{l+1}} h^{3} \, \| v'' \|^2_{T_l} \,\mathrm{d}x
= h^4 \, \| v'' \|^2_{T_l} ,
$$

and, as before, summing over all elements shows that

$$
\| e \|^2_{\Omega} \leqslant h^4 \, \| v'' \|^2_{\Omega} ,
$$

which proves {eq}`eq:interpolation-estimate-l2` with $C_0 = 1$.
```

```{prf:remark} Non-uniform meshes
:label: rem:interpolation-nonuniform

The proof used the uniformity of the mesh only when replacing the local mesh
size $h_{T_l}$ by $h$. Carrying $h_{T_l}$ along instead gives the local estimates

$$
\| e' \|_{T_l} \leqslant h_{T_l} \| v'' \|_{T_l},
\qquad
\| e \|_{T_l} \leqslant h_{T_l}^2 \| v'' \|_{T_l} ,
$$

so that {eq}`eq:interpolation-estimate-h1` and
{eq}`eq:interpolation-estimate-l2` hold on an arbitrary mesh with
$h = \max_{T \in \mathcal{T}_h} h_T$.
```
