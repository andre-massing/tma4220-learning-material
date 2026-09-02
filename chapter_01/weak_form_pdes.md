(sec:weak-form-pdes)=
# Weak formulation of elliptic PDEs

## Weak formulation of the Poisson problem
In this chapter, we briefly discuss how the functional
analysis and function space apparatus
can be employed to analyze the well-posedness
of a certain class of PDEs when given in a so-called
"weak" formulation.

Throughout this chapter, $\Omega \subset \mathbb{R}^n$ is a bounded
$C^1$-polyhedron (or, more generally, a bounded Lipschitz domain) with boundary
$\Gamma = \partial \Omega$. This guarantees that the Gauß theorem, the
{prf:ref}`trace theorem<thm:trace-spaces>` and the
{prf:ref}`Poincaré inequality<thm:poincare>` are all at our disposal.

We start by considering the Poisson problem

$$
- \nabla \cdot \nabla u =  -\Delta u  = f \quad \text{in } \Omega
$$ (eq:poisson-wo-bc)

supplemented with some suitable boundary conditions 
which $u$ should satisfy on 
the boundary $\Gamma = \partial \Omega$ of $\Omega$.

### Homogeneous Neumann problem for $-\Delta + \mathrm{Id}$ operator
Let us consider the homogeneous Neumann problem

$$
\left\{
\begin{alignedat}{2}
- \Delta u + u &= f & &\quad \text{in } \Omega \\
         \partial_n u &= 0 & &\quad \text{on } \Gamma
\end{alignedat}
\right.
$$ (eq:neumann-problem-I)

Here, we used the slightly simplified notation $\partial_n u = \mathbf{n} \cdot \nabla u$.
The idea to derive a so-called weak formulation of a PDE is very similar
to the idea behind the introduction of weak derivatives:
We multiply with a suitable test function $v$, integrate over $\Omega$ 
and perform integration by parts to transfer a number of derivatives 
to the test function $v$. 

What kind of test function space we choose is often dictated by two considerations:

1. What kind of smoothness do we require to make the derived formulation work?

2. How do we take into account the boundary conditions?

For the Neumann boundary problem, let's assume for the
moment that $u$, our boundary $\Gamma$ and our
test functions $v$ are smooth enough so that we can use Green's theorem,
e.g., $u \in C^2(\overline{\Omega})$, $\Gamma$ is a $C^1$ boundary, and
$v \in C^{\infty}(\overline{\Omega})$. Then multiplying 
the PDE in {eq}`eq:neumann-problem-I` with $v$ and integrating over $\Omega$
and applying Green's theorem leads to

$$
\begin{aligned}
\int_{\Omega} f v
&= 
- \int_{\Omega}  \nabla \cdot (\nabla u) v 
+\int_{\Omega} uv
\\
&=
- \int_{\Gamma} \underbrace{(\mathbf{n} \cdot \nabla u)}_{=0} v
+ \int_{\Omega} \nabla u \cdot \nabla  v
+\int_{\Omega} uv
\end{aligned}
$$ (eq:weak-form-neumann-derivation)

Note that the Neumann boundary condition $\mathbf{n} \cdot \nabla u = 0$
makes the boundary integrals vanish.
Also observe that the right-hand side of {eq}`eq:weak-form-neumann-derivation`
can be interpreted as taking 
the inner product associated with $H^1(\Omega)$ between $u$ and $v$.
In fact, the expression makes perfect sense even if we only
assume that both $u, v \in H^1(\Omega) =: V$.
With this assumption, we can define the bilinear form 

$$
a(v, w) := 
\int_{\Omega} \nabla v \cdot \nabla  w
+\int_{\Omega} vw
$$ (poisson-bilinear-form-a)

on $V\times V$, and it is straightforward to
show that $a(\cdot, \cdot)$ (being the $H^1$ inner product itself)
satisfies the required assumptions of {prf:ref}`the Lax-Milgram theorem<thm-lax-milgram>`:
Here is some other text.

$$
\begin{aligned}
\text{Boundedness: }  |a(v,w)| &=
\Bigl| \int_{\Omega} \nabla v \cdot \nabla  w
+\int_{\Omega} vw \Bigr| = |(v, w)_{H^1(\Omega)}| \leqslant \|v\|_{H^1(\Omega)} \|w\|_{H^1(\Omega)}
\\
\text{Coercivity: }  a(v,v) &= 
\int_{\Omega} |\nabla v|^2 + \int_{\Omega} |v|^2
= (v, v)_{H^1(\Omega)} = \|v\|_{H^1(\Omega)}^2 
\end{aligned}
$$ (eq-lax-milgram-assump)

so that the assumptions of {prf:ref}`the Lax-Milgram theorem<thm-lax-milgram>`
hold with $C_a = 1$ and $\alpha = 1$.

Next, we define the linear form $l : V \to \mathbb{R}$

$$
l(v) := \int_{\Omega} f v = (f, v)_{L^2(\Omega)}
$$ (poisson-linear-form-l)

If we assume that $f \in L^2(\Omega)$, then thanks to the Cauchy-Schwarz inequality,

$$
|l(v)| 
= |(f, v)_{L^2(\Omega)}| 
\leqslant 
\|f\|_{L^2(\Omega)} \|v\|_{L^2(\Omega)}
\leqslant
\|f\|_{L^2(\Omega)} \|v\|_{H^1(\Omega)},
$$

we can immediately conclude that $l$ is a continuous linear form with $C_l = \|f\|_{L^2(\Omega)}$.
Thus the {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>`
lets us conclude that the problem: find $u \in H^1(\Omega)=: V$ such that $\forall v \in V$

$$
a(u,v) = l(v)
$$

has a unique solution for every $f \in L^2(\Omega)$ with $\|u\|_{H^1(\Omega)} \leqslant \|f\|_{L^2(\Omega)}$.


### Homogeneous Dirichlet problem for $-\Delta + \mathrm{Id}$ operator
Next, we consider 

$$
\left\{
\begin{alignedat}{2}
- \Delta u + u &= f & &\quad \text{in } \Omega \\
         u &= 0 & &\quad \text{on } \Gamma
\end{alignedat}
\right.
$$ (eq:dirichlet-problem-I)

We proceed as for the Neumann problem: we multiply with suitable test functions $v$
and integrate by parts,
but this time, the boundary integral does not vanish since we don't have natural boundary conditions
to incorporate. 
To compensate, we only consider test functions $v \in C^{\infty}_c(\Omega)$
which vanish at the boundary. 
Then again, we obtain

$$
\begin{aligned}
\int_{\Omega} f v
&= 
- \int_{\Omega}  \nabla \cdot (\nabla u) v 
+\int_{\Omega} uv
\\
&=
- \int_{\Gamma} (\mathbf{n} \cdot \nabla u)\underbrace{v}_{=0}
+ \int_{\Omega} \nabla u \cdot \nabla  v
+\int_{\Omega} uv.
\end{aligned}
$$ (eq:weak-form-dirichlet-derivation)

Intuitively speaking, we know what the solution $u$ is going to look like
on the boundary, namely $u=0$, so we don't need test functions which
test for how the equation "behaves" at the boundary.
Also, we now require that our function $u$ comes from a function space
where the boundary condition $u=0$ is already incorporated.
This is exactly what the $H^1_0(\Omega)$ space is made for!
So the weak formulation for {eq}`eq:dirichlet-problem-I` is:

Find $u \in V := H^1_0(\Omega)$ such that

$$
a(u,v) = l(v) \quad \forall v \in V,
$$ (eq:dirichlet-problem-I-weak)

where $a(\cdot, \cdot)$ and $l(\cdot)$
are defined as in {eq}`poisson-bilinear-form-a` and
{eq}`poisson-linear-form-l`, respectively.
As in the case for the homogeneous Neumann problem {eq}`eq:neumann-problem-I`,
we can show that $a$ and $l$ satisfy the assumptions of the
{prf:ref}`Lax-Milgram theorem<thm-lax-milgram>`, and therefore
we can conclude there is a unique solution $u$ to the
weak formulation of the homogeneous Poisson problem which depends
continuously on the data $f$.

```{important}
The only but very important difference between the weak formulation
of the 
*homogeneous Neumann* problem {eq}`eq:neumann-problem-I`
and the 
*homogeneous Dirichlet* problem {eq}`eq:dirichlet-problem-I`
is the Hilbert space on which they are posed.
```

### Homogeneous Dirichlet problem for $-\Delta$ operator
Now, we consider a slightly modified Poisson problem where the
low order term $u$ is left out:

$$
\left\{
\begin{alignedat}{2}
- \Delta u &= f & &\quad \text{in } \Omega \\
         u &= 0 & &\quad \text{on } \Gamma
\end{alignedat}
\right.
$$ (eq:dirichlet-problem-II)

Repeating the steps from the previous section, we arrive at the problem:
Find $u \in V := H^1_0(\Omega)$ such that

$$
a(u,v) = l(v) \quad \forall v \in V,
$$ (eq:dirichlet-problem-II-weak)

with the only distinction that $a(\cdot, \cdot)$ is now given by

$$
a(v, w) = \int_{\Omega} \nabla v \cdot \nabla w.
$$

The boundedness of $a(\cdot, \cdot)$ and $l(\cdot)$ can be shown (almost) exactly as before.
But let's have a look at the coercivity/ellipticity: Setting $u=v$, we obtain 

$$
a(v,v) = \int_{\Omega} |\nabla v|^2
$$

But thanks to the {prf:ref}`Poincaré inequality<thm:poincare>`
and {prf:ref}`cor:poincare` we not only know that $|\nabla\ \cdot |$ 
defines a norm on the closed subspace $H^1_{0}(\Omega)$ but that this norm
is equivalent to the usual $H^1$-norm. Thanks to the proof of {prf:ref}`cor:poincare` we see that

$$
a(v,v) = \int_{\Omega} |\nabla v|^2 \geqslant (1 + C_P^2)^{-1} \|v\|_{1,\Omega}^2.
$$

Thus the {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>` applies with
$C_a = 1$ and $\alpha = (1+C_P^2)^{-1}$, and the resulting stability estimate
reads $\|u\|_{1,\Omega} \leqslant (1+C_P^2)\, \|f\|_{L^2(\Omega)}$.


### Inhomogeneous Dirichlet problem for $-\Delta$ operator
Next, we consider 

$$
\left\{
\begin{alignedat}{2}
- \Delta u &= f & &\quad \text{in } \Omega \\
         u &= g_D & &\quad \text{on } \Gamma
\end{alignedat}
\right.
$$ (eq:dirichlet-problem-inhomog)

Compared to our previous weak formulation for the homogeneous problem,
the main question is now: how can we incorporate the inhomogeneous
Dirichlet b.c. $u = g_D$? First, we realize that
the trial space $H^1_0(\Omega)$ for the solution does not 
make sense anymore. So let's start from $H^1(\Omega)$. Then
we also observe that the data $g_D$ should be in 
$H^{1/2}(\Gamma)$, see {prf:ref}`def:honehalf` to ensure
that we can satisfy the equation $u = g_D$,
and only $u$ satisfying this b.c. should be viable solution candidates
for our weak formulation. Thus we set

$$
H^1_{g_D}(\Omega) := 
\{ v \in H^1(\Omega) \;|\; \gamma(v) = g_D \}.
$$

Since $g_D \in H^{1/2}(\Gamma)$, this set is not empty.
Note that 
$H^1_{g_D}(\Omega)$ is not really a vector space whenever 
$g_D$ is not $0$ everywhere since the addition of 
two functions $u_1, u_2$ with the same non-vanishing boundary data $g_D$ will result in
a function $u$ satisfying $u = 2 g_D$!
In that sense, $H^1_{g_D}(\Omega)$ should rather be considered an **affine** subspace:
For any $u_g$ satisfying $\gamma(u_g) = g_D$, it holds that

$$
H^1_{g_D}(\Omega) = u_g + H^1_0(\Omega)
                  := \{ u_g + v \;|\; v \in H^1_0(\Omega) \}
                   = \gamma^{-1}(g_D).
$$

So the resulting weak formulation is
Find $u \in V := H^1_{g_D}(\Omega)$ such that
for all $v \in \widehat{V} := H^1_0(\Omega)$,

$$
\underbrace{\int_{\Omega} \nabla u \cdot\nabla v}_{=:a(u,v)} = \underbrace{\int_{\Omega} f v}_{=:l(v)}.
$$

Note how in this case the trial function space and test function space
are not identical any more!  How can we prove the well-posedness of
this weak formulation? Lax-Milgram usually requires that the first and
second slots of $a(\cdot, \cdot)$ invoke elements from the same
(vector) space!  The common trick here is to "**lift**" the boundary
condition, i.e. we know that by the definition of $H^{1/2}(\Gamma)$
there must be a $u_g \in H^1(\Omega)$ such that $\gamma(u_g) = g_D$. 
Then we make the ansatz $u = u_0 + u_g$ with $u_0 \in
H^1_0(\Omega)$, leading to the following weak formulation: find $u_0
\in H^1_0(\Omega) =: V$ such that

$$
a(u_0,v) = l(v) - a(u_g, v) =: \widetilde{l}(v) \quad \forall\, v \in V.
$$ (eq:weak-form-dirichlet-inhomog-lifted)

```{exercise} Well-posedness of the lifted Dirichlet formulation
Show that 
{eq}`eq:weak-form-dirichlet-inhomog-lifted` is well-posed and that the solution $u_0$ depends continuously on the data $f$ and $g_D$.
Which of the assumptions of the Lax-Milgram theorem
requires a little extra checking? 
```

### Homogeneous Neumann problem for $-\Delta$ operator
Now we drop the low order term $u$ in the PDE and consider the homogeneous
Neumann problem for the Laplace operator:

$$
\left\{
\begin{alignedat}{2}
- \Delta u  &= f & &\quad \text{in } \Omega \\
         \partial_n u &= 0 & &\quad \text{on } \Gamma
\end{alignedat}
\right.
$$ (eq:neumann-problem-II)

For this section, we assume in addition that $\Omega$ is **connected**, so that
the {prf:ref}`Poincaré--Wirtinger inequality<thm:poincare-wirtinger>` is at our
disposal.

Repeating the derivation {eq}`eq:weak-form-neumann-derivation` --- the boundary
term still vanishes, since $\partial_n u = 0$ is a *natural* boundary condition
--- but now without the low order term, we arrive at

$$
\int_{\Omega} f v = \int_{\Omega} \nabla u \cdot \nabla v,
$$

which suggests the bilinear and linear forms

$$
a(v, w) := \int_{\Omega} \nabla v \cdot \nabla w,
\qquad
l(v) := \int_{\Omega} f v.
$$ (eq:neumann-problem-II-forms)

At first sight we can now proceed exactly as for the homogeneous Neumann problem
{eq}`eq:neumann-problem-I` and simply pose the weak formulation on
$V = H^1(\Omega)$. But this time the resulting problem is **not** well-posed, and
it is instructive to see what goes wrong.

#### The problem is not uniquely solvable on $H^1(\Omega)$
Neither the differential equation nor the boundary condition in
{eq}`eq:neumann-problem-II` sees an additive constant: if $u$ is a solution and
$c \in \mathbb{R}$, then $-\Delta (u + c) = -\Delta u = f$ and
$\partial_n (u+c) = \partial_n u = 0$, so $u + c$ is a solution as well. Since
$\Omega$ is bounded, the constants belong to $H^1(\Omega)$, and the same
degeneracy is visible in the weak formulation:

$$
a(v + c, w) = a(v,w)
\quad \forall\, v, w \in H^1(\Omega), \; c \in \mathbb{R},
$$ (eq:neumann-II-constants)

because $\nabla c = 0$. Consequently the solution of
$a(u,v) = l(v) \; \forall v \in H^1(\Omega)$ can never be unique.

On the level of the {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>` this shows
up as a failure of the coercivity assumption: boundedness still holds with
$C_a = 1$, since

$$
|a(v,w)| =
\Bigl| \int_{\Omega} \nabla v \cdot \nabla w \Bigr|
\leqslant \| \nabla v\|_{L^2(\Omega)} \| \nabla w \|_{L^2(\Omega)}
\leqslant \|v\|_{1,\Omega} \|w\|_{1,\Omega},
$$

but testing the coercivity requirement $a(v,v) \geqslant \alpha \|v\|_{1,\Omega}^2$
with the constant function $v \equiv 1$ gives

$$
a(v,v) = \int_{\Omega} |\nabla v|^2 = 0
\quad \text{while} \quad
\|v\|_{1,\Omega}^2 = |\Omega| > 0,
$$

so that no $\alpha > 0$ can work. Note that this is *not* a deficiency of the
Lax-Milgram theorem: it faithfully reflects the genuine non-uniqueness
{eq}`eq:neumann-II-constants` of the problem. Note also the contrast with the
Neumann problem {eq}`eq:neumann-problem-I`, where the low order term
$\int_{\Omega} uv$ contributed the missing $\|v\|_{L^2(\Omega)}^2$ and thus
rendered $a(\cdot,\cdot)$ coercive on all of $H^1(\Omega)$.

#### Removing the constants: the space $H^1_{\#}(\Omega)$
The remedy is the same as for the Dirichlet problem: we choose a smaller
space. Since $\Omega$ is connected, $a(v,v) = \|\nabla v\|_{L^2(\Omega)}^2 = 0$
implies that $v$ is constant, so the constants are *exactly* the obstruction we
have to eliminate. In other words, a solution is determined only up to its mean
value, and we can single out one representative by prescribing that mean value
--- most conveniently by setting it to $0$. This is precisely what the space
$H^1_{\#}(\Omega)$ from {prf:ref}`def:mean-zero-lp` and
{eq}`eq:mean-zero-sobolev` was introduced for, and we take

$$
V := H^1_{\#}(\Omega) = H^1(\Omega) \cap L^2_{\#}(\Omega)
= \Bigl\{ v \in H^1(\Omega) \mid \overline{v} = 0 \Bigr\}.
$$

Note that the boundary condition $u=0$ is what excludes the constants from
$H^1_0(\Omega)$, while here it is the normalization $\overline{u} = 0$ that does
the job.
By {prf:ref}`lem:mean-zero-decomposition`, every $v \in H^1(\Omega)$ splits
uniquely as $v = \overline{v} + v_{\#}$ with $v_{\#} \in H^1_{\#}(\Omega)$, so
$H^1_{\#}(\Omega)$ contains exactly one representative of each class
$v + \mathbb{R}$ of functions which {eq}`eq:neumann-II-constants` does not
distinguish. Moreover, $H^1_{\#}(\Omega)$ is the kernel of the continuous linear
functional $v \mapsto \int_{\Omega} v$ on $H^1(\Omega)$ and thus a *closed*
subspace, hence itself a Hilbert space with respect to
$(\cdot, \cdot)_{H^1(\Omega)}$ --- exactly what the Lax-Milgram theorem requires.

#### Coercivity via the Poincaré--Wirtinger inequality
On this space the coercivity we lost is restored, and the mechanism is the exact
analogue of the one used for the homogeneous Dirichlet problem
{eq}`eq:dirichlet-problem-II`: there the
{prf:ref}`Poincaré inequality<thm:poincare>` controlled $\|v\|_{L^2(\Omega)}$ by
$\|\nabla v\|_{L^2(\Omega)}$ on $H^1_0(\Omega)$, here the
{prf:ref}`Poincaré--Wirtinger inequality<thm:poincare-wirtinger>` does so on
$H^1_{\#}(\Omega)$. In particular, {prf:ref}`cor:poincare-wirtinger` tells us
that $\|\nabla \cdot\|_{L^2(\Omega)}$ is a norm on $H^1_{\#}(\Omega)$ which is
equivalent to the usual $H^1$-norm, and its proof shows that

$$
a(v,v) = \int_{\Omega} |\nabla v|^2
\geqslant (1 + C_{PW}^2)^{-1} \| v \|_{1,\Omega}^2
\quad \forall\, v \in H^1_{\#}(\Omega).
$$ (eq:neumann-II-coercivity)

Thus $a(\cdot,\cdot)$ is bounded and coercive on $V \times V$ with $C_a = 1$ and
$\alpha = (1+C_{PW}^2)^{-1}$.

#### A compatibility condition on the right-hand side
There is one more difference to the Dirichlet case, and this time it concerns the
**data** rather than the space. Assume for a moment that $u$ is a classical
solution of {eq}`eq:neumann-problem-II`. Integrating the PDE over $\Omega$ and
applying the Gauß theorem together with the boundary condition yields

$$
\int_{\Omega} f
= - \int_{\Omega} \Delta u
= - \int_{\Gamma} \partial_n u
= 0.
$$ (eq:neumann-compatibility)

So {eq}`eq:neumann-problem-II` has *no* solution at all unless the right-hand
side satisfies the **compatibility condition** {eq}`eq:neumann-compatibility`,
that is, unless $f \in L^2_{\#}(\Omega)$. Equivalently, and more in the spirit of
the weak formulation, {eq}`eq:neumann-compatibility` is what makes testing with
$H^1_{\#}(\Omega)$ instead of the full space $H^1(\Omega)$ harmless.

````{prf:lemma} Compatibility condition
:label: lem:neumann-compatibility
Let $f \in L^2(\Omega)$ and $u \in H^1_{\#}(\Omega)$, and let $a(\cdot,\cdot)$
and $l(\cdot)$ be given by {eq}`eq:neumann-problem-II-forms`. Then

$$
a(u,v) = l(v) \quad \forall\, v \in H^1(\Omega)
$$

holds if and only if

$$
a(u,v) = l(v) \quad \forall\, v \in H^1_{\#}(\Omega)
\qquad \text{and} \qquad
\int_{\Omega} f = 0.
$$

````

````{prf:proof}
Let $v \in H^1(\Omega)$ and split it as $v = \overline{v} + v_{\#}$ according to
{prf:ref}`lem:mean-zero-decomposition`. Since the first summand is a constant, we
have $\nabla v = \nabla v_{\#}$ and therefore $a(u,v) = a(u, v_{\#})$, while

$$
l(v) = \int_{\Omega} f v_{\#} + \overline{v} \int_{\Omega} f
= l(v_{\#}) + \overline{v} \int_{\Omega} f.
$$

Consequently,

$$
a(u,v) - l(v)
= \bigl( a(u,v_{\#}) - l(v_{\#}) \bigr) - \overline{v} \int_{\Omega} f
\quad \forall\, v \in H^1(\Omega).
$$ (eq:neumann-II-splitting)

If the left-hand side vanishes for all $v \in H^1(\Omega)$, then choosing
$v \in H^1_{\#}(\Omega)$ (so that $\overline{v} = 0$) gives the first assertion,
and choosing $v \equiv 1$ (so that $a(u,v) = l(v) - \int_{\Omega} f \cdot 1$
reduces to $0 = \int_{\Omega} f$) gives the second. Conversely, if both
assertions hold, then both terms on the right-hand side of
{eq}`eq:neumann-II-splitting` vanish for every $v \in H^1(\Omega)$.
````

#### The weak formulation and its well-posedness
Collecting everything, the weak formulation of {eq}`eq:neumann-problem-II` reads:
given $f \in L^2_{\#}(\Omega)$, find $u \in V := H^1_{\#}(\Omega)$ such that

$$
a(u,v) = l(v) \quad \forall\, v \in V,
$$ (eq:neumann-problem-II-weak)

with $a(\cdot,\cdot)$ and $l(\cdot)$ from {eq}`eq:neumann-problem-II-forms`.
Boundedness of $l$ follows as before from the Cauchy--Schwarz inequality,
$|l(v)| \leqslant \|f\|_{L^2(\Omega)} \|v\|_{1,\Omega}$, so with
{eq}`eq:neumann-II-coercivity` the {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>`
applies with $C_a = 1$ and $\alpha = (1+C_{PW}^2)^{-1}$: the problem
{eq}`eq:neumann-problem-II-weak` has a unique solution, and it satisfies the
stability estimate

$$
\|u\|_{1,\Omega} \leqslant (1+C_{PW}^2) \, \|f\|_{L^2(\Omega)}.
$$

By {prf:ref}`lem:neumann-compatibility`, this $u$ also satisfies
$a(u,v) = l(v)$ for **all** $v \in H^1(\Omega)$, and by
{eq}`eq:neumann-II-constants` the complete set of solutions in $H^1(\Omega)$ is
the affine subspace $u + \mathbb{R}$. Uniqueness therefore holds "up to
constants" only, and $H^1_{\#}(\Omega)$ is nothing but a convenient way of
picking one member of that family.

```{important}
Compared with the homogeneous Dirichlet problem {eq}`eq:dirichlet-problem-II`,
the pure Neumann problem {eq}`eq:neumann-problem-II` differs in **two**
respects, and both stem from the fact that the constants solve the homogeneous
problem:
* the solution space is $H^1_{\#}(\Omega)$ instead of $H^1_0(\Omega)$, and
  coercivity comes from the Poincaré--Wirtinger inequality instead of the
  Poincaré inequality;
* the data must satisfy the compatibility condition
  $\int_{\Omega} f = 0$, i.e. $f \in L^2_{\#}(\Omega)$, otherwise no solution
  exists.
```

```{exercise} Uniqueness up to constants
:label: ex:neumann-uniqueness-up-to-constants
Let $f \in L^2(\Omega)$ and assume $u_1, u_2 \in H^1(\Omega)$ both satisfy
$a(u_i, v) = l(v)$ for all $v \in H^1(\Omega)$. Show that $u_1 - u_2$ is
constant, and deduce that exactly one function in $u_1 + \mathbb{R}$ belongs to
$H^1_{\#}(\Omega)$. Where did you use that $\Omega$ is connected?
```


```{exercise} Inhomogeneous Neumann problems
Discuss Neumann problems when $g_N \neq 0$:
$$
\left\{
\begin{alignedat}{2}
- \Delta u  &= f & &\quad \text{in } \Omega \\
         \partial_n u &= g_N & &\quad \text{on } \Gamma
\end{alignedat}
\right.
$$
```



## Robin problems
Finally, we consider the Robin problem

$$
\left\{
\begin{alignedat}{2}
- \Delta u  &= f & &\quad \text{in } \Omega, \\
         \partial_n u &= \sigma(g_R-u) & &\quad \text{on } \Gamma,
\end{alignedat}
\right.
$$

where $\sigma \in L^{\infty}(\Gamma)$ satisfies
$\sigma \geqslant \sigma_0 > 0$ almost everywhere on $\Gamma$.

```{exercise} Weak formulation of the Robin problem
:label: ex:robin-weak-formulation
Derive a weak formulation for the Robin problem and show that it is well-posed.
```

```{hint}
:class: dropdown
Multiplying with a test function $v \in H^1(\Omega)$ and inserting the Robin
condition produces a bilinear form which carries an additional boundary term
$\int_{\Gamma} \sigma u v$. Coercivity on $H^1(\Omega)$ is the assumption
requiring care here: it does not follow from the
{prf:ref}`Poincaré inequality<thm:poincare>` alone, since the constants do not
belong to $H^1_0(\Omega)$, but from a Poincaré-type inequality involving the
boundary term. This is where $\sigma \geqslant \sigma_0 > 0$ enters.
```

(sec:weak-form-general-pdes)=
## Weak formulation of general second order elliptic PDEs

The PDE {eq}`eq:poisson-wo-bc` is the prototype example
of a 2nd order elliptic operator.
More generally and without any significant complications, 
we can consider a more general PDE of the form

$$
\mathcal{A} u := - \nabla \cdot ( A \nabla u) = f
$$

where $A = (a_{ij}(x))_{i,j=1}^n$ is a pointwise defined matrix.
Note that

$$
\mathcal{A} u = - \nabla \cdot ( A(x) \nabla u(x))
=
-\sum_{i,j=1}^n \partial_{i} (a_{ij}(x) \partial_{j} u(x))
$$ (eq:def-A-operator)

We say that $\mathcal{A}$ is a **second order operator in divergence form**.
Note that the expression {eq}`eq:def-A-operator` 
  * is a generalization of the Laplace operator;
  * *does not make sense in a strong/pointwise sense* if $a_ij$ are not smooth enough, but as we see below, it can be interpreted very easily in a weak sense.

For most of the remaining lectures, we will require $A(x)$
to satisfy the following definition.

```{prf:definition} Ellipticity of $\mathcal{A}$ 
:label: def:ellipticity

The partial differential operator $\mathcal{A}$ given
by {eq}`eq:def-A-operator` with coefficients
$A = (a_{ij})_{i,j=1}^n \in (L^{\infty}(\Omega))^{n\times n}$
is called **(uniformly) elliptic** if there exists a constant $\alpha > 0$ such that
* $ \lambda \cdot  A(x) \lambda \geqslant \alpha |\lambda|^2$

for all $\lambda \in \mathbb{R}^n$ and almost every $x \in \Omega$.
```

```{prf:remark}
Note that $A \in (L^{\infty}(\Omega))^{n\times n}$ also implies that
there exists a $\beta \geqslant 0 $ such that also
* $|A(x) \lambda| \leqslant \beta  |\lambda|$

holds for all $\lambda \in \mathbb{R}^n$ and almost every $x \in \Omega$, and by
ellipticity, we can
conclude that in fact $\beta \geqslant \alpha > 0$.
```
```{exercise} Boundedness of the coefficient matrix
Prove the statements made in the previous remark
```

<!-- ```{admonition} TODO
:class: :danger :dropdown
* Relate $\mathcal{A}$ to classical Poisson problem
* Explain why general $A(x)$ is useful, e.g. anisotropic heat conduction problems
``` -->

We are now prepared to investigate the well-posedness of a number of boundary value
problems where we supplement the partial differential operator 
$\mathcal{A}$ with one of the following boundary conditions

* **Dirichlet boundary conditions** 
   Given a function  $g_D: \Gamma  \to \mathbb{R}$, we require that

  $$
  u = g_D \quad \text{on } \Gamma
  $$

* **Neumann boundary conditions** 
   Given a function  $g_N: \Gamma  \to \mathbb{R}$, we require that

  $$
  \mathbf{n} \cdot A \nabla u = g_N \quad \text{on } \Gamma 
  $$

* **Robin boundary conditions** 
   Given $g_R \in L^2(\Gamma)$ and $\sigma \in L^{\infty}(\Gamma)$ with
   $\sigma \geqslant \sigma_0 > 0$ almost everywhere on $\Gamma$, we require that

  $$
  \mathbf{n} \cdot A \nabla u = \sigma(g_R - u) \quad \text{on } \Gamma 
  $$

These boundary conditions are called *homogeneous* if $g_D$ (respectively $g_N$, $g_R$)
is zero, otherwise we deal with *inhomogeneous* boundary data.
We start by looking at the Poisson problem supplemented with Neumann boundary conditions

```{exercise} Well-posedness of 2nd order elliptic PDEs
1. Derive a weak formulation for the PDE $\mathcal{A} u = f$ with
   for each of the three types of boundary conditions discussed above.
2. Show that the weak formulation is well-posed.
```
