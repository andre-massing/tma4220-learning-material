(sec:weak-pdes-motivation)=
# Weak formulations of partial differential equations

## Towards a weak formulation of the Poisson problem
The concept of a *weak formulation* of a PDE and the corresponding
analysis of the well-posedness of this formulation play a central role
in the design and analysis of finite element methods.

The point of departure for this theory is the insight that a ("regular enough")
function $u$ on an open set $\Omega \subseteq \mathbb{R}^n$ can be identified by considering
all its weighted integrals $\int_{\Omega} u v$ if the **weight or test function**
$v$ comes from a rich enough space:

````{prf:lemma} Fundamental lemma of calculus of variations
:label: lem-fundamental-calc-var
For $u \in C(\Omega)$ we have that

$$
u = 0
\quad \Leftrightarrow \quad
\int_{\Omega} u v = 0 \quad \forall v \in C^{\infty}_c(\Omega).
$$ (eq:fundamental-lemma)

````

````{prf:proof}
"$\implies$" is trivial.

"$\Leftarrow$": Proof by contradiction, i.e. we show that
$u \neq 0 \Rightarrow \exists\, v \in C^{\infty}_c(\Omega)$
such that (s.t.) $\int_{\Omega} u v \neq 0$.

Assume without loss of generality (w.l.o.g.) that there is a $x_0 \in \Omega$
s.t. $u(x_0) = c_0 > 0$. Since $u$ is continuous and $\Omega$ is open,
there must be an open ball $B_{\varepsilon}(x_0) \subseteq \Omega$ of radius
$\varepsilon > 0$ and center $x_0$ s.t.

$$
u(x) \geqslant \frac{c_0}{2} \quad \forall x \in B_{\varepsilon}(x_0).
\qquad (+)
$$

Next, we take $v \in C^{\infty}_c(\Omega)$ to be a so-called *bump function*
satisfying

$$
v(x_0) = 1, \quad v(x) \geqslant 0 \;\; \forall x \in \Omega,
\quad \operatorname{supp} v \subseteq B_{\varepsilon}(x_0).
\qquad (*)
$$

(The existence of such bump functions is discussed in the exercise below.)
Then

$$
\int_{\Omega} u v
\overset{(*)}{=} \int_{B_{\varepsilon}(x_0)} u v
\overset{(+)}{\geqslant} \frac{c_0}{2}
\underbrace{\int_{B_{\varepsilon}(x_0)} v}_{>\, 0 \;\text{(why?)}}
> 0.
$$

````

`````{exercise} Construction of bump functions in $\mathbb{R}^n$
:label: exer-bump-functions
In this exercise, we construct the bump functions used in the proof of the
fundamental lemma.

1. Consider the function $f : \mathbb{R} \to \mathbb{R}$ defined by

   $$
   f(t) :=
   \begin{cases}
   e^{-1/t} & \text{for } t > 0, \\
   0 & \text{for } t \leqslant 0.
   \end{cases}
   $$

   Show that $f \in C^{\infty}(\mathbb{R})$.

2. Define $\varphi : \mathbb{R}^n \to \mathbb{R}$ by $\varphi(x) := f(1 - \|x\|^2)$.
   Show that $\varphi \in C^{\infty}_c(\mathbb{R}^n)$ with
   $\varphi(x) > 0$ for $\|x\| < 1$ and
   $\operatorname{supp} \varphi = \overline{B_1(0)}$.

3. For a given ball $B_{\varepsilon}(x_0) \subseteq \Omega$, define

   $$
   v(x) := \frac{1}{\varphi(0)}\,
   \varphi\Bigl( \frac{2 (x - x_0)}{\varepsilon} \Bigr).
   $$

   Show that $v \in C^{\infty}_c(\Omega)$ satisfies the properties $(*)$
   required in the proof of the fundamental lemma, namely
   $v(x_0) = 1$, $v(x) \geqslant 0$ for all $x \in \Omega$, and
   $\operatorname{supp} v \subseteq \overline{B_{\varepsilon/2}(x_0)}
   \subseteq B_{\varepsilon}(x_0)$.
`````

````{hint}
:class: hint dropdown
For 1., the only critical point is $t = 0$. Show first by induction that for
every $k \in \mathbb{N}$ there is a polynomial $p_k$ such that

$$
f^{(k)}(t) = p_k(1/t) \, e^{-1/t} \quad \text{for } t > 0,
$$

and then use that exponential decay beats polynomial growth, i.e.
$\lim_{s \to \infty} p(s) e^{-s} = 0$ for every polynomial $p$, to conclude
that $f^{(k)}(0) = 0$ exists for all $k$.

For 2., note that $x \mapsto 1 - \|x\|^2$ is a polynomial and recall that
compositions of $C^{\infty}$ functions are $C^{\infty}$.

For 3., every property transfers directly from $\varphi$ under translation
and scaling.
````

````{prf:corollary}
:label: cor-weighted-integrals
If $u_1, u_2 \in C(\Omega)$ satisfy

$$
\int_{\Omega} u_1 v = \int_{\Omega} u_2 v \quad \forall v \in C^{\infty}_c(\Omega),
$$

then $u_1 = u_2$.
````

```{prf:proof}
Exercise.
```

With the previous lemma at our disposal, we consider the
"Hello world" example of a 2nd-order elliptic PDE,
namely the **Poisson problem**: Find a function $u : \Omega \to \mathbb{R}$ s.t.

$$
\begin{alignedat}{2}
-\nabla \cdot (\nabla u(x)) = - \Delta u(x) &= f(x) &&\quad \forall x \in \Omega,
\\
u(x) &= 0 &&\quad \forall x \in \partial\Omega,
\end{alignedat}
$$ (eq:poisson-strong-form)

for a given function $f : \Omega \to \mathbb{R}$. This is called the **strong formulation**
of the Poisson problem, as the equations are supposed to hold *pointwise*
$\forall x \in \Omega$ (or $\in \partial\Omega$).

At the moment we are purposely vague about the precise assumptions
on $u$, $f$ and $\Omega$.

The idea is now to use the previous Corollary to reformulate the
strong formulation into a **weighted integral** version: If
$u \in C^2(\Omega) \cap C(\overline{\Omega})$ and $f \in C(\Omega)$, then

$$
-\Delta u(x) = f(x) \;\; \forall x \in \Omega
\quad \Leftrightarrow \quad
- \int_{\Omega} \Delta u \, v = \int_{\Omega} f v
\quad \forall v \in C^{\infty}_c(\Omega).
$$ (eq:poisson-weighted)

```{note}
We usually don't write down integration measures such as $\,\mathrm{d}x$, $\,\mathrm{d}t$,
as they can be deduced from the domain of integration using the
canonical measures.
```

## The divergence theorem and Green's formulas

Now we need to recall some analysis tools to reformulate {eq}`eq:poisson-weighted`.
In particular, we need to recall the Gauß theorem and its relatives.
The Gauß theorem holds only on certain domain types, and throughout this course
we require $\Omega$ to be at least a $C^1$-polyhedron:

````{prf:definition} $C^1$-polyhedron
:label: def-c1-polyhedron

* Let $\Omega \subseteq \mathbb{R}^n$ be open. A point $\overline{x} \in \partial\Omega$
  is called a **regular boundary point** of $\Omega$ if there is a neighborhood
  $\mathcal{N}(\overline{x}) \subseteq \mathbb{R}^n$ of $\overline{x}$ and a $C^1$ function
  $q : \mathcal{N}(\overline{x}) \to \mathbb{R}$ s.t.
  1. $\nabla q(x) \neq 0 \quad \forall x \in \mathcal{N}(\overline{x})$,
  2. $\Omega \cap \mathcal{N}(\overline{x}) = \{ x \in \mathcal{N}(\overline{x}) \;|\; q(x) < 0 \}$.

* **Regular boundary** $\partial_r \Omega := \{ x \in \partial\Omega \;|\; x \text{ is a regular boundary point} \}$

* **Singular boundary** $\partial_s \Omega := \partial\Omega \setminus \partial_r \Omega$

* A set $M$ is called a **$d$-dimensional (Hausdorff) null set** if $\forall \varepsilon > 0$
  there exists a countable collection of balls $\{B^{\varepsilon}_k\}_{k=1}^{\infty}$
  with radius $r_k$ s.t.
  1. $M \subseteq \bigcup_{k=1}^{\infty} B^{\varepsilon}_k$,
  2. $\sum_{k=1}^{\infty} r_k^d < \varepsilon$.

* An open domain $\Omega \subseteq \mathbb{R}^n$ is called a **$C^1$-polyhedron** if
  $\partial_s \Omega$ is an $(n-1)$-dimensional null set.
````

For a regular boundary point $\overline{x} \in \partial_r \Omega$ with a
defining function $q$ as in {prf:ref}`def-c1-polyhedron`, we set

$$
\mathbf{n}(\overline{x}) := \dfrac{\nabla q(\overline{x})}{|\nabla q(\overline{x})|}.
$$ (eq:outer-normal)

We call $\mathbf{n}(\overline{x})$ the **outer** (or outward) unit normal of $\Omega$
at $\overline{x}$. The name is justified as follows: $q$ is negative on
$\Omega \cap \mathcal{N}(\overline{x})$ and vanishes on the boundary piece
$\partial\Omega \cap \mathcal{N}(\overline{x})$, so $\nabla q(\overline{x})$ points in
the direction in which $q$ increases, that is, away from $\Omega$.

Of course, a regular boundary point may admit many different defining
functions $q$, and we should convince ourselves that they all lead to the same
normal vector.

````{prf:lemma} The outer unit normal is well defined
:label: lem-outer-normal
Let $\Omega \subseteq \mathbb{R}^n$ be open and let
$\overline{x} \in \partial_r \Omega$ be a regular boundary point. Then the
vector $\mathbf{n}(\overline{x})$ defined in {eq}`eq:outer-normal` satisfies
$|\mathbf{n}(\overline{x})| = 1$, and it does not depend on the choice of the
defining function: if $q_1$ and $q_2$ both satisfy the conditions of
{prf:ref}`def-c1-polyhedron` at $\overline{x}$, then they define one and the
same vector $\mathbf{n}(\overline{x})$.
````

```{exercise} Well-definedness of the outer unit normal
:label: exer-outer-normal
Prove {prf:ref}`lem-outer-normal`.
```

```{hint}
:class: dropdown
Near $\overline{x}$, both defining functions vanish exactly on the same piece
of $\partial\Omega$, so their zero level sets coincide there. Recall that the
gradient of a $C^1$ function with non-vanishing gradient is orthogonal to its
own zero level set. Hence $\nabla q_1(\overline{x})$ and
$\nabla q_2(\overline{x})$ are both orthogonal to the same $(n-1)$-dimensional
tangent space, and are therefore parallel. It remains to exclude opposite
signs, which follows since $q_1$ and $q_2$ are negative on the same side of
the boundary, namely on $\Omega$.
```

Since the outer unit normal is only available on the regular part
$\partial_r \Omega$, integration over the boundary of a $C^1$-polyhedron is
always understood as integration over the regular part of the boundary, that
is, $\int_{\partial\Omega} := \int_{\partial_r \Omega}$. This is no
restriction: by {prf:ref}`def-c1-polyhedron`, the singular boundary
$\partial_s \Omega$ is an $(n-1)$-dimensional null set and therefore does not
contribute to such integrals.

```{figure} figures/c1-polyhedron-sketch.svg
:label: fig-c1-polyhedron
:alt: A domain Omega whose boundary has both smooth parts and corners, with a regular boundary point x-bar on a smooth part, its circular neighborhood N(x-bar), the shaded region where q is negative, and a corner marked as a singular boundary point.
:width: 55%

A $C^1$-polyhedron $\Omega$: a regular boundary point $\overline{x} \in \partial_r\Omega$
on a smooth part of the boundary, with its neighborhood
$\mathcal{N}(\overline{x})$. The shaded region is
$\Omega \cap \mathcal{N}(\overline{x}) = \{ x \in \mathcal{N}(\overline{x}) \;|\; q(x) < 0 \}$,
and the highlighted piece of $\partial\Omega$ is the zero level set of $q$.
The corners of $\partial\Omega$ are singular boundary points and form the
singular boundary $\partial_s\Omega$ — a finite set of points, and thus an
$(n-1)$-dimensional null set for $n = 2$.
```

````{prf:theorem} Gauß theorem
:label: thm-gauss
Let $\Omega \subseteq \mathbb{R}^n$ be a bounded $C^1$-polyhedron and
$\mathbf{F} : \Omega \to \mathbb{R}^n$ be a vector field s.t.
1. $\mathbf{F} \in C^1(\Omega, \mathbb{R}^n) \cap C(\overline{\Omega}, \mathbb{R}^n)$,
2. $\nabla \cdot \mathbf{F}$ is integrable over $\Omega$,
3. $\mathbf{F} \cdot \mathbf{n}$ is integrable over $\partial\Omega$.

Then

$$
\int_{\Omega} \nabla \cdot \mathbf{F} = \int_{\partial\Omega} \mathbf{F} \cdot \mathbf{n}.
$$ (eq:gauss-theorem)

````

The Gauß theorem has some immediate and important consequences.

````{prf:corollary} Green's formula and partial integration
:label: cor-greens-formulas
Let $\Omega$ be a bounded $C^1$-polyhedron. Then we have

1. for $u \in C^1(\overline{\Omega})$, $\mathbf{F} \in C^1(\overline{\Omega}, \mathbb{R}^n)$:

   $$
   \int_{\Omega} \nabla u \cdot \mathbf{F}
   = - \int_{\Omega} u \, \nabla \cdot \mathbf{F}
   + \int_{\partial\Omega} u \, \mathbf{F} \cdot \mathbf{n}
   $$ (eq:green-a)

2. for $u \in C^2(\overline{\Omega})$, $v \in C^1(\overline{\Omega})$:

   $$
   \int_{\Omega} \nabla u \cdot \nabla v
   = - \int_{\Omega} \Delta u \, v
   + \int_{\partial\Omega} \partial_n u \, v,
   \qquad \text{where } \partial_n u := \nabla u \cdot \mathbf{n},
   $$ (eq:green-b)

3. for $u, v \in C^2(\overline{\Omega})$:

   $$
   \int_{\Omega} (\Delta u \, v - u \, \Delta v)
   = \int_{\partial\Omega} (\partial_n u \, v - u \, \partial_n v)
   $$ (eq:green-c)

4. for $u, v \in C^1(\overline{\Omega})$:

   $$
   \int_{\Omega} \partial_{x_i} u \, v
   = - \int_{\Omega} u \, \partial_{x_i} v
   + \int_{\partial\Omega} u v \, n_i.
   $$ (eq:green-d)

````

```{exercise} Proving Green's formulas
Deduce {eq}`eq:green-a`--{eq}`eq:green-d` from the Gauß theorem.
```

Now we can return to the integral equation in {eq}`eq:poisson-weighted` and
rewrite it via {eq}`eq:green-b`. Note that {eq}`eq:green-b` is stated for
$u \in C^2(\overline{\Omega})$, while we only assumed
$u \in C^2(\Omega) \cap C(\overline{\Omega})$. This is not a problem here:
since $v \in C^{\infty}_c(\Omega)$ vanishes outside a compact subset of
$\Omega$, we may apply {eq}`eq:green-b` on a subdomain
$\Omega' \Subset \Omega$ containing $\operatorname{supp} v$, so that no
regularity of $u$ up to $\partial\Omega$ is required. We thus see that
$-\Delta u(x) = f(x) \;\; \forall x \in \Omega \Leftrightarrow$

$$
\int_{\Omega} f v
= - \int_{\Omega} \Delta u \, v
= \int_{\Omega} \nabla u \cdot \nabla v
- \underbrace{\int_{\partial\Omega} \partial_n u \, v}_{= \, 0
\text{ since } v \text{ vanishes near } \partial\Omega}
= \int_{\Omega} \nabla u \cdot \nabla v
\quad \forall v \in C^{\infty}_c(\Omega).
$$ (eq:poisson-weak-derivation)

So instead of the strong form, we could now be tempted to
solve the following problem: Find $u \in C^2(\Omega) \cap C_0(\overline{\Omega})$ s.t.

$$
\underbrace{\int_{\Omega} \nabla u \cdot \nabla v}_{:= a(u,v)} = \underbrace{\int_{\Omega} f v}_{:= l(v)}
\quad \forall v \in C^{\infty}_c(\Omega).
$$ (eq:poisson-almost-weak)

Here

$$
C_0(\overline{\Omega}) = \{ u \in C(\overline{\Omega}) \;|\; u|_{\partial\Omega} = 0 \}
$$

is used to enforce the boundary condition $u = 0$ on $\partial\Omega$,
which is **not** directly imposed in our integral rewrite of the PDE part,
as $C^{\infty}_c(\Omega)$ test functions "don't see" the boundary.

* Note that the integral identity in {eq}`eq:poisson-almost-weak` alone does
  not determine $u$: if $u$ satisfies it, so does $u + c$ for every constant
  $c$, since $\nabla (u + c) = \nabla u$. It is precisely the membership
  $u \in C_0(\overline{\Omega})$ which removes this non-uniqueness.

* Note that {eq}`eq:poisson-almost-weak` uses different function spaces for
  $u$ and $v$, which differ in their differentiability requirements, although the expression $\int_{\Omega} \nabla u \cdot \nabla v$
  requires precisely $1$ derivative.

So to make sense of $\int_{\Omega} \nabla u \cdot \nabla v$, we need two things:
1. $\nabla u$, $\nabla v$ must be properly defined,
2. the integral of the product $\nabla u \cdot \nabla v$ is well-defined.

Point 1. says something about **differentiability**, point 2. about **integrability**
properties of $u$ and $v$.

So the main idea of the weak formulation is to
   * to cast the strong PDE into a so-called weak formulation of the form {eq}`eq:poisson-almost-weak`
   which uses the weighted integral form 
   * to use suitable function spaces with relaxed differentiability and certain integrability properties 
   so that the weak formulation still makes sense, 
   * to use the tools of functional analysis to prove well-posedness of the weak formulation.

The well-posedness of the weak formulation will then typically only
guarantee solutions in a weaker sense than the strong formulation.
To obtain "classical" solutions with the desired differentiability properties, one typically then needs to invoke **regularity theory**,
which we will only discuss very briefly in this course.

In the next two sections, we will introduce the necessary tools from [functional analysis](#sec:functional-analysis)
 and [Sobolev spaces](#ssec:sobolev-spaces) to make the above ideas precise. 
