(sec:weak-pdes-motivation)=
# Motivation: towards weak formulations of PDEs

The concept of a *weak formulation* of a PDE and the corresponding
analysis of the well-posedness of this formulation plays a central role
in the design and analysis of finite element methods.

The point of departure for this theory is the insight that a ("regular enough")
function $u$ on an open set $\Omega \subseteq \RR^n$ can be identified by considering
all its weighted integrals $\int_{\Omega} u v$ if the **weight or test function**
$v$ is coming from a rich enough space:

````{prf:lemma} Fundamental lemma of calculus of variations
:label: lem-fundamental-calc-var
For $u \in C(\Omega)$ we have that
```{math}
:label: eq:fundamental-lemma
u = 0
\quad \Leftrightarrow \quad
\int_{\Omega} u v = 0 \quad \forall v \in C^{\infty}_c(\Omega).
```
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

```{math}
u(x) \geqslant \frac{c_0}{2} \quad \forall x \in B_{\varepsilon}(x_0).
\qquad (+)
```

Next, we take $v \in C^{\infty}_c(\Omega)$ to be a so-called *bump function*
satisfying

```{math}
v(x_0) = 1, \quad v(x) \geqslant 0 \;\; \forall x \in \Omega,
\quad \operatorname{supp} v \subseteq B_{\varepsilon}(x_0).
\qquad (*)
```

(The existence of such bump functions is discussed in the exercise below.)
Then

```{math}
\int_{\Omega} u v
\overset{(*)}{=} \int_{B_{\varepsilon}(x_0)} u v
\overset{(+)}{\geqslant} \frac{c_0}{2}
\underbrace{\int_{B_{\varepsilon}(x_0)} v}_{>\, 0 \;\text{(why?)}}
> 0.
```
````

`````{exercise} Construction of bump functions in $\RR^n$
:label: exer-bump-functions
In this exercise, we construct the bump functions used in the proof of the
fundamental lemma.

1. Consider the function $f : \RR \to \RR$ defined by
   ```{math}
   f(t) :=
   \begin{cases}
   e^{-1/t} & \text{for } t > 0, \\
   0 & \text{for } t \leqslant 0.
   \end{cases}
   ```
   Show that $f \in C^{\infty}(\RR)$.

2. Define $\varphi : \RR^n \to \RR$ by $\varphi(x) := f(1 - \|x\|^2)$.
   Show that $\varphi \in C^{\infty}_c(\RR^n)$ with
   $\varphi(x) > 0$ for $\|x\| < 1$ and
   $\operatorname{supp} \varphi = \overline{B_1(0)}$.

3. For a given ball $B_{\varepsilon}(x_0) \subseteq \Omega$, define
   ```{math}
   v(x) := \frac{1}{\varphi(0)}\,
   \varphi\Bigl( \frac{2 (x - x_0)}{\varepsilon} \Bigr).
   ```
   Show that $v \in C^{\infty}_c(\Omega)$ satisfies the properties $(*)$
   required in the proof of the fundamental lemma, namely
   $v(x_0) = 1$, $v(x) \geqslant 0$ for all $x \in \Omega$, and
   $\operatorname{supp} v \subseteq \overline{B_{\varepsilon/2}(x_0)}
   \subseteq B_{\varepsilon}(x_0)$.
`````

````{hint}
:class: hint dropdown
For 1., the only critical point is $t = 0$. Show first by induction that for
every $k \in \NN$ there is a polynomial $p_k$ such that
```{math}
f^{(k)}(t) = p_k(1/t) \, e^{-1/t} \quad \text{for } t > 0,
```
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
```{math}
\int_{\Omega} u_1 v = \int_{\Omega} u_2 v \quad \forall v \in C^{\infty}_c(\Omega),
```
then $u_1 = u_2$.
````

```{prf:proof}
Exercise.
```

With the previous lemma at our disposal, we consider the
"Hello world" example of a 2nd order, elliptic PDE,
namely the **Poisson problem**: Find a function $u : \Omega \to \RR$ s.t.

```{math}
:label: eq:poisson-strong-form
\begin{alignat}{2}
-\nabla \cdot (\nabla u(x)) = - \Delta u(x) &= f(x) &&\quad \forall x \in \Omega,
\\
u(x) &= 0 &&\quad \forall x \in \partial\Omega,
\end{alignat}
```

for a given function $f : \Omega \to \RR$. This is called the **strong formulation**
of the Poisson problem, as the equations are supposed to hold *pointwise*
$\forall x \in \Omega$ (or $\in \partial\Omega$).

At the moment we are purposely vague about the precise assumptions
on $u$, $f$ and $\Omega$.

The idea is now to use the previous Corollary to reformulate the
strong formulation into a **weighted integral** version: If
$u \in C^2(\Omega) \cap C(\overline{\Omega})$ and $f \in C(\Omega)$, then

```{math}
:label: eq:poisson-weighted
-\Delta u(x) = f(x) \;\; \forall x \in \Omega
\quad \Leftrightarrow \quad
- \int_{\Omega} \Delta u \, v = \int_{\Omega} f v
\quad \forall v \in C^{\infty}_c(\Omega).
```

```{note}
We usually don't write down integration measures such as $\dx$, $\dt$,
as they can be deduced from the domain of integration using the
canonical measures.
```

Now we need to recall some analysis tools to reformulate {eq}`eq:poisson-weighted`.
In particular, we need to recall the Gauß theorem and its relatives.
Gauß' theorem holds only on certain domain types, and throughout this course
we require $\Omega$ to be at least a $C^1$-polyhedron:

````{prf:definition} $C^1$-polyhedron
:label: def-c1-polyhedron

* Let $\Omega \subseteq \RR^n$ be open. A point $\overline{x} \in \partial\Omega$
  is called a **regular boundary point** of $\Omega$ if there is a neighborhood
  $\mcN(\overline{x}) \subseteq \RR^n$ of $\overline{x}$ and a $C^1$ function
  $q : \mcN(\overline{x}) \to \RR$ s.t.
  1. $\nabla q(x) \neq 0 \quad \forall x \in \mcN(\overline{x})$,
  2. $\Omega \cap \mcN(\overline{x}) = \{ x \in \mcN(\overline{x}) \st q(x) < 0 \}$.

* **Regular boundary** $\partial_r \Omega := \{ x \in \partial\Omega \st x \text{ is a regular boundary point} \}$

* **Singular boundary** $\partial_s \Omega := \partial\Omega \setminus \partial_r \Omega$

* A set $M$ is called a **$d$-dimensional (Hausdorff) null set** if $\forall \varepsilon > 0$
  there exists a countable collection of balls $\{B^{\varepsilon}_k\}_{k=1}^{\infty}$
  with radius $r_k$ s.t.
  1. $M \subseteq \bigcup_{k=1}^{\infty} B^{\varepsilon}_k$,
  2. $\sum_{k=1}^{\infty} r_k^d < \varepsilon$.

* An open domain $\Omega \subseteq \RR^n$ is called a **$C^1$-polyhedron** if
  $\partial_s \Omega$ is an $(n-1)$-dimensional null set.
````

```{figure} figures/c1-polyhedron-sketch.svg
:label: fig-c1-polyhedron
:alt: A domain Omega whose boundary has both smooth parts and corners, with a regular boundary point x-bar on a smooth part, its circular neighborhood N(x-bar), the shaded region where q is negative, and a corner marked as a singular boundary point.
:width: 55%

A $C^1$-polyhedron $\Omega$: a regular boundary point $\overline{x} \in \partial_r\Omega$
on a smooth part of the boundary, with its neighborhood
$\mcN(\overline{x})$. The shaded region is
$\Omega \cap \mcN(\overline{x}) = \{ x \in \mcN(\overline{x}) \st q(x) < 0 \}$,
and the highlighted piece of $\partial\Omega$ is the zero level set of $q$.
The corners of $\partial\Omega$ are singular boundary points and form the
singular boundary $\partial_s\Omega$ — a finite set of points, and thus an
$(n-1)$-dimensional null set for $n = 2$.
```

````{prf:theorem} Gauß theorem
:label: thm-gauss
Let $\Omega \subseteq \RR^n$ be a bounded $C^1$-polyhedron and
$\mathbf{F} : \Omega \to \RR^n$ be a vector field s.t.
1. $\mathbf{F} \in C^1(\Omega, \RR^n) \cap C(\overline{\Omega}, \RR^n)$,
2. $\nabla \cdot \mathbf{F}$ is integrable over $\Omega$,
3. $\mathbf{F} \cdot \bfn$ is integrable over $\partial\Omega$.

Then
```{math}
:label: eq:gauss-theorem
\int_{\Omega} \nabla \cdot \mathbf{F} = \int_{\partial\Omega} \mathbf{F} \cdot \bfn.
```
````

The Gauß theorem has some immediate and important consequences.

````{prf:corollary} Green's formula and partial integration
:label: cor-greens-formulas
Let $\Omega$ be a bounded $C^1$-polyhedron. Then we have

1. for $u \in C^1(\overline{\Omega})$, $\mathbf{F} \in C^1(\overline{\Omega}, \RR^n)$:
   ```{math}
   :label: eq:green-a
   \int_{\Omega} \nabla u \cdot \mathbf{F}
   = - \int_{\Omega} u \, \nabla \cdot \mathbf{F}
   + \int_{\partial\Omega} u \, \mathbf{F} \cdot \bfn
   ```
2. for $u \in C^2(\overline{\Omega})$, $v \in C^1(\overline{\Omega})$:
   ```{math}
   :label: eq:green-b
   \int_{\Omega} \nabla u \cdot \nabla v
   = - \int_{\Omega} \Delta u \, v
   + \int_{\partial\Omega} \partial_n u \, v,
   \qquad \text{where } \partial_n u := \nabla u \cdot \bfn,
   ```
3. for $u, v \in C^2(\overline{\Omega})$:
   ```{math}
   :label: eq:green-c
   \int_{\Omega} (\Delta u \, v - u \, \Delta v)
   = \int_{\partial\Omega} (\partial_n u \, v - u \, \partial_n v)
   ```
4. for $u, v \in C^1(\overline{\Omega})$:
   ```{math}
   :label: eq:green-d
   \int_{\Omega} \partial_{x_i} u \, v
   = - \int_{\Omega} u \, \partial_{x_i} v
   + \int_{\partial\Omega} u v \, n_i.
   ```
````

```{exercise} Proving Green's formulas
Deduce {eq}`eq:green-a`--{eq}`eq:green-d` from the Gauß theorem.
```

Now we can return to the integral equation in {eq}`eq:poisson-weighted` and
rewrite it via {eq}`eq:green-b` to see that
$-\Delta u(x) = f(x) \;\; \forall x \in \Omega \Leftrightarrow$

```{math}
:label: eq:poisson-weak-derivation
\int_{\Omega} f v
= - \int_{\Omega} \Delta u \, v
= \int_{\Omega} \nabla u \cdot \nabla v
- \underbrace{\int_{\partial\Omega} \partial_n u \, v}_{= \, 0
\text{ since } v|_{\partial\Omega} = 0}
= \int_{\Omega} \nabla u \cdot \nabla v
\quad \forall v \in C^{\infty}_c(\Omega).
```

So instead of the strong form, we could now be tempted to
solve the following problem: Find $u \in C^2(\Omega) \cap C_0(\overline{\Omega})$ s.t.

```{math}
:label: eq:poisson-almost-weak
\underbrace{\int_{\Omega} \nabla u \cdot \nabla v}_{:= a(u,v)} = \underbrace{\int_{\Omega} f v}_{:= l(v)}
\quad \forall v \in C^{\infty}_c(\Omega).
```

Here
```{math}
C_0(\overline{\Omega}) = \{ u \in C(\overline{\Omega}) \st u|_{\partial\Omega} = 0 \}
```
is used to enforce the boundary condition $u = 0$ on $\partial\Omega$,
which is **not** directly imposed in our integral rewrite of the PDE part,
as $C^{\infty}_c(\Omega)$ test functions "don't see" the boundary.

* Also note that for each $u$ that solves {eq}`eq:poisson-almost-weak`,
  $u + \text{constant}$ also solves {eq}`eq:poisson-almost-weak`!

* Note that {eq}`eq:poisson-almost-weak` uses different function spaces for
  $u$ and $v$, which differ in their differentiability requirements, despite
  the fact that the expression $\int_{\Omega} \nabla u \cdot \nabla v$
  requires precisely $1$ derivative.

So to make sense of $\int_{\Omega} \nabla u \cdot \nabla v$, we need two things:
1. $\nabla u$, $\nabla v$ must be properly defined,
2. the integral of the product $\nabla u \cdot \nabla v$ is well-defined.

Point 1. says something about **differentiability**, point 2. about **integrability**
properties of $u$ and $v$.

So the main idea of the weak formulation is to
   * cast the strong PDE into a so-called weak formulation of the form {eq}`eq:poisson-almost-weak`
   which uses the weighted integral form 
   * to use suitable function spaces with relaxed differentiability and certain integrability properties 
   so that weak formulation still makes sense, 
   * to use the tools of functional analysis to prove well-posedness of the weak formulation.

The well-posedness of the weak formulation will then typically only
guarantee solutions in a weaker sense than the strong formulation.
To obtain "classical" solutions with the desired differentiability properties, one typically then needs to invoke **regularity theory**,
which we only will discuss very briefly in this course.

In the next two sections, we will introduce the necessary tools from [functional analysis](sec:functional-analysis)
 and [Sobolev spaces](ssec:sobolev-spaces) to make the above ideas precise. 