(sec:function-spaces)=
# A brief review on function spaces

In this chapter we collect some results on various function spaces we will
use throughout the book. One essential property of many function
spaces we will consider is that they are *complete*, i.e. they are either a
Banach or a Hilbert space, see Section {ref}`sec:functional-analysis`.


```{note}
This section barely scratches the surface of the topic; we will
only summarize (and not even prove) the most essential results we need
later on in this course.

Also, this chapter will be a work in progress during the entire course,
as we will add relevant results here whenever we need them elsewhere.
```

## Measure and integration theory, Lebesgue spaces
Lebesgue integration theory provides a powerful generalization
of the Riemann integral which makes sure that the set of so-called
Lebesgue-integrable functions turns into a Banach space when
endowed with a suitable norm. Nowadays, in most standard textbooks,
Lebesgue integration theory is presented as part of the
curriculum on *Measure and Integration theory*, see Chapters 9-10 in 
{cite}`Browder2012` for a quick introduction. 
To this end, 

````{prf:definition} Lebesgue spaces
Let $\Omega \subset \RR^n$ be an open domain. 

Then
the Lebesgue spaces $L^p(\Omega)$ are defined by
```{math}
:label: eq:lebesque-spaces
L^p(\Omega) = \{ f : \Omega \to \RR \text{ is measurable and } 
\|f\|_{L^p(\Omega)} 
< \infty
 \}.
```
Here, the $L^p$-norm $\|\cdot\|_{L^p(\Omega)}$ is defined by
```{math}
:label: eq:lebesque-norm
\|f\|_{L^p(\Omega)} 
:=  
\begin{cases}
\Bigl(
    \int_{\Omega} |f(x)|^p
\Bigr)^{1/p}  
&\quad 1 \leqslant p < \infty
\\
\mathrm{ess}\sup_{\Omega} |f|
&\quad p=\infty
\end{cases}
```

Sometimes we write $\|f\|_{p,\Omega}$ instead of $\|f\|_{L^p(\Omega)}$.
A function $f \in L^p(\Omega)$ is often called $L^p$-integrable. 

We also introduce the space of *locally* $L^p$-integrable functions on $\Omega$; that is,
functions that are $L^p$ integrable on every compact subset $K \Subset \Omega$,
```{math}
:label: eq:loc-lebesque-spaces
L^p_{\mathrm{loc}}(\Omega) = \{ f: \Omega \to \RR | f \in L^p(K)  \; \forall K \Subset \Omega \}.
```
````

```{prf:remark}
As usual, we identify functions which agree almost everywhere; strictly
speaking, $L^p(\Omega)$ denotes the resulting space of equivalence classes.
This is what turns $\|\cdot\|_{L^p(\Omega)}$ into a genuine norm rather than
merely a seminorm: a function vanishing outside a set of measure zero has
$\|f\|_{L^p(\Omega)} = 0$ without being the zero function.
```

````{prf:lemma} Hölder's inequality
:label: lem:hoelder
Let $1 \leqslant p, q \leqslant \infty$ with $\frac{1}{p} + \frac{1}{q} = 1$
(with the convention $\frac{1}{\infty} := 0$).
If $f \in L^p(\Omega)$ and $g \in L^q(\Omega)$, then $f g \in L^1(\Omega)$ and
```{math}
:label: eq:hoelder
\| f g \|_{L^1(\Omega)} = \int_{\Omega} |f g|
\leqslant \|f\|_{L^p(\Omega)} \|g\|_{L^q(\Omega)}.
```
````

```{prf:remark}
For an open domain with finite measure $|\Omega| = \int_{\Omega} 1 < \infty$, 
Hölder's inequality implies that 
$L^p(\Omega) \subset L^q(\Omega)$ for $1 \leqslant q \leqslant p \leqslant \infty$.
```

```{exercise}
Proof this.
```


Among the Lebesgue spaces, $L^2(\Omega)$ plays a special role, as it can be
endowed with an inner product.

````{prf:definition} Inner product on $L^2(\Omega)$
:label: def:l2-inner-product
For $f, g \in L^2(\Omega)$, we define
```{math}
:label: eq:l2-inner-product
(f, g)_{L^2(\Omega)} := \int_{\Omega} f g.
```
Sometimes we write $(f,g)_{\Omega}$ instead of $(f, g)_{L^2(\Omega)}$.
````

Note that the inner product is *well-defined*: applying
{prf:ref}`Hölder's inequality<lem:hoelder>` with $p = q = 2$ shows that
$f g \in L^1(\Omega)$, so the integral in {eq}`eq:l2-inner-product` is finite. It is easy
to check that $(\cdot, \cdot)_{L^2(\Omega)}$ is bilinear, symmetric and positive
definite, so it is indeed an inner product on $L^2(\Omega)$, and it
induces exactly the $L^2$-norm from {eq}`eq:lebesque-norm`,
```{math}
(f, f)_{L^2(\Omega)}^{1/2} = \Bigl( \int_{\Omega} |f|^2 \Bigr)^{1/2}
= \|f\|_{L^2(\Omega)}.
```
Moreover, rewriting {prf:ref}`Hölder's inequality<lem:hoelder>` for $p = q = 2$
in terms of the inner product and the $L^2$-norm yields
```{math}
:label: eq:l2-cauchy-schwarz
| (f, g)_{L^2(\Omega)} | \leqslant \int_{\Omega} |f g|
\leqslant \|f\|_{L^2(\Omega)} \|g\|_{L^2(\Omega)},
```
that is, the *Cauchy-Schwarz inequality* on $L^2(\Omega)$ arises as the
special case $p = q = 2$ of Hölder's inequality.

````{prf:lemma} Determining uniqueness through testing
:label: lem:uniqueness-by-testing
Let $u_1, u_2 \in L^1_{\mathrm{loc}}(\Omega)$ and assume that
```{math}
\int_{\Omega} (u_1 - u_2) \phi = 0  \quad \forall \phi \in C^{\infty}_c(\Omega).
```
Then $u_1 = u_2$ almost everywhere in $\Omega$; that is, up to a set of measure $0$.
````

```{prf:proof}
Note that this lemma generalizes the
{prf:ref}`fundamental lemma of calculus of variations<lem-fundamental-calc-var>`
from continuous functions to $L^1_{\mathrm{loc}}(\Omega)$; the almost
everywhere conclusion requires a mollification argument rather than the bump
function construction used there. For a proof, we refer to
{cite}`Brezis2011` (Corollary 4.24).
```

```{prf:remark}
In this setting, $\phi$ is typically called a *test function*. 
When determining 
whether two functions are equal,
the previous lemma roughly states that you can do this
by comparing their "actions" on suitable test functions $\phi$
instead of comparing their values at (almost) every point.

Here, the "action" is simply the resulting number computed from multiplying the functions in question with the test function $\phi$
and integrating over $\Omega$. 
```

### Lebesgue spaces with vanishing mean value

For domains of finite measure it is often convenient to single out those
$L^p$ functions whose average vanishes. Note that $|\Omega| < \infty$
implies $L^p(\Omega) \subset L^1(\Omega)$, so that the integral in the
following definition is well-defined for every $v \in L^p(\Omega)$.

````{prf:definition} Mean value and Lebesgue spaces with vanishing mean
:label: def:mean-zero-lp
Let $\Omega \subset \RR^n$ be an open domain of finite measure,
$|\Omega| = \int_{\Omega} 1 < \infty$.
For $v \in L^1(\Omega)$, we call
```{math}
:label: eq:mean-value
\overline{v} := \dfrac{1}{|\Omega|} \int_{\Omega} v
```
the *mean value* (or *average*) of $v$ over $\Omega$.
For $1 \leqslant p \leqslant \infty$, we then define
```{math}
:label: eq:mean-zero-lp
L^p_{\#}(\Omega) :=
\Bigl\{ v \in L^p(\Omega) \mid \int_{\Omega} v = 0 \Bigr\}
= \{ v \in L^p(\Omega) \mid \overline{v} = 0 \},
```
the subspace of $L^p$ functions with *vanishing mean value*.
````

Since $|\Omega| > 0$, the two descriptions in {eq}`eq:mean-zero-lp` indeed
agree. Every $L^p$ function can be split into its mean value and a
fluctuation around it.

````{prf:lemma} Decomposition into mean value and fluctuation
:label: lem:mean-zero-decomposition
Let $\Omega$ be of finite measure and $1 \leqslant p \leqslant \infty$.
Then every $v \in L^p(\Omega)$ can be written uniquely as
```{math}
:label: eq:mean-zero-decomposition
v = \overline{v} + v_{\#},
\qquad
v_{\#} := v - \overline{v} \in L^p_{\#}(\Omega),
```
where the first summand is understood as a constant function on $\Omega$.
In other words, $L^p(\Omega) = \RR \oplus L^p_{\#}(\Omega)$, and
$L^p_{\#}(\Omega)$ is a closed subspace of $L^p(\Omega)$.
````

````{prf:proof}
Constants belong to $L^p(\Omega)$ because $|\Omega| < \infty$, so
$v_{\#} = v - \overline{v} \in L^p(\Omega)$, and
```{math}
\int_{\Omega} v_{\#}
= \int_{\Omega} v - \overline{v} \, |\Omega|
= \int_{\Omega} v - \int_{\Omega} v = 0,
```
so $v_{\#} \in L^p_{\#}(\Omega)$. For uniqueness, assume
$v = c_1 + w_1 = c_2 + w_2$ with constants $c_i$ and
$w_i \in L^p_{\#}(\Omega)$. Integrating over $\Omega$ and using
$\int_{\Omega} w_i = 0$ gives $c_1 |\Omega| = c_2 |\Omega|$, hence
$c_1 = c_2$ and therefore $w_1 = w_2$.

Closedness follows since $v \mapsto \int_{\Omega} v$ is a continuous
linear functional on $L^p(\Omega)$ --- by
{prf:ref}`Hölder's inequality<lem:hoelder>` applied with $g \equiv 1$ ---
and $L^p_{\#}(\Omega)$ is precisely its kernel.
````

```{prf:remark}
For $p = 2$, the space $L^2_{\#}(\Omega)$ is a closed subspace of the
Hilbert space $L^2(\Omega)$ and hence itself a Hilbert space with respect
to the inner product from {prf:ref}`def:l2-inner-product`. Spaces of this
type appear naturally whenever a boundary value problem determines its
solution only up to an additive constant --- the pure Neumann problem being
the prototypical example --- since normalizing the mean value to zero
singles out exactly one representative.
```

```{exercise} Mean value as an orthogonal projection
Show that the map $P : L^2(\Omega) \to L^2_{\#}(\Omega)$,
$P v := v - \overline{v}$, is linear, satisfies $P^2 = P$, and is the
orthogonal projection onto $L^2_{\#}(\Omega)$ with respect to
$(\cdot, \cdot)_{L^2(\Omega)}$. Deduce that
$\|v - \overline{v}\|_{L^2(\Omega)} \leqslant \|v\|_{L^2(\Omega)}$.
```


(ssec:sobolev-spaces)=
## Sobolev spaces

### Weak derivatives 

Let us start with a motivating example. Let $u \in C^k(\Omega)$ and $\phi \in C^{\infty}_c(\Omega)$.
Using Green's theorem and taking into account that $\phi = 0$ on an open neighborhood of the boundary of $\Omega$, 
we see that
```{math}
:label: eq:weak-deriv-first
\int_{\Omega} \partial_{x_i} u \phi = 
- \int_{\Omega} u \partial_{x_i} \phi,
```
and iterating this formula, we observe that for any multiindex
$\alpha \in \NN_0^n$, where $\NN_0 := \NN \cup \{0\}$, with
$|\alpha| \leqslant k$,
```{math}
:label: eq:weak-deriv-alpha
\int_{\Omega} \partial^{\alpha} u \phi = 
(-1)^{|\alpha|} \int_{\Omega} u \partial^{\alpha} \phi,
```

where $|\alpha| = \alpha_1 + \cdots + \alpha_n$.
Note that the integral expression on the right-hand side of {eq}`eq:weak-deriv-alpha` makes perfect
sense even for $u\in L^1_{\mathrm{loc}}$ and not only $u\in C^k(\Omega)$.
This leads to a possibility to generalize or weaken the classical definition of derivatives. 

````{prf:definition} Weak derivative
Let $\alpha \in \NN_0^n$ be a multiindex and $u, u_{\alpha} \in L^1_{\mathrm{loc}}(\Omega)$.
We say that $u_{\alpha}$ is *$\alpha$-th weak derivative* of $u$ if
```{math}
\int_{\Omega}  u_{\alpha} \phi = 
(-1)^{|\alpha|} \int_{\Omega} u \partial^{\alpha} \phi
```
holds for all $\phi \in C^{\infty}_c(\Omega)$.
````

```{prf:lemma} Uniqueness of weak derivatives
If $u \in L^1_{\mathrm{loc}}(\Omega)$ possesses an $\alpha$-th weak derivative, it is uniquely defined
in $L^1_{\mathrm{loc}}(\Omega)$.
```

````{prf:proof}
For
two weak derivatives $u_{\alpha}$ and $\tilde{u}_{\alpha}$ we have that
```{math}
\begin{aligned}
\int_{\Omega} u_{\alpha} \phi &= (-1)^{|\alpha|} \int_{\Omega} u \partial^{\alpha} \phi
\\
\int_{\Omega} \tilde{u}_{\alpha} \phi &= (-1)^{|\alpha|} \int_{\Omega} u \partial^{\alpha} \phi
\end{aligned}
```
and by subtracting the second from the first equation, we obtain that 
```{math}
\int_{\Omega}
(u_{\alpha} - \tilde{u}_{\alpha} ) \phi =  0 \quad \forall \phi \in C_c^{\infty}(\Omega),
```
and thus $u_{\alpha} = \tilde{u}_{\alpha}$ almost everywhere
by {prf:ref}`lem:uniqueness-by-testing`.
````

````{exercise} Relation between the modulus function and the sign function
Let $\Omega = (-1,1)$ and set
```{math}
\begin{aligned}
u(x) &= |x| \\
\mathrm{sgn}(x) &= \begin{cases} -1 &\quad x \in (-1,0) \\ 
                      1  &\quad x \in [0, 1)
        \end{cases}
\end{aligned}
```
By simply using the definition of the weak derivative, show that
$\mathrm{sgn}(x)$ is the weak derivative of $u$.
````

````{prf:definition} Sobolev spaces
* $W^{k,p}(\Omega) := 
  \{
  u \in L^p(\Omega) |\, \partial^{\alpha}u \text{ exists as a weak derivative and belongs to } L^p(\Omega) 
  \, \forall \alpha \text{ with } |\alpha| \leqslant k
  \}
  $
* For $p=2$, we usually write
  ```{math}
  H^k(\Omega) := W^{k,2}(\Omega)
  ```
  Note that the $\| \cdot \|_{H^k(\Omega)}$ is induced by the inner product
  ```{math}
  (v,w)_{H^k(\Omega)} := 
    \sum_{|\alpha| \leqslant k} (\partial^{\alpha} v, \partial^{\alpha} w)_{L^2(\Omega)}
  ```
* For $u \in W^{k,p}(\Omega)$, we set 
  ```{math}
  \| u \|_{W^{k,p}({\Omega})} := \|u\|_{k,p,\Omega} 
  :=
  \begin{cases}
  \Bigl( 
    \sum_{|\alpha| \leqslant k} \| \partial^{\alpha} u \|_{L^p(\Omega)}^p
  \Bigr)^{1/p} 
  & 1\leqslant p <  \infty,
  \\
    \sum_{|\alpha| \leqslant k} \| \partial^{\alpha} u \|_{L^{\infty}(\Omega)}
  & p =  \infty.
  \end{cases}
  ```
* We set
  ```{math}
  W_0^{k,p}(\Omega) := \overline{C_c^{\infty}(\Omega)}^{\|\cdot\|_{k,p,\Omega}},
  ```
  that is, the topological closure of $C_c^{\infty}(\Omega)$ in $W^{k,p}(\Omega)$.
* Finally, we introduce the common notation for the dual space of $H^1_0(\Omega)$, 
  ```{math}
  H^{-1}(\Omega) := (H^1_0(\Omega))'. 
  ``` 
````

```{prf:remark}
$W_0^{k,p}(\Omega)$ can be understood as the closed subspace 
consisting of those functions $\phi$ in $W^{k,p}(\Omega)$ which are limits
of sequences $\{\phi_n\}_{n=1}^\infty \subset C_c^{\infty}(\Omega)$.
```
For first order Sobolev spaces, we abbreviate the norm of the gradient by
```{math}
:label: eq:gradient-norm
\| \nabla u \|_{L^p(\Omega)}^p := \sum_{i=1}^n \| \partial_{x_i} u \|_{L^p(\Omega)}^p,
```
so that $\| u \|_{W^{1,p}(\Omega)}^p = \| u \|_{L^p(\Omega)}^p + \| \nabla u \|_{L^p(\Omega)}^p$.
For $p = 2$, this agrees with $\bigl( \int_{\Omega} |\nabla u|^2 \bigr)^{1/2}$,
where $|\cdot|$ denotes the Euclidean norm on $\RR^n$.

Later we will need the following important result known as Poincaré inequality.

````{prf:theorem} Poincaré inequality
:label: thm:poincare
Let $\Omega$ be an open and bounded subset of $\RR^n$ and let
$1 \leqslant p < \infty$. Then there is a constant $C_P = C_P(p,n,\Omega)$ such that
```{math}
:label: eq:poincare
\|u \|_{L^p(\Omega)} \leqslant C_P \|\nabla u \|_{L^p(\Omega)}.
```
for any $u \in W^{1,p}_0(\Omega)$.
````
```{prf:proof}
For a proof we refer to {cite}`Evans2010` (p. 279).
```

````{prf:corollary}
:label: cor:poincare
On $W^{1,p}_0(\Omega)$, the norm $\| \cdot\|_{W^{1,p}(\Omega)}$ is equivalent to
```{math}
\| u \|_{W^{1,p}_0(\Omega)} := \| \nabla u \|_{L^{p}(\Omega)}.
```
More precisely,
```{math}
:label: eq:poincare-norm-equivalence
\| \nabla u \|_{L^p(\Omega)}
\leqslant \| u \|_{W^{1,p}(\Omega)}
\leqslant (1+C_P^p)^{1/p} \| \nabla u \|_{L^p(\Omega)}
\quad \forall u \in W^{1,p}_0(\Omega).
```
````

````{prf:proof}
A simple application of the Poincaré inequality yields
```{math}
\|\nabla u\|_{L^p(\Omega)}^p
\leqslant
\| u \|_{L^p(\Omega)}^p +
\|\nabla u\|_{L^p(\Omega)}^p
= \| u \|_{W^{1,p}(\Omega)}^p
\leqslant
(1+C_P^p) \|\nabla u\|_{L^p(\Omega)}^p,
```
and taking $p$-th roots gives {eq}`eq:poincare-norm-equivalence`.
````

<!-- ### Approximation results -->

<!-- ### Poincaré inequalties  -->

### Trace operators
Next, we very briefly discuss whether and how functions of certain Sobolev spaces defined
on the domain $\Omega$ can be restricted to the boundary $\partial \Omega$. This plays
an important role in the well-posedness of boundary value problems, as we need to determine 
the correct spaces for the boundary data in, e.g., a Dirichlet or Neumann boundary-value problem
when the data is **non-homogeneous**.

For the remaining part of this Chapter, we assume that $\Omega$ is a bounded domain with 
a "well-behaved" boundary, that is, it is either a $C^1$-polyhedron, a Lipschitz domain, or --- if this doesn't tell  you much ---
simply a $C^{\infty}$ domain.
Integrals over $\Gamma = \partial \Omega$ are always understood with respect
to the $(n-1)$-dimensional surface measure on $\Gamma$, and $L^2(\Gamma)$
denotes the corresponding Lebesgue space.

````{prf:theorem} Traces of $H^1(\Omega)$ spaces
:label: thm:trace-spaces
For a bounded domain $\Omega$ with Lipschitz (or $C^{\infty}$) boundary
$ \Gamma = \partial \Omega$, there
exists a bounded operator $\gamma : H^1(\Omega) \to L^2(\Gamma)$ (the so-called *trace operator*) such that 
$\gamma(u) = u|_{\Gamma}$ for all $u \in H^1(\Omega) \cap C(\overline{\Omega})$.

Moreover, under the same assumptions on $\Omega$, one can show that
```{math}
H^1_0(\Omega) = \mathrm{ker} \gamma = \{v \in H^1(\Omega) \mid \gamma(v) = 0 \}.
```
````

It turns out that the trace operator $\gamma$ **is not onto $L^2(\Gamma)$**. Thus, when we later want
to find certain weak formulations and solutions $u \in H^1(\Omega)$ which also need to satisfy certain
inhomogeneous boundary conditions such as $u = u_D$ on $\Gamma$, we need to be careful about
the choice of function space from which we take the boundary data $u_D$.
That motivates the following

````{prf:definition} $H^{1/2}(\Gamma)$
:label: def:honehalf
We set
```{math}
H^{1/2}(\Gamma) := \{ v \in L^2(\Gamma) \mid v = \gamma(w) \text{ for some } w \in H^1(\Omega) \}
= \mathrm{ran}\, \gamma
```
and define a corresponding norm by
```{math}
  \|v \|_{H^{1/2}(\Gamma)} := \|v\|_{1/2, \Gamma} := \inf \{ \|w\|_{1,\Omega} \mid w \in H^1(\Omega), \; \gamma(w) = v\}.
```
Consequently,
```{math}
\|\gamma(w)\|_{1/2,\Gamma} \leqslant \| w\|_{1, \Omega} \quad \forall w \in H^1(\Omega),
```
that is, $\gamma : H^1(\Omega) \to H^{1/2}(\Gamma)$ is bounded with norm at most $1$.
````
