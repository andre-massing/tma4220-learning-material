(sec:function-spaces)=
# A brief review on function spaces

In this chapter we collect some results on various function space we will
use throughout the book. One essential property of many function
space we will consider is that they are *complete*, i.e. they are either 
Banach or Hilbert space, see Section {ref}`sec:functional-analysis`.


```{note}
This section barely scratches at the surface of the topic, we will
only summarize (and not even prove) the most essential results we need
later one in this course.

Also, this chapter will be a work in progress during the entire course,
as we will add relevant results here whenever we need them elsewhere.
```

## Measure and integration theory, Lebesque spaces
Lebesque integration theory provides a powerful generalization
of the Riemann integral which makes sure that the set of so-called
Lebesque-integrable functions turns into a Banach space when
endowed with a suitable norm. Nowdays, in most standard text books,
Lebesque integration theory is presented as part of the
curriculum on *Measure and Integration theory*, see Chapter 9-10 in 
{cite}`Browder2012` for a quick introduction. 
To this end, 

````{prf:definition} Lebesque spaces
Let $\Omega \subset \RR^n$ be a open domain. 

Then
the Lebesque spaces $L^p(\Omega)$ are defined by
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
    \int_{\Omega} |f(x)|^p \dx
\Bigr)^{1/p}  
&\quad 1 \leqslant p < \infty
\\
\mathrm{ess}\sup_{\Omega} |u|
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

````{prf:lemma} Hölder's inequality
:label: lem:hoelder
Let $1 \leqslant p, q \leqslant \infty$ with $\frac{1}{p} + \frac{1}{q} = 1$
(with the convention $\frac{1}{\infty} := 0$).
If $f \in L^p(\Omega)$ and $g \in L^q(\Omega)$, then $f g \in L^1(\Omega)$ and
```{math}
:label: eq:hoelder
\| f g \|_{L^1(\Omega)} = \int_{\Omega} |f g| \dx
\leqslant \|f\|_{L^p(\Omega)} \|g\|_{L^q(\Omega)}.
```
````

Among the Lebesgue spaces, $L^2(\Omega)$ plays a special role, as it can be
endowed with an inner product.

````{prf:definition} Inner product on $L^2(\Omega)$
:label: def:l2-inner-product
For $f, g \in L^2(\Omega)$, we define
```{math}
:label: eq:l2-inner-product
(f, g)_{L^2(\Omega)} := \int_{\Omega} f g \dx.
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
(f, f)_{L^2(\Omega)}^{1/2} = \Bigl( \int_{\Omega} |f|^2 \dx \Bigr)^{1/2}
= \|f\|_{L^2(\Omega)}.
```
Moreover, rewriting {prf:ref}`Hölder's inequality<lem:hoelder>` for $p = q = 2$
in terms of the inner product and the $L^2$-norm yields
```{math}
:label: eq:l2-cauchy-schwarz
| (f, g)_{L^2(\Omega)} | \leqslant \int_{\Omega} |f g| \dx
\leqslant \|f\|_{L^2(\Omega)} \|g\|_{L^2(\Omega)},
```
that is, the *Cauchy-Schwarz inequality* on $L^2(\Omega)$ arises as the
special case $p = q = 2$ of Hölder's inequality.

````{prf:lemma} Determining uniqueness through testing
:label: lem:uniqueness-by-testing
Let $u_1, u_2 \in L^1_{\mathrm{loc}}(\Omega)$ and assume that
```{math}
\int_{\Omega} (u_1 - u_2) \phi \dx = 0  \quad \forall \phi \in C^{\infty}_c(\Omega).
```
Then $u_1 = u_2$ almost every in $\Omega$; that is, up to set of measure $0$.
````

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

(ssec:sobolev-spaces)=
## Sobolev spaces

### Weak derivatives 

Let start with a motivating example. Let $u \in C^k(\Omega)$ and $\phi \in C^{\infty}_c(\Omega)$.
Using Green's theorem and taking into account that $\phi = 0$ on a open neighborhood of the boundary of $\Omega$, 
we see that
```{math}
:label: eq:weak-deriv-first
\int_{\Omega} \partial_{x_i} u \phi \dx = 
- \int_{\Omega} u \partial_{x_i} \phi \dx,
```
and iterating this formula, we observe that for any multiindex $\alpha
\in \NN^n$ with $\alpha \leqslant k$, it holds that
```{math}
:label: eq:weak-deriv-alpha
\int_{\Omega} \partial^{\alpha} u \phi \dx = 
(-1)^{|\alpha|} \int_{\Omega} u \partial^{\alpha} \phi \dx,
```

where $|\alpha| = \alpha_1 + \cdots \alpha_n$.
Note that the integral expression on the right-hand side of {eq}`eq:weak-deriv-alpha` makes perfectly
sense even for $u\in L^1_{\mathrm{loc}}$ and not only $u\in C^k(\Omega)$.
This leads to a possibility to generalize or weakened the classical definition of derivatives. 

````{prf:definition} Weak derivative
Let $\alpha \in \NN^n$ be a multiindex and $u, u_{\alpha} \in L^1_{\mathrm{loc}}(\Omega)$.
We say that $u_{\alpha}$ is *$\alpha$-th weak derivative* of $u$ if
```{math}
\int_{\Omega}  u_{\alpha} \phi \dx = 
(-1)^{|\alpha|} \int_{\Omega} u \partial^{\alpha} \phi \dx
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
\begin{align}
\int_{\Omega} u_{\alpha} \phi &= (-1)^{|\alpha|} \int_{\Omega} u \partial^{\alpha} \phi dx
\\
\int_{\Omega} \tilde{u}_{\alpha} \phi &= (-1)^{|\alpha|} \int_{\Omega} u \partial^{\alpha} \phi dx
\end{align}
and by substracting the second from the first inequality, we obtain that 
```{math}
\int_{\Omega}
(u_{\alpha} - \tilde{u}_{\alpha} ) \phi dx =  0 \quad \forall \phi \in C_c^{\infty}(\Omega),
```
and thus $u_{\alpha} = \tilde{u}_{\alpha}$ almost everywhere
by {prf:ref}`lem:uniqueness-by-testing`.
````

````{exercise} Relation between the modulus function and the Heaviside function
Let $\Omega = (-1,1)$ and set
```{math}
\begin{align}
u(x) &= |x| \\
H(x) &= \begin{cases} -1 &\quad x \in (-1,0) \\ 
                      1  &\quad x \in [0, 1)
        \end{cases}
\end{align}
```
By simply using the definition of the weak derivative, show that
$H(x)$ is the weak derivative of $u$.
````

````{prf:definition} Sobolev spaces
* $W^{k,p}(\Omega) := 
  \{
  u \in L^p(\Omega) |\, \partial^{\alpha}u \text{ exists and belongs to } L^p(\Omega) 
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
consisting of those function $\phi$ in $W^{k,p}(\Omega)$ which are limits
of sequences $\{\phi_n\}_{n=1}^\infty \subset C_c^{\infty}(\Omega)$.
```
Later we will need the following important result known as Poincaré inequality.

````{prf:theorem} Poincaré inequality
:label: thm:poincare
Let $\Omega$ be an open and bounded subset of $\RR^n$ and suppose
then there is a constant $C_P = C_P(p,n,\Omega)$ such that
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
On $W^{1,p}_0(\Omega)$, the $\| \cdot\|_{W^{1,p}(\Omega)}$ is equivalent to the norm
```{math}
\| u \|_{W^{1,p}_0(\Omega)} := \| \nabla u \|_{L^{p}(\Omega)}  
```
````

````{prf:proof}
A simple application of the Poincaré application yields
```{math}
\|\nabla u\|_{\Omega}^p
\leqslant 
\| u \|_{\Omega}^p +
\|\nabla u\|_{\Omega}^p
\leqslant
(1+C_P^p) \|\nabla u\|_{\Omega}^p.
```
````

<!-- ### Approximation results -->

<!-- ### Poincaré inequalties  -->

### Trace operators
Next, we very briefly discuss whether and how functions of certain Sobolev spaces defined
on the domain $\Omega$ can be restricted to the boundary $\partial \Omega$. This plays
an important role in the well-posedness of boundary value problems, as we need to determine 
the correct spaces for the boundary data in a e.g. Dirichlet or Neumann boundary problem
when the data is **non-homogeneous**.

For the remaining part of this Chapter, we assume that $\Omega$ is a bounded domain with 
a "well-behaving" boundary, that is, it is either a $C^1$-polyhedron, a Lipschitz domain, or --- if this doesn't tell  you much ---
simply a $C^{\infty} domain$.

````{prf:theorem} Traces of $H^1(\Omega)$ spaces
:label: thm:trace-spaces
For a bounded domain $\Omega$ with Lipschitz (or $C^{\infty}$) boundary
$ \Gamma = \partial \Omega$, there
exists a bounded operator $\gamma : H^1(\Omega) \to L^2(\Gamma)$ (the so-called *Trace Operator*) such that 
$\gamma(u) = u|_{\Gamma}$ whenever $u \in C(\overline{\Omega})$.

If such a trace operator exists, then one can show that
```{math}
H^1_0(\Omega) = \mathrm{ker} \gamma = \{v \in H^1(\Omega) \mid \gamma(v) = 0 \}.
```
````

It turns out that the trace operator $\gamma$ **is not onto $L^2(\Omega)$**. Thus, when we later want
to find certain weak formulations and solutions $u \in H^1(\Omega)$ which also need to satisfy certain
inhomogeneous boundary conditions such as $u = u_D$ on $\Gamma$, we need to be careful about
the choice of function space from which we take the boundary data $u_D$.
That motivates the following

````{prf:definition} $H^{1/2}(\Gamma)$
:label: def:honehalf
We set
```{math}
H^{1/2}(\Omega) = \{ v \in L^2(\Omega) \mid \gamma(\overline{v}) = v \text{ for some } \overline{v} \in H^1(\Omega) \}
```
and define a corresponding norm by
```{math}
  \|v \|_{H^{1/2}(\Gamma)} := \|v\|_{1/2, \Gamma} := \inf \{ \|\overline{v}\|_{1,\Omega} \mid \gamma(\overline{v}) = v\}.
```
Consequently, 
```{math}
\|v\|_{1/2,\Gamma} \leqslant \| v\|_{1, \Omega}.
```
````
