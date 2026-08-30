(sec:weak-form-pdes)=
# Weak formulation of elliptic PDEs

## Weak formulation of the Poisson problem
In this chapter, we briefly discuss how the functional
analysis and function space apparatus
can be employed to analyse the well-posedness
of certain class of PDEs when given in a so-called
"weak" formulation. 

We start by considering the Poisson problem
```{math}
:label: eq:poisson-wo-bc
\nabla \cdot \nabla u =  -\Delta u  = f \quad \text{in } \Omega
```
supplemented with some suitable boundary conditions 
which $u$ should satisfy on 
the boundary $\Gamma = \partial \Omega$ of $\Omega$.

## Neumann problems
Let us consider the homogenous Neumann problem
```{math}
:label: eq:neumann-problem-I
\left\{
\begin{alignedat}{2}
- \Delta u + u &= f & &\quad \text{in } \Omega \\
         \partial_n u &= 0 & &\quad \text{on } \Gamma
\end{alignedat}
\right.
```
Here, we used the slightly simplified notation $\partial_n u = \bfn \cdot \nabla u$.
The idea to derive a so-called weak formulation of an PDE is very similar
to the idea behind the introduction of weak derivatives:
We multiply with a suitable test function $v$, integrate over $\Omega$ 
and perform integration by parts to transfer a number of derivatives 
to the test function $v$. 

What kind of test function space we choose is often dictated by 2 considerations:

1. What kind of smoothness do we require to make the derived formulation work?

2. How do we take into account the boundary conditions?

For the Neumann boundary problem, let's assume for the
moment that $u$, our boundary $\Gamma$ and our
test functions $v$ are smooth enough so that we can use Green's theorem,
e.g, $u \in C^2(\overline{\Omega})$, $\Gamma$ is a $C^1$ boundary, and
$v \in C^{\infty}(\overline{\Omega})$. Then multiplying 
the PDE in {eq}`eq:neumann-problem-I` with $v$ and integrating over $\Omega$
and applying Green's theorem leads to
```{math}
:label: eq:weak-form-neumann-derivation
\begin{aligned}
\int_{\Omega} f v \dx
&= 
- \int_{\Omega}  \nabla \cdot (\nabla u) v \dx 
+\int_{\Omega} uv \dx
\\
&=
- \int_{\Gamma} \underbrace{(\bfn \cdot \nabla u)}_{=0} v \dx
+ \int_{\Omega} \nabla u \cdot \nabla  v \dx
+\int_{\Omega} uv \dx
\end{aligned}
```

Note that the Neumann boundary condition $\bfn \cdot \nabla u = 0$
makes the boundary integrals vanish.
Also observe that the right-hand side of {eq}`eq:weak-form-neumann-derivation`
can be interpreted as taking 
the inner product associated with $H^1(\Omega)$ between $u$ and $v$.
In fact, the expression makes perfectly sense even if we assume only
assume that both $u, v \in H^1(\Omega) =:V$!.
With this assumption, we can define the bilinear form 
```{math}
:label: poisson-bilinear-form-a
a(v, w) := 
\int_{\Omega} \nabla v \cdot \nabla  w \dx
+\int_{\Omega} vw \dx
```
on $V\times V$, and it is straightforward to
show that $a(\cdot, \cdot)$ (being the $H^1$ inner product itself)
satisfies the required assumptions of {prf:ref}`the Lax-Milgram theorem<thm-lax-milgram>`:
```{math}
:label: eq-lax-milgram-assump
\begin{aligned}
\text{Boundedness: }  a(v,w) &:= 
\int_{\Omega} \nabla v \cdot \nabla  w \dx
+\int_{\Omega} vw \dx = (v, w)_{H^1(\Omega)} \leqslant \|v\|_{H^1(\Omega)} \|w\|_{H^1(\Omega)}
\\
\text{Coercivity: }  a(v,v) &= 
\int_{\Omega} |\nabla v|^2 \dx + \int_{\Omega} |v|^2 \dx
= (v, v)_{H^1(\Omega)} = \|v\|_{H^1(\Omega)}^2 
\end{aligned}
```

Next, we define the linear form $l : V \to \RR$
```{math}
:label: poisson-linear-form-l
l(v) := \int_{\Omega} f v \dx = (f, v)_{L^2(\Omega)}
```
If we assume that $f \in L^2(\Omega)$, then thanks to the Cauchy-Schwarz inequality,
```{math}
|l(v)| 
= |(f, v)_{L^2(\Omega)}| 
\leqslant 
f\|_{L^2(\Omega)} \|v\|_{L^2(\Omega)}
\leqslant
f\|_{L^2(\Omega)} \|v\|_{H^1(\Omega)},
```
we can immediately conclude that $l$ is a continuous bilinear form with $C_l = \|f\|_{L^2(\Omega)}$.
Thus the {prf:ref}`the Lax-Milgram theorem<thm-lax-milgram>`
let us conclude that the problem: find $u \in H^1(\Omega)=: V$ such that $\forall v \in V$
```{math}
a(u,v) = l(v)
```
has a unique solution for every $f \in L^2(\Omega)$ with $\|u\|_{H^1(\Omega)} \leqslant \|f\|_{L^2(\Omega)}$.

````{admonition} (variants of Neumann problems)
:class: :danger :dropdown
Discuss Neumann problems when a) lower term is not present b) $g_N \neq 0$:
```{math}
\left\{
\begin{alignedat}{2}
- \Delta u  &= f & &\quad \text{in } \Omega \\
         \partial_n u &= 0 & &\quad \text{on } \Gamma
\end{alignedat}
\right.
```
and
```{math}
\left\{
\begin{alignedat}{2}
- \Delta u  &= f & &\quad \text{in } \Omega \\
         \partial_n u &= g_N & &\quad \text{on } \Gamma
\end{alignedat}
\right.
```
````

## Dirichlet problems

### Homogeneous Dirichlet problem for $-\Delta + \mathrm{Id}$ operator
Next, we consider 
```{math}
:label: eq:dirichlet-problem-I
\left\{
\begin{alignedat}{2}
- \Delta u + u &= f & &\quad \text{in } \Omega \\
         u &= 0 & &\quad \text{on } \Gamma
\end{alignedat}
\right.
```
We proceed as for the Neumann problem: we multiply with suitable test functions $v$
and integrate by part,
but this time, the boundary integral does not vanish since we don't have natural boundary conditions
to incorporate. 
To compensate, we only consider test functions $v \in C^{\infty}_c(\Omega)$
which vanish at the boundary. 
Then again, we obtain
```{math}
:label: eq:weak-form-dirichlet-derivation
\begin{aligned}
\int_{\Omega} f v \dx
&= 
- \int_{\Omega}  \nabla \cdot (\nabla u) v \dx 
+\int_{\Omega} uv \dx
\\
&=
- \int_{\Gamma} (\bfn \cdot \nabla u)\underbrace{v}_{=0} \dx
+ \int_{\Omega} \nabla u \cdot \nabla  v \dx
+\int_{\Omega} uv \dx.
\end{aligned}
```
Intuitively speaking we know how the solution $u$ is going to look like
on the boundary, namely $u=0$, so we don't need test functions which
test for how the equation "behaves" at the boundary.
Also, we now require that our function $u$ comes from a function space
where the boundary condition $u=0$ is already incorporated.
This is exactly what the $H^1_0(\Omega)$ space is made for!
So the weak formulation for {eq}`eq:dirichlet-problem-I` is:

Find $u \in V := H^1_0(\Omega)$ such that
```{math}
:label: eq:dirichlet-problem-I-weak
a(u,v) = l(v) \quad \forall v \in V,
```
where $a(\cdot, \cdot)$ and $l(\cdot)$
are defined as in {eq}`poisson-bilinear-form-a` and
{eq}`poisson-linear-form-l`, respectively.
As in the case for the homogeneous Neumann problem {eq}`eq:neumann-problem-I`,
we can show that $a$ and $l$ satisfy the assumption of the
{prf:ref}`the Lax-Milgram theorem<thm-lax-milgram>`, and therefore
we can conclude there there is a unique solution $u$ to the
weak formulation of the homogeneous Poisson problem which depends
continuously on the data $f$.

```{important}
The only but very important difference between the weak formulation
of the 
*homogeneous Neumann* problem {eq}`eq:neumann-problem-I`
and the 
*homogeneous Dirichlet* problem {eq}`eq:dirichlet-problem-I`
is the Hilbert space on which they are posed on.
```

### Homogeneous Dirichlet problem for $-\Delta$ operator
Now, we consider a slightly modified problem Poisson problem where the
low order term $u$ is left out:
```{math}
:label: eq:dirichlet-problem-II
\left\{
\begin{alignedat}{2}
- \Delta u &= f & &\quad \text{in } \Omega \\
         u &= 0 & &\quad \text{on } \Gamma
\end{alignedat}
\right.
```
Refer back to @eq:dirichlet-problem-II

Repeating the steps from the previous section, we arrive at the problem:
Find $u \in V := H^1_0(\Omega)$ such that
```{math}
:label: eq:dirichlet-problem-II-weak
a(u,v) = l(v) \quad \forall v \in V,
```
with the only distinction that $a(\cdot, \cdot)$ is now given by
```{math}
a(v, w) = \int_{\Omega} \nabla v \cdot \nabla w dx.
```
The boundedness of $a(\cdot, \cdot)$ and $l(\cdot)$ can be shown (almost) exactly as before.
But let's have a look at the coercivity/ellipticity: Setting $u=v$, we obtain 
```{math}
a(v,v) = \int_{\Omega} |\nabla v|^2
```
But thanks to the {prf:ref}`Poincaré inequality<thm:poincare>`
and {prf:ref}`cor:poincare` we not only know that $|\nabla\ \cdot |$ 
defines norm on the closed subspace $H^1_{0}(\Omega)$ but that this norm
is equivalent to the usual $H^1$-norm. Thanks to the proof of {prf:ref}`cor:poincare` we see that
```{math}
a(v,v) = \int_{\Omega} |\nabla v|^2 \geqslant (1 + C_p^2)^{-1/2} \|u\|_{1,\Omega}^2.
```


### Inhomogeneous Dirichlet problem for $-\Delta$ operator
Next, we consider 
```{math}
:label: eq:dirichlet-problem-inhomog
\left\{
\begin{alignedat}{2}
- \Delta u + u &= f & &\quad \text{in } \Omega \\
         u &= g_D & &\quad \text{on } \Gamma
\end{alignedat}
\right.
```
Compared to our previous weak formulation for the homogenous,
the main question is now: how can we incorporate the non-homogenous
Dirichlet b.c. $u = g_D$? First, we realize that
the trial function $H^1_0(\Omega)$ for the solution does not 
make sense anymore. So let's start from $H^1(\Omega)$. Then
we also observe that the data $g_D$ should be in 
$H^{1/2}(\Gamma)$, see {prf:ref}`def:honehalf` to ensure
that we can satisfy the equation $u = g_D$,
and only $u$ satisfying this b.c. should be vialable solution candidates
for our weak formulation. Thus we set
```{math}
H^1_{g_D}(\Omega) := 
\{ v \in H^1(\Omega) \st \gamma(v) = g_D \}.
```
Since $g_D \in H^{1/2}(\Gamma)$, this set is not empty.
Note that 
$H^1_{g_D}(\Omega)$ is not really a vector space whenever 
$g_D$ is not $0$ everywhere since the addition of 
two functions $u_1, u_2$ with the same non-vanishing boundary data $g_D$ will result in
a function $u$ satisfying $u = 2 g_D$!
Is that sense, $H^1_{g_D}(\Omega)$ should rather be considered as **affine** subspace:
For any $u_{g_D}$ satfying $\gamma(u_{g_D}) = g_D$, it holds that
```{math}
H^1_{g_D}(\Omega) = u_{g_D} + H^1_0(\Omega)
                  := \{ u_{g_D} + v \st v \in H^1_0(\Omega) \}
                   = \gamma^{-1}(g_D).
```
So the resulting weak formulation is
Find $u \in V := H^1_{g_D}(\Omega)$ such that
for all $v \in \widehat{V} := H^1_0(\Omega)$,
```{math}
\underbrace{\int_{\Omega} \nabla u \cdot\nabla v\dx}_{=:a(u,v)} = \underbrace{\int_{\Omega} f v \dx}_{=:l(v)}.
```
Note how in this case the trial function space and test function space
are not identical any more!  How can we prove the well-posedness of
this weak formulation? Lax-Milgram usually requires that the first and
second slot of $a(\cdot, \cdot)$ invokes elements from the same
(vector) space!  The common trick here is to "**lift**" the boundary
condition, i.e. we know that by the definition of $H^{1/2}(\Gamma)$
there must a $u_g \in H^1(\Omega)$ such that $\gamma_{u_g} = g_D$. 
Then we make the ansatz $u = u_0 + u_g$ and with $u_0 \in
H^1_0(\Omega)$, leading to the following weak formulation: find $u_0
\in H^1_0(\Omega) =: V$ such that
```{math}
:label: eq:weak-form-dirichlet-inhomog-lifted
a(u_0,v) = l(v) - a(u_g, v) =: \widetilde{l}(v) \quad \foralls v \in V.
``` 

```{exercise} Well-posedness of the lifted Dirichlet formulation
Show that 
{eq}`eq:weak-form-dirichlet-inhomog-lifted` is well-posed and that the solution $u_0$ depends continuously on the data $f$ and $g_D$.
Which of the assumptions of the Lax-Milgram theorem
requires a little extra checking? 
```

## Robin problems
Finally, we consider the Robin problem
```{math}
\left\{
\begin{alignedat}{2}
- \Delta u  &= f & &\quad \text{in } \Omega, \\
         \partial_n u &= \sigma(g_R-u) & &\quad \text{on } \Gamma,
\end{alignedat}
\right.
```
for some constant $\sigma > 0$.

```{exercise} Weak formulation of the Robin problem
:label: ex:robin-weak-formulation
Derive a weak formulation for the Robin problem and show that it is well-posed.
```

(sec:weak-form-general-pdes)=
## Weak formulation of general second order elliptic PDEs

The PDE {eq}`eq:poisson-wo-bc` is the prototype example
of a 2nd order elliptic operator.
More generally and without any significant complications, 
we can consider a more general PDE of the form
```{math}
\mcA u := - \nabla \cdot ( A \nabla u) = f
```
where $A = (a_{ij}(x))_{i,j=1}^n$ is a pointwise defined matrix.
Note that
```{math}
:label: eq:def-A-operator
\nabla \cdot ( A(x) \nabla u(x))
= 
-\sum_{i,j=1}^n \partial_{i} (a_{ij}(x) \partial_{j} u(x))
```

For most part of the remaining lectures, we will require $A(x)$
to satisfy the following definition.

```{prf:definition} Ellipticity of $\mcA$ 
:label: def:ellipticity

The partial differential operator $\mcA$ given
by {eq}`eq:def-A-operator` with coefficients
$A = (a_{ij})_{i,j=1}^N \in (L^{\infty}(\Omega))^{n\times n}$
is called elliptic  constant $\alpha > 0$ such that
* $ \lambda \cdot  A(x) \lambda \geqslant \alpha |\lambda|^2$

for any $\lambda \in \RR^n$.
```

```{prf:remark}
Note that $A \in (L^{\infty}(\Omega))^{n\times n}$ also implies that
that there exists an $\beta \geqslant 0 $ such that also
* $|A(x) \lambda| \leqslant \beta  |\lambda|$

holds for any $\lambda \in \RR^n$, and by ellipticity, we can
conclude that in fact $\beta \geqslant \alpha > 0$.
```
```{exercise} Boundedness of the coefficient matrix
Prove the statements made in the previous remark
```

```{admonition} TODO
:class: :danger :dropdown
* Relate $\mcA$ to classical Poisson problem
* Explain why general $A(x)$ is useful, e.g. anisotropic heat conduction problems
```

We now prepared to investigate the well-posedness of a number of boundary value
problems where we supplement the partial differential operator 
$\mcA$ with one of the following boundary conditions

* **Dirichlet boundary conditions** 
   Given function  $g_D: \Gamma  \to \RR$, we require that

  ```{math}
  u = u_D \quad \text{on } \Gamma
  ```
* **Neumann boundary conditions** 
   Given function  $g_N: \Gamma  \to \RR$, we require that

  ```{math}
  \bfn \cdot A \nabla u = g_N \quad \text{on } \Gamma 
  ```
* **Robin boundary conditions** 
   Given functions $g_R, \sigma: \Gamma  \to \RR$, we require that

  ```{math}
  \bfn \cdot A \nabla u = \sigma(g_R - u) \quad \text{on } \Gamma 
  ```

These boundary conditions are called *homogeneous* if $g_D$ (respectively $g_N$, $g_R$)
is zero, otherwise we deal with *inhomogeneous* boundary data.
We start by looking at the Poisson supplemented with Neumann boundary conditions

```{admonition} TODO
:class: :danger :dropdown
Later, mention possible impact on test and trial spaces.
```

