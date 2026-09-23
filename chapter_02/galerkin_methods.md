# A primer on finite element methods

## Galerkin's method: A recipe to discretize partial differential equations

### The general recipe

Galerkin's method is an general approach to solve partial differential equation numerically 
by transforming them into a system of discrete equations. 
The computed solution to the discrete equations can then be thought of as an approximation
to the solution of the original PDE. Figure {numref}`fem-steps` summarizes the four stages
of the discretization approach which we describe next. As usual, we will use
the Poisson problem as a guiding prototype example.

```{note}
The following recipe is deliberately kept vague. Details such as boundary conditions,
suitable function spaces etc. need to be taken into account accordingly. 
```

% TODO Refer back to PP

```{figure} ./fem_steps_stand_alone.png
---
alt: 4 steps in the FEM recipe
align: center
name: fem-steps
---

The four stages of Galerkin's method to discretize partial differential equations.
```

Stage 1 
: Strong formulation of the PDE

Starting point is a partial differential equation

$$
\mcA u = f 
$$ (strong-form)

in *strong form*: for given function $f:\Omega \subset \RR^n \to \RR$ and partial differential operator $\mcA$,
we assume that the function $u: \Omega \to  \RR$ satisfies
the relation {eq}`strong-form` *pointwise* so that $\mcA u(x) = f(x)\: \forall x \in \Omega$.

Our favorite example is of course the Poisson problem with *homogeneous Dirichlet boundary condition*:
\begin{alignat}{2}
- \Delta u  &= f & &\quad \text{in } \Omega, \\
         u  &= 0 & &\quad \text{on } \Gamma.
\end{alignat} (poisson-strong-form)

Stage 2
: Continuous weak formulation of the PDE

Find $u \in V$ such that

$$
a(u,v) = l(v)  \quad \forall v \in V.
$$ (cont-weak-form)

The standard approach to obtain a weak formulation is to multiply
with the strong form of the PDE with appropriate test functions
which satisfy appropriate smoothness assumptions and ---if required--- boundary conditions.

For the Poisson problem with homogeneous Dirichlet b.c. we saw previously that we arrived 
at the weak formulation:
Find $u \in V := H^1_0(\Omega)$ s.t. $\forall v \in V$

$$
\underbrace{\int_{\Omega} \nabla u \cdot \nabla v \dx}_{a(u,v)} = \underbrace{\int_{\Omega} f v \dx}_{l(v)}.
$$

Stage 3
: Discrete weak formulation of the PDE

By choosing a suitable approximation space $V_h \subset V$
with finite dimension $N = \dim V_h$,
we obtain the discrete weak formulation:

Find $u_h \in V_h$ such that

$$
a(u_h,v_h) = l(v_h)  \quad \forall v_h \in V_h.
$$ (disc-weak-form)

Stage 4
: Formulation as system of discrete equations

To translate the now finite-dimensional problem {eq}`disc-weak-form` into a discrete
system of equations, we choose a basis for the discrete function space

$$
V_h = \spann \{\phi_i\}_{i=1}^N.
$$

First observe that the problem {eq}`disc-weak-form` is then equivalent to
seek a $u_h \in V_h$ such that

$$
a(u_h,\phi_h) = l(\phi_i)  \quad i = 1, \ldots N,
$$ (disc-weak-form-with-basis-function)

since $a(\cdot, \cdot)$ and $l(\cdot)$
are linear in their respectively second and first slot.

The second step is now to rewrite $u_h = \sum_{j=1}^N U_j \phi_j$ with the help of the basis functions
and to insert this representation into {eq}`disc-weak-form-with-basis-function` to obtain
a discrete system of equations for the coefficient vector $(U_i)_i^N \in \RR^N$.
If $a(\cdot,\cdot)$ is *not* linear in the first slot ---think of nonlinear PDEs such as the Navier-Stokes equations---,
then the resulting system is truly nonlinear. But for weak formulation one *linear* PDEs, 
the linearity in the first slot of $a(\cdot, \cdot)$ allows us to to cast
{eq}`disc-weak-form-with-basis-function` into a linear system

$$
A U = b
$$ (discrete-system)

with $A \in \RR^{N\times N}$ and $b \in \RR^N$
since

$$
a(u_h,\phi_h)
= 
a\Bigl(\sum_{j=1}^N U_j \phi_j, \phi_h \Bigr)
=
\sum_{j=1}^N \underbrace{a(\phi_j, \phi_i)}_{=: A_{ij}} U_j
= \underbrace{l(\phi_i)}_{=: b_i}  \quad i = 1, \ldots N.
$$ (discrete-system-derivation)

### Abstract error theory

````{prf:lemma} Galerkin orthogonality
:label: lem-galerkin-ortho

Assume that $V_h \subset V$ and that
$u \in V$ and 
$u_h \in V_h$
solve the continuous weak formulation {eq}`cont-weak-form`
and the discrete weak formulation {eq}`disc-weak-form`, respectively.
Then the error $u-u_h$ satisfies the orthogonality relation

$$
a(u-u_h, v_h) = 0 \quad \forall\, v_h \in V_h.
$$ (eq:galerkin-ortho)

````

````{prf:proof}
If $v_h \in V_h \subset V$, then the continuous $u$ and the discrete solution $u_h$ satisfy

$$
\begin{aligned}
a(u, v_h)   &= l(v_h), \\
a(u_h, v_h) &= l(v_h).
\end{aligned}
$$

Subtracting the second equality from the first yields {eq}`eq:galerkin-ortho`.
````

````{prf:lemma} Cea's lemma
:label: ceas-lemma

Assume that 
* $V_h \subset V$ 
* the continuous weak formulation {eq}`cont-weak-form` satisfies the assumptions of the Lax-Milgram theorem.
* $u \in V$   is solution to the continuous weak formulation {eq}`cont-weak-form`
* $u_h \in V_h$ is solution to the discrete weak formulation {eq}`disc-weak-form`

Then $u_h$ satisfies a quasi best approximation property in the sense that

$$
\| u - u_h \| \leqslant \dfrac{C_a}{\alpha} \inf_{v \in V_h} \| u - v_h\|.
$$ (eq:ceas-lemma)

holds for the error $u-u_h$.
Here $C_a$ and $\alpha$ are the boundedness and ellipticity constants for $a(\cdot, \cdot)$
appearing the assumptions for the Lax-Milgram theorem.

````

````{prf:proof}

Let $v_h \in V_h$ be fixed but arbitrary, then we wish to show that

$$
\| u - u_h \| \leqslant \dfrac{C_a}{\alpha} \| u - v_h\|.
$$ (eq:ceas-lemma-alt)

The proof of this inequality is rather short. Its main essence consists of three estimates 
where we first use coercivity of $a(\cdot, \cdot)$ to relate $\|u-u_h\|$ to $a(\cdot, \cdot)$,
then add and substract the given $v_h$ and apply Galerkin orthogonality, and finally 
boundedness of $a(\cdot, \cdot)$ is exploited to estimate the resulting expression.
To this end, we see that

$$
\begin{aligned}
\alpha \| u - u_h \|^2 &\leqslant a(u - u_h, u - u_h)  \\
                       &= a(u-u_h, u - v_h + v_h - u_h)  \\ 
                       &= a(u-u_h,  u - v_h) + \underbrace{a(u-u_h, v_h - u_h)}_{= 0} \\
                       &\leqslant C_a \| u - u_h \| \| u - v_h\|
\end{aligned}
$$ (eq:ceas-lemma-proof)

Assuming that $\| u - u_h \| \neq 0$ (otherwise 
{eq}`eq:ceas-lemma-alt` is trivially satisfied), we can divide {eq}`eq:ceas-lemma-proof` by $\| u - u_h \|$ and 
$\alpha$ to see that

$$
\| u - u_h \| \leqslant \dfrac{C}{\alpha}\| u - v_h\|.
$$

````
