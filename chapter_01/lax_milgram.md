(sec:lax-milgram)=
# The Lax-Milgram theorem

<!-- Later, when we have introduced the concept of weak formulation of partial differential equations, we will make heavily use of the  Lax-Milgram theorem. -->
Throughout this section, we rely on the concepts and results collected in
{ref}`sec:functional-analysis`, in particular the
{prf:ref}`Riesz representation theorem<thm-riesz-representation>`, the
{prf:ref}`orthogonal decomposition theorem<thm-orthogonal-decomposition>` and the
{prf:ref}`Banach fixed-point theorem<thm-banach-fixed-point>`.

````{prf:theorem} Lax-Milgram
:label: thm-lax-milgram
Given a real Hilbert space $\bigl(V, (\cdot,\cdot)\bigr)$ with induced norm
$\| \cdot \| := (\cdot,\cdot)^{1/2}$, a bilinear form
$a(\cdot, \cdot): V \times V \to \mathbb{R}$, and a linear form
$l(\cdot): V \to \mathbb{R}$, the problem: Find $u \in V$ such 
that 
```{math}
:label: eq:lax-milgram-problem
a(u, v) = l(v) \quad \forall v\in V
```
possesses a unique solution $u \in V$ if the following assumptions are satisfied.
1. The linear form $l$ is bounded, i.e. there exists a constant $C_l \geqslant 0$ such that
    ```{math}
    :label: eq:lax-milgram-bounded-l
    | l(v) | \leqslant C_l \| v\| \quad \forall v \in V.
    ```
2. The bilinear form $a$ is bounded, i.e. there exists a constant $C_a \geqslant 0$ such that
    ```{math}
    :label: eq:lax-milgram-bounded-a
    | a(v, w) | \leqslant C_a \| v\| \|w\| \quad \forall v,w \in V.
    ```
3. The bilinear form $a$ is coercive, i.e. there is a constant $\alpha > 0$ such that
    ```{math}
    :label: eq:lax-milgram-coerc
    a(v, v)  \geqslant \alpha \|v\|^2 \quad \forall v \in V.
    ```
Moreover, the solution $u$ satisfies the stability (or a priori) estimate
```{math}
:label: eq:lax-milgram-stab
\|u\| \leqslant \dfrac{C_l}{\alpha}.
```
````

````{prf:proof}

If $V = \{0\}$, the assertion is trivial, so we assume $V \neq \{0\}$ from now on.
In that case, testing {eq}`eq:lax-milgram-bounded-a` and {eq}`eq:lax-milgram-coerc`
with $v = w \neq 0$ shows that $0 < \alpha \leqslant C_a$; in particular $C_a > 0$.

**Step 1: Reformulation as an operator equation.**

Fix $u \in V$. Then $v \mapsto a(u, v)$ is a linear form on $V$, which by the
boundedness {eq}`eq:lax-milgram-bounded-a` of $a$ is continuous with
$\| a(u, \cdot) \|_{V'} \leqslant C_a \|u\|$. The
{prf:ref}`Riesz representation theorem<thm-riesz-representation>` therefore provides a **unique** element
of $V$, which we denote by $A u$, such that
```{math}
:label: eq:lax-milgram-def-A
a(u, v) = (A u, v) \quad \forall v \in V,
\qquad \text{and} \qquad
\|A u\| = \| a(u, \cdot) \|_{V'} \leqslant C_a \|u\|.
```
This defines a map $A: V \to V$, which is linear: for $u_1, u_2 \in V$ and
$\lambda \in \mathbb{R}$, the bilinearity of $a$ gives
```{math}
(A(u_1 + \lambda u_2), v)
= a(u_1 + \lambda u_2, v)
= a(u_1, v) + \lambda \, a(u_2, v)
= (A u_1 + \lambda A u_2, v)
\quad \forall v \in V,
```
and since the only element of $V$ that is orthogonal to *all* $v \in V$ is $0$, we
conclude that $A(u_1 + \lambda u_2) = A u_1 + \lambda A u_2$. Together with
{eq}`eq:lax-milgram-def-A`, the operator $A$ is thus linear and bounded with
$\|A\|_{V \to V} \leqslant C_a$.

In exactly the same way, the {prf:ref}`Riesz representation theorem<thm-riesz-representation>` applied to the
bounded linear form $l$ yields a unique $f \in V$ with
```{math}
:label: eq:lax-milgram-def-f
l(v) = (f, v) \quad \forall v \in V,
\qquad \text{and} \qquad
\|f\| = \|l\|_{V'} \leqslant C_l.
```

In total, we thus have
```{math}
(A u, v) = a(u, v) = l(v) = (f, v) \quad \forall v \in V,
```
and therefore, for the element $f \in V$ representing $l$ in the sense of
{eq}`eq:lax-milgram-def-f`, problem {eq}`eq:lax-milgram-problem`
is equivalent to the operator equation
```{math}
:label: eq:lax-milgram-operator-eq
A u = f.
```
Unique solvability of {eq}`eq:lax-milgram-operator-eq` for *every* $f \in V$ is
in turn equivalent to the statement that $A$ is bijective,
while the stability estimate {eq}`eq:lax-milgram-stab` is equivalent to
the boundedness of the inverse operator $A^{-1}$.

**Step 2: Analysis of the operator $A$.**
For a non-zero $u \in V$, 
```{math}
\alpha \|u\|^2 \leqslant a(u, u) = (A u, u) \leqslant \|A u\| \, \|u\|
```
thanks to coercivity {eq}`eq:lax-milgram-coerc` and the Cauchy-Schwarz
inequality, and dividing by $\|u\|$ and $\alpha>0$ gives
the important lower bound for $Au$
```{math}
:label: eq:lax-milgram-A-lower-bound
\|u\| \leqslant \frac{1}{\alpha} \|A u\|.
```
(For $u=0$ the inequality is trivially satisfied.) This lower bound is crucial;
we will use it to show that $A$ is injective and that its range
```{math}
:label: eq:lax-milgram-range
\mathrm{ran}(A) := \{ A u : u \in V \}
```
is closed.

**Injectivity of $A$.** If $A u = 0$ for some $u \in V$, then {eq}`eq:lax-milgram-A-lower-bound` implies that $\|u\| = 0$, and hence $u = 0$.

**Closedness of $\mathrm{ran}(A)$.** Since $A$ is linear, its range
is a linear subspace of $V$. To show that it is closed, let
$(w_n)_{n \in \mathbb{N}} \subset \mathrm{ran}(A)$ be a sequence converging to some
$w \in V$, 
and pick $u_n \in V$ with $A u_n = w_n$. Applying the lower bound
{eq}`eq:lax-milgram-A-lower-bound` to $u_n - u_m$ and using the linearity of $A$,
we obtain
```{math}
:label: eq:lax-milgram-cauchy
\|u_n - u_m\| \leqslant \dfrac{1}{\alpha} \| A(u_n - u_m) \|
= \dfrac{1}{\alpha} \| w_n - w_m \|
\quad \forall n, m \in \mathbb{N}.
```
Being convergent, $(w_n)_n$ is a Cauchy sequence in $V$, and hence so is $(u_n)_n$
by {eq}`eq:lax-milgram-cauchy`. Since the Hilbert space $V$ is complete, $(u_n)_n$
converges to some $u \in V$. As $A$ is bounded and therefore continuous, we may pass
to the limit and conclude
```{math}
w = \lim_{n \to \infty} w_n = \lim_{n \to \infty} A u_n = A u \in \mathrm{ran}(A).
```
Thus $\mathrm{ran}(A)$ contains all its limit points, i.e. $\mathrm{ran}(A)$ is closed.

Note how the lower bound {eq}`eq:lax-milgram-A-lower-bound` is exactly what allows us
to transfer the Cauchy property from the image sequence $(w_n)_n$ back to the
preimage sequence $(u_n)_n$; without coercivity this step would fail.

**Step 3: Surjectivity of $A$.**

Since $\mathrm{ran}(A)$ is a *closed* subspace of the Hilbert space $V$, the
{prf:ref}`orthogonal decomposition theorem<thm-orthogonal-decomposition>` gives
```{math}
:label: eq:lax-milgram-orth-decomp
V = \mathrm{ran}(A) \oplus \mathrm{ran}(A)^{\perp}.
```
By {eq}`eq:orthogonal-decomposition-trivial`, it therefore suffices to show that
$\mathrm{ran}(A)^{\perp} = \{0\}$.

So let $v \in \mathrm{ran}(A)^{\perp}$, that is,
```{math}
(A u, v) = 0 \quad \forall u \in V.
```
Choosing the particular test element $u = v$ and employing the coercivity
{eq}`eq:lax-milgram-coerc` together with the definition {eq}`eq:lax-milgram-def-A`
of $A$, we find
```{math}
\alpha \|v\|^2 \leqslant a(v, v) = (A v, v) = 0,
```
and hence $\|v\| = 0$, i.e. $v = 0$, since $\alpha > 0$. Consequently
$\mathrm{ran}(A)^{\perp} = \{0\}$ and thus $\mathrm{ran}(A) = V$ by
{eq}`eq:lax-milgram-orth-decomp`, which means that $A$ is surjective.

Combining this with the injectivity established in Step 2, the operator
$A : V \to V$ is a **bijection**. Hence, for the element $f \in V$ representing $l$
in the sense of {eq}`eq:lax-milgram-def-f`, the vector $u := A^{-1} f$ is the unique
solution of the operator equation {eq}`eq:lax-milgram-operator-eq`, and by the
equivalence established in Step 1, the unique solution of problem
{eq}`eq:lax-milgram-problem`. This proves both **existence** and **uniqueness**.

**Step 4: Stability estimate.**

The lower bound {eq}`eq:lax-milgram-A-lower-bound` states precisely that the inverse
operator $A^{-1} : V \to V$ is bounded with $\|A^{-1}\|_{V \to V} \leqslant 1/\alpha$.
Together with the norm identity in {eq}`eq:lax-milgram-def-f`, the solution
$u = A^{-1} f$ therefore satisfies
```{math}
\|u\| \leqslant \dfrac{1}{\alpha} \| A u \| = \dfrac{1}{\alpha} \| f \|
= \dfrac{1}{\alpha} \| l \|_{V'} \leqslant \dfrac{C_l}{\alpha},
```
which is the desired stability estimate {eq}`eq:lax-milgram-stab`.

Alternatively, {eq}`eq:lax-milgram-stab` can be derived directly from the variational
problem, without any reference to the operator $A$. Assume $u$ solves
{eq}`eq:lax-milgram-problem`. Then set $v=u$ and
successively employ the coercivity of $a$ and boundedness of $l$ to see that
```{math}
\alpha \|u\|^2 \leqslant a(u, u) = l(u) \leqslant C_l \|u\|.
```
Dividing the previous chain of inequalities by $\alpha$ and $\|u\|$ if $\| u \| \neq 0$
yields {eq}`eq:lax-milgram-stab`. For $\|u\| = 0$ the stability estimate is trivially satisfied.

In the same spirit, uniqueness --- which we already obtained from the injectivity of
$A$ --- also follows directly from the stability estimate: if $u_1$ and $u_2$
both satisfy problem {eq}`eq:lax-milgram-problem`, then
thanks to linearity of $a$ in the first slot,
the difference $u_1-u_2$ satisfies problem {eq}`eq:lax-milgram-problem` but with $l =  0$
instead. In that case $C_l = 0$ and thus $0\leqslant\|u_1 - u_2 \| \leqslant \tfrac{0}{\alpha} = 0$, and thus $u_1 = u_2$.
````

The existence part of the previous proof was non-constructive: it relied on the
{prf:ref}`orthogonal decomposition theorem<thm-orthogonal-decomposition>` to conclude
that $\mathrm{ran}(A) = V$. The following alternative argument replaces Step 3
by a fixed-point iteration, which has the advantage of being constructive;
the stability estimate is then obtained as in the second, direct argument of
Step 4.

````{admonition} Alternative existence proof via a contraction argument
:class: dropdown

The idea is to turn {eq}`eq:lax-milgram-operator-eq` into a fixed-point problem.
For a parameter $\rho > 0$, which will be chosen below, define
```{math}
:label: eq:lax-milgram-fixed-point-map
T_{\rho} : V \to V, \qquad T_{\rho} v := v - \rho \, (A v - f).
```
Since $\rho > 0$, an element $u \in V$ satisfies $T_{\rho} u = u$ if and only if
$A u = f$, that is, if and only if $u$ solves {eq}`eq:lax-milgram-problem`. We now
show that $\rho$ can be chosen such that $T_{\rho}$ is a contraction.

Let $v, w \in V$ and set $z := v - w$. By the linearity of $A$ we have
$T_{\rho} v - T_{\rho} w = z - \rho A z$, and therefore, using the coercivity
{eq}`eq:lax-milgram-coerc` in the form $(A z, z) = a(z,z) \geqslant \alpha \|z\|^2$
together with the bound $\|A z\| \leqslant C_a \|z\|$ from
{eq}`eq:lax-milgram-def-A`,
```{math}
\begin{aligned}
\| T_{\rho} v - T_{\rho} w \|^2
&= \| z - \rho A z \|^2 \\
&= \|z\|^2 - 2 \rho \, (A z, z) + \rho^2 \|A z\|^2 \\
&\leqslant \|z\|^2 - 2 \rho \alpha \|z\|^2 + \rho^2 C_a^2 \|z\|^2 \\
&= \underbrace{\bigl( 1 - 2 \rho \alpha + \rho^2 C_a^2 \bigr)}_{=: \, k(\rho)^2}
   \, \| v - w \|^2 .
\end{aligned}
```
The quadratic $\rho \mapsto 1 - 2\rho\alpha + \rho^2 C_a^2$ is $< 1$ precisely for
$0 < \rho < 2\alpha / C_a^2$. Choosing for instance $\rho := \alpha / C_a^2$ gives
```{math}
:label: eq:lax-milgram-contraction
k(\rho)^2 = 1 - \dfrac{\alpha^2}{C_a^2} \in [0, 1),
```
so $T_{\rho}$ is a contraction with contraction constant $k(\rho) < 1$. As a Hilbert
space, $V$ is complete with respect to the metric induced by $\|\cdot\|$, so the
{prf:ref}`Banach fixed-point theorem<thm-banach-fixed-point>` applies and yields a
**unique** $u \in V$ with $T_{\rho} u = u$. By the equivalence observed above, this
$u$ is the unique solution of {eq}`eq:lax-milgram-problem`.
````

```{prf:remark}
Note that the proof does **not** require the bilinear form $a$ to be symmetric.
For symmetric $a$, the Lax-Milgram theorem reduces to the
{prf:ref}`Riesz representation theorem<thm-riesz-representation>`,
as you are asked to show in following exercise.
```

```{exercise} Lax-Milgram for symmetric bilinear forms
:label: exer-lax-milgram-symmetric
Assume in addition that $a$ is *symmetric*, i.e. $a(v, w) = a(w, v)$ for all
$v, w \in V$. 
Prove the {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>` in this simplified setting,
including the stability estimate {eq}`eq:lax-milgram-stab`.
```

```{hint}
:class: dropdown
First show that $a(\cdot, \cdot)$ defines an inner product on $V$ whose
induced *energy norm* $\|v\|_a := \sqrt{a(v,v)}$ is equivalent to $\|\cdot\|$,
thanks to {eq}`eq:lax-milgram-bounded-a` and {eq}`eq:lax-milgram-coerc`.
Then apply the
{prf:ref}`Riesz representation theorem<thm-riesz-representation>` to the Hilbert
space $\bigl(V, a(\cdot, \cdot)\bigr)$.
```

```{prf:remark}
The alternative contraction argument is constructive: the Picard iteration
$u_{n+1} = u_n - \rho\,(A u_n - f)$ associated with {eq}`eq:lax-milgram-fixed-point-map`
is nothing but a *preconditioned Richardson iteration* for the operator equation
{eq}`eq:lax-milgram-operator-eq`, and {eq}`eq:lax-milgram-contraction` quantifies its
convergence rate in terms of $\alpha$ and $C_a$. 
```

```{prf:remark}
The {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>` ensures that 
problem {eq}`eq:lax-milgram-problem`
is well-posed, i.e.,
* **Existence** of a solution
* **Uniqueness** of the solution
* **Continuous dependence on the data** (or **Stability**) of the solution.
  In the particular case of the {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>`,
  stability is guaranteed through {eq}`eq:lax-milgram-stab`: applying it to the
  difference of the solutions $u_1, u_2$ associated with two right-hand sides
  $l_1, l_2$ yields
  $\|u_1 - u_2\| \leqslant \alpha^{-1} \| l_1 - l_2 \|_{V'}$, so that "small
  changes" in the data $l$ will only lead to small changes in the solution $u$.
  Perturbations of the bilinear form $a$ require a separate argument, which we
  will encounter later in the form of Strang's lemma.
```
