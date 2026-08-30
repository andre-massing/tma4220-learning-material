(sec:functional-analysis)=
# Relevant concepts from functional analysis

Vector space
: A definition can be found on [wiki](https://en.wikipedia.org/wiki/Vector_space).

[Metric space](https://en.wikipedia.org/wiki/Metric_space)
: A metric space is a set $X$ which is equipped with a [distance function (or metric)](https://en.wikipedia.org/wiki/Metric_space#Definition)

```{math}
d(x,y): X \times X \to \RR.
```

[Complete metric space](https://en.wikipedia.org/wiki/Complete_metric_space)
: A metric space is called *complete* if every [Cauchy sequence](https://en.wikipedia.org/wiki/Cauchy_sequence#In_a_metric_space)
converges to some $x \in X$.

[Normed vector space](https://en.wikipedia.org/wiki/Metric_space)
: A vector space $(V, \|\cdot\|_V)$ consist of a vector space $V$ which is 
equipped with a [norm](https://en.wikipedia.org/wiki/Norm_(mathematics)) 

```{math}
\| \cdot \|_V : V \to \RR
``` 

Note that every norm induces a natural metric $d(x, y) := \|x-y\|_V$.
Typically we do not use the verbose notation $(V, \|\cdot\|_V)$, instead
we simply speak of a normed vector space $V$, and we omit the subscript $_V$
in the norm symbol when the norm is clear from
the context.

[Banach space](https://en.wikipedia.org/wiki/Banach_space)
: A normed vector space which is complete with respect to the induced metric.

[Inner product space](https://en.wikipedia.org/wiki/Inner_product_space)
: An inner product space $\bigl(V, (\cdot, \cdot)\bigr)_V$ is a real 
(or complex) vector space $V$ equipped with a [inner product](https://en.wikipedia.org/wiki/Inner_product_space#Basic_properties)
    
```{math}
(\cdot, \cdot)_V : V \times V \to \RR \quad (\text{or } \CC )
``` 

Every inner product induces a natural norm $\| \cdot \| := \sqrt{(\cdot, \cdot)}$, and thereby a metric. 
Again, we typically do not use the verbose notation $\bigl(V, (\cdot,
\cdot)\bigr)$, instead we simply speak of a inner product space $V$,
and we often omit the subscript $_V$
in $(\cdot, \cdot)_V$ symbol when the inner product is clear from
the context.

Inner products satisfy the *Cauchy-Schwarz inequality:*
```{math}
(u,v)_V \leqslant \|u\|_V \|v\|_V.
```
<!-- and inequality only holds -->


[Bounded linear operator](https://en.wikipedia.org/wiki/Operator_(mathematics)#Bounded_operators)
: A linear operator $L: V \to W$ between two normed vector spaces 
$(V, \|\cdot\|_V)$ and $(W, \|\cdot\|_W)$
is call bounded if there is a constant $C \in \RR^+_0$ such that
```{math}
\| L v \|_W \leqslant C \|v\|_V.
```

The *operator norm* $\|L\|_{V\to W}$ of $T$ is then the smallest such constant given by 
```{math}
\begin{align}
\|L\|
&= \inf \{C \in \RR^+_0 : \|L v \|_W \leqslant \|v\|_V \, \forall v \in V\} \\
& = \sup_{v \in V \setminus \{0\}} \dfrac{\|L v \|_W}{\|v\|_V} \\
& = \sup_{v \in V, \|v\|_V = 1} \|L v \|_W.
\end{align}
```
It can be shown that the the following statements are equivalent for **linear operators**:
* $L: V \to W$ is bounded
* $L: V \to W$ is continuous

<!-- See {cite}`Brezis2011` for a proof. -->

```{exercise} Equivalence of boundedness and continuity
Before you look up the proof, try to prove the previous claim yourself.
```

A linear operator $l : V \to \RR \; (\text{or } \CC )$ is often called
a *linear functional* or a *linear form* on $V$.

Dual space
: The dual space $V^*$ for a normed vector space $(V, \|\cdot\|)$ consist
of all **continuous** linear functionals defined on $V$.

Note that for inner product spaces $V$, every $u \in V$ give rise to a 
continuous linear functional $l_u$ defined by
```{math}
l_u(v) := (v, u)_V \quad \forall v \in V.
```

For Hilbert space $H$, that is in essence all the continuous linear functionals
you can construct on $H$ thanks to the following theorem.

Riesz representation theorem 
````{prf:theorem} Riesz representation theorem
:label: thm-riesz-representation

Let $H$ be a Hilbert space with a inner product $(\cdot, \cdot)$. Then for
every continuous functional $l:H \to \RR$, there is a unique vector $u_l \in H$
such that
```{math}
l(v) = (v, u_l) \quad \forall v \in H,
```
and we have that
```{math}
\| l_u \|_{V^*} = \| u \|_{V}.
```
````

```{prf:proof}
For a proof, we refer to Section 5.2 in {cite}`Brezis2011`. 
```

The proof of the {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>` in
{ref}`sec:lax-milgram` will also rely on the orthogonal
decomposition of a Hilbert space with respect to a *closed* subspace.

````{prf:theorem} Orthogonal decomposition
:label: thm-orthogonal-decomposition
Let $H$ be a Hilbert space and let $M \subset H$ be a **closed** subspace. Then
```{math}
H = M \oplus M^{\perp},
\qquad \text{where} \quad
M^{\perp} := \{ w \in H : (v, w) = 0 \;\; \forall v \in M \},
```
that is, every $x \in H$ can be written uniquely as $x = m + m^{\perp}$ with
$m \in M$ and $m^{\perp} \in M^{\perp}$. In particular,
```{math}
:label: eq:orthogonal-decomposition-trivial
M^{\perp} = \{0\} \quad \implies \quad M = H.
```
````

```{prf:proof}
For a proof, we refer to Section 5.1 in {cite}`Brezis2011`.
```

Finally, the Banach fixed-point theorem turns a contraction on a complete metric
space into a solvable equation. We will use it to give a second, constructive
existence proof for the {prf:ref}`Lax-Milgram theorem<thm-lax-milgram>`.

````{prf:theorem} Banach fixed-point theorem
:label: thm-banach-fixed-point
Let $(X, d)$ be a non-empty complete metric space and let $T: X \to X$ be a
*contraction*, that is, there exists a constant $k \in [0, 1)$ such that
```{math}
:label: eq:banach-contraction
d(T x, T y) \leqslant k \, d(x, y) \quad \forall x, y \in X.
```
Then $T$ possesses exactly one fixed point $x^{\ast} \in X$, i.e. $T x^{\ast} = x^{\ast}$.
````


