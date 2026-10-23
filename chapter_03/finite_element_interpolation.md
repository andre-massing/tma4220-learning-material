#  Finite Element Interpolation - Error Estimates

In this section, we will lay out the abstract framework for deriving error estimates
for finite element interpolation we considered in the previous section. 
This framework is built on several key ideas which we introduce in a series of lemma,
where we adapt the presentation from [@ErnGuermond2021, Chapter 11], [@GrossmannRoosStynes2007, Chapter 4.4]

:::{prf:lemma} Bramble-Hilbert
Given a bounded and connected $C^1$ polyhedron $B \subset \mathbb{R}^d$, and
let $q$ be a bounded subfunctional of $H^{k+1}(B)$, i.e.,


* $ |q(v + w)| \leqslant 
|q(v)| + |q(w)| $ for all $v,w \in H^{k+1}(B)$, and
* $|q(v)| \leqslant C \|v\|_{k+1,B}$ for all $v \in H^{k+1}(B)$,

with some constant $C\geqslant 0 $.
If $q$ vanishes on the space of polynomials $\mathbb{P}_k$, i.e.,
$q(p) = 0$ for all $p \in \mathbb{P}_k$, then there exists a constant $C>0$ such that

$$
|q (v)| \leqslant C |v|_{k+1,B} \quad \forall v \in H^{k+1}(B).
$$

:::

:::{prf:proof}

**Step 1.**
First we show that for any $v \in H^{k+1}(B)$, there exists a polynomial $p \in \mathbb{P}_k$ such that

$$
\int_{\Omega} \partial^{\alpha} (v - p) \mathrm{d}x = 0 \quad \forall |\alpha| \leqslant k.
$$ (eq:polynomial-approximation-condition)

Making the ansatz $p(x) = \sum_{|\beta| \leqslant k} c_{\beta} x^{\beta}$,
where $x^{\beta} = x_1^{\beta_1} \cdots x_d^{\beta_d}$, we see that the above condition is equivalent to solving the linear system for
the coefficients $c_{\beta}$:

$$
\sum_{|\beta| \leqslant k} c_{\beta} \int_{\Omega} \partial^{\alpha} x^{\beta} \mathrm{d}x = \int_{\Omega} \partial^{\alpha} v \mathrm{d}x, \quad |\alpha| \leqslant k.
$$ (eq:linear-system-polynomial-approximation)

We will now show that this corresponds to a trianguar linear system, which allows us to determine the coefficients $c_{\beta}$ recursively.
The crucial observation is that 
$\partial^{\alpha} x^{\beta} = 0$ 
for all $\alpha$ and $\beta$ with $\alpha_i > \beta_i$ for some $i$.
Now starting with some multi-index $\alpha$ with $|\alpha| = k$, we thus have that $\partial^{\alpha} x^{\beta} = 0$ for all $|\beta| < k.$
Moreover, for $|\beta| = k$, we have that 
$$
\partial^{\alpha} x^{\beta} = 
\begin{cases}
0 & \text{if } \alpha \neq \beta, \\
\alpha! & \text{if } \alpha = \beta,
\end{cases}
$$
where $\alpha! = \alpha_1! \cdots \alpha_d!$.

Consequently, for given $\alpha$ with $|\alpha| = k$, the equation @eq:linear-system-polynomial-approximation 
reduces to

$$
\sum_{|\beta| \leqslant k} c_{\beta} \int_{B} \partial^{\alpha} x^{\beta} \mathrm{d}x = 
c_{\alpha} \int_{B} \alpha! \, \mathrm{d}x
= \int_{B} \partial^{\alpha} v \mathrm{d}x
$$

and therefore we can solve for the coefficient $c_{\alpha}$ as

$$
c_{\alpha} = \frac{1}{\alpha! |B|}
\int_{B} \partial^{\alpha} v \mathrm{d}x.
$$

This allows us to determine all coefficients $c_{\beta}$ with $|\beta| = k$.

Next, we consider some multi-index $\alpha$ with $|\alpha| = k-1$. 
Then for any $\beta$ with $|\beta| \leqslant k-1,$ we have 
we have again that $\partial^{\alpha} x^{\beta} =
\delta_{\alpha \beta} \cdot \alpha!$. Of course, for any $\beta$ with $|\beta| = k$, $\partial^{\alpha} x^{\beta} \neq 0$ in general, but those coefficients $c_{\beta}$ have already been determined in the previous step.
Consequently, we deduce that in this case, the equation @eq:linear-system-polynomial-approximation reduces to

$$
\begin{aligned}
\int_{B} \partial^{\alpha} v \mathrm{d}x
&=
\sum_{|\beta| \leqslant k} c_{\beta} \int_{B} \partial^{\alpha} x^{\beta} \mathrm{d}x 
\\ 
&= 
\sum_{|\beta| = k - 1 } c_{\beta} \int_{B} \partial^{\alpha} x^{\beta} \mathrm{d}x 
+ 
\sum_{|\beta| = k} c_{\beta} \int_{B} \partial^{\alpha} x^{\beta} \mathrm{d}x
\\ 
&= 
 c_{\alpha} \int_{B} \alpha! \, \mathrm{d}x 
+ 
\sum_{|\beta| = k} c_{\beta} \int_{B} \partial^{\alpha} x^{\beta} \mathrm{d}x
\end{aligned}
$$

and we can solve for the coefficient $c_{\alpha}$ as

$$
c_{\alpha} = \frac{1}{\alpha! |B|}
\left(
\int_{B} \partial^{\alpha} v \mathrm{d}x
-
\sum_{|\beta| = k} c_{\beta} \int_{B} \partial^{\alpha} x^{\beta} \mathrm{d}x
\right).
$$

Proceeding recursively, we can determine all coefficients $c_{\beta}$ with $|\beta| = k-2, k-3, \ldots, 0$ in the same way, showing the existence of a polynomial $p \in \mathbb{P}_k$ satisfying @eq:polynomial-approximation-condition.

**Step 2.**
Next, we use a Poincaré inequality of the form 

$$
\| v \|_{k+1,B}^2 \leqslant C
\bigg(
|v|_{k+1,B}^2 + \sum_{|\alpha| \leqslant k} 
\Big|\int_{B} \partial^{\alpha} v \mathrm{d}x\Big|^2
\bigg),
$$ (eq:poincare-inequality-sobolev-norm)

which implies that the rhs of @eq:poincare-inequality-sobolev-norm defines a norm on $H^{k+1}(B)$ equivalent to the standard Sobolev norm $\| \cdot \|_{k+1,B}$.
For a given $v$, we combine this with 
the assumption that $q$ vanishes on $\mathbb{P}_k$ and the polynomial $p$ constructed for $v$ in  Step 1 to deduce that 

$$
\begin{aligned}
|q(v)| &= |q(v - p + p)| \leqslant |q(v - p)| + |q(p)| 
\\
&= |q(v - p)| 
\\
&\leqslant C \|v-p\|_{k+1,B} 
\\
&\leqslant \widetilde{C} \bigg( |v-p|_{k+1,B}^2 + \sum_{|\alpha| \leqslant k} 
\Big|\underbrace{\int_{B} \partial^{\alpha} (v-p) \mathrm{d}x}_{= 0}\Big|^2\bigg)^{1/2}
= \widetilde{C} |v-p|_{k+1,B},
\end{aligned}
$$
where in the last step we used $p$ is a polynomial of degree $k$ and thus $\partial^{\alpha} p = 0$ for all $|\alpha| = k+1$.

:::

:::{exercise} Equivalent norms on Sobolev spaces
Based on the Poincaré inequality from {prf:ref}`lem:poincare-mean-value`,

$$
\|u\|_{0,\Omega}^2 \leqslant C \big( |u|_{1,\Omega}^2 + |\int_{\Omega} u \mathrm{d}x|^2 \big)
$$

show that @eq:poincare-inequality-sobolev-norm holds by an induction argument on $k$.

:::




:::{exercise} Exercise: Transformation formula for Sobolev seminorms for a simple affine mapping

Let $\widehat{T}$ be a non-degenerate simplex in $\mathbb{R}^d$ and set
$ T = h \widehat{T} 
= {\{ x \in \mathbb{R}^d : x = h \hat{x}, \hat{x} \in \widehat{T} \}}$ for some $h>0$.
Show that in this case for any multi-index $\alpha$ with $|\alpha| = l$, we have that
$$
\int_{\widehat{T}} |\partial^{\alpha}_{\widehat{x}} \hat{v}(\widehat{x})|^2 \mathrm{d}\widehat{x} = h^{2l-d} \int_{T} |\partial^{\alpha}_{x} v(x)|^2 \mathrm{d}x
$$
and consequently that
$$
|\hat{v}|_{l,\hat{T}} = h^{l-d/2
} |v|_{l,T}.
$$
:::

:::{hint}
Consider the affine mapping $F_T: \hat{T} \to T$ with $F_T(\hat{x}) = h \hat{x}$ and use the chain rule and change of variables formula for integrals.
:::

For a general affine mapping $F_T: \hat{T} \to T$ with $F_T(\hat{x}) = B_T\hat{x} + b_T$, we 
can generalize the previous exercise to the following lemma, which is a key ingredient in the derivation of error estimates for finite element interpolation.


:::{prf:lemma} Transformation formula for Sobolev seminorms
Let $F_T: \hat{T} \to T$ with $F_T(\hat{x}) = B_T\hat{x} + b_T$ be an affine mapping from a simplex element $\hat{T}$ to another simplex element $T$.
Assume that $v \in H^k(T)$. Set $\hat{v} := v \circ F_T \in H^k(\hat{T})$. 
Then $v \in H^k(T)$ if and only if $\hat{v} \in H^k(\hat{T})$,
and there exists a constant $C>0$ independent of $T$ such that $0 \leqslant l \leqslant k$,
the following inequalities hold:

$$
\begin{aligned}
|\hat{v}|_{l,\hat{T}} &\leqslant C \|B_T\|^l |\det B_T|^{-1/2} |v|_{l,T}, 
\\
|v|_{l,T} &\leqslant C \|B_T^{-1}\|^k |\det B_T|^{1/2} |\hat{v}|_{l,\hat{T}}.
\end{aligned}
$$ (eq:transformation-formula-sobolev-seminorms)

:::

:::{prf:proof}
We only sketch the proof.
Let $F_T(\widehat{x}) = B_T \widehat{x} + b_T$ be an affine mapping from a non-degenerate simplex $\hat{T}$ to another simplex $T$.
Set $\widehat{v}(\widehat{x}) := v(F_T(\widehat{x}))$, we also 
write $x = x(\widehat{x}) = F_T(\widehat{x})$ and $\widehat{x} = \widehat{x}(x) = F_T^{-1}(x)$ to denote the change of variables.
Now the chain rule gives

$$
\begin{aligned}
\partial_{\widehat{x}_i} \widehat{v}(\widehat{x}) 
&= \partial_{\widehat{x}_i} v(F_T(\widehat{x})) 
\\
&= 
\sum_{k=1}^d \partial_{x_k} v(x) \cdot \partial_{\widehat{x}_i} x_k(\widehat{x})
= \sum_{k=1}^d \partial_{x_k} v(x) \cdot (B_T)_{ki}
\end{aligned}
$$ (eq:chain-rule-first-derivative)

For second order derivatives, we have that

$$
\begin{aligned}
\partial_{\widehat{x}_j} \partial_{\widehat{x}_i} \hat{v}(\widehat{x}) 
&= \sum_{k=1}^d \partial_{\widehat{x}_j} (\partial_{x_k} v(x(\widehat{x}))) \cdot (B_T)_{ki})
\\
&= \sum_{k=1}^d \sum_{l=1}^d 
 \partial_{x_l} \partial_{x_k} 
 v(x)
 \partial_{\widehat{x}_j} x_l(\widehat{x})
 \cdot (B_T)_{ki}
\\
&= \sum_{k=1}^d \sum_{l=1}^d 
 \partial_{x_l} \partial_{x_k} 
 v(x)
 \cdot
 (B_T)_{lj}
 \cdot (B_T)_{ki},
\end{aligned}
$$ (eq:chain-rule-second-derivative)

where we used the fact that $B_T$ is constant, and thus $\partial_{\widehat{x}_j} (B_T)_{ki} = 0$.
For higher order derivatives, we can proceed similarly.

Based on @eq:chain-rule-first-derivative, we see that any first order partial derivative of $\hat{v}$ can be bounded as follows:

$$
\begin{aligned}
|\partial_{\widehat{x}_i} \widehat{v}(\widehat{x})|
&\leqslant 
d 
\max_{1 \leqslant k \leqslant d} |\partial_{x_k} v(x)|
\cdot 
\underbrace{\max_{1 \leqslant i,k \leqslant d} |(B_T)_{ki}| }_{=: \| B_T \|_{\max}}
\\
&\leqslant 
C(d)
\max_{1 \leqslant k \leqslant d} |\partial_{x_k} v(x)|
\| B_T \|_2
\end{aligned}
$$

where we used the fact that the maximum norm of a matrix is bounded by its spectral norm, up to a constant depending on the dimension $d$, which is denoted by $C(d)$.

Similarly, based on @eq:chain-rule-second-derivative,
we infer that any second order partial derivative of $\hat{v}$ can be bounded as follows:

$$
\begin{aligned}
|\partial_{\widehat{x}_j} \partial_{\widehat{x}_i} \hat{v}(\widehat{x})|
&\leqslant
C(d)^2
\max_{1 \leqslant k,l \leqslant d} |\partial_{x_l} \partial_{x_k} v(x)|
\| B_T \|_2^2
\end{aligned}
$$

By induction, one can show that for any $l$-th order partial derivative of $\partial^{\alpha}_{\widehat{x}} \hat{v}$ with $|\alpha| = l$, the inequality

$$
\begin{aligned}
|\partial^{\alpha}_{\widehat{x}} \hat{v}(\widehat{x})|
&\leqslant
C(d)^l
\max_{|\beta| = l} |\partial^{\beta}_{x} v(x)|
\| B_T \|_2^{l}
\end{aligned}
$$ (eq:chain-rule-l-th-derivative-est)

must hold.

Thus for any multi-index $\alpha$ with $|\alpha| = l$, we estimate

$$
\begin{aligned}
\int_{\widehat{T}} |\partial^{\alpha}_{\widehat{x}} \hat{v}(\widehat{x})|^2 \mathrm{d}\widehat{x}
&\leqslant
C(d)^{2l} \| B_T \|_2^{2l}
\int_{\widehat{T}} \max_{|\beta| = l} |\partial^{\beta}_{x} v(x(\widehat{x}))|^2 \mathrm{d}\widehat{x}
\\
&\leqslant
C(d)^{2l} \| B_T \|_2^{2l}
\int_{\widehat{T}} \sum_{|\beta| = l} |\partial^{\beta}_{x} v(x(\widehat{x}))|^2 \mathrm{d}\widehat{x}
\\
&=
C(d)^{2l} \| B_T \|_2^{2l}
\int_{{T}} \sum_{|\beta| = l} |\partial^{\beta}_{x} v(x)|^2 
|\det DF_T^{-1}(x)| \mathrm{d} x
\\
&=
C(d)^{2l} \| B_T \|_2^{2l} |\det B_T|^{-1}
\int_{{T}} \sum_{|\beta| = l} |\partial^{\beta}_{x} v(x)|^2 
\\
&=
C(d)^{2l} \| B_T \|_2^{2l} |\det B_T|^{-1}
||v|_{l,T}^2
\end{aligned}
$$

Summing over all multi-indices $\alpha$ with $|\alpha| = l$ to compute the complete Sobolev seminorm $|\hat{v}|_{l,\hat{T}}$
and taking the square root,
 we obtain the first inequality in @eq:transformation-formula-sobolev-seminorms.
 The second inequality can be derived similarly by considering the inverse mapping $F_T^{-1}: T \to \hat{T}$.

:::

Next, we wish to estimate the norm of the matrix $B_T$ 
and $B_T^{-1}$. For that, we recall the notion of the diameter $h_T = \sup_{x,y \in T} \|x-y\|$
of an element $T$ and denote by $\rho_T$ of the inscribed ball of the element $T$.

```{figure} figures/element-diameter-inradius.svg
:label: fig-element-diameter-inradius
:alt: Left, a triangle with its longest edge highlighted and labelled h_T, and its inscribed circle shaded with a horizontal diameter labelled rho_T. Right, a tetrahedron in a three-dimensional view with its longest edge highlighted and labelled h_T, and its inscribed sphere shaded with a diameter labelled rho_T.
:width: 80%

The diameter $h_T$ (magenta) and the diameter $\rho_T$ of the inscribed ball (cyan) of a triangle
and of a tetrahedron. For a simplex, $h_T$ is the length of the longest edge, and
$\rho_T$ is the diameter of the inscribed circle or sphere.
```

:::{prf:lemma} Norm of the affine mapping
Let $F_T(\widehat{x}) = B_T \widehat{x} + b_T$ be an affine mapping from a non-degenerate simplex $\hat{T}$ to another simplex $T$.
Then

$$
\|B_T\|_2 \leqslant \dfrac{h_T}{\rho_{\hat{T}}}\
\quad \text{and} \quad
\|B_T^{-1}\|_2 \leqslant \dfrac{h_{\hat{T}}}{\rho_T}.
$$

:::

:::{prf:proof}

We only show the first inequality, the second one can be derived in the exact same way by considering the inverse mapping $F_T^{-1}: T \to \hat{T}$.
Recall that the norm of a matrix $B_T$ is defined as

$$
\| B_T \|_2 = \sup_{\widehat{\boldsymbol{e}} \in \mathbb{R}^d, \|\widehat{\boldsymbol{e}}\|_2 = 1} \| B_T \boldsymbol{e} \|_2.
$$

For any unit vector $\widehat{\boldsymbol{e}} \in \mathbb{R}^d$, we   can find
points $\widehat{x}$ and $\widehat{y}$ on the inscribed
circle of $\widehat{T}$ which satisfy $\widehat{x} - \widehat{y} = \rho_T \,
\widehat{\boldsymbol{e}}$,
see @fig-affine-map-inscribed-ball.
Their images
$F_T(\widehat{x})$ and $F_T(\widehat{y})$ lie in $T$, and therefore

$$
\| B_T \rho_T \widehat{\boldsymbol{e}} \| = \| F_T(\widehat{x}) - F_T(\widehat{y}) \| \leqslant h_T. 
$$
or in other words
$
\| B_T \widehat{\boldsymbol{e}} \| \leqslant \frac{h_T}{\rho_{\hat{T}}} 
$ for any unit vector $\widehat{\boldsymbol{e}} \in \mathbb{R}^d$.
:::

```{figure} figures/affine-map-inscribed-ball.svg
:label: fig-affine-map-inscribed-ball
:alt: Left, the reference triangle with its inscribed circle and two diametrically opposite points x-hat and y-hat on the circle, joined by a segment labelled rho e-hat; its longest edge is labelled h_T-hat. An arrow labelled F_T points to the right, where the image triangle T shows the points F_T(x-hat) and F_T(y-hat) joined by the image segment, the image of the inscribed circle as a dashed ellipse, the inscribed circle of T with diameter rho_T, and the longest edge h_T.
:width: 90%

The affine map $F_T$ from the reference triangle $\widehat{T}$ to $T$. For a unit vector
$\widehat{\boldsymbol{e}}$, the points $\widehat{x}$ and $\widehat{y}$ on the inscribed
circle of $\widehat{T}$ (cyan) satisfy $\widehat{x} - \widehat{y} = \rho \,
\widehat{\boldsymbol{e}}$ with $\rho = \rho_{\widehat{T}}$ (orange). Their images
$F_T(\widehat{x})$ and $F_T(\widehat{y})$ lie in $T$, and
$F_T(\widehat{x}) - F_T(\widehat{y}) = \rho \, B_T \widehat{\boldsymbol{e}}$, so
$\rho \, \| B_T \widehat{\boldsymbol{e}} \| \leqslant h_T$. The dashed ellipse is the
image of the inscribed circle of $\widehat{T}$. Magenta: the longest edges, of length
$h_{\widehat{T}}$ and $h_T$. Cyan: the inscribed circles; $\rho_T$ marks the diameter of
the one of $T$.
```

The previous lemma motivates the
notion of **shape-regularity** of a family of meshes, which is a key assumption in the derivation of error estimates for finite element interpolation.

Let us consider a series of affine meshes $\{\mathcal{T}_h\}_{h>0}$ of a domain $\Omega \subset \mathbb{R}^d$ with $d=1,2,3$.
Such a family of meshes is said to be **shape-regular** if there exists a constant $\sigma>0$ such that for all $h>0$ and all $T \in \mathcal{T}_h$, we have 

$$
\sigma_T := \frac{h_T}{\rho_T} \leqslant \sigma, \forall T \in \mathcal{T}_h, \forall h>0,
$$

where $h_T$ is the diameter of the element $T$ and $\rho_T$ is the diameter of the largest ball contained in $T$.

```{figure} figures/shape-regularity.svg
:label: fig-shape-regularity
:alt: Three triangles with their longest edge highlighted and their inscribed circle shaded. An equilateral triangle with sigma_T about 1.7, a long thin needle-shaped triangle with sigma_T about 8, and a flat triangle with one angle close to 180 degrees with sigma_T about 10.
:width: 90%

Shape regularity. Left: an equilateral triangle, for which $\sigma_T = h_T/\rho_T = \sqrt{3}$
is as small as possible. Middle: a very acute triangle, with one angle of about $7^\circ$.
Right: a very obtuse triangle, with one angle of about $156^\circ$. In both degenerate
cases the inscribed circle is small compared with $h_T$, so $\sigma_T$ is large. A
shape-regular family of meshes must avoid both, as $\sigma_T$ grows without bound when the
small angle tends to $0^\circ$ or the large angle to $180^\circ$.
```