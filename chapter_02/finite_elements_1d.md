---
jupytext:
  cell_metadata_filter: -jupyter
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.5
kernelspec:
  display_name: tma4220-finite-elements (3.14.7.final.0)
  language: python
  name: python3
---

# The finite element method

Finite element methods are specific realizations of
the Galerkin method.
In a nutshell, the basis functions of the discrete
function spaces are constructed from in three steps:

1. Generation of mesh $\mathcal{T}_h = \{T\}$ which is a partition of the domain $\Omega$ into smaller subunits $T$, the mesh elements.
2. On each mesh element, a (local) finite dimensional function space $V_T$ is constructed which has a fixed dimension.
3. The global discrete function space $V_h$ is patched together suitably from the local, 
element-wise define function spaces. How the local function spaces 
are patched together will determine the resulting conformity of the
global function space $V_h$, e.g. whether $V_h \subset L^2(\Omega)$
or $V_h \subset H^1(\Omega)$. 

To illustrate some of the core concepts when realizing a finite
element discretization, we start with the simplest 1D example, before we
formalize and extend those to cover the plethora of finite element
families which exist out there today.

(ssec:lagrange-interpolation)=
## A quick recap of Lagrange interpolation

Before we construct finite element spaces, we recall how a polynomial on a
single interval is determined by its values at a set of nodes. Throughout
this section, $T = [a, b]$ is a fixed interval of length $h = b - a$, and
$k \geqslant 1$ denotes the polynomial degree.

```{prf:definition} Lagrange basis polynomials
:label: def:lagrange-basis-1d

Let $a = \xi_0 < \xi_1 < \ldots < \xi_k = b$ be $k+1$ distinct
**interpolation nodes** in $T$. Unless stated otherwise, we use
equispaced nodes $\xi_j = a + j h/k$. The **Lagrange basis polynomials**
$\{\phi_j\}_{j=0}^k$ associated with these nodes are

$$
\phi_j(x) = \prod_{\substack{m=0 \\ m \neq j}}^{k} \frac{x - \xi_m}{\xi_j - \xi_m},
\qquad j = 0, 1, \ldots, k.
$$ (eq:lagrange-basis-poly)
```

Each $\phi_j$ is a product of $k$ linear factors, so $\phi_j \in \mathbb{P}_k(T)$.
The numerator vanishes at every node except $\xi_j$, where the fraction
equals $1$. Hence

$$
\phi_j(\xi_m) = \delta_{jm} \quad \text{for } j, m = 0, 1, \ldots, k.
$$ (eq:lagrange-delta)

This property makes $\{\phi_j\}_{j=0}^k$ a basis of $\mathbb{P}_k(T)$. If
$\sum_{j=0}^k c_j \phi_j = 0$, then evaluating at $\xi_m$ gives $c_m = 0$
for every $m$. The $k+1$ polynomials are therefore linearly independent,
and $\dim \mathbb{P}_k(T) = k+1$.

The widget below shows the Lagrange basis polynomials on $T = [0, 1]$ for
equispaced nodes. Each $\phi_j$ equals $1$ at its own node and $0$ at all
other nodes, as stated in [](#eq:lagrange-delta). Expand the code cells to
see how the plots are generated.

```{code-cell} ipython3
:tags: [hide-input]

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np


def lagrange_basis(xis, j):
    """Return the Lagrange basis polynomial phi_j for the nodes xis."""
    def phi_j(x):
        x = np.asarray(x, dtype=float)
        result = np.ones_like(x)
        for m, xi_m in enumerate(xis):
            if m != j:
                result *= (x - xi_m) / (xis[j] - xi_m)
        return result
    return phi_j


def generate_lagrange_polynomial(a, b, order):
    num_points = order + 1
    xis = np.linspace(a, b, num_points)
    lagrange_poly_fs = [np.eye(num_points)[i] for i in range(num_points)]
    lagrange_polys = [lagrange_basis(xis, i) for i in range(num_points)]
    return xis, lagrange_poly_fs, lagrange_polys
```

```{code-cell} ipython3
:tags: [hide-input]

def plot_lagrange_basis(order):
    a, b = 0, 1
    xs = np.linspace(a, b, 200)
    xis, fs, polys = generate_lagrange_polynomial(a, b, order)

    fig, ax = plt.subplots(figsize=(6, 4))
    for i, p in enumerate(polys):
        line, = ax.plot(xs, p(xs), label=fr"$\phi_{{{i}}}$")
        ax.scatter(xis, fs[i], color=line.get_color())

    # Draw the interval as a thick black line, with ticks at its endpoints
    ax.plot([a, b], [0, 0], 'k-', linewidth=2, zorder=1)
    ax.set_xticks([a, b])
    ax.set_xticklabels(["$a$", "$b$"])
    ax.tick_params(axis='x', length=0)

    ax.set_title(f"Lagrange basis functions on $T=[{a}, {b}]$, degree {order}")
    ax.set_xlabel("x")
    ax.legend()
    ax.grid(True)
    plt.show()

widgets.interact(
    plot_lagrange_basis,
    order=widgets.IntSlider(
        value=2,
        min=1,
        max=5,
        step=1,
        description='Degree:',
        continuous_update=False
    )
);
```

```{prf:definition} Lagrange interpolation operator
:label: def:lagrange-interpolant-1d

Let $\{\phi_j\}_{j=0}^k$ be the Lagrange basis polynomials for the nodes
$\{\xi_j\}_{j=0}^k$ in $T$. The **Lagrange interpolation operator**
$I_k : C(T) \to \mathbb{P}_k(T)$ is defined by

$$
I_k f = \sum_{j=0}^k f(\xi_j)\, \phi_j .
$$ (eq:lagrange-interpolant)

We call $I_k f$ the **Lagrange interpolant** of $f$.
```

By [](#eq:lagrange-delta), we have
$(I_k f)(\xi_m) = \sum_{j=0}^k f(\xi_j)\, \phi_j(\xi_m) = f(\xi_m)$,
so $I_k f$ agrees with $f$ at every node. It is the only polynomial in
$\mathbb{P}_k(T)$ with this property. If $p, q \in \mathbb{P}_k(T)$ both agree with $f$ at
the $k+1$ nodes, then $p - q \in \mathbb{P}_k(T)$ has $k+1$ distinct roots and must
vanish identically. Clearly, $I_k$ is also linear.

```{prf:remark} $I_k$ is a projection
:label: rem:lagrange-projection

Let $p \in \mathbb{P}_k(T)$. Both $p$ and $I_k p$ are polynomials of degree at most
$k$ which agree with $p$ at the nodes, so by uniqueness $I_k p = p$.
Consequently, $I_k (I_k f) = I_k f$ for all $f \in C(T)$, and $I_k$ is a
linear projection of $C(T)$ onto $\mathbb{P}_k(T)$.
```

```{prf:remark} Nodal values as functionals
:label: rem:nodal-functionals

The nodal evaluations can be viewed as linear functionals

$$
\sigma_j : C(T) \ni f \mapsto \sigma_j(f) = f(\xi_j) \in \mathbb{R},
\qquad j = 0, 1, \ldots, k.
$$

In this language, [](#eq:lagrange-delta) reads $\sigma_m(\phi_j) = \delta_{jm}$.
The Lagrange basis $\{\phi_j\}_{j=0}^k$ is thus the basis of $\mathbb{P}_k(T)$ that
is *dual* to the functionals $\{\sigma_j\}_{j=0}^k$, and the interpolant takes
the form $I_k f = \sum_{j=0}^k \sigma_j(f)\, \phi_j$. We will use exactly
this viewpoint later to arrive at a very general definition of finite elements: a local polynomial space, together
with a set of functionals which uniquely determine its elements.
```

```{prf:theorem} Interpolation error
:label: thm:lagrange-interpolation-error

Let $f \in C^{k+1}(T)$. Then for every $x \in T$ there is an
$\eta = \eta(x) \in (a, b)$ such that

$$
f(x) - (I_k f)(x) = \frac{f^{(k+1)}(\eta)}{(k+1)!} \prod_{j=0}^k (x - \xi_j).
$$ (eq:lagrange-interpolation-error)

In particular, since $|x - \xi_j| \leqslant h$ for all $x \in T$,

$$
\| f - I_k f \|_{L^\infty(T)}
\leqslant \frac{h^{k+1}}{(k+1)!} \, \| f^{(k+1)} \|_{L^\infty(T)}.
$$ (eq:lagrange-interpolation-bound)
```

This is a standard result from numerical analysis. Its proof applies
Rolle's theorem $k+1$ times to a suitable auxiliary function.

The bound [](#eq:lagrange-interpolation-bound) suggests two ways to make the
interpolation error small: increase the degree $k$, or shrink the interval
length $h$. The widget below lets you explore the first option on
$T = [0, 1]$. Choose a function and vary the degree. For the polynomial
functions, observe that the error vanishes (up to round-off) once $k$
reaches the degree of the polynomial, as predicted by
[](#rem:lagrange-projection).

```{code-cell} ipython3
:tags: [hide-input]

def f_1(x):
    return x*np.sin(3*np.pi*x)

def f_poly1(x):
    return 2 - 10*x

def f_poly2(x):
    return x*(1 - x)

def f_poly3(x):
    return x*(1 - x)*(x - 0.5)

def f_runge(x):
    return 1/(1 + 25*(2*x - 1)**2)

functions = {
    "x sin(3πx)": f_1,
    "Polynomial, degree 1: 2 - 10x": f_poly1,
    "Polynomial, degree 2: x(1-x)": f_poly2,
    "Polynomial, degree 3: x(1-x)(x-0.5)": f_poly3,
    "Runge function: 1/(1 + 25(2x-1)²)": f_runge,
}


def plot_lagrange_interpolation(order, function_name):
    a, b = 0, 1
    f = functions[function_name]
    xs = np.linspace(a, b, 400)
    xis, fs, polys = generate_lagrange_polynomial(a, b, order)

    # I_k f = sum_j f(xi_j) phi_j
    interpolant = sum(f(xi) * phi(xs) for xi, phi in zip(xis, polys))
    max_error = np.max(np.abs(f(xs) - interpolant))

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(xs, f(xs), 'k-', linewidth=1.5, alpha=0.6, label="$f$")
    ax.plot(xs, interpolant, color='C1', linewidth=2, label=f"$I_{{{order}}} f$")
    ax.scatter(xis, f(xis), color='C1', zorder=3, label="nodes")

    ax.set_xticks([a, b])
    ax.set_xticklabels(["$a$", "$b$"])

    ax.set_title(f"Lagrange interpolation of degree {order}\n"
                 f"max error ≈ {max_error:.2e}")
    ax.set_xlabel("x")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.show()


widgets.interact(
    plot_lagrange_interpolation,
    order=widgets.IntSlider(
        value=2, min=1, max=12, step=1,
        description='Degree:',
        continuous_update=False
    ),
    function_name=widgets.Dropdown(
        options=list(functions.keys()),
        value=list(functions.keys())[0],
        description='Function:'
    )
);
```

```{prf:remark} Runge's phenomenon
:label: rem:runge-phenomenon

Select the Runge function in the widget above and increase the degree.
Although this function is infinitely differentiable, the interpolation error
for equispaced nodes does not decrease. Beyond moderate degrees it *grows*
with $k$, and large oscillations appear near the endpoints of the interval.
This does not contradict [](#thm:lagrange-interpolation-error): the
derivatives $f^{(k+1)}$ of the Runge function grow so quickly with $k$ that
the right-hand side of [](#eq:lagrange-interpolation-bound) does not tend to
zero.

Raising the polynomial degree on a fixed interval is therefore not a reliable
way to approximate a function. One remedy is to cluster the nodes towards the
endpoints, e.g. by using Chebyshev nodes. The remedy central to this course is
the second option offered by [](#eq:lagrange-interpolation-bound): keep the
degree $k$ fixed and small, divide the domain into many small elements, and
interpolate on each of them, so that $h$ becomes small. This is precisely the
idea behind the finite element spaces constructed next.
```

(ssec:finite-element-1d-method)=
## Finite element in 1D
Consider a the finite interval $\Omega = (a, b)$. To create a mesh, we introduce $N+2$
**nodes** or **vertices** $\{x_i\}_{i=0}^{N+1}$ such that
$a = x_0 < x_1 < \ldots < x_N < x_{N+1} = b$. 
Then the corresponding **mesh** $\mathcal{T}_h = \{T_i\}_{i=1}^{N+1}$
is simply the finite collection of the $N+1$ subelements 
$T_i = [x_{i-1}, x_i]$ which are the **mesh elements**.
We set mesh size $h_T$
\begin{align}
\text{Local mesh size}\; &h_T = \mathrm{diam}(T) \text{ for } T \in \mathcal{T}_h 
\\
\text{Global mesh size}\;  &h = \max_{T \in \mathcal{T}_h} h_T
\end{align}
As usual, for a set $U \subset \mathbb{R}^d$, we set $\mathrm{diam}(U) := \sup \{\| \mathbf{x} - \mathbf{y} \| : \mathbf{x}, \mathbf{y} \in U  \}$, 
which for our particular case means that
$\mathrm{diam}(T_i) = |x_i - x_{i-1}|$.

### Continuous piecewise linear elements
Next, we construct continuous, piecewise linear finite elements defined on $\mathcal{T}_h$

$$
\mathbb{P}^{\mathrm{c}}_1(\mathcal{T}_h)= \{ v \in C(\Omega) : v|_T \in \mathbb{P}_1(T)\; \forall T \in \mathcal{T}_h \}
$$

In the widget below, you can explore typical functions
$u_h \in \mathbb{P}^{\mathrm{c}}_1(\mathcal{T}_h)$ which interpolate a given
function on a uniform mesh of $[0, 1]$. On each element, the function is interpolated with
the linear Lagrange basis from the previous section. Choose a function and vary
the number of elements. For the Runge function, compare the result with the
high-degree interpolation on a single interval above.

```{code-cell} ipython3
:tags: [hide-input]

def plot_fe_interpolation_p1(num_elements, function_name):
    a, b = 0, 1
    f = functions[function_name]
    mesh = np.linspace(a, b, num_elements + 1)
    x_fine = np.linspace(a, b, 400)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x_fine, f(x_fine), 'k-', linewidth=1.5, alpha=0.6, label="$f$")

    max_error = 0.0
    for e in range(num_elements):
        # Interpolate f on element [x_e, x_{e+1}] with the linear Lagrange basis
        xis, fs, polys = generate_lagrange_polynomial(mesh[e], mesh[e + 1], 1)
        x_local = np.linspace(mesh[e], mesh[e + 1], 100)
        interpolant = sum(f(xi) * phi(x_local) for xi, phi in zip(xis, polys))
        max_error = max(max_error, np.max(np.abs(f(x_local) - interpolant)))
        ax.plot(x_local, interpolant, color='C1', linewidth=2,
                label="FE interpolant" if e == 0 else None)
        ax.scatter(xis, f(xis), color='C1', zorder=3)

    ax.axhline(0, color='black', linewidth=1, zorder=0)
    for x in mesh:
        ax.axvline(x, color='gray', linestyle='--', linewidth=0.5, zorder=0)
    # Label the nodes only while the labels still fit
    if num_elements <= 10:
        ax.set_xticks(mesh)
        ax.set_xticklabels([f"$x_{{{k}}}$" for k in range(len(mesh))])
    else:
        ax.set_xticks([a, b])
        ax.set_xticklabels(["$a$", "$b$"])

    ax.set_title(f"Piecewise linear interpolation of {function_name}\n"
                 f"({num_elements} elements, max error ≈ {max_error:.2e})")
    ax.set_xlabel("x")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.show()


widgets.interact(
    plot_fe_interpolation_p1,
    num_elements=widgets.IntSlider(
        value=4, min=1, max=20, step=1,
        description='#Elements:',
        continuous_update=False
    ),
    function_name=widgets.Dropdown(
        options=list(functions.keys()),
        value=list(functions.keys())[0],
        description='Function:'
    )
);
```

This global discrete function space can be constructed by defining suitably chosen basis functions
for the space $\mathbb{P}_1(T_i)$
on each mesh element $T_i$ and then patching them together to form a global basis for $V_h$.
Let's focus on some element $ T = T_i = [x_{i-1}, x_i]$ for the moment and
introduce the interpolation nodes $\xi_0 = x_{i-1}$ and $\xi_1 = x_i$, 
then our beloved Lagrange basis functions $\{\phi_j^T\}_{j=0}^1$ are defined as
the unique linear functions which satisfy

$$
\phi_j^T(\xi_k) = \delta_{jk} \quad \text{for } j, k = 0, 1.
$$ (lagrange_basis_1d)

Here, we already note that the function evaluation $\phi_j^T(\xi_k)$
can also be interpreted as evaluation of certain **functionals**
$\{\sigma_k \}_{k=0}^1$ if we define those functionals as follows

$$
\sigma_k: C(T_i)  \ni f \mapsto \sigma_k(f) = f(\xi_k) \in \mathbb{R}
$$

If we repeat this construction on each element $T_i$, we see that that elementwise defined
basis functions are patched together to global, so-called **hat functions** $\{\phi_i\}_{i=0}^{N+1}$ which satisfy $\phi_j(x_i) = \delta_{ji}$ for a given mesh node $x_i$ and which is linear on each element.

The widget below shows the hat functions on a uniform mesh of $[0, 1]$. Each
color marks one global basis function $\phi_i$. It is glued together from the
local linear basis functions on the (at most two) elements which share the
node $x_i$. Vary the number of elements and observe how the hat functions
become narrower as the mesh is refined.

```{code-cell} ipython3
:tags: [hide-input]

def plot_hat_functions(num_elements):
    a, b = 0, 1
    mesh = np.linspace(a, b, num_elements + 1)
    cmap = plt.get_cmap('tab10')

    fig, ax = plt.subplots(figsize=(8, 5))
    for e in range(num_elements):
        xis, fs, polys = generate_lagrange_polynomial(mesh[e], mesh[e + 1], 1)
        xs_local = np.linspace(mesh[e], mesh[e + 1], 200)
        for j, p in enumerate(polys):
            # Local basis function j on element e belongs to global node e + j
            color = cmap((e + j) % 10)
            ax.plot(xs_local, p(xs_local), color=color)
            ax.scatter(xis[j], 1, color=color, zorder=3)

    for x in mesh:
        ax.axvline(x, color='gray', linestyle='--', linewidth=0.5, zorder=1)

    # Draw the elements as a black line, with node markers on top
    ax.plot(mesh, np.zeros_like(mesh), 'k-|', linewidth=2, markersize=16, markeredgewidth=1.5, zorder=4)
    ax.set_xticks(mesh)
    ax.set_xticklabels([f"$x_{{{k}}}$" for k in range(len(mesh))])
    ax.tick_params(axis='x', length=0)

    ax.set_title(f"Hat functions on a mesh with {num_elements} elements")
    ax.set_xlabel("x")
    ax.grid(True, alpha=0.3)
    plt.show()

widgets.interact(
    plot_hat_functions,
    num_elements=widgets.IntSlider(
        value=4, min=1, max=10, step=1,
        description='#Elements:',
        continuous_update=False
    )
);
```

### Continuous piecewise quadratic elements
We can repeat the entire procedure to construct piecewise *quadratic* elements.
This time, we simply need to construct for each $T \in \mathcal{T}_h$ the three Lagrange basis polynomials $\{\phi^{T_i}_j\}_{j=0}^2$ which form a basis for the space $\mathbb{P}_2(T)$. Recall that after introducing the local interpolation nodes $\xi_0 = x_{i-1}, \xi_{1} = x_{i-1} + h_i/2$, and $\xi_{2} = x_{i}$, the local Lagrange basis polynomials
are again characterized by the property

$$
\sigma_{k}(\phi_j^T) :=  \phi_j^T(\xi_k) = \delta_{jk} \quad \text{for } j, k = 0, 1, 2,
$$

Again, the three functionals $\{\sigma_k\}_{k=0}^2$ are defined as functionals
which evaluates a given function in on of the respective interpolation points $\xi_k$.

Repeating the construction on each $T_i$ and patching those local basis functions together, we get elementwise quadratic functions of the following form:

The widget below shows the global basis functions for continuous piecewise
linear, quadratic, and cubic elements on a uniform mesh of $[0, 1]$. Each color
marks one global basis function. A basis function associated with a mesh node
$x_i$ is glued together from the two elements which share that node, while a
basis function associated with a node in the interior of an element vanishes
outside that element. Choose the degree and vary the number of elements.

```{code-cell} ipython3
:tags: [hide-input]

degree_names = {1: "linear", 2: "quadratic", 3: "cubic"}


def plot_global_basis(num_elements, order):
    a, b = 0, 1
    mesh = np.linspace(a, b, num_elements + 1)
    cmap = plt.get_cmap('tab10')

    fig, ax = plt.subplots(figsize=(8, 5))
    for e in range(num_elements):
        xis, fs, polys = generate_lagrange_polynomial(mesh[e], mesh[e + 1], order)
        xs_local = np.linspace(mesh[e], mesh[e + 1], 200)
        for j, p in enumerate(polys):
            # Local basis function j on element e belongs to global node e*order + j,
            # so neighbouring elements share the node at their common endpoint
            color = cmap((e * order + j) % 10)
            ax.plot(xs_local, p(xs_local), color=color)
            ax.scatter(xis[j], 1, color=color, zorder=3)

    for x in mesh:
        ax.axvline(x, color='gray', linestyle='--', linewidth=0.5, zorder=1)

    # Draw the elements as a black line, with node markers on top
    ax.plot(mesh, np.zeros_like(mesh), 'k-|', linewidth=2, markersize=16, markeredgewidth=1.5, zorder=4)
    ax.set_xticks(mesh)
    ax.set_xticklabels([f"$x_{{{k}}}$" for k in range(len(mesh))])
    ax.tick_params(axis='x', length=0)

    ax.set_title(f"Global basis functions: {num_elements} elements, "
                 f"piecewise {degree_names[order]}")
    ax.set_xlabel("x")
    ax.grid(True, alpha=0.3)
    plt.show()

widgets.interact(
    plot_global_basis,
    num_elements=widgets.IntSlider(
        value=4, min=1, max=10, step=1,
        description='#Elements:',
        continuous_update=False
    ),
    order=widgets.Dropdown(
        options=[(name, k) for k, name in degree_names.items()],
        value=2,
        description='Degree:'
    )
);
```

As for piecewise linear elements, the widget below shows the continuous
piecewise polynomial interpolant of a given function on a uniform mesh of
$[0, 1]$. Choose the number of elements, the polynomial degree, and the
function. Doubling the number of elements lets you compare how quickly the
error decreases for the different degrees.

```{code-cell} ipython3
:tags: [hide-input]

def plot_fe_interpolation(num_elements, order, function_name):
    a, b = 0, 1
    f = functions[function_name]
    mesh = np.linspace(a, b, num_elements + 1)
    x_fine = np.linspace(a, b, 400)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x_fine, f(x_fine), 'k-', linewidth=1.5, alpha=0.6, label="$f$")

    max_error = 0.0
    for e in range(num_elements):
        # Interpolate f on element [x_e, x_{e+1}] with the local Lagrange basis
        xis, fs, polys = generate_lagrange_polynomial(mesh[e], mesh[e + 1], order)
        x_local = np.linspace(mesh[e], mesh[e + 1], 100)
        interpolant = sum(f(xi) * phi(x_local) for xi, phi in zip(xis, polys))
        max_error = max(max_error, np.max(np.abs(f(x_local) - interpolant)))
        ax.plot(x_local, interpolant, color='C1', linewidth=2,
                label="FE interpolant" if e == 0 else None)
        ax.scatter(xis, f(xis), color='C1', zorder=3)

    ax.axhline(0, color='black', linewidth=1, zorder=0)
    for x in mesh:
        ax.axvline(x, color='gray', linestyle='--', linewidth=0.5, zorder=0)
    # Label the mesh nodes only while the labels still fit
    if num_elements <= 10:
        ax.set_xticks(mesh)
        ax.set_xticklabels([f"$x_{{{k}}}$" for k in range(len(mesh))])
    else:
        ax.set_xticks([a, b])
        ax.set_xticklabels(["$a$", "$b$"])

    ax.set_title(f"Piecewise {degree_names[order]} interpolation of {function_name}\n"
                 f"({num_elements} elements, max error ≈ {max_error:.2e})")
    ax.set_xlabel("x")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.show()


widgets.interact(
    plot_fe_interpolation,
    num_elements=widgets.Dropdown(
        options=[1, 2, 4, 8, 16],
        value=4,
        description='#Elements:'
    ),
    order=widgets.Dropdown(
        options=[(name, k) for k, name in degree_names.items()],
        value=2,
        description='Degree:'
    ),
    function_name=widgets.Dropdown(
        options=list(functions.keys()),
        value=list(functions.keys())[0],
        description='Function:'
    )
);
```
