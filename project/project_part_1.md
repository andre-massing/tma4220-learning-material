---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.5
kernelspec:
  display_name: tma4220-finite-elements (3.14.10)
  language: python
  name: python3
---

# Implementation of a finite element solver, part 1: L²-projection

```{code-cell} ipython3
# Make plotted figures interactive
%matplotlib widget

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import sympy as sy
import meshtools as mt
import plottools as pt
import tracemalloc
```

## Problem 1. Assembly of the mass matrix

### The local mass matrix on the reference triangle

+++

The nodal basis/shape functions for the first-order Lagrange element $\mathbb{P}_1$ on the unit triangle $\widehat{T}$ are

$$
\begin{aligned}
\lambda_0(\widehat{x},\widehat{y}) &= 1 - \widehat{x} -\widehat{y}
\\
\lambda_1(\widehat{x},\widehat{y}) &= \widehat{x}
\\
\lambda_2(\widehat{x},\widehat{y}) &= \widehat{y}
\end{aligned}
$$

+++

The local mass matrix $\widehat{M}$ on the reference triangle $\widehat{T}$
is defined with entries   

$$
\widehat{M}_{\alpha\beta} = \dfrac{1}{|\widehat{T}|}\int_{\widehat{T}} \lambda_\beta(\widehat{x},\widehat{y}) \lambda_\alpha(\widehat{x},\widehat{y}) \, d\widehat{x} d\widehat{y}, \quad \alpha,\beta = 0,1,2.
$$

so that the local mass matrix on an element $T$ can be computed as

$$
M_T = |T| \cdot \widehat{M}
$$

where $|T|$ is the area of the element $T$. Now compute $\widehat{M}$
using symbolic integration with the `sympy` package. 
<!-- **b)** a **suitable** quadrature rule for the unit triangle $\widehat{T}$ from the list given in the lectures. Justify your choice of quadrature rule and make sure you get the same in both cases (up to machine precision, i.e $\leqslant 10^{-13}$). -->

### Assembly of the global mass matrix 

Next, you are asked to implement three different variants of the assembly of the global mass matrix $M$ in Python, and to benchmark their memory usage and execution time. 

The first version should use a dense matrix representation using standard NumPy arrays, the second version should use a sparse matrix representation in DOK format, and the third version should use a sparse matrix representation in COO format.

For all three variants, use the function signature

```python
def assemble_mass_matrix_MATRIX_CLASS(P, C):
    """Assemble the P1 mass matrix on the triangulation (P, C).

    P : (N_p, 2) array, coordinates of the mesh vertices (point matrix).
    C : (N_t, 3) integer array, vertex indices of each triangle (connectivity matrix).

    Returns the (N_p, N_p) mass matrix M, stored as MATRIX_CLASS.
    """
    ...
    return M
```

where `MATRIX_CLASS` is either `dense`, `dok`, or `coo`.

To measure the memory usage, you can use the following function as a starting point:

```{code-cell} ipython3
def memit(f, *args, **kwargs):
    """Call f(*args, **kwargs) and print the memory it allocated.

    Uses tracemalloc to track the memory allocated by Python (including NumPy
    and SciPy arrays) during the call and prints two numbers in MB:

    kept : memory still allocated when f returns (essentially the result).
    peak : maximum memory allocated at any time during the call,
           including temporary arrays that were freed again.

    Returns the return value of f, so that memit(f, ...) can replace f(...).
    """
    tracemalloc.start()
    result = f(*args, **kwargs)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"kept: {current/1e6:.3f} MB, peak: {peak/1e6:.3f} MB")
    return result
```

To print the memory usage of the assembly functions, you can simply type

```python
M = memit(assemble_mass_matrix_dense, P, C)
```

Feel free to modify the function to suit your needs, e.g. if you want to return the measurement as a tuple of (kept, peak) values.

+++

For the execution time, you can use the `%timeit` and `%%timeit` magic functions in IPython,
see [corresponding documentation](https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-timeit).

In a nutshell, `%%timeit` measures the execution time of an entire cell, while `%timeit`
only measures the execution time of a single line, e.g. as in

```python
%timeit my_function()
```
__Regarding the usage of `timeit`:__
To obtain reliable timings, `timeit` does not perform a single run, but
rather a number of runs `R`, and in each run, the given statement is executed
`N` times in a loop. This can sometimes lead to long waiting times, so you can set `R` and `N` yourself:
```python
%timeit -r <R> -n <N> my_function()
```
You can also store the timing results by passing the option `-o`:

```python
timings_data = %timeit -o my_function()
```
    
which stores the data from the timing experiment. You can access the best time measured
in seconds by
```python
timings_data.best
```

With these preparations, let us get started!

+++

**a)** Implement the assembly of the mass matrix using a dense, DOK and COO matrix representation. 

For the first one you can use the following code outline:

```{code-cell} ipython3
def assemble_mass_matrix_dense(P, C):
    """Assemble the P1 mass matrix as a dense NumPy array.

    P : (N_p, 2) array, coordinates of the mesh vertices (point matrix).
    C : (N_t, 3) integer array, vertex indices of each triangle (connectivity matrix).

    Returns the (N_p, N_p) mass matrix M as a dense np.ndarray, assembled
    element by element by adding |T_k| * M_ref into the rows and columns C[k].
    """
    
    # Define constant part M_ref here
    # M_ref = ...
    
    # Deduce number of unknowns from dimensions/shape of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of C
    n_t = C.shape[0]
    
    # Create full (dense) matrix M
    M = np.zeros((n_p, n_p))
    
    for T in range(n_t):
        # Remove the pass statement and implement the loop body here
        pass
        # Get local to global map from C
        # ...
        # Get triangle nodes from P
        # ...
        # Unpack nodes into N0,N1,N2
        # ...
        # Compute area of T
        # ...
        # area = ...
        
        # Use area and M_ref to compute M_T
        # M_T = ...
        
        # Instead of 2 loops we can slice out the blocks which
        # correspond to the entries in l2g by using the funny
        # function ix_ in numpy
        # M[np.ix_(l2g, l2g)] += M_T
    
    return M
```

**Debugging/testing hints**: You can compute the mass matrix manually for a small mesh,
e.g. with $N=1$ using the `unitsquare_mesh` function from the `meshtools.py` file, which gives you a mesh with 2 elements, or use the more flexible
`rectangle_mesh` function, where you can specify the number of elements in the $x$ and $y$ directions, e.g. with $N_x=2$ and $N_y=1$, which will produce a mesh with 4 elements.

+++

Now implement the `dok` matrix variant. 

**Hint.** In principle, you only have to change one line compared to the dense matrix variant :)

```{code-cell} ipython3
def assemble_mass_matrix_dok(P, C):
    # Remove the pass statement and implement function body here
    pass 
```

Finally, implement the `coo` matrix variant. Here, try to
avoid the use of a Python for-loop, and rather use vectorized operations from NumPy. 

**Debugging tip.** Since the implementation of the `coo` matrix variant deviates quite a bit from the other two variants, it might
be a good idea to check whether you get the same matrix for all three variants, up to machine precision.

```{code-cell} ipython3
def assemble_mass_matrix_coo(P, C):
    
    # Define local mass matrix on reference element
    M_ref = 1.0/12.0*np.array([[2, 1, 1],[1, 2, 1],[1, 1, 2]])
    
    # Compute e1 and e2 and area for all triangles at once
    # e1 = 
    # e2 = 
    # area =
    
    # Create rows, cols, data for the sparse matrix in COO format 
    
    # You can create rows and cols by repeating and tiling the triangle indices in C
    # rows = ...
    # cols = 
    
    # Compute the values for the sparse matrix in COO format by multiplying the area of each triangle with M_ref.
    # vals = 
    
    # Create sparse matrix M in COO format and return it
    # n_p = P.shape[0]
    # M = ...
    
    # return M
```

**b)** Benchmark the **memory usage** for each of the three implementations
by assembling the mass matrix for a unit square mesh generated
with `unitsquare_mesh` function for $N=5, 10, 20, 40, 80, 160$.
(If you don't have enough memory, you can skip the last value of $N$.)

Report and compare the results, and explain the results as well as you can.

**c)** Benchmark the **execution time** for each of the three implementations
by assembling the mass matrix for a unit square mesh generated
with `unitsquare_mesh` function for $N=5, 10, 20, 40, 80, 160$.
(If you don't have enough memory, you can skip the last one or two values of $N$.)

Report and compare the results, and explain the results as well as you can.

## Problem 2. Quadrature rules

For each degree $d \in \{1,2,3,4\}$, implement one quadrature rule for the unit triangle $\widehat{T}$ from the list given in the lectures.
If several quadrature rules are available for a given degree, give a short justification of your choice.

Check your implementation by verifying that the quadrature rule integrates polynomials of degree $d$ exactly (i.e. up to machine precision error of $10^{-13}$),
and that it does not integrate polynomials of degree $d+1$ exactly.

Throughout this project, we use the convention that the quadrature weights sum to 1, so that
$\int_{\widehat{T}} g \, d\widehat{x} \approx |\widehat{T}| \sum_{q=1}^{n_q} w_q \, g(\widehat{\boldsymbol{x}}_q)$.
If you take a rule from a table where the weights sum to $|\widehat{T}| = 1/2$, rescale the weights accordingly.

## Problem 3. Assembly of the load vector

Recall that the load vector $\boldsymbol{b} = (b_i)_{i=1}^{n_{\mathrm{sh}}}$ is given by  

$$
b_i = \int_\Omega f(x,y) \phi_i(x,y)\, dx dy
$$

where $\{\phi_i(x,y)\}_{i=1}^{n_{\mathrm{sh}}} $ are the basis functions of the finite element space $V_h = \mathbb{P}_1^c(\mathcal{T}_h)$.

Implement a function that assembles the load vector $\boldsymbol{b}$
by completing the code snippet below.

```{code-cell} ipython3
def assemble_load_vector(f, P, C, xq, wq):
    """Assemble the P1 load vector b_i = int f phi_i dx using quadrature.

    f      : function f(x, y), vectorized over NumPy arrays x and y.
    P      : (N_p, 2) array, coordinates of the mesh vertices (point matrix).
    C      : (N_t, 3) integer array, vertex indices of each triangle (connectivity matrix).
    xq, wq : quadrature points (n_q, 2) and weights (n_q,) on the reference
             triangle, with the weights summing to 1.

    Returns the load vector b as a NumPy array of length N_p.
    """
    # Define constant part here: values of the barycentric coordinates
    # lambda_0 = 1 - x - y, lambda_1 = x, lambda_2 = y at the quadrature points
    # Lambda_ref = 
    # Multiply with the weights: L_ref[alpha, q] = w_q * lambda_alpha(xq_q)
    # L_ref = 

    # Deduce number of unknowns from dimensions/shape of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of C
    n_t = C.shape[0]
    
    # Create load vector b
    b = np.zeros(n_p)

    for T in range(n_t):
        pass
        # Get local to global map from C
        
        # Get triangle nodes from P
        
        # Unpack nodes into N0,N1,N2
        
        # Compute area of T
        
        # Map the quadrature points to T: row q is sum_beta lambda_beta(xq_q) * N_beta
        # xq_T = ...
        # Evaluate f at the mapped quadrature points
        # f_q = ... 
        # Local load vector b_T[alpha] = |T| sum_q w_q lambda_alpha(xq_q) f(F_T(xq_q))
        # b_T = 

        # Add the local contributions to the global load vector
        # ...

    return b
```

**Testing/Debugging tips:**

You can test your implementation as follows:
  * If $f = 1$, then the sum of the entries of the load vector $\boldsymbol{b}$ should be equal to the area of the domain $\Omega$, up to a machine precision error of $10^{-13}$. Why?
  * If $f \in \mathbb{P}_1^c(\mathcal{T}_h)$, and you choose a suitable quadrature rule, then $\boldsymbol{b} = M \boldsymbol{f}_I$, where $\boldsymbol{f}_I$ is the vector of nodal values of $f$ at the mesh vertices. Why is this the case, and which order of quadrature rule is suitable?

## Problem 4. The $L^2$ projection

Recall that the $L^2$ projection $P_h f \in V_h = \mathbb{P}_1^c(\mathcal{T}_h)$ of a function $f$ is defined by

$$
\int_\Omega P_h f \, v \, dx = \int_\Omega f \, v \, dx \quad \forall\, v \in V_h .
$$

Writing $P_h f = \sum_{j=1}^{n_{\mathrm{sh}}} u_j \phi_j$ and testing with $v = \phi_i$ leads to the linear system

$$
M \boldsymbol{u} = \boldsymbol{b},
$$

where $M$ is the mass matrix from Problem 1 and $\boldsymbol{b}$ is the load vector from Problem 3.
The solution vector $\boldsymbol{u}$ contains the values of $P_h f$ at the mesh vertices.

### Task 1. Computing the $L^2$ projection

Write a function that computes the coefficient vector $\boldsymbol{u}$ of $P_h f$ for a given function $f$ and a given mesh $(P, C)$.
Reuse your assembly routines for the mass matrix and the load vector, and solve the linear system with `spsolve` from `scipy.sparse.linalg`.
Note that `spsolve` expects the matrix in CSR or CSC format, so convert your mass matrix with `.tocsr()` if necessary.

To assemble the load vector, use two different quadrature rules from Problem 2:
one that is exact for polynomials of degree 1, and one that is exact for polynomials of degree 2.

**a)** Compute $P_h f$ for the functions

1. $f(x,y) = 1$,
2. $f(x,y) = 1 + x + y$,

on the unit square mesh generated by `unitsquare_mesh` with $N = 10$, using both quadrature rules.
Compare the vector $\boldsymbol{u}$ with the vector $\boldsymbol{f}_I$ of nodal values of $f$ by computing $\max_i |u_i - f(\boldsymbol{x}_i)|$.

* What do you expect $P_h f$ to be for functions $f \in V_h$, and why?
* Which of the four computations reproduce $f$ up to machine precision? Explain your observations
  using the degree of the integrand $f \phi_i$ in the load vector.

**b)** Now consider the function

$$
f(x,y) = \sin(2 \pi x) \cos(\pi y).
$$

Compute $P_h f$ on the unit square meshes generated by `unitsquare_mesh` with $N = 10, 20, 40, 80$, again with both quadrature rules.
For each value of $N$, plot $f$, $P_h f$ and their difference $f - P_h f$.
You can use the function `plot_comparison_2D` from `plottools.py`, which shows all three functions in one figure.
Evaluate $f$ at the mesh vertices to plot it.

* Do you observe convergence in the "eyeball norm"?
* By roughly which factor does the size of the difference $f - P_h f$ decrease each time $N$ is doubled?
  Which convergence rate in $h$ does this suggest?
* Does the choice of quadrature rule for the load vector affect what you see?

### Task 2. Convergence study in the $L^2$ norm

To quantify the observations from Task 1 b), we now compute the error $\|f - P_h f\|_{L^2(\Omega)}$.
Splitting the integral into contributions from the individual elements gives

$$
\|f - P_h f\|_{L^2(\Omega)}^2 = \int_\Omega |f - P_h f|^2 \, dx = \sum_{T \in \mathcal{T}_h} \int_T |f - P_h f|^2 \, dx .
$$

On each element $T$ we approximate the integral with a quadrature rule $\{\widehat{\boldsymbol{x}}_q, w_q\}_{q=1}^{n_q}$ on the reference triangle $\widehat{T}$ that is exact for polynomials of degree 4.
With our convention that the weights sum to 1, and with $\boldsymbol{x}_q^T = F_T(\widehat{\boldsymbol{x}}_q)$ denoting the quadrature points mapped to $T$, this reads

$$
\int_T |f - P_h f|^2 \, dx \approx |T| \sum_{q=1}^{n_q} w_q \, \bigl|f(\boldsymbol{x}_q^T) - P_h f(\boldsymbol{x}_q^T)\bigr|^2 .
$$

**Hint:** You do not need to evaluate the global basis functions to compute $P_h f(\boldsymbol{x}_q^T)$.
On $T$, the function $P_h f$ is the linear interpolant of its three vertex values, so

$$
P_h f(\boldsymbol{x}_q^T) = \sum_{\alpha=0}^{2} u_{\iota(T,\alpha)} \, \lambda_\alpha(\widehat{\boldsymbol{x}}_q),
$$

where $\lambda_\alpha$ are the barycentric coordinates on $\widehat{T}$ and $\iota(T, \alpha)$ is the local-to-global map given by the connectivity matrix $C$.
The values $\lambda_\alpha(\widehat{\boldsymbol{x}}_q)$ are the same for all elements, just as in your load vector assembly.

**a)** Implement a function `compute_L2_error(f, u, P, C, xq, wq)` that returns the $L^2$ error $\|f - P_h f\|_{L^2(\Omega)}$
for a given function $f$ and the coefficient vector $\boldsymbol{u}$ of $P_h f$ on the mesh $(P, C)$.

*Testing tip:* For $f \in V_h$ and $\boldsymbol{u} = \boldsymbol{f}_I$, the error must vanish up to machine precision.
For $\boldsymbol{u} = \boldsymbol{0}$ you obtain an approximation of $\|f\|_{L^2(\Omega)}$, which should converge to the exact value
as $N$ grows (for $f(x,y) = \sin(2\pi x)\cos(\pi y)$ the exact value is $1/2$).

**b)** For the function $f$ from Task 1 b), compute $P_h f$ on the unit square meshes with $N = 10, 20, 40, 80$
and the error $e_N = \|f - P_h f\|_{L^2(\Omega)}$ for each $N$.
Do this for both quadrature rules used for the load vector, but always use the degree-4 rule to compute the error.

* Plot $e_N$ as a function of $h = 1/N$ in a log-log plot. Add a reference line proportional to $h^2$.
* Present the results in a table with the columns $N$, $h$, $e_N$ and the experimental order of convergence

  $$
  \mathrm{EOC}_N = \frac{\log(e_{N/2} / e_N)}{\log 2},
  $$

  computed for each $N$ except the coarsest one.
* Which convergence rate do you observe for each of the two quadrature rules?
  Does this agree with your observations from Task 1 b)?
* Why does the error computation need a more accurate quadrature rule than the load vector assembly?
