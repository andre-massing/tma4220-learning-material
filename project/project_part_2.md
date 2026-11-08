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

# Implementation of a finite element solver, part 2: Solving the Poisson problem

+++

The second part of the project will guide you through the implementation of a finite element solver for the Poisson problem in 2D with mixed boundary conditions (i.e. different types of boundary conditions on different parts of the boundary). 
The Poisson problem we consider is of the form

\begin{alignat}{3}
-\Delta u &= f && \quad \text{in } \Omega 
\\
u &= g_D  && \quad \text{on } \partial \Omega_D
\\
\nabla u \cdot n &= g_{N,1} && \quad \text{on } \partial \Omega_{N,1}
\\
\nabla u \cdot n &= g_{N,2} && \quad \text{on } \partial \Omega_{N,2}
\\
\nabla u \cdot n &= \kappa(g_R - u) && \quad \text{on } \partial \Omega_R
\end{alignat}

on the following domain, which is a simplified version of the face of our friendly robot:

+++

![The robot domain with its Dirichlet, Neumann and Robin boundary parts](complex-domain.png)

```{code-cell} ipython3
# Make plotted figures interactive
%matplotlib widget

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from numpy.linalg import matrix_rank
import sympy as sy
import meshtools as mt
import plottools as pt
import tracemalloc

import meshtools as mt
import plottools as pt
```

## Problem 1 (Green formulas and friends)

+++

Derive a weak formulation of the Poisson problem above with mixed boundary conditions
and prove the well-posedness of it using the Lax-Milgram theorem.
You may assume that the boundary of your domain is split up into 4 disjoint parts as above,
 $\partial \Omega = \partial \Omega_D \cup \partial \Omega_{N,1} \cup \partial \Omega_{N,2} \cup \partial \Omega_R$ with $\partial \Omega_D$ having positive measure.

**Hint.** You might find it useful to

+++

**b)**
Choosing $u = uv$ in a), we get the partial integration formula
$$
\int_{\partial \Omega}{uv n_i} \,\mathrm{d}S
= \int_{\Omega} \partial_{x_i} (uv) \,\mathrm{d}x 
= \int_{\Omega} (\partial_{x_i} u)v \,\mathrm{d}x+\int_{\Omega} u \partial_{x_i} v \,\mathrm{d}x,\quad i = 1,\ldots,d
$$
We have
$$
\int_{\Omega} (\partial_{x_i} u) v \,\mathrm{d}x =
\int_{\partial \Omega}{u n_i} v \,\mathrm{d}S
-\int_{\Omega}  u \partial_{x_i} v \,\mathrm{d}x, \quad i = 1,\ldots,d
$$

+++

**c**)
Choosing $u=\beta_i$, the compents of a vector field $\beta$ on $\Omega$ in b), we get
$$
\int_{\Omega} (\partial_{x_i} \beta_i) v \,\mathrm{d}x =
\int_{\partial \Omega}{\beta_i n_i} v \,\mathrm{d}S
-\int_{\Omega}  \beta_i \partial_{x_i} v \,\mathrm{d}x, \quad i = 1,\ldots,d
$$
Taking the sum over $i = 1,\ldots,d$, we have
$$
\int_{\Omega} (\partial_{x_1} \beta_1) v \,\mathrm{d}x =
\int_{\partial \Omega}{\beta_1 n_1} v \,\mathrm{d}S -\int_{\Omega} \beta_1 \partial_{x_1} v \,\mathrm{d}x \\
\int_{\Omega} (\partial_{x_2} \beta_2) v \,\mathrm{d}x =
\int_{\partial \Omega}{\beta_2 n_2} v \,\mathrm{d}S -\int_{\Omega} \beta_2 \partial_{x_2} v \,\mathrm{d}x \\
\ldots \\
\ldots \\
\ldots \\
\int_{\Omega} (\partial_{x_i} \beta_i) v \,\mathrm{d}x =
\int_{\partial \Omega}{\beta_i n_i} v \,\mathrm{d}S -\int_{\Omega} \beta_i \partial_{x_i} v \,\mathrm{d}x 
$$
the sum is 
$$
\int_{\Omega} (\partial_{x_1} \beta_1 +\partial_{x_2} \beta_2+\ldots+\partial_{x_i}
\beta_i) v \,\mathrm{d}x \\
= \int_{\partial \Omega}({\beta_1 n_1+\beta_2 n_2+\ldots+\beta_i n_i}) v \,\mathrm{d}S -\int_{\Omega} (\beta_1 \partial_{x_1} v+ \beta_2 \partial_{x_2} v+\ldots+\beta_i \partial_{x_i} v)\,\mathrm{d}x
$$
For $\beta=\beta_1 x_1+\beta_2 x_2+\ldots+\beta_i x_i$, there are
$$
\nabla \cdot \beta 
=(\partial_{x_1},\partial_{x_2},\ldots,\partial_{x_i})\cdot (\beta_1,\beta_2,\ldots,\beta_i)
=\partial_{x_1} \beta_1 +\partial_{x_2} \beta_2+\ldots+\partial_{x_i}\beta_i \\
\beta \cdot n
= (\beta_1,\beta_2,\ldots,\beta_i) \cdot (n_1,n_2,\ldots,n_i)
= \beta_1 n_1+\beta_2 n_2+\ldots+\beta_i n_i \\
\beta \cdot \nabla
= (\beta_1,\beta_2,\ldots,\beta_i) \cdot (\partial_{x_1}v,\partial_{x_2}v,\ldots,\partial_{x_i}v)
= (\beta_1 \partial_{x_1} v+ \beta_2 \partial_{x_2} v+\ldots+\beta_i \partial_{x_i} v)
$$
Therefore, we have
$$
\int_{\Omega} (\nabla \cdot \beta) v \,\mathrm{d}x =
\int_{\partial \Omega}\beta \cdot n v \,\mathrm{d}S
-\int_{\Omega} \beta \cdot \nabla v \,\mathrm{d}x
$$

+++

**d**)
Choosing $\beta=-A \nabla u$ in c), we get
$$
-\int_{\Omega} \nabla \cdot ( A \nabla u) v \,\mathrm{d}x =
-\int_{\partial \Omega}{  A \nabla u \cdot n} v \,\mathrm{d}S
+\int_{\Omega}  A \nabla u \cdot \nabla v \,\mathrm{d}x,
$$
where $A(x) \in \mathbb{R}^{d\times d}$ is a pointwise defined matrix.

+++

## Problem 2 (Your first own 2D FEM solver)

+++

## Step 1: Assembly of the stiffness matrix (volume contributions)

+++

**a)** Start by implementing a first version of the  ```assemble_stiffness_matrix(P, T, a)``` which assemble the volume contributions corresponding to $\int_K a \nabla u \nabla v dx$.

```{code-cell} ipython3
# 2-a1 Define a function which computes the area of K and the cofficient vector b and c
def compute_hat_gradients(tri):
    
    # Compute area
    N0, N1, N2 = tri
    # area=abs(0.5*np.cross(N1-N0, N2-N0))
    e1, e2 = N1 - N0, N2 - N0
    area = 0.5 * np.abs((e1[0] * e2[1] - e1[1] * e2[0]))
    
    # Compute b and c 
    ones = np.ones(3)
    b = np.cross(tri[:,1], ones)/(2*area)
    c = np.cross(ones, tri[:,0])/(2*area)

    return (area, b, c)
```

Follow closely the implementation of the mass matrix assembly with only minor modifications.

```{code-cell} ipython3
# 2-a2
def assemble_stiffness_matrix(P, T, a):
    # Create matrix as before As in the mass matrix assembly:
    # Deduce number of unkowns from dimensions of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of T and sparse matrix A
    n_t = T.shape[0]
    # Create sparse matrix A
    A = sp.dok_matrix((n_p, n_p))

    for  K in range(n_t):
        # Get local to global map
        l2g = T[K]
        # Get triangle coordinates and compute area
        tri = P[l2g]
        # Unpack nodes into N1,N2,N3
        N0,N1,N2 = tri 
        
        # Compute abar = a((N0 + N1 + N2)/3) 
        abar = a((N0+N1+N2)/3)
        #print ("abar=",abar)
        
        # Compute the area and the coefficient for the hat gradients
        area, b, c = compute_hat_gradients(tri)

        # Numpy Arrays does not behave exactly like n x 1 matrices
        # To compute the outer product b*b.T or  b*b' in Matlab notation
        # we need the np.outer function
        A_K = abar*(np.outer(b, b) + np.outer(c, c))*area
        
        # Add local element matrix to global matrix as before
        A[np.ix_(l2g, l2g)] += A_K

    return A
```

**b**)

```{code-cell} ipython3
# 2-b1 #part1 3-e
def assemble_load_vector_modified(P, T, f, qr="midpoint_2d"):
    # Deduce number of unkowns from dimensions/shape of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of T
    n_t = T.shape[0]
    
    # Create sparse matrix M
    b = np.zeros(n_p)
    
    for  K in range(n_t):
        # Get local to global map from T
        l2g = T[K]   
        # Get triangle nodes from P
        tri = P[l2g]
        # Unpack nodes into N1,N2,N3
        N0,N1,N2 = tri 
        # Compute area of K. 
        # area=abs(0.5*np.cross(N1-N0,N2-N0))
        e1, e2 = N1 - N0, N2 - N0
        area = 0.5 * np.abs((e1[0] * e2[1] - e1[1] * e2[0]))
        
        # Use area and b_ref to compute b_K
        if qr == "midpoint_2d":
            mp =  0.5*np.array([N0+N1, N1+N2, N2+N0])
            b_K = area/3*0.5*np.array([f(mp[i-1,0], mp[i-1,1]) + f(mp[i,0], mp[i,1]) for i in range(0,3)])
        else:
            b_K = area/3*np.array([f(tri[i,0], tri[i,1]) for i in range(0,3)])
        b[np.ix_(l2g)] += b_K
    return b
```

```{code-cell} ipython3
# 2-b2

# Create unitsquare mesh
N = 20
P, T = mt.unitsquare_mesh(N)

x = P[:,0]
y = P[:,1]
# Numpy way to define constant functions such that an
# array of constants is returned if x,y were arrays.
def f(x,y):
    return 1*np.ones_like(x)  # f(x,y)=1
def a(x=None):
    return 1                  # a(x)=1
def u(x,y):
    return 1*np.ones_like(x) 
qr = "midpoint_2d"

A = assemble_stiffness_matrix(P, T, a)
b = assemble_load_vector_modified(P, T, f, qr)
#print (b)
A_csr = A.tocsr() 
#print (A_csr)
U = spsolve(A_csr, b)
#print (U)
print("N =",N)
pt.plot_comparison_2D(x, y, u(x,y), U, triangles=T,
                      title_f1='U_ex', title_f2='U')
A_mat = np.matrix(A_csr.todense())
rank = matrix_rank(A_mat)
print ("rank =",rank)
print("Number of unknowns:", len(U))
```

## Step 2: Incorporation of Dirichlet boundary conditions

+++

**c**) Next, we start with implementing the incorporation of Dirchlet boundary conditions.

```{code-cell} ipython3
# Define a function which checks 
# whether a point x=(x0,x1)is on the boundary of unitsquare
def on_dirichlet_boundary(x):
    eps = 1e-12 # Accounts for finite precision
    return  (x[0] < eps or x[1] < eps or x[0] > 1-eps or x[1] > 1-eps  )
```

```{code-cell} ipython3
# Define a function which extract the boundary point 
# and return a list of all node indices i
def extract_nodes(P, on_dirichlet_boundary):
    return [ i for i in range(P.shape[0]) if on_dirichlet_boundary(P[i]) ]
```

Now define unitsquare meshes with 1, 2 and 3 subdivision in each direction and 
check whether you extract all Dirichlet boundary nodes correctly by
using the improved `plot_mesh_2d` function now residing in the `plottools` module. Use the new optional 3 function argument "dirichlet_nodes" in `plot_mesh_2d` to plot the extracted nodes.

```{code-cell} ipython3
# Create unitsquare mesh
for N in range(1,4):
    P, T = mt.unitsquare_mesh(N)
    dirichlet_nodes = extract_nodes(P, on_dirichlet_boundary)
    pt.plot_mesh_2d(P, T,dirichlet_nodes)
```

**d**) Next you need a piece of code to evaluate boundary functions at the nodes you extracted and to incorporate these values into the matrix $A$ and $b$.

+++

Example of evaluation of $g_D$ at Dirichlet nodes

```{code-cell} ipython3
# Generat a mesh
P, T = mt.unitsquare_mesh(5)

# Extract mesh coordinates
x = P[:, 0]
y = P[:, 1]

# Define/use your on_dirichlet_boundary and extract_nodes functions
# to compute the dirichlet nodes
dirichlet_nodes = extract_nodes(P, on_dirichlet_boundary)

# Extract coordinates associated with Dirichlet nodes 
# using an array slicing operation
X_dc = x[dirichlet_nodes]
Y_dc = y[dirichlet_nodes]

# Evaluate g at boundary nodes and 
def g_D(x,y):
    return (x+y)              # g_D=x+y
g_D_values = g_D(X_dc, Y_dc)
```

Incorporation of Dirichlet BC.

```{code-cell} ipython3
def apply_bcs(A, b, dirichlet_nodes, g_D_values):
    # Incorporate boundary condition in vector b
    # by slicing out entries of b associated with dirichlet nodes
    # and setting them to the corresponding value
    b[dirichlet_nodes] = g_D_values

    # Incorporate boundary conditions in matrix A
    # 1. Set all rows corresponding to Dirichlet nodes to 0
    A[dirichlet_nodes, :] = 0 
    
    # 2. Set diagonal in Dirichlet nodes to one 
    # (Opposed to Matlab this code line practically
    # extracts a part of the diagonal corresponding to passes indices in  
    # dirichlet_nodes and set them to 1
    A[dirichlet_nodes, dirichlet_nodes] = 1
```

**e**) Next, do a first check of the correctness of your implementation by solving

+++

\begin{alignat*}{3}
-\Delta u &= 0 && \quad \text{in } \Omega \\
 u(x,y) &= x + y &&\quad \text{on } \partial \Omega
\end{alignat*}

+++

and plotting $u$, $u_h$ and their difference. 

What is the exact solution $u$? And how does the finite element solution $u_h$ look like?
*Hint*: Best approximation property.

```{code-cell} ipython3
for N in [1,5,15]:
    P, T = mt.unitsquare_mesh(N)
    x = P[:,0]
    y = P[:,1]
    def g_D(x,y):
        return (x+y)              # g_D=x+y
    def f(x,y):
        return 0*np.ones_like(x)  # f(x,y)=0
    def a(x=None):
        return 1                  # a(x)=1
    def u(x,y):
        return x+y 
    qr = "midpoint_2d"
    
    dirichlet_nodes = extract_nodes(P, on_dirichlet_boundary)
    X_dc = x[dirichlet_nodes]
    Y_dc = y[dirichlet_nodes]
    g_D_values = g_D(X_dc, Y_dc)

    A = assemble_stiffness_matrix(P, T, a)
    b = assemble_load_vector_modified(P, T, f, qr)
    apply_bcs(A, b, dirichlet_nodes, g_D_values)
    A_csr = A.tocsr() 
    U = spsolve(A_csr, b)
    print("N =",N)
    pt.plot_comparison_2D(x, y, u(x,y), U, triangles=T,
                          title_f1='exact', title_f2='U')
```

**f**) Next, solve the problem
\begin{alignat*}{3}
-\Delta u &= 1 && \quad \text{in } \Omega \\
 u &= 0 &&\quad \text{on } \partial \Omega.
\end{alignat*}
on unit square meshes with mesh size $h=1/N$ for $N = 5, 10, 20, 40$
and plot the FEM solution.

```{code-cell} ipython3
for N in [5,10,20,40]:
    P, T = mt.unitsquare_mesh(N)
    x = P[:,0]
    y = P[:,1]
    def g_D(x,y):
        return 0*np.ones_like(x)  # g_D=0
    def f(x,y):
        return 1*np.ones_like(x)  # f(x,y)=1
    def a(x=None):
        return 1                  # a(x)=1
    qr = "midpoint_2d"

    dirichlet_nodes = extract_nodes(P, on_dirichlet_boundary)
    X_dc = x[dirichlet_nodes]
    Y_dc = y[dirichlet_nodes]
    g_D_values = g_D(X_dc, Y_dc)

    A = assemble_stiffness_matrix(P, T, a)
    b = assemble_load_vector_modified(P, T, f, qr)
    apply_bcs(A, b, dirichlet_nodes, g_D_values)
    A_csr = A.tocsr() 
    U = spsolve(A_csr, b)
    print("N =",N)
    pt.plot2D(x, y, U, title='U', triangles=T)
```

**g**) Repeat the experiment, this time for the boundary value problem 
\begin{alignat*}{3}
-\Delta u &= 1 && \quad \text{in } \Omega \\
 u &= 0 &&\quad \text{on } \partial \Omega_D,
 \\
 \nabla u \cdot n &= 0 &&\quad \text{on } \partial \Omega_N,
\end{alignat*}
where $\partial \Omega_N = \{1\} \times [0,1]$ is right edge of the the unit square
and $\partial \Omega_D = \partial \Omega \setminus \partial \Omega_N$.

How can you interpret the Neuman condition $\nabla u \cdot n = 0$ on $\partial \Omega_N$ geometrically?

```{code-cell} ipython3
for N in [5,10,20,40]:
    P, T = mt.unitsquare_mesh(N)
    x = P[:,0]
    y = P[:,1]
    def on_neuman_boundary(x):
        eps = 1e-12 # Accounts for finite precision
        return  (x[0] > 1-eps)
    def on_dirichlet_boundary(x):
        eps = 1e-12 # Accounts for finite precision
        return  (x[0] < eps or x[1] < eps or x[1] > 1-eps)
    def g_D(x,y):
        return 0*np.ones_like(x)  # g_D=0
    def g_N(x,y):
        return 0*np.ones_like(x)  # g_N=0
    def f(x,y):
        return 1*np.ones_like(x)  # f(x,y)=1
    def a(x=None):
        return 1                  # a(x)=1
    qr = "midpoint_2d"

    dirichlet_nodes = extract_nodes(P, on_dirichlet_boundary)
    X_dc = x[dirichlet_nodes]
    Y_dc = y[dirichlet_nodes]
    g_D_values = g_D(X_dc, Y_dc)

    A = assemble_stiffness_matrix(P, T, a)
    b = assemble_load_vector_modified(P, T, f, qr)
    apply_bcs(A, b, dirichlet_nodes, g_D_values)
    A_csr = A.tocsr() 
    U = spsolve(A_csr, b)
    print("N =",N)
    pt.plot2D(x, y, U, title='U', triangles=T)
```

## Step 3: Incorporation of inhomogeneous Neumann boundary conditions

+++

**h**) The next task is to incorporate inhomogeneous Neumann conditions
\begin{align*}
\nabla u \cdot n &= g_N \quad \text{on } \partial \Omega_N
\end{align*}
for some $g_N \neq 0$. Thus to obtain the load vector $b$, you now need to compute integrals of the form   
\begin{align*}
\int_{\Omega} f \varphi_i
\,\mathrm{d}x + \int_{\partial \Omega_N} g_N \varphi_i \quad i = 0,\ldots,n_p-1.
\end{align*}

```{code-cell} ipython3
def extract_edges(P, T, inside_domain):
    """Short function which extract all edge indices i
    for which the inside_domain(P[i]) evaluates to True
    for both edge vertices.
    """
    edges = []
    for t in T:
        # Iterate over edges
        for i in range(3):
            # Check whether end points of edges are inside domain
            # Note that we use the index notation t[-1] to access the last node
            if inside_domain(P[t[i-1]]) and inside_domain(P[t[i]]):
                edges.append([t[i-1], t[i]])
    # Return the edge list as numpy array to treat it similar
    # to T. Note that we have to explicitly tell the array constructor
    # that the element of the array should be treat as integers!
    # Otherwise you will get indexing problems later one.
    return np.array(edges, dtype=np.int32)
```

```{code-cell} ipython3
# Define proper 'inside_domain' function 
def inhomog_neumann_boundary(x):
    eps = 1e-12 # Accounts for finite precision
    return  (x[0] < eps) 
```

```{code-cell} ipython3
# Use 'plot_mesh_2d' function to visually check 'inhomog_neumann_boundary'
# from matplotlib.path import Path
# from matplotlib.patches import PathPatch

for N in range(1,4):
    P, T = mt.unitsquare_mesh(N)
    inhomog_neumann_nodes = extract_nodes(P, inhomog_neumann_boundary)
    inhomog_neuman_edges = extract_edges(P, T, inhomog_neumann_boundary)
    print (inhomog_neuman_edges)
    pt.plot_mesh_2d(P, T, inhomog_neumann_nodes, inhomog_neuman_edges)
```

```{code-cell} ipython3
# The load veator with boundary condition.
def assemble_load_vector_modified_withbc(P, T, f, qr, E_N, g_N):
    # Deduce number of unkowns from dimensions/shape of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of T
    n_t = T.shape[0]
    n_e = E_N.shape[0]
    
    # Create sparse matrix M
    b = np.zeros(n_p)
    
    for  K in range (n_t):
        # Get local to global map from T
        l2g = T[K]   
        # Get triangle nodes from P
        tri = P[l2g]
        # Unpack nodes into N1,N2,N3
        N0,N1,N2 = tri 
        # Compute area of K. 
        # area=abs(0.5*np.cross(N1-N0,N2-N0))
        e1, e2 = N1 - N0, N2 - N0
        area = 0.5 * np.abs((e1[0] * e2[1] - e1[1] * e2[0]))
        
        # Use area and b_ref to compute b_K
        if qr == "midpoint_2d":
            mp =  0.5*np.array([N0+N1, N1+N2, N2+N0])
            b_K = area/3*0.5*np.array([f(mp[i-1,0], mp[i-1,1]) + f(mp[i,0], mp[i,1]) 
                                       for i in range(0,3)])
            # Here the value of 0.5 comes from the hat function at the mid-point.
        else:
            b_K = area/3*np.array([f(tri[i,0], tri[i,1]) for i in range(0,3)])
        b[np.ix_(l2g)] += b_K
        
    for l in range (n_e):
        l2g = E_N[l]
        tri = P[l2g]
        N0,N1 = tri
        E = np.sqrt((tri[0,0]-tri[1,0])**2+(tri[0,1]-tri[1,1])**2)
        # the length of E : \E\=((x(0)-x(1))^2+(y(0)-y(1))^2)^(1/2)
        #b_l = E/2*np.array([g_N(tri[i,0], tri[i,1])for i in range(0,2)])
        #if qr == "midpoint_2d":
            #mp =  0.5*(N0+N1)
            #print("mp")
            #print(mp)
            #print("mp stop")
            #b_l = 0.5*E*np.array([g_N(mp[0], mp[1]), g_N(mp[0], mp[1])])
            #b_l = 0.5*E*np.array([g_N(mp[0], mp[1])],[ g_N(mp[0], mp[1])])
                    
        #else:
        b_l = E/2*np.array([g_N(tri[i,0], tri[i,1]) for i in range(0,2)])
        b[np.ix_(l2g)] += b_l
    return b
```

**i**) 
Next, to appreciate your craftsmanship, validate your new functionality using the method
of manufactured solutions. Start from the exact solution
$$
u_{ex}(x,y) = \sin(\pi x) \cos(\pi y)
$$
to construct a manufactured solution for the boundary value problem
\begin{alignat*}{3}
-\Delta u &= f && \quad \text{in } \Omega \\
 u &= g_D &&\quad \text{on } \partial \Omega_D,
 \\
 \nabla u \cdot n &= g_N &&\quad \text{on } \partial \Omega_N,
\end{alignat*}
by computing the proper problem data $f$, $g_D$ and $g_N$. Here,
$\partial \Omega_N = \{0\} \times[0,1]$ is the left edge of the unit
square and
$\partial \Omega_D = \partial \Omega \setminus \partial \Omega_N$.

+++

**Computing the error in the $H^1$ semi-norm.**
The error $u - u_h$ in the $H^1$ semi-norm is given by

$$
|u - u_h|_{H^1(\Omega)}^2 = \int_\Omega |\nabla u - \nabla u_h|^2 \,\mathrm{d}x
= \sum_{K \in \mathcal{T}_h} \int_K |\nabla u - \nabla u_h|^2 \,\mathrm{d}x .
$$

On each element $K$ with vertices $N_0, N_1, N_2$ and global indices $\iota(K, \alpha)$, the discrete solution is $u_h|_K = \sum_{\alpha=0}^{2} U_{\iota(K,\alpha)} \varphi_\alpha$.
The local basis functions $\varphi_\alpha$ are linear, so $\nabla u_h$ is *constant* on $K$:

$$
\nabla u_h|_K = \sum_{\alpha=0}^{2} U_{\iota(K,\alpha)} \nabla \varphi_\alpha
= \begin{pmatrix} \sum_\alpha b_\alpha U_{\iota(K,\alpha)} \\ \sum_\alpha c_\alpha U_{\iota(K,\alpha)} \end{pmatrix},
$$

where $\nabla \varphi_\alpha = (b_\alpha, c_\alpha)$ are exactly the hat function gradients returned by `compute_hat_gradients`.
In the code, this is `grad_uh = np.array([b @ U[l2g], c @ U[l2g]])`.
The exact gradient $\nabla u$ is not constant, so we approximate the element integral with a quadrature rule
$\{\widehat{\boldsymbol{x}}_q, w_q\}_{q=1}^{n_q}$ on the reference triangle (weights summing to 1), mapped to $K$:

$$
\int_K |\nabla u - \nabla u_h|^2 \,\mathrm{d}x \approx |K| \sum_{q=1}^{n_q} w_q \, \bigl|\nabla u(\boldsymbol{x}_q^K) - \nabla u_h|_K\bigr|^2,
\qquad \boldsymbol{x}_q^K = \sum_{\alpha=0}^{2} \lambda_\alpha(\widehat{\boldsymbol{x}}_q) N_\alpha .
$$

In the code, the mapped points are `xq_K = Lambda_ref.T @ tri`, and the sum over $q$ is `area * np.sum(wq * ...)`.
Since $\nabla u - \nabla u_h$ behaves locally like a linear polynomial, its square is essentially quadratic.
We use the six-point rule of degree 4, which integrates this leading part exactly with a safety margin.

Note that, unlike the energy-norm formula above, this computation does not rely on Galerkin orthogonality,
so it is also valid for inhomogeneous Dirichlet data.

```{code-cell} ipython3
def six_point_rule():
    """Six-point quadrature rule on the reference triangle, exact for polynomials
    of degree 4. Weights sum to 1."""
    a1, w1 = 0.445948490915965, 0.223381589678011
    a2, w2 = 0.091576213509771, 0.109951743655322
    xq = np.array([[a1, a1], [1 - 2*a1, a1], [a1, 1 - 2*a1],
                   [a2, a2], [1 - 2*a2, a2], [a2, 1 - 2*a2]])
    wq = np.array([w1, w1, w1, w2, w2, w2])
    return xq, wq


def compute_H1_semi_error(grad_u, U, P, T, xq=None, wq=None):
    """Compute the H1 semi-norm error |u - u_h|_{H^1} of a P1 solution using quadrature.

    grad_u : function grad_u(x, y) returning the exact gradient as a tuple
             (du/dx, du/dy), vectorized over NumPy arrays x and y.
    U      : (n_p,) array, nodal values of the P1 solution u_h.
    P, T   : point matrix (n_p, 2) and connectivity matrix (n_t, 3).
    xq, wq : quadrature rule on the reference triangle with weights summing to 1
             (default: six-point rule of degree 4).

    Returns the approximate error |u - u_h|_{H^1(Omega)} as a float.
    """
    if xq is None:
        xq, wq = six_point_rule()
    # Barycentric coordinates at the quadrature points, shape (3, n_q)
    Lambda_ref = np.array([1 - xq[:, 0] - xq[:, 1], xq[:, 0], xq[:, 1]])

    err2 = 0.0
    for K in range(T.shape[0]):
        # Get local to global map and triangle nodes
        l2g = T[K]
        tri = P[l2g]
        # Area and hat function gradients (b_alpha, c_alpha)
        area, b, c = compute_hat_gradients(tri)
        # Gradient of u_h on K (constant)
        grad_uh = np.array([b @ U[l2g], c @ U[l2g]])                  # (2,)
        # Map quadrature points to K and evaluate the exact gradient there
        xq_K = Lambda_ref.T @ tri                                      # (n_q, 2)
        ux, uy = grad_u(xq_K[:, 0], xq_K[:, 1])                        # (n_q,), (n_q,)
        # Quadrature approximation of int_K |grad u - grad u_h|^2 dx
        err2 += area * np.sum(wq * ((ux - grad_uh[0])**2 + (uy - grad_uh[1])**2))

    return np.sqrt(err2)

def compute_L2_error(f, u, P, C, xq, wq):
    """Compute the L2 error ||f - u_h|| of the P1 function u_h = sum_i u_i phi_i using quadrature.

    f      : function f(x, y), vectorized over NumPy arrays x and y.
    u      : (N_p,) array, values of u_h at the mesh vertices.
    P, C   : point matrix (N_p, 2) and connectivity matrix (N_t, 3).
    xq, wq : quadrature rule on the reference triangle, with the weights summing to 1.

    Returns the approximate L2 error as a float.
    """
    # Values of the barycentric coordinates at the quadrature points, shape (3, n_q)
    Lambda_ref = np.array([1 - xq[:, 0] - xq[:, 1], xq[:, 0], xq[:, 1]])

    n_t = C.shape[0]
    err2 = 0.0
    for T in range(n_t):
        # Get local to global map from C
        l2g = C[T]
        # Get triangle nodes from P
        tri = P[l2g]
        # Unpack nodes into N0,N1,N2
        N0, N1, N2 = tri
        # Compute area of T
        e1, e2 = N1 - N0, N2 - N0
        area = 0.5 * np.abs((e1[0] * e2[1] - e1[1] * e2[0]))
        # Map the quadrature points to T
        xq_T = Lambda_ref.T @ tri                                  # (n_q, 2)
        # Evaluate f and u_h at the mapped quadrature points
        f_q = f(xq_T[:, 0], xq_T[:, 1])                            # (n_q,)
        uh_q = Lambda_ref.T @ u[l2g]                               # (n_q,)
        # Add the quadrature approximation of int_T |f - u_h|^2 dx
        err2 += area * np.sum(wq * (f_q - uh_q)**2)

    return np.sqrt(err2)
```

We test the function first.
For a linear function, $u_h$ reproduces $u$ exactly, so the error must vanish up to machine precision.
For $U = 0$ we obtain $|u_{ex}|_{H^1(\Omega)} = \pi/\sqrt{2}$.
Then we repeat the convergence study from **i)** and compare with the energy-norm formula used above.

```{code-cell} ipython3
def u_ex(x,y):
    return np.sin(np.pi*x)*np.cos(np.pi*y)
def grad_u_ex(x, y):
    return (np.pi*np.cos(np.pi*x)*np.cos(np.pi*y),
            -np.pi*np.sin(np.pi*x)*np.sin(np.pi*y))
    
def f(x,y):
    return 2*np.pi**2*np.sin(np.pi*x)*np.cos(np.pi*y)

def on_dirichlet_boundary(x):
    eps = 1e-12 # Accounts for finite precision
    return  (x[0] > 1-eps or x[1] < eps or x[1] > 1-eps)
def g_D(x,y):
    return np.sin(np.pi*x)*np.cos(np.pi*y)  
    
def inhomog_neumann_boundary(x):
    eps = 1e-12 # Accounts for finite precision
    return  (x[0] < eps) 
def g_N(x,y):
    return -np.pi*np.cos(np.pi*x)*np.cos(np.pi*y)

def a(x=None):
    return 1  
```

```{code-cell} ipython3
# Sanity checks
P, T = mt.unitsquare_mesh(10)
U_lin = 1 + 2*P[:, 0] - 3*P[:, 1]
print("linear u, U = u_I:  error =", compute_H1_semi_error(lambda x, y: (2 + 0*x, -3 + 0*y), U_lin, P, T))
print("U = 0:  |u_ex|_H1 =", compute_H1_semi_error(grad_u_ex, np.zeros(P.shape[0]), P, T),
      " exact:", np.pi/np.sqrt(2))

# Convergence study, same setup as in i)
refs = np.array([ 5*2**j for j in range(0,4)])
errors_H1 = np.zeros(refs.size)
errors_L2 = np.zeros(refs.size)
# Quadrature rule for error computation, exact for polynomials of degree 4
xq_err, wq_err = six_point_rule()

for j, N in enumerate(refs):
    P, T = mt.unitsquare_mesh(N)
    x, y = P[:, 0], P[:, 1]
    E_N = extract_edges(P, T, inhomog_neumann_boundary)
    dirichlet_nodes = extract_nodes(P, on_dirichlet_boundary)
    A = assemble_stiffness_matrix(P, T, a)
    b = assemble_load_vector_modified_withbc(P, T, f, qr, E_N, g_N)
    apply_bcs(A, b, dirichlet_nodes, g_D(x[dirichlet_nodes], y[dirichlet_nodes]))
    U = spsolve(A.tocsr(), b)
    errors_H1[j] = compute_H1_semi_error(grad_u_ex, U, P, T, xq_err, wq_err)
    errors_L2[j] = compute_L2_error(u_ex, U, P, T, xq_err, wq_err)
    if j in [0,2,3]:
        pt.plot_comparison_2D(x, y, u_ex(x,y), U, triangles=T,
                              title_f1='u_ex', title_f2='U')
```

```{code-cell} ipython3
# Plot convergence in H1 semi-norm and L2 norm
eocs_H1 = np.log(errors_H1[:-1] / errors_H1[1:]) / np.log(2)
eoc_L2 = np.log(errors_L2[:-1] / errors_L2[1:]) / np.log(2)
hs = 1.0 / refs
```

```{code-cell} ipython3
import pandas as pd

# Convergence table: the EOC is undefined on the coarsest mesh (NaN, shown as "–")
conv_table = pd.DataFrame({
    "N": refs,
    "h": hs,
    "L2 error": errors_L2,
    "EOC (L2)": np.concatenate([[np.nan], eoc_L2]),
    "H1 error": errors_H1,
    "EOC (H1)": np.concatenate([[np.nan], eocs_H1]),
})
display(conv_table.style
        .format({"h": "{:.4f}", "L2 error": "{:.3e}", "H1 error": "{:.3e}",
                 "EOC (L2)": "{:.2f}", "EOC (H1)": "{:.2f}"}, na_rep="–")
        .hide(axis="index"))

# Log-log plot of the errors against h, with reference slopes h and h^2
fig, ax = plt.subplots()
ax.loglog(hs, errors_L2, "o-", label=r"$\|u - u_h\|_{L^2(\Omega)}$")
ax.loglog(hs, errors_H1, "s-", label=r"$|u - u_h|_{H^1(\Omega)}$")
ax.loglog(hs, errors_L2[0] * (hs / hs[0])**2, "k--", label=r"$\mathcal{O}(h^2)$")
ax.loglog(hs, errors_H1[0] * (hs / hs[0]), "k:", label=r"$\mathcal{O}(h)$")
ax.set_xlabel("$h = 1/N$")
ax.set_ylabel("error")
ax.set_title("Convergence of the P1 solution")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
plt.show()
```

**j)** Finally we want to solve 

\begin{alignat}{3}
- \Delta u &= f && \quad \text{in } \Omega 
\\
u &= g_D  && \quad \text{on } \partial \Omega_D
\\
\nabla u \cdot n &= g_{N,1} && \quad \text{on } \partial \Omega_{N,1}
\\
\nabla u \cdot n &= g_{N,2} && \quad \text{on } \partial \Omega_{N,1}
\\
 \nabla u \cdot n &= \kappa(g_R - u) && \quad \text{on } \partial \Omega_R
\end{alignat}

posed on our friendly robot face domain given as described in the following sketch.

```{code-cell} ipython3
# The load veator with boundary condition.
def assemble_load_vector_modified_withbc_all(P, T, f, qr, E_N, g_N, E_R, g_R):
    # Deduce number of unkowns from dimensions/shape of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of T
    n_t = T.shape[0]
    n_e = E_N.shape[0]
    r_e = E_R.shape[0]
    
    # Create sparse matrix M
    b = np.zeros(n_p)
    
    for  K in range (n_t):
        # Get local to global map from T
        l2g = T[K]   
        # Get triangle nodes from P
        tri = P[l2g]
        # Unpack nodes into N1,N2,N3
        N0,N1,N2 = tri 
        # Compute area of K. 
        area=abs(0.5*np.cross(N1-N0,N2-N0))
        # Use area and b_ref to compute b_K
        if qr == "midpoint_2d":
            mp =  0.5*np.array([N0+N1, N1+N2, N2+N0])
            b_K = area/3*0.5*np.array([f(mp[i-1,0], mp[i-1,1]) + f(mp[i,0], mp[i,1]) for i in range(0,3)])
        else:
            b_K = area/3*np.array([f(tri[i,0], tri[i,1]) for i in range(0,3)])
        b[np.ix_(l2g)] += b_K
    for l in range (n_e):
        l2g = E_N[l]
        tri = P[l2g]
        E = np.sqrt((tri[0,0]-tri[1,0])**2+(tri[0,1]-tri[1,1])**2)
        # the length of E : \E\=((x(0)-x(1))^2+(y(0)-y(1))^2)^(1/2)
        if qr == "midpoint_2d":
            mp =  0.5*np.array([N0+N1])
            b_l = E/2*0.5*np.array([g_N(mp[i,0], mp[i,1])for i in range(0,1)])
        else:
            b_l = E/2*np.array([g_N(tri[i,0], tri[i,1]) for i in range(0,1)])
        b[np.ix_(l2g)] += b_l
    for l in range (r_e):
        l2g = E_R[l]
        tri = P[l2g]
        E = np.sqrt((tri[0,0]-tri[1,0])**2+(tri[0,1]-tri[1,1])**2)
        # the length of E : \E\=((x(0)-x(1))^2+(y(0)-y(1))^2)^(1/2)
        if qr == "midpoint_2d":
            mp =  0.5*np.array([N0+N1])
            #b_l = kappa*E/2*0.5*np.array([g_R(mp[0], mp[1]), g_R(mp[0], mp[1])])
            b_l = kappa*E/2*0.5*np.array([g_R(mp[i,0], mp[i,1])for i in range(0,1)])
        else:
            b_l = kappa*E/2*np.array([g_R(tri[i,0], tri[i,1]) for i in range(0,1)])
        b[np.ix_(l2g)] += b_l
    return b
```

```{code-cell} ipython3
#  The stiffness_matrix with all boundary condition modified
def assemble_stiffness_matrix_modified(P, T, a, E_R):
    # Create matrix as before As in the mass matrix assembly:
    # Deduce number of unkowns from dimensions of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of T and sparse matrix A
    n_t = T.shape[0]
    # Create sparse matrix A
    A = sp.dok_matrix((n_p, n_p))
    r_e = E_R.shape[0]

    for  K in range(n_t):
        # Get local to global map
        l2g = T[K]
        # Get triangle coordinates and compute area
        tri = P[l2g]
        # Unpack nodes into N1,N2,N3
        N0,N1,N2 = tri 
        
        # Compute abar = a((N0 + N1 + N2)/3) 
        abar = a((N0+N1+N2)/3)
        #print ("abar=",abar)
        
        # Compute the area and the coefficient for the hat gradients
        area, b, c = compute_hat_gradients(tri)

        # Numpy Arrays does not behave exactly like n x 1 matrices
        # To compute the outer product b*b.T or  b*b' in Matlab notation
        # we need the np.outer function
        A_K = abar*(np.outer(b, b) + np.outer(c, c))*area
        
        # Add local element matrix to global matrix as before
        A[np.ix_(l2g, l2g)] += A_K
    for k in range (r_e):
        l2g = E_R[k]
        tri = P[l2g]
        E = np.sqrt((tri[0,0]-tri[1,0])**2+(tri[0,1]-tri[1,1])**2)
        # the length of E : \E\=((x(0)-x(1))^2+(y(0)-y(1))^2)^(1/2)
        A_K_R = 1/6*kappa*E*np.array([[2, 1],[1, 2]])
        A[np.ix_(l2g, l2g)] += A_K_R
    return A
```

```{code-cell} ipython3
def on_dirichlet_boundary(x):
        eps = 1e-12 
        return  (x[0] > 5-eps or x[0] < eps or x[1] > 5-eps
                or (x[1] < eps and x[1] > -eps and x[0] < 1.8+eps and x[0] > 3.2-eps)
                or  x[1] < -2+eps 
                or (x[0] < 1.8+eps and x[1] < eps) 
                or (x[0] > 3.2-eps and x[1] < eps))
ref = np.array([ 5*2**j for j in range(1,4)])
print (ref)
for j in range(1,4):
    P = np.loadtxt("meshes/mesh_points_ref_%d.txt" % ref[j-1]) 
    T = np.loadtxt("meshes/mesh_triangles_ref_%d.txt" % ref[j-1], dtype=np.int32) 
    neuman_nodes = extract_nodes(P, on_dirichlet_boundary)
    print ("N =",ref[j-1])
    pt.plot_mesh_2d(P, T, neuman_nodes)
```

```{code-cell} ipython3
def on_neuman_boundary_1(x):
        eps = 1e-12 # Accounts for finite precision
        return  (
                ((x[0]-1)**2+(x[1]-3.9)**2) > (0.5-eps)**2
                 and ((x[0]-1)**2+(x[1]-3.9)**2) < (0.5+eps)**2
                 or 
                 ((x[0]-4)**2+(x[1]-3.9)**2)> (0.5-eps) **(2)
                 and ((x[0]-4)**2+(x[1]-3.9)**2) < 0.5+eps**(2)
                )
ref = 10
P = np.loadtxt("meshes/mesh_points_ref_%d.txt" % ref) 
T = np.loadtxt("meshes/mesh_triangles_ref_%d.txt" % ref, dtype=np.int32) 
neuman_nodes_1 = extract_nodes(P, on_neuman_boundary_1)
print ("N =",ref)
pt.plot_mesh_2d(P, T, neuman_nodes_1)
```

```{code-cell} ipython3
def on_neuman_boundary_2(x):
        eps = 1e-12 
        return  (
                 (x[0] > 2.25-eps and x[0] < 2.75+eps and 
                  x[1] > 2-eps and x[1] < 3.5+eps) 
                 )
ref=10
P = np.loadtxt("meshes/mesh_points_ref_%d.txt" % ref) 
T = np.loadtxt("meshes/mesh_triangles_ref_%d.txt" % ref, dtype=np.int32) 
neuman_nodes_2 = extract_nodes(P, on_neuman_boundary_2)
print ("N =",ref)
pt.plot_mesh_2d(P, T, neuman_nodes_2)
```

```{code-cell} ipython3
def on_robin_boundary(x):
        eps = 1e-12 
        return  (
                (x[0] > 1.8-eps and x[0] < 3.2+eps and x[1] > 0.7-eps and x[1] < 1.3+eps)
                )
ref=10
P = np.loadtxt("meshes/mesh_points_ref_%d.txt" % ref) 
T = np.loadtxt("meshes/mesh_triangles_ref_%d.txt" % ref, dtype=np.int32) 
robin_nodes = extract_nodes(P, on_robin_boundary)
print ("N =",ref)
pt.plot_mesh_2d(P, T, robin_nodes)
```

```{code-cell} ipython3
# j 
refs = np.array([ 5*2**j for j in range(1,4)])
for j in range(1,4):
    N = refs[j-1]
    h[j-1] = 1/refs[j-1]
    P = np.loadtxt("meshes/mesh_points_ref_%d.txt" % refs[j-1]) 
    T = np.loadtxt("meshes/mesh_triangles_ref_%d.txt" % refs[j-1], dtype=np.int32) 
    x = P[:,0]
    y = P[:,1]
    
    def on_neuman_boundary_1(x):
        eps = 1e-12 # Accounts for finite precision
        return  (
                ((x[0]-1)**2+(x[1]-3.9)**2) > (0.5-eps)**2
                 and ((x[0]-1)**2+(x[1]-3.9)**2) < (0.5+eps)**2
                 or 
                 ((x[0]-4)**2+(x[1]-3.9)**2)> (0.5-eps) **(2)
                 and ((x[0]-4)**2+(x[1]-3.9)**2) < 0.5+eps**(2)
                )
    def on_neuman_boundary_2(x):
        eps = 1e-12 
        return  ((x[0] > 2.25-eps and x[0] < 2.75+eps and x[1] > 2-eps 
                  and x[1] < 3.5+eps) )
    def on_dirichlet_boundary(x):
        eps = 1e-12 
        return  (x[0] > 5-eps or x[0] < eps or x[1] > 5-eps
                or (x[1] < eps and x[1] > -eps and x[0] < 1.8+eps and x[0] > 3.2-eps)
                or  x[1] < -2+eps 
                or (x[0] < 1.8+eps and x[1] < eps) 
                or (x[0] > 3.2-eps and x[1] < eps))
    def on_robin_boundary(x):
        eps = 1e-12 
        return  ((x[0] > 1.8-eps and x[0] < 3.2+eps and x[1] > 0.7-eps 
                  and x[1] < 1.3+eps))
    
    def g_D(x,y):
        return np.cos(2*np.pi*x/5)*np.sin(2*np.pi*y/5)  
    def g_N_1(x,y):
        return 0*np.ones_like(x)
    def g_N_2(x,y):
        return np.cos(2*np.pi*x/5)*np.sin(2*np.pi*y/5)
    def g_R(x,y):
        return np.cos(2*np.pi*x/5)*np.sin(2*np.pi*y/5)
    def f(x,y):
        return np.cos(2*np.pi*x/5)*np.sin(2*np.pi*y/5)
    def a(x=None):
        return 1  
    # def kappa(x=None):
     #   return 1
    kappa =1
    qr = "midpoint_2d"
    E_N = extract_edges(P, T, inhomog_neumann_boundary)
    E_R = extract_edges(P, T, on_robin_boundary)

    dirichlet_nodes = extract_nodes(P, on_dirichlet_boundary)
    X_dc = x[dirichlet_nodes]
    Y_dc = y[dirichlet_nodes]
    g_D_values = g_D(X_dc, Y_dc)
    
    robin_nodes = extract_nodes(P,on_robin_boundary)
    X_rb = x[robin_nodes]
    Y_rb = y[robin_nodes]
    g_R_values = g_R(X_rb, Y_rb)
    
    A = assemble_stiffness_matrix_modified(P, T, a, E_R)
    b = assemble_load_vector_modified_withbc_all(P, T, f, qr, E_N, g_N, E_R, g_R)
    apply_bcs(A, b, dirichlet_nodes, g_D_values)
    A_csr = A.tocsr() 
    U = spsolve(A_csr, b)
    print("N =",N)
    pt.plot2D(x, y, U, title='U', triangles=T, block=True)
```
