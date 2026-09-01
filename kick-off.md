---
jupytext:
  formats: md:myst,ipynb
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.5
kernelspec:
  display_name: tma4220-finite-elements (3.14.7)
  language: python
  name: python3
---

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

# TMA4220 - Numerical Solution of Partial Differential Equations using Finite Element Methods
(sec:kick-off-meeting)=
## Kick-off Meeting

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

## Course content and objectives

* The course gives an introduction to **finite element methods** for the numerical solution of partial differential equations (PDEs). 

* We will discuss the **strong** and **weak formulation** of prototypical PDEs including a brief review of the functional analysis related tools like Hilbert spaces, Sobolev spaces, Lax-Milgram.  

* For prototypical PDE problems, we will discuss the typical **4-step recipe** to go from PDE to a numerically solvable linear system.

* During this course we will have an in-depth discussion of **theoretical analysis** of the finite element methods we develop for our prototype PDEs including stability and error estimates.

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

* We will also discuss the **efficient implementation and solution** of the finite element based solution schemes.

* For time-dependent problems, we discuss how to combine the finite element with various **time discretization** approaches such as one-step or multi-step methods.

* The **treatment of general boundary conditions** and **interface conditions** will be discussed.

* If time permits, we will consider some examples of **nonlinear problems** as well.

To this end we will develop, theoretically analyze, and practically implement numerical methods for the following *prototype problems*:

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

### Poisson problem 
has a wide range of applications in, e.g. thermodynamics, Newtonian gravity and electrostatics.
In its simplest form, it takes the form

$$
\begin{alignedat}{2}
-\Delta u &= f &&\quad \text{in } \Omega,
\\
u &= g &&\quad \text{on } \partial\Omega,
\end{alignedat} 
$$

consisting of **PDE** part and a **boundary condition** part. Throughout the course we will discuss different types of boundary conditions
and how they can be incorporated into the weak formulation and numerical FEM-schemes.

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

### Parabolic equation
is a time-dependent variant of the Poisson problem

$$
\begin{alignedat}{2}
\partial_t u(x,t) -\Delta u(x,t) &= f(x,t) &&\quad \text{for } (x,t) \in \Omega \times (0,T],
\\
u(x,t) &= g(x,t) &&\quad \text{for } (x,t) \in \partial\Omega \times (0,T],
\\
u(x,0) &= u_0(x) &&\quad \text{for } x \in \Omega
\end{alignedat} 
$$

consisting of the PDE, boundary conditions and **initial conditions**

If time permits, we will also have a look at how to numerically solve nonlinear time-dependent problems such as

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

### Gray-Scott model
describes an autocatalytic chemical reaction
between two species $u$ and $v$ diffusing in a domain $\Omega$. It is
a classical example of a reaction-diffusion system that produces
Turing patterns (spots, stripes, self-replicating spots) purely from
diffusion and local reaction kinetics, and serves as a good nonlinear,
time-dependent test problem for the finite element method

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "fragment"}}

Find $u, v : \Omega \times (0, T] \to \mathbb{R}$ such that

$$
\begin{alignedat}{2}
\partial_t u &= D_u \Delta u - u v^2 + F(1 - u) &&\quad \text{in } \Omega \times (0, T],
\\
\partial_t v &= D_v \Delta v + u v^2 - (F + k) v &&\quad \text{in } \Omega \times (0, T],
\end{alignedat}
$$ (eq:gray-scott)

which must be supplemented with boundary conditions (b.c.) and initial conditions (i.c.)

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "fragment"}}

Here $D_u, D_v > 0$ are diffusion coefficients for $u$ and $v$, $F$ is
the feed rate at which $u$ is replenished, and $k$ is the rate
constant for the removal of $v$. The term $u v^2$ models the
autocatalytic reaction $u + 2v \to 3v$, in which $v$ converts $u$ into
more of itself. 

Tuning $F$ and $k$ results in different patterns (dots, stripes, coral).

```{code-cell} ipython3
---
deletable: true
editable: true
slideshow:
  slide_type: slide
---
from IPython.display import YouTubeVideo, HTML
YouTubeVideo('nw2bPnhtxN8', width=800, height=500)
```

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

## Exercises

* As part of (mostly handwritten) lectures notes and Jupyter-based learning material, I will distribute a number of problem sets containing both  **theoretical and implementation tasks** to help you digest the material. You don't have to submit solutions for those.

* **We won't have a weekly exercise session**, but we will schedule regular exercise sessions when we are close to the conclusion of topical unit.

* **Group work** is very much encouraged during the exercise sessions!

+++ {"deletable": true, "editable": true, "jp-MarkdownHeadingCollapsed": true, "slideshow": {"slide_type": "slide"}}

* For **programming exercises**, implementation guides will be given when necessary, and you can use the exercise sessions to work on the implementation while I am hanging around :)

* For a typical **theoretical exercise** meet-up, the idea is that each group picks a problem to work on the first 45 min. Then each group presents their ideas for the solution after the break in the second 45 min (even if you won't be able to completely solve the exercise problem!). Presentation is supposed to be very informal at the black/whiteboard, and open discussion is very much encouraged!

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

## Passing and Grading

* One final grade composed of the exam grade (75%) and the project grade (25%)
* Project **is not mandatory** but you lose 25% of your points -> you can at most get a $C$ (if you get an $A$ in the exam!).

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

## Project
* Group work for the project is highly encouraged and project reports can be submitted **in groups up to 3 persons**.

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

* Project work will consist of 3 phases
  * Phase 1: Initial implementation of the project and report writing (~2 weeks) with a first deadline.

  * Phase 2: Peer review: Each group will assess and review another group's project, write a short assessment/feedback summary (~1 page) and present it to/discuss it with that group. This will be done during the exercise sessions. I will provide some assessment guidelines.
    

  * Phase 3: After receiving the feedback, each group has time to improve their reports and incorporate the given feedback before final submission. If you don't agree with certain points in the received feedback, you can account for this in the report.

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

* Each project submission consists of
  1. the final project report (90 %)
  2. the assessment summary provided to another group (10%)
  3. the assessment summary received from another group

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

## Timeslots for lectures and tutorials

* Tuesday 8:15-10:00
* Wednesday 14:15-16:00
* Thursday 8:15-10:00

We will settle for the lecture times (or possible alternative times) during a menti.

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

## Reading material

A digital reading list consisting of book chapters from multiple sources, handwritten and jupyter-based lectures
will be provided at
<a href="https://wiki.math.ntnu.no/tma4220/2026h/start">wiki page</a>.

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

## Software

* In the beginning, we will use the traditional Python software stack for scientific computing, i.e. numpy, scipy, matplotlib etc.
* When we have understood the algorithmic realization of the finite element methods, we will switch to the open source finite element library
* <a href="https://ngsolve.org/">ngsolve</a> which provides a high-level Python interface to solve complex PDE problems

+++ {"deletable": true, "editable": true, "slideshow": {"slide_type": "slide"}}

**Any questions?**

## Mentimeter time!
