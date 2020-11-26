# geo-hpc-course
Parallel CPU and GPU high-performance computing - crash course

## Description
This crash course aims at providing an interactive and applied approach in an hands-on format to parallel and high-performance computing in Julia. This crash course covers trendy areas in modern geocomputing. Seeking at solutions of differential equations requires efficient numerical schemes optimally leveraging modern hardware. These solvers permit to resolve scientific problems that where technically not possible a decade ago.

The goal of this crash course is to offer an interactive and tutorial-like hands-on to solve systems of differential equations in parallel on many-core hardware accelerators such as GPUs using the Julia language. Julia combines high-level language simplicity to low-level language performance. The resulting codes and applications are fast, short and readable.

## Objectives
We will design and implement an iterative numerical algorithm that resolves (non-linear) diffusion in 2D for two applications:

1. The diffusion of heat:
```julia
∂🔥/∂t	= 1/ρCp*(-∂qx/∂x -∂qx/∂x)
qx     	= -λ*∂🔥/∂x
qy     	= -λ*∂🔥/∂y
```
For an initial Gaussian distribution, the heat diffusion code produces following output:

![heat diffusion 2D](/docs/heat_2D.gif)

2. The non-linear diffusion of ice topography (simplified shallow-ice):
```julia
∂❄/∂t	= -∂qx/∂x -∂qx/∂x + b
qx     	= -❄^n*∂❄/∂x
qy     	= -❄^n*∂❄/∂y
```
For an initial Gaussian distribution of ice and a circular and centred source/sink term, the simplified shallow-ice code produces following output:

![heat diffusion 2D](/docs/sia_2D_ss.png)

These two examples will enable to address the technical objectives of this course.

We will use (1) as playground to address:
- vectorised plain Julia implementation _CPU_ (idem as python, Matlab, Octave)
- vectorised versus loop versions on _CPU_
- "kernel"-style loops and multi-threading _multi-core CPU_
- vectorised plain Julia with GPU backend _GPU_ (idem as abstract GPU functions for e.g. python, Matlab)
- explicit "kernel"-style _GPU_ (Julia's power: C CUDA low-level handles)
- using ParallelStencil enabling both _multi-core CPU_ and _GPU_

We will use (2) as playground to address:
- explicit time stepping
- transient, steady-state solutions
- explicit vs implicit solutions

## Pre-requisites


