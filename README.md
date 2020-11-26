# geo-hpc-course
Parallel CPU and GPU high-performance computing - crash course

## Description
This crash course aims at providing an interactive and applied approach in an hands-on format to parallel and high-performance computing in Julia. This crash course covers trendy areas in modern geocomputing. Seeking at solutions of differential equations requires efficient numerical schemes optimally leveraging modern hardware. These solvers permit to resolve scientific problems that where technically not possible a decade ago.

The goal of this crash course is to offer an interactive and tutorial-like hands-on to solve systems of differential equations in parallel on many-core hardware accelerators such as GPUs using the Julia language. Julia combines high-level language simplicity to low-level language performance. The resulting codes and applications are fast, short and readable \[[1][JuliaCon20a], [2][JuliaCon20b], [3][JuliaCon19]\].

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

## Pre-requisite
_... work in progress ..._
The hands-on format prioritises the _learning-by-doing_ thus not much preliminary knowledge is required. Basic programming skills won't hurt though. The course will build upon the use of the [Julia] programming language. 

#### Performance metric to compare the various code implementations
Majority of stencil based codes as in this course are memory bounded, meaning the limiting factor in performance is the rate at which memory is transferred from and back between the memory and the arithmetic units.

#### Technical side
On the CPU, multi-threading is made accessible via [Base.Threads] and the environment variable [JULIA_NUM_THREADS] can be used to define the number of cores to use on the CPU, e.g. `export JULIA_NUM_THREADS=2` to enable 2 threads (2 CPU cores). The [CUDA.jl] module permits to launch compute kernels on Nvidia GPUs within Julia. [JuliaGPU] provides further reading and introductory material about GPU ecosystem within [Julia].
- 

## Material
_... work in progress ..._

## Get started
_... work in progress ..._


## Further reading
\[1\] [Omlin, S., Räss, L., Kwasniewski, G., Malvoisin, B., & Podladchikov, Y. Y. (2020). Solving Nonlinear Multi-Physics on GPU Supercomputers with Julia. JuliaCon Conference, virtual.][JuliaCon20a]

\[2\] [Räss, L., Reuber, G., Omlin, S. (2020). Multi-Physics 3-D Inversion on GPU Supercomputers with Julia. JuliaCon Conference, virtual.][JuliaCon20b]

\[3\] [Räss, L., Omlin, S., & Podladchikov, Y. Y. (2019). Porting a Massively Parallel Multi-GPU Application to Julia: a 3-D Nonlinear Multi-Physics Flow Solver. JuliaCon Conference, Baltimore, USA.][JuliaCon19]


### Contact
Ludovic Räss (ludovic.rass@gmail.com)

[JuliaCon20a]: https://www.youtube.com/watch?v=vPsfZUqI4_0
[JuliaCon20b]: https://www.youtube.com/watch?v=1t1AKnnGRqA
[JuliaCon19]: https://www.youtube.com/watch?v=b90qqbYJ58Q
[Julia]: https://julialang.org
[Base.Threads]: https://docs.julialang.org/en/v1/base/multi-threading/
[JULIA_NUM_THREADS]:https://docs.julialang.org/en/v1.0.0/manual/environment-variables/#JULIA_NUM_THREADS-1
[CUDA.jl]: https://github.com/JuliaGPU/CUDA.jl
[Julia REPL]: https://docs.julialang.org/en/v1/stdlib/REPL/
[JuliaGPU]: https://juliagpu.org
