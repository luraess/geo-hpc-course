# geo-hpc-course
Parallel CPU and GPU high-performance computing - crash course
----
This crash course aims at providing an interactive and applied approach in an hands-on format to parallel and high-performance computing in Julia. This crash course covers trendy areas in modern geocomputing. Seeking at solutions of differential equations requires efficient numerical schemes optimally leveraging modern hardware. These solvers permit to resolve scientific problems that where technically not possible a decade ago.

The goal of this crash course is to offer an interactive and tutorial-like hands-on to solve systems of differential equations in parallel on many-core hardware accelerators such as GPUs using the Julia language. Julia combines high-level language simplicity to low-level language performance. The resulting codes and applications are fast, short and readable \[[1][JuliaCon20a], [2][JuliaCon20b], [3][JuliaCon19]\].


## Content
* [Objectives](#objectives)
* [Pre-requisite](#pre-requisite)
* [Material](#material)
* [Getting started](#getting-started)
* [Course outline](#course-outline)
* [Advanced start](#advanced-start)
* [Further reading](#further-reading)


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

#### These two examples will enable to address the technical objectives of this course.

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
The hands-on format prioritises the _learning-by-doing_ thus not much preliminary knowledge is required. Basic programming skills won't hurt though. The course will build upon the use of the [Julia] programming language. 

#### Performance metric to compare the various code implementations
Majority of stencil based codes as in this course are memory bounded, meaning the limiting factor in performance is the rate at which memory is transferred from and back between the memory and the arithmetic units. The maximal rate at which the memory transfers occur is the memory copy rate, in the order of 50 GB/s for CPUs and about 1 TB/s for modern GPUs. The effective memory throughput metric (T_eff) measures how good an iterative stencil-based algorithm performs in terms of memory throughput, to be compared to the memory copy rate. The T_eff formula reads: `T_eff = (nIO/1e9*nxy*PRECIS)/(time_s/nt)`, where `nIO` is the number of read/write operations (2 for an update rule), `nxy` is the numerical grid resolution, `PRECIS` is the arithmetic precision (8 or 4 bytes per number), `time_s` is the execution time in second to perform `nt` iterations \[[1][JuliaCon20a]].

#### Programming in Julia
On the CPU, multi-threading is made accessible via [Base.Threads] and the environment variable [JULIA_NUM_THREADS] can be used to define the number of cores to use on the CPU, e.g. `export JULIA_NUM_THREADS=2` to enable 2 threads (2 CPU cores). The [CUDA.jl] module permits to launch compute kernels on Nvidia GPUs within Julia. [JuliaGPU] provides further reading and introductory material about GPU ecosystem within [Julia].


## Material
The course material contains some ready-to-run example scripts, _tmp_ scripts that we will complete together and the solution scripts for the _tmp_ scripts.

#### example scripts
The active working directory during the course will be [/scripts/](/scripts/), that contains the example scripts and the _tmp_ scripts we will work on.

#### example scripts
All _tmp_ scripts have the corresponding solution scripts located in [/solutions/](/solutions/)


## Getting started
Clone this repo (or download it otherwise) to run the example [/scripts/](/scripts/) and access the draft [/scripts/](/scripts/) to be completed during the course. Solution or "cheat-sheets" can be found here [/solutions/](/solutions/).

There are two ways of executing a Julia script, from the Julia command window known as the [Julia REPL], or from the terminal shell directly.

To run Julia interactively, start Julia typing
```sh
julia
```
in the shell (or Terminal). Then, in the [Julia REPL] type
```julia-repl
include("<my_script>.jl")
```
to execute the script. Note that typing `;` in the [Julia REPL] permits you to execute shell commands (like `cd ..`).

For optimal performance (like measuring T_eff), it is more optimal to run Julia as executable from the shell directly, using the optimisation flag `-O3` and disabling bound checking `--check-bounds=no` as following
```sh
julia -O3 --check-bounds=no <my_script>.jl
```
Note that plotting may fail then.


## Course outline
During the course, we will go through the following steps:

_part 1_
1. Intro 1
2. **TODO** Finalise the 1D heat diffusion code [/scripts/heat_1D_tmp.jl](/scripts/heat_1D_tmp.jl).
3. See how the diffusion looks like in 2D [/scripts/heat_2D.jl](/scripts/heat_2D.jl).
4. **TODO** Finalise the 2D loop version of the heat code [/scripts/heat_2D_loop_tmp.jl](/scripts/heat_2D_loop_tmp.jl).
5. **TODO** Import the flux and 🔥 loop calculations in the heat code using external "kernel"-like compute functions [/scripts/heat_2D_loop_fun_tmp.jl](/scripts/heat_2D_loop_fun_tmp.jl).
6. See how one can simply use the GPU to perform the 2D heat diffusion calculations [/scripts/heat_2D_gpu.jl](/scripts/heat_2D_gpu.jl).
7. **TODO** The performance "magic"; update the script [/scripts/heat_2D_gpu_fun_tmp.jl](/scripts/heat_2D_gpu_fun_tmp.jl) based on previous knowledge and step (5.).
8. See how steps 5. and 7. can be combined into a single code using ParallelStencil.jl in [/scripts/heat_2D_xpu.jl](/scripts/heat_2D_xpu.jl)
9. Discussion on CPU vs GPU architectures and performance concerns (T_eff). Q&A.

 _part 2_
10. Intro 2
11. **TODO** Based on your acquired experience, finalise the [/scripts/sia_2D_tmp.jl](/scripts/sia_2D_tmp.jl) script to convert the heat diffusion into an ice cap evolution over time.
12. **TODO** Modify the the script from (11.) to have an implicit steady-state solution.
13. Discussion about pseudo-transient solvers, damping and convergence. Q&A.


## Advanced start
Steps already done on the GPU server you are running on (CentOS 8 linux)

Starting in the shell:
```sh
wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.3-linux-x86_64.tar.gz
tar -xzf julia-1.5.3-linux-x86_64.tar.gz
vim ~/.bashrc # Add the line: PATH=~/julia-1.5.3/bin/:$PATH
export JULIA_CUDA_USE_BINARYBUILDER=false
julia
```
Then in Julia:
```julia-repl
julia> ]
(@v1.5) pkg> add Plots
(@v1.5) pkg> add PyPlot
(@v1.5) pkg> add CUDA
julia> using Plots
julia> using CUDA 
```
_Note: ParallelStencil.jl is about to be publicly released and will then be listed here as well._

_Note: [emoji] keyboard can be added to CentOS as following:_
```sh
sudo dnf install google-noto-emoji-color-fonts.noarch
```

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
[emoji]: https://opensource.com/article/19/10/how-type-emoji-linux
