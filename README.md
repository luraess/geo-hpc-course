# geo-HPC course

**Parallel CPU and GPU high-performance computing course**

This short course aims at providing an interactive and applied approach in an hands-on format to parallel and high-performance computing in Julia. This short course covers trendy areas in modern geocomputing. Seeking at solutions of differential equations requires efficient numerical schemes optimally leveraging modern hardware. These solvers permit to resolve scientific problems that where technically not possible a decade ago.

The goal of this short course is to offer an interactive and tutorial-like hands-on to solve systems of differential equations in parallel on many-core hardware accelerators such as GPUs using the Julia language. Julia combines high-level language simplicity to low-level language performance. The resulting codes and applications are fast, short and readable \[[1][JuliaCon20a], [2][JuliaCon20b], [3][JuliaCon19]\].

The iterative algorithms can be converted into efficient linear and non-linear solvers relying on a second order Richardson type of iterations strategy \[[4][Frankel50]\].


# Content
* [Objectives](#objectives)
* [Pre-requisite](#pre-requisite)
* [Material](#material)
* [Getting started](#getting-started)
* [Course outline](#course-outline)
* [Advanced start](#advanced-start)
* [Extras](#extras)
* [Further reading](#further-reading)


# Objectives
We will design and implement an iterative numerical algorithm that resolves (non-linear) diffusion in 2D for two applications:

1. The diffusion of heat:
```julia
dT/dt	= 1/ρCp*(-dqx/dx -dqy/dy)
qx     	= -λ*dT/dx
qy     	= -λ*dT/dy
```
For an initial Gaussian distribution, the heat diffusion code produces following output:

![heat diffusion 2D](/docs/heat_2D.gif)

2. The non-linear diffusion of ice topography (simplified shallow-ice):
```julia
dH/dt	= -dqx/dx -dqy/dy + b
qx     	= -H^n*dH/dx
qy     	= -H^n*dH/dy
```
For an initial Gaussian distribution of ice and a circular and centred source/sink term, the simplified shallow-ice code produces following output while reaching a steady state:

![sia non-linear diffusion 2D](/docs/sia_2D.gif)

**These two examples will enable to address the technical objectives of this course (3 parts).**

## Technical objectives of this course (3 parts)

**_Part 1_** | We will use (1) as playground to address:
- vectorised plain Julia implementation _CPU_ (idem as python, Matlab, Octave)
- vectorised versus loop versions on _CPU_
- "kernel"-style loops and multi-threading _multi-core CPU_
- vectorised plain Julia with _GPU_ backend (similar to abstract GPU functions for e.g. python, Matlab)
- explicit "kernel"-style _GPU_ (Julia's power: C CUDA low-level handles)
- enabling single/multi-XPU (both _multi-core CPU_ and _GPU_) using [ParallelStencil.jl] 

**_Part 2_** | We will use (2) as playground to address:
- explicit time stepping
- transient, steady-state solutions
- explicit vs implicit solutions

**_Part 3_** | We will use (1) and (2) to address:
- distributed memory parallelisation ("fake" and "real" parallelisation)
- local and global domain, internal and global boundaries, initial condition, boundary update, synchronisation
- visualisation and I/O
- message passing interface (MPI), MPI + GPU (CUDA-aware MPI)
- communication/computation overlap (hide communication)
- using [ImplicitGlobalGrid.jl] for high-level implementation along with [ParallelStencil.jl]


# Pre-requisite
The hands-on format prioritises the _learning-by-doing_ thus not much preliminary knowledge is required. Basic programming skills won't hurt though. The course will build upon the use of the [Julia] programming language. 

## Performance metric to compare the various code implementations
Majority of stencil based codes as in this course are memory bounded, meaning the limiting factor in performance is the rate at which memory is transferred from and back between the memory and the arithmetic units. The maximal rate at which the memory transfers occur is the memory copy rate, in the order of 50 GB/s for CPUs and about 1 TB/s for modern GPUs. The effective memory throughput metric (T_eff) measures how good an iterative stencil-based algorithm performs in terms of memory throughput, to be compared to the memory copy rate. The T_eff formula reads: `T_eff = (nIO/1e9*nxy*PRECIS)/(time_s/nt)`, where `nIO` is the number of read/write operations (2 for an update rule), `nxy` is the numerical grid resolution, `PRECIS` is the arithmetic precision (8 or 4 bytes per number), `time_s` is the execution time in second to perform `nt` iterations \[[1][JuliaCon20a]].

## Programming in Julia
On the CPU, multi-threading is made accessible via [Base.Threads] and the environment variable [JULIA_NUM_THREADS] can be used to define the number of cores to use on the CPU, e.g. `export JULIA_NUM_THREADS=2` to enable 2 threads (2 CPU cores). The [CUDA.jl] module permits to launch compute kernels on Nvidia GPUs within Julia. [JuliaGPU] provides further reading and introductory material about GPU ecosystem within [Julia].


# Material
The course material contains some ready-to-run _example_ scripts, draft _tmp_ scripts to be complete as tasks during the course and their corresponding _solution_ scripts.

- **_example scripts_ |** The active working directory for the course will be [/scripts/](/scripts/), that contains the example scripts and the _tmp_ scripts to work on.
- **_solution scripts_ |** All _tmp_ scripts have their corresponding solution scripts located in [/solutions/](/solutions/)


# Getting started
If it applies, follow the instructions provided on the course's private channel. 

## Julia quick start
In general, clone this repo (or download it otherwise) to run the example [/scripts/](/scripts/) and access the draft [/scripts/](/scripts/) to be completed during the course. Solutions or "cheat-sheets" can be found in the [/solutions/](/solutions/) folder. The examples rely on 3 main Julia modules, `Plots.jl` (and `PyPlot.jl`) and `CUDA.jl`. The XPU examples require [ParallelStencil.jl] to be installed. The MPI examples require `MPI.jl` to be installed. The multi-XPU scripts require [ImplicitGlobalGrid.jl] to be installed.

There are two ways of executing a Julia script, from the Julia command window known as the [Julia REPL], or from the terminal shell directly. The MPI and multi-XPU examples need to be executed from the terminal shell.

To run Julia interactively, start Julia from the shell (or Terminal). Go to the `geo-hpc-course` folder. Then start Julia appending the `--project` flag to gain access to the required modules:
```sh
$ cd <path-to>/geo-hpc-course/
$ julia --project
```
Then, in the [Julia REPL], you can execute the script as following:
```julia-repl
julia> include("<my_script>.jl")
```
> 💡 Note that typing `;` in the [Julia REPL] permits you to execute shell commands (like `cd ..`).

For optimal performance (like measuring T_eff) and for running Julia MPI, run Julia as executable from the shell directly, using the optimisation flag `-O3` and disabling bound checking `--check-bounds=no` as following:
```sh
$ julia --project -O3 --check-bounds=no <my_script>.jl
```
Note that interactive plotting may fail then.

Set the default `viz = false` flag to `true` if you want to plot output in all codes beyond step 2.

## Running Julia MPI
This section is about launching a Julia MPI script. For [MPI.jl] install notes, refer to the [Advanced start - Julia MPI](#julia-mpi) section and the [MPI.jl] doc. In the proposed approach, each MPI process will handle one CPU thread. In the MPI GPU case (multi-GPUs), each MPI process handles one GPU.

Assuming a working Julia MPI installation, a Julia MPI program can be launched using the Julia MPI wrapper `mpiexecjl` (located in `~/.julia/bin`).

Running the Julia MPI [/scripts/hello_mpi.jl](/scripts/hello_mpi.jl) script on 4 processes can be achieved following:
```sh
$ mpiexecjl -n 4 julia --project scripts/hello_mpi.jl
$ Hello world, I am 0 of 3
$ Hello world, I am 1 of 3
$ Hello world, I am 2 of 3
$ Hello world, I am 3 of 3
```

The 2D Julia MPI diffusion script [/solutions/heat_2D_mpi.jl](/solutions/heat_2D_mpi.jl) executed on 4 MPI processes (global grid: 2x2) produces the following output (see [Extras](#extras) for infos about the gif-making scripts).

![heat diffusion 2D](/docs/heat_2D_mpi_4procs.gif)

_Note: The presented concise Julia MPI scripts are inspired from [this 2D python script](https://github.com/omlins/adios2-tutorial/blob/main/example/mpi_diffusion2D.py)._

Advanced documentation on running the multi-XPU codes can be found in the [ParallelStencil.jl module documentation](https://github.com/omlins/ParallelStencil.jl#miniapp-content).


# Course outline
During the course, we will go through the following steps:

1. **Intro _part 1_**
2. **TODO** Finalise the 1D heat diffusion code [/scripts/heat_1D_tmp.jl](/scripts/heat_1D_tmp.jl), implementing the equations from [Objectives 1.](#objectives)
3. See how the diffusion looks like in 2D [/scripts/heat_2D.jl](/scripts/heat_2D.jl).
4. **TODO** Finalise the 2D loop version of the heat code [/scripts/heat_2D_loop_tmp.jl](/scripts/heat_2D_loop_tmp.jl).
5. **TODO** Import the flux and `T` loop calculations in the heat code using external "kernel"-like compute functions [/scripts/heat_2D_loop_fun_tmp.jl](/scripts/heat_2D_loop_fun_tmp.jl).
6. See how the function-based loop version (5.) can be further extended to checking bounds with `if` statement in the `ix` and `iy` loops, including "loop-fusion" for the flux computations [/scripts/heat_2D_loop_fun_gpustyle.jl](/scripts/heat_2D_loop_fun_gpustyle.jl).
7. See how one can simply use the **GPU** to perform the 2D heat diffusion calculations [/scripts/heat_2D_gpu.jl](/scripts/heat_2D_gpu.jl).
8. **TODO** The performance "magic"; update the script [/scripts/heat_2D_gpu_fun_tmp.jl](/scripts/heat_2D_gpu_fun_tmp.jl) based on previous knowledge and step (5.).
9. See how steps 5., 6. and 7. can be combined into a single code using [ParallelStencil.jl] in [/scripts/heat_2D_xpu.jl](/scripts/heat_2D_xpu.jl)
10. Discussion on CPU vs GPU architectures and performance evaluation (T_eff). Q&A.

---

11. **Intro _part 2_**
12. **TODO** Based on your acquired experience, finalise the [/scripts/sia_2D_tmp.jl](/scripts/sia_2D_tmp.jl) script to convert the heat diffusion `T` into an ice cap thickness `H` evolution over time.
13. **TODO** Modify the the script from (12.) to have an implicit solver while reaching a steady-state solution.
14. **TODO** (NEW!) You demystified GPU computing with completing step 9. Update now the script [/scripts/sia_2D_xpu_tmp.jl](/scripts/sia_2D_xpu_tmp.jl) to have an XPU (CPU or GPU) code ready!
15. Discussion about pseudo-transient solvers, damping and convergence. Q&A.

---

16. **Intro _part 3_ (NEW!)**
17. **TODO** Finalise the 1D heat diffusion code [/scripts/heat_1D_2procs_tmp.jl](/scripts/heat_1D_2procs_tmp.jl), splitting the calculation of temperature evolution on one left and one right domain. This "fake-parallelisation" requires left and right temperature arrays, `TL` and `TR`, respectively.
18. **TODO** Generalise the "fake-parallel" approach on 2 processes to `n` processes by modifying the code [/scripts/heat_1D_nprocs_tmp.jl](/scripts/heat_1D_nprocs_tmp.jl), taking care of implementing the initial condition, heat diffusion physics and the boundary update.
19. The script [/scripts/hello_mpi.jl](/scripts/hello_mpi.jl) shows a "Hello World" example implementing [MPI.jl]. Use this script to test your [MPI.jl] install (see the [Running Julia MPI](#running-julia-mpi) section for more infos on installing and running Julia MPI).
20. Discover a concise MPI 1D heat diffusion example [/scripts/heat_1D_mpi.jl](/scripts/heat_1D_mpi.jl). Learn about the minimal requirements to initialise a Cartesian MPI topology and how to code the boundary update functions (here using blocking messages). Use the [/scripts/vizme1D_mpi.jl](/scripts/vizme1D_mpi.jl) script to visualise the results (each MPI process saving it's local output).
21. **TODO** Yay, you have your MPI 1D Julia script running! Finalise the MPI 2D heat diffusion script [/scripts/heat_2D_mpi_tmp.jl](/scripts/heat_2D_mpi_tmp.jl) to solve the 2D diffusion equation using MPI. Use the [/scripts/vizme2D_mpi.jl](/scripts/vizme2D_mpi.jl) script to visualise the results (each MPI process saving it's local output).
22. Now that you demystified distributed memory parallelisation, see how using [ImplicitGlobalGrid.jl] along with [ParallelStencil.jl] leads to concise and efficient distributed memory parallelisation on multiple _XPUs_ in 2D [/scripts/heat_2D_multixpu.jl](/scripts/heat_2D_multixpu.jl). Also, take a closer look at the [@hide_communication](https://github.com/luraess/geo-hpc-course/blob/0a722ac5f6da47779dfceadfec79b92c95e9e40e/scripts/heat_2D_multixpu.jl#L61) feature. Further infos can be found [here](https://github.com/omlins/ParallelStencil.jl#seamless-interoperability-with-communication-packages-and-hiding-communication).
23. **TODO** Instrument the 2D shallow ice code sia_2D_xpu.jl (task 14.) to enable distributed memory parallelisation using [ImplicitGlobalGrid.jl] along with [ParallelStencil.jl].
> 💡 Use [/solutions/sia_2D_xpu.jl](/solutions/sia_2D_xpu.jl) for a quick start, and [/solutions/sia_2D_multixpu.jl](/solutions/sia_2D_multixpu.jl) for a solution.
24. Yay, you made it - you demystified running Julia codes in parallel on multi-XPU :-) Q&A.


# Advanced start
Steps already done on the GPU server you are running on (CentOS 8 linux)

## Julia
Starting in the shell:
```sh
$ wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.4-linux-x86_64.tar.gz
$ tar -xzf julia-1.5.4-linux-x86_64.tar.gz
$ vim ~/.bashrc # Add the line: PATH=~/julia-1.5.4/bin/:$PATH
$ export JULIA_CUDA_USE_BINARYBUILDER=false
$ cd <path-to>/geo-hpc-course/
$ julia --project
```
Once Julia launched, enter the package mode `]` and `instantiate` the project. This will download all the project's dependencies and install them:
```julia-repl
julia> ]

(geo-hpc-course) pkg> instantiate

julia>
```

If your GPU system contains more than one GPU, you can add following at the beginning of each `gpu` named code to target a specific device identified by its unique `ID` (default being `0`):
```julia
GPU_ID = ID
CUDA.device!(GPU_ID)
```

## Julia MPI
The following steps permit you to install [MPI.jl] on your machine:
1. Add `MPI.jl`:
```julia-repl
(project) pkg> add MPI

julia> using MPI
[ Info: Precompiling MPI [da04e1cc-30fd-572f-bb4f-1f8673147195]

julia> MPI.install_mpiexecjl()
[ Info: Installing `mpiexecjl` to `HOME/.julia/bin`...
[ Info: Done!
```
2. Then, one should add `HOME/.julia/bin` to PATH in order to launch the Julia MPI wrapper `mpiexecjl`.

3. Running the Julia MPI code on 3 processes:
```sh
$ HOME/.julia/bin/mpiexecjl -n 4 julia --project scripts/hello_mpi.jl
```
> 💡 Note: On MacOS, there seems to be an issue (https://github.com/JuliaParallel/MPI.jl/issues/407). To fix it, define following `ENV` variable:
```sh
$ export MPICH_INTERFACE_HOSTNAME=localhost
```
> and add `-host localhost` to the execution script:
```sh
$ HOME/.julia/bin/mpiexecjl -n 3 -host localhost julia --project scripts/hello_mpi.jl
```


# Extras
[Julia] supports UTF-8 (Unicode) characters. Also, the plotting package [Plots.jl] permits to create gif animation out-of-the-box. The [/extras/heat_2D_gif_unicode.jl](/extras/heat_2D_gif_unicode.jl) exemplifies these two fun capabilities.

The code [/extras/sia_2D_ss_gif.jl](/extras/sia_2D_ss_gif.jl) uses the out-of-the-box gif-making capabilities to produce the SIA non-linear diffusion gif.

The code [/extras/heat_2D_mpi_gif.jl](/extras/heat_2D_mpi_gif.jl) produces time-dependent output to be visualised as a gif using the [/extras/vizme2D_mpi_gif.jl](/extras/vizme2D_mpi_gif.jl) script.

> 💡 Note: On Linux machines, [emoji] keyboard may need to be installed in order to display the Unicode emoticons.
```sh
$ sudo dnf install google-noto-emoji-color-fonts.noarch
```


# Further reading
\[1\] [Omlin, S., Räss, L., Kwasniewski, G., Malvoisin, B., & Podladchikov, Y. Y. (2020). Solving Nonlinear Multi-Physics on GPU Supercomputers with Julia. JuliaCon Conference, virtual.][JuliaCon20a]

\[2\] [Räss, L., Reuber, G., Omlin, S. (2020). Multi-Physics 3-D Inversion on GPU Supercomputers with Julia. JuliaCon Conference, virtual.][JuliaCon20b]

\[3\] [Räss, L., Omlin, S., & Podladchikov, Y. Y. (2019). Porting a Massively Parallel Multi-GPU Application to Julia: a 3-D Nonlinear Multi-Physics Flow Solver. JuliaCon Conference, Baltimore, USA.][JuliaCon19]

\[4\] [Frankel, S. P. (1950). Convergence rates of iterative treatments of partial differential equations, Mathe. Tables Other Aids Comput., 4, 65–75.][Frankel50]


[JuliaCon20a]: https://www.youtube.com/watch?v=vPsfZUqI4_0
[JuliaCon20b]: https://www.youtube.com/watch?v=1t1AKnnGRqA
[JuliaCon19]: https://www.youtube.com/watch?v=b90qqbYJ58Q
[Frankel50]: /docs/frankel_1950.pdf
[Julia]: https://julialang.org
[Base.Threads]: https://docs.julialang.org/en/v1/base/multi-threading/
[JULIA_NUM_THREADS]:https://docs.julialang.org/en/v1.0.0/manual/environment-variables/#JULIA_NUM_THREADS-1
[CUDA.jl]: https://github.com/JuliaGPU/CUDA.jl
[Plots.jl]: https://github.com/JuliaPlots/Plots.jl
[ParallelStencil.jl]: https://github.com/omlins/ParallelStencil.jl
[ImplicitGlobalGrid.jl]: https://github.com/eth-cscs/ImplicitGlobalGrid.jl
[MPI.jl]: https://juliaparallel.github.io/MPI.jl/stable/examples/01-hello/
[Julia REPL]: https://docs.julialang.org/en/v1/stdlib/REPL/
[JuliaGPU]: https://juliagpu.org
[emoji]: https://opensource.com/article/19/10/how-type-emoji-linux
