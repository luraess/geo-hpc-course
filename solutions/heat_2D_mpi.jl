# run: ~/.julia/bin/mpiexecjl -n 4 julia --project solutions/heat_2D_mpi.jl
using Plots, Printf, MAT
import MPI

do_save = true

# MPI functions
@views function update_halo(A, neighbors_x, neighbors_y, comm)
    if neighbors_x[1] >= 0 # MPI_PROC_NULL?
        sendbuf = A[2,:]
        recvbuf = zeros(size(A[1,:]))
        MPI.Send(sendbuf,  neighbors_x[1], 0, comm)
        MPI.Recv!(recvbuf, neighbors_x[1], 1, comm)
        A[1,:] = recvbuf
    end
    if neighbors_x[2] >= 0 # MPI_PROC_NULL?
        sendbuf = A[end-1,:]
        recvbuf = zeros(size(A[end,:]))
        MPI.Send(sendbuf,  neighbors_x[2], 1, comm)
        MPI.Recv!(recvbuf, neighbors_x[2], 0, comm)
        A[end,:] = recvbuf
    end
    if neighbors_y[1] >= 0 # MPI_PROC_NULL?
        sendbuf = A[:,2]
        recvbuf = zeros(size(A[:,1]))
        MPI.Send(sendbuf,  neighbors_y[1], 2, comm)
        MPI.Recv!(recvbuf, neighbors_y[1], 3, comm)
        A[:,1] = recvbuf
    end
    if neighbors_y[2] >= 0 # MPI_PROC_NULL?
        sendbuf = A[:,end-1]
        recvbuf = zeros(size(A[:,end]))
        MPI.Send(sendbuf,  neighbors_y[2], 3, comm)
        MPI.Recv!(recvbuf, neighbors_y[2], 2, comm)
        A[:,end] = recvbuf
    end
    return
end

@views function heat_2D_mpi()
    # MPI
    MPI.Init()
    dims        = [0,0]
    comm        = MPI.COMM_WORLD
    nprocs      = MPI.Comm_size(comm)
    MPI.Dims_create!(nprocs, dims)
    comm_cart   = MPI.Cart_create(comm, dims, [0,0], 1)
    me          = MPI.Comm_rank(comm_cart)
    coords      = MPI.Cart_coords(comm_cart)
    neighbors_x = MPI.Cart_shift(comm_cart, 0, 1)
    neighbors_y = MPI.Cart_shift(comm_cart, 1, 1)
    if (me==0) println("nprocs=$(nprocs), dims[1]=$(dims[1]), dims[2]=$(dims[2])") end
    # Physics
    lx, ly     = 10.0, 10.0
    λ          = 1.0
    ρCp        = 1.0
    nt         = 100
    # Numerics
    nx, ny     = 32, 32                             # local
    nx_g, ny_g = dims[1]*(nx-2)+2, dims[2]*(ny-2)+2 # global
    # Derived numerics
    dx, dy     = lx/nx_g, ly/ny_g                   # global
    dt         = min(dx,dy)^2/ρCp/λ/4.1
    # Array allocation
    qx         = zeros(nx-1,ny-2)
    qy         = zeros(nx-2,ny-1)
    # Initial condition
    x0, y0     = coords[1]*(nx-2)*dx, coords[2]*(ny-2)*dy
    xc         = [x0 + ix*dx - dx/2 - 0.5*lx  for ix=1:nx]
    yc         = [y0 + iy*dy - dy/2 - 0.5*ly  for iy=1:ny]
    T          = exp.(.-xc.^2 .-yc'.^2)
    # Time loop
    for it = 1:nt
        if (it==11) global t0 = Base.time() end
        qx .= .-λ.*diff(T[:,2:end-1], dims=1)./dx
        qy .= .-λ.*diff(T[2:end-1,:], dims=2)./dy
        T[2:end-1,2:end-1] .= T[2:end-1,2:end-1] .- dt./ρCp.*(diff(qx, dims=1)./dx + diff(qy, dims=2)./dy)
        update_halo(T, neighbors_x, neighbors_y, comm_cart)
    end
    time_s = (Base.time()-t0)
    if (me==0) @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/(nt-10)), sigdigits=2)) end
    # Save to visualise
    if do_save file = matopen("$(@__DIR__)/T_$(me).mat", "w"); write(file, "T", Array(T)); close(file) end
    MPI.Finalize()
    return
end

heat_2D_mpi()
