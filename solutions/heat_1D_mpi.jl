using Plots, Printf, MAT
import MPI
do_save = true
# MPI functions
@views function update_halo(A, neighbors_x, comm)
    if neighbors_x[1] >= 0 # MPI_PROC_NULL?
        sendbuf = A[2]
        recvbuf = zeros(1)
        MPI.Send(sendbuf,  neighbors_x[1], 0, comm)
        MPI.Recv!(recvbuf, neighbors_x[1], 1, comm)
        A[1] = recvbuf[1]
    end
    if neighbors_x[2] >= 0 # MPI_PROC_NULL?
        sendbuf = A[end-1]
        recvbuf = zeros(1)
        MPI.Send(sendbuf,  neighbors_x[2], 1, comm)
        MPI.Recv!(recvbuf, neighbors_x[2], 0, comm)
        A[end] = recvbuf[1]
    end
    return
end

@views function heat_1D_mpi()
    # MPI
    dims        = [0]
    MPI.Init()
    comm        = MPI.COMM_WORLD
    nprocs      = MPI.Comm_size(comm)
    MPI.Dims_create!(nprocs, dims)
    comm_cart   = MPI.Cart_create(comm, dims, [0], 1)
    me          = MPI.Comm_rank(comm_cart)
    coords      = MPI.Cart_coords(comm_cart)
    neighbors_x = MPI.Cart_shift(comm_cart, 0, 1)
    if (me==0) println("nprocs=$(nprocs), dims[1]=$(dims[1])") end
    # physics
    lx   = 10.0
    λ    = 1.0
    ρCp  = 1.0
    nt   = 20
    # numerics
    nx   = 31
    nx_g = dims[1]*(nx-2) + 2
    dx   = lx/nx_g # global
    dt   = dx^2/ρCp/λ/2.1
    # array initialisation
    xc   = zeros(nx  )  
    qx   = zeros(nx-1)
    T    = zeros(nx  )
    # initial condition
    x0   = coords[1]*(nx-2)*dx + dx/2
    xc  .= [x0 + ix*dx - 0.5*lx for ix=1:nx]
    T   .= exp.(.-xc.^2)
    # action
    t0  = Base.time()
    for it = 1:nt
        qx         .= .-λ.*diff(T)./dx
        T[2:end-1] .= T[2:end-1] .- dt./ρCp.*diff(qx)./dx
        update_halo(T, neighbors_x, comm_cart)
     end
    time_s = (Base.time()-t0)
    if (me==0) @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*sizeof(lx))/(time_s/nt), sigdigits=2)) end
    if do_save file = matopen("T_$me.mat", "w"); write(file, "T", Array(T)); close(file) end
    MPI.Finalize()
    return
end

heat_1D_mpi()
