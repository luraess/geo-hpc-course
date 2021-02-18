const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using ImplicitGlobalGrid, Plots, Printf, Statistics
import MPI
# pyplot()
do_viz = false

@parallel function compute_flux!(qx::Data.Array, qy::Data.Array, T::Data.Array, λ::Data.Number, dx::Data.Number, dy::Data.Number)
    @all(qx) = -λ*@d_xi(T)/dx
    @all(qy) = -λ*@d_yi(T)/dy
    return
end

@parallel function update_T!(T::Data.Array, qx::Data.Array, qy::Data.Array, dt::Data.Number, ρCp::Data.Number, dx::Data.Number, dy::Data.Number)
    @inn(T) = @inn(T) - dt/ρCp*(@d_xa(qx)/dx + @d_ya(qy)/dy)
    return
end

@views function heat_2D_xpu()
    # Physics
    lx, ly = 10.0, 10.0
    λ      = 1.0
    ρCp    = 1.0
    nt     = 200
    # Numerics
    nx, ny = 128-1, 128-1
    nout   = 10
    # Derived numerics
    me, dims, nprocs, coords, comm = init_global_grid(nx, ny) # MPI initialisation
    select_device()                                           # select one GPU per MPI local rank (if >1 GPU per node)
    dx, dy = lx/nx_g(), ly/ny_g()                             # cell sizes
    dt     = min(dx^2,dy^2)/ρCp/λ/4.1
    # Array allocation
    qx     = @zeros(nx-1,ny-2)
    qy     = @zeros(nx-2,ny-1)
    # Initial condition
    T      = Data.Array([exp(-(x_g(ix,dx,T)+dx/2 -0.5*lx)^2 -(y_g(iy,dy,T)+dy/2 -0.5*ly)^2) for ix=1:size(T,1), iy=1:size(T,2)])
    # Preparation of visualisation
    if do_viz
        nx_v, ny_v = (nx-2)*dims[1], (ny-2)*dims[2]
        if (nx_v*ny_v*sizeof(Data.Number) > 0.8*Sys.free_memory()) error("Not enough memory for visualization.") end
        T_v   = zeros(nx_v, ny_v) # global array for visu
        T_inn = zeros(nx-2, ny-2) # no halo local array for visu
        Xi_g, Yi_g = LinRange(dx+dx/2, lx-dx-dx/2, nx_v), LinRange(dy+dy/2, ly-dy-dy/2, ny_v) # inner points only
    end
    # Time loop
    for it = 1:nt
        if (it==11) global t0 = Base.time() end
        @parallel compute_flux!(qx, qy, T, λ, dx, dy)
        @hide_communication b_width begin # communication/computation overlap
            @parallel update_T!(T, qx, qy, dt, ρCp, dx, dy)
            update_halo!(T)
        end
        if mod(it,nout)==0 && do_viz
            T_inn .= T[2:end-1,2:end-1]; gather!(T_inn, T_v)
            if (me==0) display(heatmap(Xi_g, Yi_g, Array(T)', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.))) end
            # sleep(.01)
        end
    end
    time_s = (Base.time()-t0)
    if (me==0) @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/(nt-10)), sigdigits=2)) end
    finalize_global_grid()
    return
end

@time heat_2D_xpu()
