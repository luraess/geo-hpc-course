const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using Plots, Printf, Statistics
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
    nx, ny = 128, 128
    nout   = 10
    # Derived numerics
    dx, dy = lx/nx, ly/ny
    dt     = min(dx^2,dy^2)/ρCp/λ/4.1
    xc     = LinRange(dx/2, lx-dx/2, nx)
    yc     = LinRange(dy/2, ly-dy/2, ny)
    # Array allocation
    qx     = @zeros(nx-1,ny-2)
    qy     = @zeros(nx-2,ny-1)
    # Initial condition
    T      = Data.Array( exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2) )
    # Time loop
    for it = 1:nt
        if (it==11) global t0 = Base.time() end
        @parallel compute_flux!(qx, qy, T, λ, dx, dy)
        @parallel update_T!(T, qx, qy, dt, ρCp, dx, dy)
        if mod(it,nout)==0 && do_viz
            display(heatmap(xc, yc, Array(T)', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
            # sleep(.01)
        end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/(nt-10)), sigdigits=2))
    return
end

@time heat_2D_xpu()
