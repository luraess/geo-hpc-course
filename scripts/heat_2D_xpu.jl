const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using Plots, Printf, Statistics
pyplot()
viz = false

@parallel function compute_flux!(qx::Data.Array, qy::Data.Array, T::Data.Array, λ::Data.Number, dx::Data.Number, dy::Data.Number)
    @inn_x(qx) = -λ*@d_xa(T)/dx
    @inn_y(qy) = -λ*@d_ya(T)/dy
    return
end

@parallel function update_T!(T::Data.Array, qx::Data.Array, qy::Data.Array, dt::Data.Number, ρCp::Data.Number, dx::Data.Number, dy::Data.Number)
    @all(T) = @all(T) - dt/ρCp*(@d_xa(qx)/dx + @d_ya(qy)/dy)
    return
end

@views function heat_2D_xpu()
    # physics
    lx     = 10.0
    ly     = 10.0
    λ      = 1.0
    ρCp    = 1.0
    nt     = 200
    # numerics
    nx     = 128-1
    ny     = 128-1
    nout   = 10
    dx     = lx/nx
    dy     = ly/ny
    xc     = LinRange(dx/2, lx-dx/2, nx)
    yc     = LinRange(dy/2, ly-dy/2, ny)
    T      = @zeros(nx  ,ny  )
    qx     = @zeros(nx+1,ny  )
    qy     = @zeros(nx  ,ny+1)
    T      = Data.Array( exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2) )
    dt     = min(dx^2,dy^2)/ρCp/λ/4.1
    # action
    t0     = Base.time()
    for it = 1:nt
        @parallel compute_flux!(qx, qy, T, λ, dx, dy)
        @parallel update_T!(T, qx, qy, dt, ρCp, dx, dy)
        if mod(it,nout)==0 && viz
            display(heatmap(xc, yc, Array(T)', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
            # sleep(.01)
        end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/nt), sigdigits=2))
    return
end

@time heat_2D_xpu()
