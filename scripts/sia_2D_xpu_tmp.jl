const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
    macro pow(args...)  esc(:(CUDA.pow($(args...)))) end
else
    @init_parallel_stencil(Threads, Float64, 2)
    pow(x,y) = x^y
    macro pow(args...)  esc(:(pow($(args...)))) end
end
using Plots, Printf, Statistics
# pyplot()
do_viz = false

@parallel function check_Err!(ErrH::Data.Array, H::Data.Array)
    @all(ErrH) = @all(ErrH) - @all(H)
    return
end

@parallel function compute_Err!(ErrH::Data.Array, H::Data.Array)
    @all(ErrH) = @all(H)
    return
end

@parallel function compute_flux!(# TODO add args)
    # TODO flux calculations qx, qy
    return
end

@parallel function compute_ResH!(# TODO add args)
    # TODO dtau and ResH calculations
    return
end

@parallel function update_H!(# TODO add args)
    # TODO H calculation
    return
end

@views function sia_2D_xpu()
    # Physics
    lx, ly = 10.0, 10.0
    n      = 3.0
    # Numerics
    nx, ny = 128, 128
    niter  = 10000
    nout   = 100
    dmp    = 0.96
    ε      = 1e-8
    # Derived numerics
    dx, dy = lx/nx, ly/ny
    xc     = LinRange(dx/2, lx-dx/2, nx)
    yc     = LinRange(dy/2, ly-dy/2, ny)
    # Array allocation
    qx     = @zeros(nx-1,ny-2)
    qy     = @zeros(nx-2,ny-1)
    dtau   = @zeros(nx-2,ny-2)
    ResH   = @zeros(nx-2,ny-2)
    ErrH   = @zeros(nx  ,ny  )
    # Initial condition
    b      = 0.5*ones(nx,ny)
    H      = Data.Array( exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2) )
    rad    = (xc.-lx./2.0).^2 .+(yc.-ly./2.0)'.^2; b[rad.>lx/4] .= -0.5
    b      = Data.Array(b) # move data back to XPU
    # Time loop
    for iter = 1:niter
        if (iter==11) global t0 = Base.time() end
        @parallel compute_Err!(ErrH, H)
        # TODO complete de function calls
        @parallel compute_flux!(# TODO add args)
        @parallel compute_ResH!(# TODO add args)
        @parallel update_H!(# TODO add args)
        @parallel check_Err!(ErrH, H)
        # Visualise
        if mod(iter,nout)==0 && do_viz 
            p1 = heatmap(xc, yc, Array(H)', xlabel="lx", ylabel="ly", title="shallow ice, iter=$iter", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:viridis)
            p2 = plot(xc, Array(H)[:,Int(round(ny/2))], xlabel="lx", ylabel="height", xlims=(xc[1], xc[end]), ylims=(0., 1.2), legend=false, framestyle=:box, aspect_ratio=2.5)
            l = @layout [ a{0.8h} ; b{0.2h} ]
            display(plot(p1, p2, layout=l ))
        end
        if mod(iter,nout)==0 maxErr=maximum(abs.(ErrH)); @printf("iter=%d, max(err)=%1.2e \n", iter, maxErr); if maxErr<ε global itg=iter; break; end end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/(itg-10)), sigdigits=2))
    return
end

@time sia_2D_xpu()
