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
pyplot()
viz = false

@parallel function check_Err!(ErrH::Data.Array, H::Data.Array)
    @all(ErrH) = @all(ErrH) - @all(H)
    return
end

@parallel function compute_H3_Err!(H3::Data.Array, ErrH::Data.Array, H::Data.Array, n::Data.Number)
    @all(H3)   = @pow( @all(H), n)
    @all(ErrH) = @all(H)
    return
end

@parallel function compute_flux!(qx::Data.Array, qy::Data.Array, H::Data.Array, H3::Data.Array, dx::Data.Number, dy::Data.Number)
    @inn_x(qx) = -@av_xa(H3)*@d_xa(H)/dx
    @inn_y(qy) = -@av_ya(H3)*@d_ya(H)/dy
    return
end

@parallel function compute_ResH!(dtau::Data.Array, ResH::Data.Array, H3::Data.Array, qx::Data.Array, qy::Data.Array, b::Data.Array, dmp::Data.Number, dx::Data.Number, dy::Data.Number)
    @all(dtau) = min(dx*dx, dy*dy)/(1.0 + @all(H3))/4.1./4.0
    @all(ResH) = (-dx*@d_xa(qx)/dx  - @d_ya(qy)/dy + @all(b)) + dmp*@all(ResH)
    return
end

@parallel function update_H!(H::Data.Array, dtau::Data.Array, ResH::Data.Array)
    @all(H) = max(0.0, @all(H) + @all(dtau)*@all(ResH))
    return
end

@views function sia_2D_xpu()
    # physics
    lx     = 10.0
    ly     = 10.0
    n      = 3.0
    # numerics
    nx     = 128-1
    ny     = 128-1
    niter  = 10000
    nout   = 100
    dmp    = 0.96
    ε      = 1e-8
    dx     = lx/nx
    dy     = ly/ny
    xc     = LinRange(dx/2, lx-dx/2, nx)
    yc     = LinRange(dy/2, ly-dy/2, ny)
    H      = @zeros(nx  ,ny  )
    H3     = @zeros(nx  ,ny  )
    qx     = @zeros(nx+1,ny  )
    qy     = @zeros(nx  ,ny+1)
    dtau   = @zeros(nx  ,ny  )
    ResH   = @zeros(nx  ,ny  )
    ErrH   = @zeros(nx  ,ny  )
    b      = 0.5*ones(nx  ,ny  )
    H      = Data.Array( exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2) )
    rad    = (xc.-lx./2.0).^2 .+(yc.-ly./2.0)'.^2; b[rad.>lx/4] .= -0.5
    b      = Data.Array(b) # move data back to XPU
    # action
    for iter = 1:niter
        if (iter==11) global t0 = Base.time() end
        @parallel compute_H3_Err!(H3, ErrH, H, n)
        @parallel compute_flux!(qx, qy, H, H3, dx, dy)
        @parallel compute_ResH!(dtau, ResH, H3, qx, qy, b, dmp, dx, dy)
        @parallel update_H!(H, dtau, ResH)
        @parallel check_Err!(ErrH, H)
        if mod(iter,nout)==0 && viz 
            p1 = heatmap(xc, yc, Array(H)', xlabel="lx", ylabel="ly", title="shallow ice, iter=$iter", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:viridis)
            p2 = plot(xc, Array(H)[:,Int(round(ny/2))], xlabel="lx", ylabel="height", xlims=(xc[1], xc[end]), ylims=(0., 1.), legend=false, framestyle=:box, aspect_ratio=2.5)
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
