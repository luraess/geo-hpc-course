using Plots, Printf
pyplot()
viz = false

@views function sia_2D()
    # physics
    lx   = 10.0
    ly   = 10.0
    n    = 3
    nt   = 5000
    # numerics
    nx   = 128-1
    ny   = 128-1
    nout = 100
    dx   = lx/nx
    dy   = ly/ny
    xc   = LinRange(dx/2, lx-dx/2, nx)
    yc   = LinRange(dy/2, ly-dy/2, ny)
    H    = zeros(nx  ,ny  )
    qx   = zeros(nx+1,ny  )
    qy   = zeros(nx  ,ny+1)
    b    = 0.5*ones(nx  ,ny  )
    H    = exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2)
    rad  = (xc.-lx./2.0).^2 .+(yc.-ly./2.0)'.^2; b[rad.>lx/4] .= -0.5
    # action
    t0   = Base.time()
    for it = 1:nt
        # TODO add flux calculations here
        dt             = min(dx^2,dy^2)/maximum(H)^n/4.1/4.0
        H             .= max.(0.0, # TODO add H updated here )
        if mod(it,nout)==0 && viz
            p1 = heatmap(xc, yc, H', xlabel="lx", ylabel="ly", title="shallow ice, it=$it", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:viridis)
            p2 = plot(xc, H[:,Int(round(ny/2))], xlabel="lx", ylabel="height", xlims=(xc[1], xc[end]), ylims=(0., 1.), legend=false)
            display(plot(p1, p2, layout=(2, 1)))
        end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/nt), sigdigits=2))
    return
end

@time sia_2D()
