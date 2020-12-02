using Plots, Printf
pyplot()
viz = true

@views function sia_2D()
    # physics
    lx   = 10.0
    ly   = 10.0
    n    = 3
    nt   = 5000
    # numerics
    nx   = 127
    ny   = 127
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
        qx[2:end-1,:] .= .-0.5.*(H[1:end-1,:].+H[2:end,:]).^n .*diff(H, dims=1)./dx
        qy[:,2:end-1] .= .-0.5.*(H[:,1:end-1].+H[:,2:end]).^n .*diff(H, dims=2)./dy
        dt             = min(dx^2,dy^2)/maximum(H)^n/4.1/4.0
        H             .= max.(0.0, H .+ dt.*(.-diff(qx, dims=1)./dx .-diff(qy, dims=2)./dy .+ b))
        if mod(it,nout)==0 && viz 
            p1 = heatmap(xc, yc, H', xlabel="lx", ylabel="ly", title="shallow ice, it=$it", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:viridis)
            p2 = plot(xc, H[:,Int(round(ny/2))], xlabel="lx", ylabel="height", xlims=(xc[1], xc[end]), ylims=(0., 1.), legend=false)
            display(plot(p1, p2, layout=(2, 1)))
        end
    end
    @printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/nt))
    return
end

@time sia_2D()
