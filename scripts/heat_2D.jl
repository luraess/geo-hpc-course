using Plots, Printf
pyplot()
viz = true

@views function heat_2D()
    # physics
    lx   = 10.0
    ly   = 10.0
    λ    = 1.0
    ρCp  = 1.0
    nt   = 200
    # numerics
    nx   = 127
    ny   = 127
    nout = 10
    dx   = lx/nx
    dy   = ly/ny
    xc   = LinRange(dx/2, lx-dx/2, nx)
    yc   = LinRange(dy/2, ly-dy/2, ny)
    T    = zeros(nx  ,ny  )
    qx   = zeros(nx+1,ny  )
    qy   = zeros(nx  ,ny+1)
    T    = exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2)
    dt   = min(dx^2,dy^2)/ρCp/λ/4.1
    # action
    t0   = Base.time()
    for it = 1:nt
        qx[2:end-1,:] .= .-λ.*diff(T, dims=1)./dx
        qy[:,2:end-1] .= .-λ.*diff(T, dims=2)./dy
        T             .= T .- dt./ρCp.*(diff(qx, dims=1)./dx + diff(qy,dims=2)./dy)
        if mod(it,nout)==0 && viz 
            display(heatmap(xc, yc, T', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
            # sleep(.01)
        end
    end
    @printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/nt))
    return
end

@time heat_2D()
