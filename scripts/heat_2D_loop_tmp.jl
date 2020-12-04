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
    nx   = 128-1
    ny   = 128-1
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
        for iy=1:ny
            for ix=2:nx
                qx[ix,iy] = -λ*(T[ix,iy]-T[ix-1,iy])/dx
            end
        end
        # TODO add qy computation in a loop fashion
        # TODO add T  computation in a loop fashion
        if mod(it,nout)==0 && viz
            display(heatmap(xc, yc, T', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
            # sleep(.01)
        end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2e GB/s \n", time_s, (2/1e9*nx*ny*sizeof(lx))/(time_s/nt))
    return
end

@time heat_2D()
