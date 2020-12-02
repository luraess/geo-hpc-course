using Plots, Printf, Statistics
pyplot()
viz = true

function compute_flux!(qx, qy, T, λ, dx, dy, nx, ny)
    # Threads.@threads for iy=1:ny+1
    for iy=1:ny+1
        for ix=1:nx+1
            if (2<=ix<=nx && iy<=ny) qx[ix,iy] = -λ*(T[ix,iy]-T[ix-1,iy])/dx; end
            if (ix<=nx && 2<=iy<=ny) qy[ix,iy] = -λ*(T[ix,iy]-T[ix,iy-1])/dy; end
        end
    end
    return
end

function update_T!(T, qx, qy, dt, ρCp, dx, dy, nx, ny)
    # Threads.@threads for iy=1:ny+1
    for iy=1:ny+1
        for ix=1:nx+1
            if (ix<=nx && iy<=ny) T[ix,iy] = T[ix,iy] - dt/ρCp*((qx[ix+1,iy]-qx[ix,iy])/dx + (qy[ix,iy+1]-qy[ix,iy])/dy); end
        end
    end
    return
end

@views function heat_2D()
    @show Threads.nthreads()
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
    compute_flux!(qx, qy, T, λ, dx, dy, nx, ny)
    update_T!(T, qx, qy, dt, ρCp, dx, dy, nx, ny)
        if mod(it,nout)==0 && viz
            display(heatmap(xc, yc, T', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
            # sleep(.01)
        end
    end
    @printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/nt))
    return
end

@time heat_2D()
