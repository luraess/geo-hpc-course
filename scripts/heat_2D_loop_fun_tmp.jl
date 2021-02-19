using Plots, Printf, Statistics
# pyplot()
do_viz = false

function compute_flux!(qx, qy, T, λ, dx, dy, nx, ny)
    # Threads.@threads for iy=2:ny-1

        # TODO add qx computation in a loop fashion here

    # Threads.@threads for iy=1:ny-1

        # TODO add qx computation in a loop fashion here

    return
end

function update_T!(T, qx, qy, dt, ρCp, dx, dy, nx, ny)
    # Threads.@threads for iy=2:ny-1
    for iy=2:ny-1
        for ix=2:nx-1
            T[ix,iy] = T[ix,iy] - dt/ρCp*((qx[ix,iy-1]-qx[ix-1,iy-1])/dx + (qy[ix-1,iy]-qy[ix-1,iy-1])/dy)
        end
    end
    return
end

@views function heat_2D()
    @show Threads.nthreads()
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
    qx     = zeros(nx-1,ny-2)
    qy     = zeros(nx-2,ny-1)
    # Initial condition
    T      = exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2)
    # Time loop
    for it = 1:nt
        if (it==11) global t0 = Base.time() end
        compute_flux!(qx, qy, T, λ, dx, dy, nx, ny)
        update_T!(T, qx, qy, dt, ρCp, dx, dy, nx, ny)
        # Visualise
        if mod(it,nout)==0 && do_viz
            display(heatmap(xc, yc, T', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
            # sleep(.01)
        end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/(nt-10)), sigdigits=2))
    return
end

@time heat_2D()
