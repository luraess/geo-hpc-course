using Plots, Printf
# pyplot()
do_viz = true

@views function sia_2D()
    # physics
    lx, ly = 10.0, 10.0
    n      = 3.0
    nt     = 5000
    # numerics
    nx, ny = 128, 128
    nout   = 100
    # Derived numerics
    dx, dy = lx/nx, ly/ny
    xc     = LinRange(dx/2, lx-dx/2, nx)
    yc     = LinRange(dy/2, ly-dy/2, ny)
    # Array allocation
    qx     = zeros(nx-1,ny-2)
    qy     = zeros(nx-2,ny-1)
    # Initial condition
    b      = 0.5*ones(nx  ,ny  )
    H      = exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2)
    rad    = (xc.-lx./2.0).^2 .+(yc.-ly./2.0)'.^2; b[rad.>lx/4] .= -0.5
    # Time loop
    for it = 1:nt
        if (it==11) global t0 = Base.time() end
        qx  .= .-(0.5.*(H[1:end-1,2:end-1].+H[2:end,2:end-1])).^n .*diff(H[:,2:end-1], dims=1)./dx
        qy  .= .-(0.5.*(H[2:end-1,1:end-1].+H[2:end-1,2:end])).^n .*diff(H[2:end-1,:], dims=2)./dy
        dt   = min(dx^2,dy^2)/maximum(H)^n/4.1/4.0
        H[2:end-1,2:end-1] .= max.(0.0, H[2:end-1,2:end-1] .+ dt.*(.-diff(qx, dims=1)./dx .-diff(qy, dims=2)./dy .+ b[2:end-1,2:end-1]))
        # Visualise
        if mod(it,nout)==0 && do_viz 
            p1 = heatmap(xc, yc, H', xlabel="lx", ylabel="ly", title="shallow ice, it=$it", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:viridis)
            p2 = plot(xc, H[:,Int(round(ny/2))], xlabel="lx", ylabel="height", xlims=(xc[1], xc[end]), ylims=(0., 1.2), legend=false, framestyle=:box)
            display(plot(p1, p2, layout=(2, 1)))
        end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/(nt-10)), sigdigits=2))
    return
end

@time sia_2D()
