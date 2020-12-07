using Plots, Printf
pyplot()
viz = false

@views function heat_1D()
    # physics
    lx  = 10.0
    λ   = 1.0
    ρCp = 1.0
    nt  = 200
    # numerics
    nx  = 128-1
    dx  = lx/nx
    xc  = LinRange(dx/2, lx-dx/2, nx)
    t   = zeros(nx  )
    qx  = zeros(nx+1)
    # TODO add initial condition
    dt  = dx^2/ρCp/λ/2.1
    # action
    t0  = Base.time()
    for it = 1:nt
        # TODO add physics
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/nt), sigdigits=2))
    if viz display(plot(xc, T, legend=false, xlabel="lx", ylabel="heat", title="diffusion")) end
    return
end

@time heat_1D()
