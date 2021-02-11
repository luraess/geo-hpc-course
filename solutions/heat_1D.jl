using Plots, Printf
pyplot()
viz = true

@views function heat_1D()
    # physics
    lx  = 10.0
    λ   = 1.0
    ρCp = 1.0
    nt  = 200
    # numerics
    nx  = 127
    dx  = lx/nx
    xc  = LinRange(dx/2, lx-dx/2, nx)
    qx  = zeros(nx+1)
    T   = exp.(.-(xc.-lx./2.0).^2)
    dt  = dx^2/ρCp/λ/2.1
    # action
    t0  = Base.time()
    for it = 1:nt
        qx[2:end-1] .= .-λ.*diff(T)./dx
        T           .= T .- dt./ρCp.*diff(qx)./dx
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*sizeof(lx))/(time_s/nt), sigdigits=2))
    if viz display(plot(xc, T, legend=false, xlabel="lx", ylabel="heat", title="diffusion")) end
    return
end

@time heat_1D()
