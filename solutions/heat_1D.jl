using Plots, Printf
# pyplot()
do_viz = true

@views function heat_1D()
    # Physics
    lx  = 10.0
    λ   = 1.0
    ρCp = 1.0
    nt  = 200
    # Numerics
    nx  = 128
    # Derived numerics
    dx  = lx/nx
    dt  = dx^2/ρCp/λ/2.1
    xc  = LinRange(dx/2, lx-dx/2, nx)
    # Array allocation
    qx  = zeros(nx-1)
    # Initial condition
    T   = exp.(.-(xc.-lx./2.0).^2)
    # Time loop
    for it = 1:nt
        if (it==11) global t0 = Base.time() end
        qx         .= .-λ.*diff(T)./dx
        T[2:end-1] .= T[2:end-1] .- dt./ρCp.*diff(qx)./dx
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*sizeof(lx))/(time_s/(nt-10)), sigdigits=2))
    # Visualise
    if do_viz display(plot(xc, T, legend=false, framestyle=:box, xlabel="lx", ylabel="heat", title="diffusion")) end
    return
end

@time heat_1D()
