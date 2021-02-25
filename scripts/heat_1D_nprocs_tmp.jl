using Plots
# pyplot()
do_viz = true

@views function heat_1D_nprocs()
    # Physics
    lx  = 10.0
    λ   = 1.0
    ρCp = 1.0
    nt  = 200
    # Numerics
    np  = 4             # number of procs
    nx  = 32            # local number of gridpoints
    # Derived numerics
    nxg = (nx-2)*np+2   # global number of grid points
    dxg = lx/nxg
    dt  = dxg^2/ρCp/λ/2.1
    # Array allocation
    x   = zeros(nx,np)  # local coord array
    T   = zeros(nx,np)  # local Temp array
    xt  = zeros(nxg)    # global coord array
    Tt  = zeros(nxg)    # global initial Temp array
    Tg  = zeros(nxg)    # global Temp array
    # Initial condition
    for ip = 1:np
        i1 = 1 + (ip-1)*(nx-2)
        for ix = 1:nx
            x[ix,ip] = # TODO: define the global coordinate vector such that x = (ix-0.5)*dx - 0.5*lx
            T[ix,ip] = exp(-x[ix,ip]^2)
        end
        xt[i1:i1+nx-2] .= x[1:end-1,ip]; if (ip==np) xt[i1+nx-1] = x[end,ip] end
        Tt[i1:i1+nx-2] .= T[1:end-1,ip]; if (ip==np) Tt[i1+nx-1] = T[end,ip] end
    end
    # Time loop
    for it = 1:nt
        for ip = 1:np # compute physics locally
            T[2:end-1,ip] .= # TODO: implement heat diffusion: ∆T/dt = λ/ρCp ∂^2(T)/dx^2
        end
        for ip = 1:np-1 # update boundaries
            # TODO: update the inner boundaries to ensure a correct solution
            
        end
        for ip = 1:np # global picture
            i1 = 1 + (ip-1)*(nx-2)
            Tg[i1:i1+nx-2] .= T[1:end-1,ip]
        end
        # Visualise
        if do_viz
            plot(xt, Tt, legend=false, linewidth=5, xlabel="lx", ylabel="heat")
            display(plot!(xt, Tg, legend=false, linewidth=5, framestyle=:box, xlabel="lx", ylabel="heat", title="diffusion it=$(it)"))
        end
    end
end

@time heat_1D_nprocs()
