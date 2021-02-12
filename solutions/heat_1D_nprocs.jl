using Plots
# pyplot()
viz = true

@views function heat_1D_nprocs()
    # physics
    lx  = 10.0
    λ   = 1.0
    ρCp = 1.0
    nt  = 200
    # numerics
    np  = 4             # number of procs
    nx  = 40            # local number of gridpoints
    nxg = (nx-2)*np+2   # global number of grid points
    dxg = lx/(nxg-1)
    dt  = dxg^2/ρCp/λ/2.1
    # initialise local vectors
    x   = zeros(nx,np)  # local coord array
    T   = zeros(nx,np)  # local Temp array
    xt  = zeros(nxg)    # global coord array
    Tt  = zeros(nxg)    # global initial Temp array
    Tg  = zeros(nxg)    # global Temp array
    for ip = 1:np
        i1 = 1 + (ip-1)*(nx-2)
        for ix = 1:nx
            x[ix,ip] = ( (ip-1)*(nx-2) + (ix-1) )*dxg - lx/2
            T[ix,ip] = exp(-x[ix,ip]^2)
        end
        xt[i1:i1+nx-2] .= x[1:end-1,ip]; if (ip==np) xt[i1+nx-1] = x[end,ip] end
        Tt[i1:i1+nx-2] .= T[1:end-1,ip]; if (ip==np) Tt[i1+nx-1] = T[end,ip] end
    end
    # action
    for it = 1:nt
        for ip = 1:np # compute physics locally
            T[2:end-1,ip] .= T[2:end-1,ip] .+ dt.*λ./ρCp.*diff(diff(T[:,ip])./dxg)./dxg
        end
        for ip = 1:np-1 # update boundaries
            T[end,ip  ] = T[    2,ip+1]
            T[  1,ip+1] = T[end-1,ip  ]
        end
        for ip = 1:np # global picture
            i1 = 1 + (ip-1)*(nx-2)
            Tg[i1:i1+nx-2] .= T[1:end-1,ip]
        end
        # visualise
        plot(xt, Tt, legend=false, linewidth=5, xlabel="lx", ylabel="heat")
        display(plot!(xt, Tg, legend=false, linewidth=5, xlabel="lx", ylabel="heat", title="diffusion it=$(it)"))
    end
end

@time heat_1D_nprocs()
