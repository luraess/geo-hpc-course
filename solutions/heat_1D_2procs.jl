using Plots
# pyplot()
do_viz = true

@views function heat_1D_2procs()
    # Physics
    Tl  = 10.0   # left Temperature
    Tr  = 1.0    # right Temperature
    λ   = 1.0
    ρCp = 1.0
    nt  = 200    # number of time steps
    # Numerics
    nx  = 32     # number of local grid points
    dx  = 1.0    # cell size
    # Derived numerics
    dt  = dx^2/ρCp/λ/2.1
    # Initial condition
    TL  = Tl*ones(nx)
    TR  = Tr*ones(nx)
    T   = [TL[1:end-1]; TR[2:end]]
    Tg  = T
    # Time loop
    for it = 1:nt
        # Compute physics
        TL[2:end-1] .= TL[2:end-1] .+ dt.*λ./ρCp.*diff(diff(TL)./dx)./dx
        TR[2:end-1] .= TR[2:end-1] .+ dt.*λ./ρCp.*diff(diff(TR)./dx)./dx
        # Update boundaries (MPI)
        TL[end] = TR[2]
        TR[1]   = TL[end-1]
        # Global picture
        T .= [TL[1:end-1]; TR[2:end]]
        # Compute physics locally
        Tg[2:end-1] .= Tg[2:end-1] .+ dt.*λ./ρCp.*diff(diff(Tg)./dx)./dx
        # Visualise
        if do_viz
            plot(Tg, legend=false, linewidth=0, markershape=:circle, markersize=5)
            display(plot!(T, legend=false, linewidth=2, framestyle=:box, xlabel="lx", ylabel="heat", title="diffusion it=$(it)"))
        end
    end
    return
end

@time heat_1D_2procs()
