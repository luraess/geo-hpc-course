using Plots
# pyplot()
viz = true

@views function heat_1D_2procs()
    # physics
    Tl  = 10.0   # left Temperature
    Tr  = 1.0    # right Temperature
    λ   = 1.0
    ρCp = 1.0
    nt  = 200    # number of time steps
    # numerics
    nx  = 40     # number of local grid points
    dx  = 1.0    # cell size
    TL  = Tl*ones(nx)
    TR  = Tr*ones(nx)
    T   = [TL[1:end-1]; TR[2:end]]
    Tg  = T
    dt  = dx^2/ρCp/λ/2.1
    # action
    for it = 1:nt
        # compute physics
        # TODO: Implement ∆TL/dt = λ/ρCp ∂^2(TL)/dx^2
        # TODO: Implement ∆TR/dt = λ/ρCp ∂^2(TR)/dx^2

        # update boundaries (MPI)
        # TODO: Implement the boundary exchange - the core of distributed memory parallelisation

        # global picture
        T .= [TL[1:end-1]; TR[2:end]]
        # compute physics locally
        Tg[2:end-1] .= # update Tg (∆Tg/dt = λ/ρCp ∂^2(Tg)/dx^2) to cross check your fake parallelisation implementation
        # visualise
        plot(Tg, legend=false, linewidth=0, markershape=:circle, markersize=5)
        display(plot!(T, legend=false, linewidth=2, xlabel="lx", ylabel="heat", title="diffusion it=$(it)"))
    end
    return
end

@time heat_1D_2procs()
