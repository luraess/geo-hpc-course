using Plots, Printf
do_viz = true
if do_viz
    ENV["GKSwstype"]="nul"; if isdir("viz2D_out")==false mkdir("viz2D_out") end; loadpath = "./viz2D_out/"; anim = Animation(loadpath,String[])
    println("Animation directory: $(anim.dir)")
end

@views function heat_2D()
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
    🔥    = exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2)
    # Time loop
    for it = 1:nt
        if (it==11) global t0 = Base.time() end
        qx .= .-λ.*diff(🔥[:,2:end-1], dims=1)./dx
        qy .= .-λ.*diff(🔥[2:end-1,:], dims=2)./dy
        🔥[2:end-1,2:end-1] .= 🔥[2:end-1,2:end-1] .- dt./ρCp.*(diff(qx, dims=1)./dx .+ diff(qy,dims=2)./dy)
        # Visualise
        if mod(it,nout)==0 && do_viz 
            heatmap(xc, yc, 🔥', xlabel="Lx", ylabel="Ly", title="heat diffusion, it=$it", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), clims=(0.,1.)); frame(anim)
        end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/(nt-10)), sigdigits=2))
    gif(anim, "heat_2D.gif", fps = 15)
    return
end

@time heat_2D()
