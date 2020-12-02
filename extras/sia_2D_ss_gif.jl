using Plots, Printf
# pyplot()
viz = true
ENV["GKSwstype"]="nul"; if isdir("viz2D_out")==false mkdir("viz2D_out") end; loadpath = "./viz2D_out/"; anim = Animation(loadpath,String[])
println("Animation directory: $(anim.dir)")

@views function sia_2D()
    # physics
    lx    = 10.0
    ly    = 10.0
    n     = 3
    niter = 10000
    # numerics
    nx    = 127
    ny    = 127
    nout  = 50
    dmp   = 0.96
    ε     = 1e-8
    dx    = lx/nx
    dy    = ly/ny
    xc    = LinRange(dx/2, lx-dx/2, nx)
    yc    = LinRange(dy/2, ly-dy/2, ny)
    H     = zeros(nx  ,ny  )
    qx    = zeros(nx+1,ny  )
    qy    = zeros(nx  ,ny+1)
    dtau  = zeros(nx  ,ny  )
    ResH  = zeros(nx  ,ny  )
    ErrH  = zeros(nx  ,ny  )
    b     = 0.5*ones(nx  ,ny  )
    H     = exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2)
    rad   = (xc.-lx./2.0).^2 .+(yc.-ly./2.0)'.^2; b[rad.>lx/4] .= -0.5
    # action
    itr=[]; mEr=[]
    t0   = Base.time()
    for iter = 1:niter
        ErrH .= H
        qx[2:end-1,:] .= .-0.5.*(H[1:end-1,:].+H[2:end,:]).^n .*diff(H, dims=1)./dx
        qy[:,2:end-1] .= .-0.5.*(H[:,1:end-1].+H[:,2:end]).^n .*diff(H, dims=2)./dy
        dtau          .= min(dx^2,dy^2)./(1.0.+H.^n)./4.1./4.0
        ResH          .= (.-diff(qx, dims=1)./dx .-diff(qy, dims=2)./dy .+ b) .+ dmp.*ResH
        H             .= max.(0.0, H .+ dtau.*ResH)
        ErrH .-= H
        if mod(iter,nout)==0 && viz 
            maxErr = maximum(abs.(ErrH)); @printf("iter=%d, max(err)=%1.2e \n", iter, maxErr); if maxErr<ε global itg=iter; break; end
            push!(itr, iter); push!(mEr, maxErr)
            p1 = heatmap(xc, yc, H', xlabel="lx", ylabel="ly", title="shallow ice", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:viridis)
            p2 = plot(xc, H[:,Int(round(ny/2))], xlabel="lx", ylabel="height", xlims=(xc[1], xc[end]), ylims=(0., 1.), legend=false, framestyle=:box, aspect_ratio=2.5, linewidth=2)
            p3 = plot(itr, mEr, xlabel="iters", ylabel="error", ylims=(ε/2., 1e-1), legend=false, framestyle=:box, yaxis=:log, linewidth=2)
            l = @layout [ a{0.5w} [ b{0.6h} ; c{0.3h} ] ]
            display(plot(p1, p2, p3, layout = l)); frame(anim)
        end
    end
    @printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/itg))
    gif(anim, "sia_2D.gif", fps = 5)
    return
end

@time sia_2D()
