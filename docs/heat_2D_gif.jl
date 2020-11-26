using Plots, Printf
viz = true
ENV["GKSwstype"]="nul"; if isdir("viz2D_out")==false mkdir("viz2D_out") end; loadpath = "./viz2D_out/"; anim = Animation(loadpath,String[])
println("Animation directory: $(anim.dir)")
@views function heat_2D()
	# physics
	lx   = 10.0
	ly   = 10.0
	λ    = 1.0
	ρCp  = 1.0
	nt   = 200
	# numerics
	nx   = 100
	ny   = 101
	nout = 10
	∂x   = lx/nx
	∂y   = ly/ny
	xc	 = LinRange(∂x/2, lx-∂x/2, nx)
	yc	 = LinRange(∂y/2, ly-∂y/2, ny)
	🔥	 = zeros(nx  ,ny  )
	qx   = zeros(nx+1,ny  )
	qy   = zeros(nx  ,ny+1)
	🔥	 = exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2)
	∂t   = min(∂x^2,∂y^2)/ρCp/λ/4.1
	# action
	t0   = Base.time()
	for it = 1:nt
		qx[2:end-1,:] .= .-λ.*diff(🔥, dims=1)./∂x
		qy[:,2:end-1] .= .-λ.*diff(🔥, dims=2)./∂y
		🔥  		  .= 🔥 .- ∂t./ρCp.*(diff(qx, dims=1)./∂x + diff(qy, dims=2)./∂y)
		if mod(it,nout)==0 && viz 
			heatmap(xc, yc, 🔥', xlabel="Lx", ylabel="Ly", title="heat diffusion, it=$it", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), clims=(0.,1.)); frame(anim)
		end
	end
	@printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/nt))
	gif(anim, "heat_2D.gif", fps = 15)
	return
end

@time heat_2D()
