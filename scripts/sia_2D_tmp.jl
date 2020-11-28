using Plots, Printf
pyplot()
viz = true

@views function sia_2D()
	# physics
	lx   = 10.0
	ly   = 10.0
	n    = 3
	nt   = 5000
	# numerics
	nx   = 127
	ny   = 127
	nout = 100
	∂x   = lx/nx
	∂y   = ly/ny
	xc   = LinRange(∂x/2, lx-∂x/2, nx)
	yc   = LinRange(∂y/2, ly-∂y/2, ny)
	❄    = zeros(nx  ,ny  )
	qx   = zeros(nx+1,ny  )
	qy   = zeros(nx  ,ny+1)
	b    = 0.5*ones(nx  ,ny  )
	❄    = exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2)
	rad  = (xc.-lx./2.0).^2 .+(yc.-ly./2.0)'.^2; b[rad.>lx/4] .= -0.5
	# action
	t0   = Base.time()
	for it = 1:nt
		# TODO add flux calculations here
		∂t             = min(∂x^2,∂y^2)/maximum(❄)^n/4.1/4.0
		❄             .= max.(0.0, # TODO add ❄ updated here )
		if mod(it,nout)==0 && viz
			p1 = heatmap(xc, yc, ❄', xlabel="lx", ylabel="ly", title="shallow ice, it=$it", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:viridis)
			p2 = plot(xc, ❄[:,Int(round(ny/2))], xlabel="lx", ylabel="height", xlims=(xc[1], xc[end]), ylims=(0., 1.), legend=false)
			display(plot(p1, p2, layout=(2, 1)))
		end
	end
	@printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/nt))
	return
end

@time sia_2D()
