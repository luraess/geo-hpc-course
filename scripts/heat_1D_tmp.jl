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
	nx  = 100
	∂x  = lx/nx
	xc  = LinRange(∂x/2, lx-∂x/2, nx)
	🔥  = zeros(nx  )
	qx  = zeros(nx+1)
	# TODO add initial condition
	∂t  = ∂x^2/ρCp/λ/2.1
	# action
	t0  = Base.time()
	for it = 1:nt
		# TODO add physics
	end
	@printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*sizeof(lx))/((Base.time()-t0)/nt))
	if viz display(plot(xc, 🔥, legend=false, xlabel="lx", ylabel="heat", title="diffusion")) end
	return
end

@time heat_1D()
