using CUDA, Plots, Printf
pyplot()
viz = true

@views function heat_2D_gpu()
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
	xc   = LinRange(∂x/2, lx-∂x/2, nx)
	yc   = LinRange(∂y/2, ly-∂y/2, ny)
	🔥   = CUDA.zeros(nx  ,ny  )
	qx   = CUDA.zeros(nx+1,ny  )
	qy   = CUDA.zeros(nx  ,ny+1)
	🔥   = CuArray( exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2) )
	∂t   = min(∂x^2,∂y^2)/ρCp/λ/4.1
	# action
	t0   = Base.time()
	for it = 1:nt
		qx[2:end-1,:] .= .-λ.*diff(🔥, dims=1)./∂x
		qy[:,2:end-1] .= .-λ.*diff(🔥, dims=2)./∂y
		🔥            .= 🔥 .- ∂t./ρCp.*(diff(qx, dims=1)./∂x + diff(qy, dims=2)./∂y)
		if mod(it,nout)==0 && viz
			display(heatmap(xc, yc, Array(🔥)', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
		end
	end
	@printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/nt))
	return
end

@time heat_2D_gpu()
