using CUDA, Plots, Printf
pyplot()
viz = true

function compute_flux!(qx, qy, 🔥, λ, ∂x, ∂y, nx, ny)
	ix = (blockIdx().x-1) * blockDim().x + threadIdx().x # thread ID, dimension x
	iy = (blockIdx().y-1) * blockDim().y + threadIdx().y # thread ID, dimension y

	if (2<ix<=nx && iy<=ny) qx[ix,iy] = -λ*(🔥[ix,iy]-🔥[ix-1,iy])/∂x end
	if (ix<=nx && 2<iy<=ny) qy[ix,iy] = -λ*(🔥[ix,iy]-🔥[ix,iy-1])/∂y end
	return
end

function update_🔥!(🔥, qx, qy, ∂t, ρCp, ∂x, ∂y, nx, ny)
	ix = (blockIdx().x-1) * blockDim().x + threadIdx().x # thread ID, dimension x
	iy = (blockIdx().y-1) * blockDim().y + threadIdx().y # thread ID, dimension y

	if (ix<=nx && iy<=ny) 🔥[ix,iy] = 🔥[ix,iy] - ∂t/ρCp*((qx[ix+1,iy]-qx[ix,iy])/∂x + (qy[ix,iy+1]-qy[ix,iy])/∂y) end
	return
end

@views function heat_2D_gpu()
	# physics
	lx     = 10.0
	ly     = 10.0
	λ      = 1.0
	ρCp    = 1.0
	nt     = 200
	# numerics
	BLOCKX = 16
	BLOCKY = 16
	GRIDX  = 8
	GRIDY  = 8
	nx     = BLOCKX*GRIDX-1
	ny     = BLOCKY*GRIDY-1
	nout   = 10
	∂x     = lx/nx
	∂y     = ly/ny
	xc	   = LinRange(∂x/2, lx-∂x/2, nx)
	yc	   = LinRange(∂y/2, ly-∂y/2, ny)
	🔥	   = CUDA.zeros(nx  ,ny  )
	qx     = CUDA.zeros(nx+1,ny  )
	qy     = CUDA.zeros(nx  ,ny+1)
	🔥	   = CuArray( exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2) )
	∂t     = min(∂x^2,∂y^2)/ρCp/λ/4.1
	cuthreads = (BLOCKX, BLOCKY, 1)
	cublocks  = (GRIDX,  GRIDY,  1)
	# action
	t0     = Base.time()
	for it = 1:nt
		@cuda blocks=cublocks threads=cuthreads compute_flux!(qx, qy, 🔥, λ, ∂x, ∂y, nx, ny)
		@cuda blocks=cublocks threads=cuthreads update_🔥!(🔥, qx, qy, ∂t, ρCp, ∂x, ∂y, nx, ny)
		if mod(it,nout)==0 && viz
			display(heatmap(xc, yc, Array(🔥)', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
		end
	end
	@printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/nt))
	return
end

@time heat_2D_gpu()
