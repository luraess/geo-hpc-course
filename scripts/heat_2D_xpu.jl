const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using Plots, Printf, Statistics
pyplot()
viz = true

@parallel function compute_flux!(qx::Data.Array, qy::Data.Array, 🔥::Data.Array, λ::Data.Number, ∂x::Data.Number, ∂y::Data.Number)
	@inn_x(qx) = -λ*@d_xa(🔥)/∂x
	@inn_y(qy) = -λ*@d_ya(🔥)/∂y
	return
end

@parallel function update_🔥!(🔥::Data.Array, qx::Data.Array, qy::Data.Array, ∂t::Data.Number, ρCp::Data.Number, ∂x::Data.Number, ∂y::Data.Number)
	@all(🔥) = @all(🔥) - ∂t/ρCp*(@d_xa(qx)/∂x + @d_ya(qy)/∂y)
	return
end

@views function heat_2D_xpu()
	# physics
	lx     = 10.0
	ly     = 10.0
	λ      = 1.0
	ρCp    = 1.0
	nt     = 200
	# numerics
	nx     = 128-1
	ny     = 128-1
	nout   = 10
	∂x     = lx/nx
	∂y     = ly/ny
	xc     = LinRange(∂x/2, lx-∂x/2, nx)
	yc     = LinRange(∂y/2, ly-∂y/2, ny)
	🔥     = @zeros(nx  ,ny  )
	qx     = @zeros(nx+1,ny  )
	qy     = @zeros(nx  ,ny+1)
	🔥     = Data.Array( exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2) )
	∂t     = min(∂x^2,∂y^2)/ρCp/λ/4.1
	# action
	t0     = Base.time()
	for it = 1:nt
		@parallel compute_flux!(qx, qy, 🔥, λ, ∂x, ∂y)
		@parallel update_🔥!(🔥, qx, qy, ∂t, ρCp, ∂x, ∂y)
		if mod(it,nout)==0 && viz
			display(heatmap(xc, yc, Array(🔥)', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
		end
	end
	@printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/nt))
	return
end

@time heat_2D_xpu()
