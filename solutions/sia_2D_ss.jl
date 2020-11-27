using Plots, Printf
pyplot()
viz = true

@views function sia_2D()
	# physics
	lx   = 10.0
	ly   = 10.0
	n    = 3
	nτ   = 10000
	# numerics
	nx   = 100
	ny   = 101
	nout = 100
	γ    = 0.96
	ε    = 1e-8
	∂x   = lx/nx
	∂y   = ly/ny
	xc	 = LinRange(∂x/2, lx-∂x/2, nx)
	yc	 = LinRange(∂y/2, ly-∂y/2, ny)
	❄	 = zeros(nx  ,ny  )
	qx   = zeros(nx+1,ny  )
	qy   = zeros(nx  ,ny+1)
	∂τ   = zeros(nx  ,ny  )
	Res❄ = zeros(nx  ,ny  )
	Err❄ = zeros(nx  ,ny  )
	b    = 0.5*ones(nx  ,ny  )
	❄	 = exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2)
	rad  = (xc.-lx./2.0).^2 .+(yc.-ly./2.0)'.^2; b[rad.>lx/4] .= -0.5
	# action
	t0   = Base.time()
	for iτ = 1:nτ
		Err❄ .= ❄
		qx[2:end-1,:] .= .-0.5.*(❄[1:end-1,:].+❄[2:end,:]).^n .*diff(❄, dims=1)./∂x
		qy[:,2:end-1] .= .-0.5.*(❄[:,1:end-1].+❄[:,2:end]).^n .*diff(❄, dims=2)./∂y
		∂τ            .= min(∂x^2,∂y^2)./(1.0.+❄.^n)./4.1./4.0
		Res❄          .= (.-diff(qx, dims=1)./∂x .-diff(qy, dims=2)./∂y .+ b) .+ γ.*Res❄
		❄  		  	  .= max.(0.0, ❄ .+ ∂τ.*Res❄)
		Err❄ .-= ❄
		if mod(iτ,nout)==0 && viz 
			p1 = heatmap(xc, yc, ❄', xlabel="lx", ylabel="ly", title="shallow ice, iτ=$iτ", aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:viridis)
			p2 = plot(xc, ❄[:,Int(round(ny/2))], xlabel="lx", ylabel="height", xlims=(xc[1], xc[end]), ylims=(0., 1.), legend=false, framestyle=:box, aspect_ratio=2.5)
			l = @layout [ a{0.8h} ; b{0.2h} ]
			display(plot(p1, p2, layout=l ))
		end
		if mod(iτ,nout)==0 maxErr=maximum(abs.(Err❄)); @printf("iτ=%d, max(err)=%1.2e \n", iτ, maxErr); if maxErr<ε global iτg=iτ; break; end end
	end
	@printf("T_eff = %1.2e GB/s \n", (2/1e9*nx*ny*sizeof(lx))/((Base.time()-t0)/iτg))
	return
end

@time sia_2D()
