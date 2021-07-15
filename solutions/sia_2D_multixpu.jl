const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using ImplicitGlobalGrid, Plots, Printf, Statistics
import MPI
# pyplot()
do_viz = false

maximum_g(A) = (max_l  = maximum(A); MPI.Allreduce(max_l,  MPI.MAX, MPI.COMM_WORLD))

@parallel function check_Err!(ErrH::Data.Array, H::Data.Array)
    @all(ErrH) = @all(ErrH) - @all(H)
    return
end

@parallel function compute_Err!(ErrH::Data.Array, H::Data.Array)
    @all(ErrH) = @all(H)
    return
end

@parallel function compute_flux!(qx::Data.Array, qy::Data.Array, H::Data.Array, n::Data.Number, dx::Data.Number, dy::Data.Number)
    @all(qx) = -@av_xi(H)^n*@d_xi(H)/dx
    @all(qy) = -@av_yi(H)^n*@d_yi(H)/dy
    return
end

@parallel function compute_ResH!(dtau::Data.Array, ResH::Data.Array, H::Data.Array, qx::Data.Array, qy::Data.Array, b::Data.Array, n::Data.Number, dmp::Data.Number, dx::Data.Number, dy::Data.Number)
    @all(dtau) = min(dx*dx, dy*dy)/(1.0 + @inn(H)^n)/4.1./4.0
    @all(ResH) = (-@d_xa(qx)/dx  - @d_ya(qy)/dy + @inn(b)) + dmp*@all(ResH)
    return
end

@parallel function update_H!(H::Data.Array, dtau::Data.Array, ResH::Data.Array)
    @inn(H) = max(0.0, @inn(H) + @all(dtau)*@all(ResH))
    return
end

@views function sia_2D_xpu()
    # Physics
    lx, ly = 10.0, 10.0
    n      = 3.0
    # Numerics
    nx, ny = 64, 64
    niter  = 10000
    nout   = 100
    dmp    = 0.96
    ε      = 1e-8
    # Derived numerics
    me, dims, nprocs, coords, comm = init_global_grid(nx, ny, 1) # MPI initialisation
    select_device()                                              # select one GPU per MPI local rank (if >1 GPU per node)
    dx, dy = lx/nx_g(), ly/ny_g()                                # cell sizes
    # Array allocation
    H      = @zeros(nx  ,ny  )
    qx     = @zeros(nx-1,ny-2)
    qy     = @zeros(nx-2,ny-1)
    dtau   = @zeros(nx-2,ny-2)
    ResH   = @zeros(nx-2,ny-2)
    ErrH   = @zeros(nx  ,ny  )
    # Initial condition
    b      = 0.5*ones(nx,ny)
    H     .= Data.Array([exp(-(x_g(ix,dx,H)+dx/2 -0.5*lx)^2 -(y_g(iy,dy,H)+dy/2 -0.5*ly)^2) for ix=1:size(H,1), iy=1:size(H,2)])
    rad    = [(x_g(ix,dx,H)+dx/2 -0.5*lx)^2 + (y_g(iy,dy,H)+dy/2 -0.5*ly)^2 for ix=1:size(H,1), iy=1:size(H,2)]; b[rad.>lx/4] .= -0.5
    b      = Data.Array(b) # move data back to XPU
    # Preparation of visualisation
    if do_viz
        if me==0
            ENV["GKSwstype"]="nul"; if isdir("viz2D_xpu_out")==false mkdir("viz2D_xpu_out") end; loadpath = "./viz2D_xpu_out/"; anim = Animation(loadpath,String[])
            println("Animation directory: $(anim.dir)")
        end
        nx_v, ny_v = (nx-2)*dims[1], (ny-2)*dims[2]
        if (nx_v*ny_v*sizeof(Data.Number) > 0.8*Sys.free_memory()) error("Not enough memory for visualization.") end
        H_v   = zeros(nx_v, ny_v) # global array for visu
        H_inn = zeros(nx-2, ny-2) # no halo local array for visu
        Xi_g, Yi_g = LinRange(dx+dx/2, lx-dx-dx/2, nx_v), LinRange(dy+dy/2, ly-dy-dy/2, ny_v) # inner points only
    end
    # Time loop
    for iter = 1:niter
        if (iter==11) global t0 = Base.time() end
        @parallel compute_Err!(ErrH, H)
        @parallel compute_flux!(qx, qy, H, n, dx, dy)
        @parallel compute_ResH!(dtau, ResH, H, qx, qy, b, n, dmp, dx, dy)
        @hide_communication (8, 4) begin # communication/computation overlap
            @parallel update_H!(H, dtau, ResH)
            update_halo!(H)
        end
        @parallel check_Err!(ErrH, H)
        # Visualise
        if mod(iter,nout)==0 && do_viz
            H_inn .= H[2:end-1,2:end-1]; gather!(H_inn, H_v)
            if me==0
                p1 = heatmap(Xi_g, Yi_g, H_v', xlabel="lx", ylabel="ly", title="shallow ice, iter=$iter", aspect_ratio=1, xlims=(Xi_g[1], Xi_g[end]), ylims=(Yi_g[1], Yi_g[end]), c=:viridis)
                p2 = plot(Xi_g, H_v[:,Int(round(ny_g()/2))], xlabel="lx", ylabel="height", xlims=(Xi_g[1], Xi_g[end]), ylims=(0., 1.2), legend=false, framestyle=:box, aspect_ratio=2.5)
                l = @layout [ a{0.8h} ; b{0.2h} ]
                plot(p1, p2, layout=l ); frame(anim)
            end
        end
        if mod(iter,nout)==0 maxErr=maximum_g(abs.(ErrH)); if (me==0) @printf("iter=%d, max(err)=%1.2e \n", iter, maxErr) end; if maxErr<ε global itg=iter; break; end end
    end
    time_s = (Base.time()-t0)
    if (me==0) @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/(itg-10)), sigdigits=2)) end
    if (do_viz && me==0) gif(anim, "sia2D_multixpu.gif", fps = 5)  end
    finalize_global_grid()
    return
end

sia_2D_xpu()
