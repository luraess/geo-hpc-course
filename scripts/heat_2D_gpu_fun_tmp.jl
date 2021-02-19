using CUDA, Plots, Printf
# pyplot()
do_viz = false

function compute_flux!(qx, qy, T, λ, dx, dy, nx, ny)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x # thread ID, dimension x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y # thread ID, dimension y

    # TODO add qx and qy computation in a vectorised fashion here based on index ix, iy
    return
end

function update_T!(T, qx, qy, dt, ρCp, dx, dy, nx, ny)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x # thread ID, dimension x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y # thread ID, dimension y

    # TODO add T computation in a vectorised fashion here based on index ix, iy
    return
end

@views function heat_2D_gpu()
    # Physics
    lx, ly = 10.0, 10.0
    λ      = 1.0
    ρCp    = 1.0
    nt     = 200
    # Numerics
    BLOCKX = 16
    BLOCKY = 16
    GRIDX  = 8
    GRIDY  = 8
    nx, ny = BLOCKX*GRIDX, BLOCKY*GRIDY
    nout   = 10
    # Derived numerics
    dx, dy = lx/nx, ly/ny
    dt     = min(dx^2,dy^2)/ρCp/λ/4.1
    xc     = LinRange(dx/2, lx-dx/2, nx)
    yc     = LinRange(dy/2, ly-dy/2, ny)
    # Array allocation
    qx     = CUDA.zeros(nx-1,ny-2)
    qy     = CUDA.zeros(nx-2,ny-1)
    # Initial condition
    T      = CuArray( exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2) )
    cuthreads = (BLOCKX, BLOCKY, 1)
    cublocks  = (GRIDX,  GRIDY,  1)
    # Time loop
    t0     = Base.time()
    for it = 1:nt
        if (it==11) global t0 = Base.time() end
        @cuda blocks=cublocks threads=cuthreads compute_flux!(qx, qy, T, λ, dx, dy, nx, ny)
        synchronize()
        @cuda blocks=cublocks threads=cuthreads update_T!(T, qx, qy, dt, ρCp, dx, dy, nx, ny)
        synchronize()
        # Visualise
        if mod(it,nout)==0 && do_viz
            display(heatmap(xc, yc, Array(T)', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
            # sleep(.01)
        end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/(nt-10)), sigdigits=2))
    return
end

@time heat_2D_gpu()
