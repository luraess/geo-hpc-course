using CUDA, Plots, Printf
pyplot()
viz = false

function compute_flux!(qx, qy, T, λ, dx, dy, nx, ny)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x # thread ID, dimension x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y # thread ID, dimension y

    if (2<=ix<=nx && iy<=ny)  qx[ix,iy] = -λ*(T[ix,iy]-T[ix-1,iy])/dx  end
    if (ix<=nx && 2<=iy<=ny)  qy[ix,iy] = -λ*(T[ix,iy]-T[ix,iy-1])/dy  end
    return
end

function update_T!(T, qx, qy, dt, ρCp, dx, dy, nx, ny)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x # thread ID, dimension x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y # thread ID, dimension y

    if (ix<=nx && iy<=ny)  T[ix,iy] = T[ix,iy] - dt/ρCp*((qx[ix+1,iy]-qx[ix,iy])/dx + (qy[ix,iy+1]-qy[ix,iy])/dy)  end
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
    dx     = lx/nx
    dy     = ly/ny
    xc     = LinRange(dx/2, lx-dx/2, nx)
    yc     = LinRange(dy/2, ly-dy/2, ny)
    T      = CUDA.zeros(nx  ,ny  )
    qx     = CUDA.zeros(nx+1,ny  )
    qy     = CUDA.zeros(nx  ,ny+1)
    T      = CuArray( exp.(.-(xc.-lx./2.0).^2 .-(yc.-ly./2.0)'.^2) )
    dt     = min(dx^2,dy^2)/ρCp/λ/4.1
    cuthreads = (BLOCKX, BLOCKY, 1)
    cublocks  = (GRIDX,  GRIDY,  1)
    # action
    t0     = Base.time()
    for it = 1:nt
        @cuda blocks=cublocks threads=cuthreads compute_flux!(qx, qy, T, λ, dx, dy, nx, ny)
        synchronize()
        @cuda blocks=cublocks threads=cuthreads update_T!(T, qx, qy, dt, ρCp, dx, dy, nx, ny)
        synchronize()
        if mod(it,nout)==0 && viz
            display(heatmap(xc, yc, Array(T)', xlabel="lx", ylabel="ly", title="heat diffusion, it=$it", clims=(0.,1.)))
            # sleep(.01)
        end
    end
    time_s = (Base.time()-t0)
    @printf("Time = %1.4e s, T_eff = %1.2f GB/s \n", time_s, round((2/1e9*nx*ny*sizeof(lx))/(time_s/nt), sigdigits=2))
    return
end

@time heat_2D_gpu()
