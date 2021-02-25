using Plots, MAT

nprocs = (2, 2) # nprocs (x, y) dim

@views function vizme2D_mpi(nprocs)
    T  = []
    ip = 1
    for ipx = 1:nprocs[1]
        for ipy = 1:nprocs[2]
            file = matopen("T_$(ip-1).mat"); T_loc = read(file, "T"); close(file)
            nx_i, ny_i = size(T_loc,1)-2, size(T_loc,2)-2
            ix1, iy1   = 1+(ipx-1)*nx_i, 1+(ipy-1)*ny_i
            if (ip==1)  T = zeros(nprocs[1]*nx_i, nprocs[2]*ny_i)  end
            T[ix1:ix1+nx_i-1,iy1:iy1+ny_i-1] .= T_loc[2:end-1,2:end-1]
            ip += 1
        end
    end
    display(heatmap(T', framestyle=:box, aspect_ratio=1, xlims=(1, size(T,1)), ylims=(1, size(T,2)), c=:hot, title="heat diffusion"))
    return
end

vizme2D_mpi(nprocs)
