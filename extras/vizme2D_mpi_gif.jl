using Plots, MAT

nprocs = (2, 2) # nprocs (x, y) dim
nt     = 100

@views function vizme2D_mpi(nprocs, nt)
    dname = "tmp_gif"; if isdir("./out2D/$(dname)")==false mkdir("./out2D/$(dname)") end; loadpath = "./out2D/$(dname)"; anim = Animation(loadpath,String[])
    for it = 1:nt
        T  = []
        ip = 1
        for ipx = 1:nprocs[1]
            for ipy = 1:nprocs[2]
                file = matopen("./out2D/T_$(ip-1)_$(it).mat"); T_loc = read(file, "T"); close(file)
                nx_i, ny_i = size(T_loc,1)-2, size(T_loc,2)-2
                ix1, iy1   = 1+(ipx-1)*nx_i, 1+(ipy-1)*ny_i
                if (ip==1)  T = zeros(nprocs[1]*nx_i, nprocs[2]*ny_i)  end
                T[ix1:ix1+nx_i-1,iy1:iy1+ny_i-1] .= T_loc[2:end-1,2:end-1]
                ip += 1
            end
        end
        heatmap(T', framestyle=:box, aspect_ratio=1, clims=(0.0,1.0), xlims=(1, size(T,1)), ylims=(1, size(T,2)), c=:hot, title="heat diffusion step=$(it)"); frame(anim)
    end
    gif(anim, "heat_2D_mpi_$(nprocs[1]*nprocs[2])procs.gif", fps = 15)
    return
end

vizme2D_mpi(nprocs, nt)
