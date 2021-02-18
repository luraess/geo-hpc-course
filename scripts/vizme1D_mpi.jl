using Plots, MAT

nprocs = 4

@views function vizme1D_mpi(nprocs)
    T = []
    for ip = 1:nprocs
        file = matopen("T_$(ip-1).mat"); T_loc = read(file, "T"); close(file)
        nx_i = length(T_loc)-2
        i1   = 1 + (ip-1)*nx_i
        if (ip==1)  T = zeros(nprocs*nx_i)  end
        T[i1:i1+nx_i-1] .= T_loc[2:end-1]
    end
    display(plot(T, legend=false, framestyle=:box, linewidth=3, xlims=(1, length(T)), ylims=(0, 1), xlabel="nx", ylabel="heat", title="diffusion"))
    return
end

vizme1D_mpi(nprocs)
