using Plots, MAT

nprocs = 4

@views function vizme1Dmpi(nprocs)
    T = []
    for ip = 1:nprocs
        file = matopen("T_$(ip-1).mat"); T_loc = read(file, "T"); close(file)
        nx_i = length(T_loc)-2
        if (ip==1)  T = zeros(nprocs*nx_i)  end
        i1   = 1 + (ip-1)*nx_i
        @show((i1,i1+nx_i-1))
        T[i1:i1+nx_i-1] .= T_loc[2:end-1]
    end
    display(plot(T, legend=false, framestyle=:box, linewidth=2, xlabel="nx", ylabel="heat", title="diffusion"))
    return
end

vizme1Dmpi(nprocs)