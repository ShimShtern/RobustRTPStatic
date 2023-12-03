phi_under = ϕ.- δ
phi_under[phi_under.<0].= 0
phi_bar = ϕ.+ δ
phi_bar[phi_bar.>1].= 1

gamma_func_str  = FileIO.load(RegFileName, "gamma_func")
print(gamma_func_str)
eval(Meta.parse(gamma_func_str*"+gamma_const*(x>0)"))

println("Now solving with min phi bar = ", minimum(phi_bar), " min phi_under = " , minimum(phi_under), " max phi_bar = ", maximum(phi_bar), " max phi_under = ", maximum(phi_under))
phi_u_n=[]
phi_b_n=[]
dists=[]
file_name_gamma=@sprintf("./Data/RS_Dists/%s_%s_Gamma_dist_new_%1.3f.jld2",Patient,Visit,gamma_const)
file_name_proj=@sprintf("./Data/Projections/%s_%s_Projection_new_%1.3f_%1.3f.jld2",Patient,Visit,gamma_const,δ)
if bLOAD_FROM_FILE_gamma
    dists = FileIO.load(file_name_gamma,"dists")
else
    dists = []
end
if bLOAD_FROM_FILE_projection*bLOAD_FROM_FILE_gamma
    phi_u_n, phi_b_n = FileIO.load(file_name_proj,"phi_u_n","phi_b_n")
else
    phi_u_n, phi_b_n, dists = maxmin_twostage_subprob.computeProjections(γ, gamma_func, phi_under, phi_bar, dists)
    if bSAVE_DISTORPROJ_FILES
        if !bLOAD_FROM_FILE_gamma
            FileIO.save(file_name_gamma,"dists",dists)
        end
        FileIO.save(file_name_proj,"phi_u_n",phi_u_n,"phi_b_n",phi_b_n)
    end
end
