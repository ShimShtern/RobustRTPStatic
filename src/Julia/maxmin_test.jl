include("maxmin_twostage_subprob.jl")
include("evaluate_solution.jl")
#using maxmin_twostage_subprob
#using maxmin_twostage_subprob.initModel, maxmin_twostage_subprob.solveModel
#import maxmin_twostage_subprob: initModel, solveModel
using MAT
using JuMP
using SparseArrays
using FileIO, JLD2
using Printf
bLOAD_FROM_FILE=false
bLOAD_FROM_FILE_gamma=false
bLOAD_FROM_FILE_projection = true
bSAVE_FILES = true
bSAVE_DISTORPROJ_FILES = false

MAX_HOMOGEN_CONSTR = 2000
#file = matopen("liverEx2.mat")
#ρ = [0.99; 1; 1]
#ρ = [1; 1; 1]
#t = [40.0 ; 40.0]
#t = [62; 54; 100]
#t=[60; 54; 100]
ρ = [1; 1; 1]
t = [62; 54; 100]
tmax = [62; 54; 100]
#β = 0.01
β = 0 #1e-8 # coeficient do deal with dose volume constraint
λ = [0; 0] # coeficient for regularizer preventing multiple solutions
#the first parameter coincides with the sum(x_i) regularizer
#the first parameter coincides with the g_nom - minimum nominal dose regularizer
δ=0.04 #0.01:0.01:0.1
μ = 1.1 #1.25 #1.1
gamma_const=0
if Sys.islinux()
    μ = parse(Float64,ARGS[1])#1.2 #1.25 #1.1
    gamma_const = parse(Float64,ARGS[2])
    δ = parse(Float64,ARGS[3])#0.05
end
#for patient 4
α=[0.0309449+gamma_const,-0.0014883, 0.0136452]
max_dist=10
max_γ=α'*[1;max_dist;log(max_dist)]+gamma_const
gamma_func(x) = (x<=max_dist)*(x>0)*(α'*[1;x;log(x)])+(x>max_dist)*max_γ #0.04
MatDataFile="Patient4_Visit1_16beams_refpoint5_notincludingdeadvoxels_20220810.mat"
file = matopen(MatDataFile) #changed from 13 since it did not cover the PTV
γ = read(file, "neighbors_Mat")
ϕ = read(file, "omf_Vec")
close(file)

D = []
firstIndices = []
dvrhs = []

#D_file="D_Liv_formatted.jld2"
D_file="Patient4_Visit1_D_formatted.jld2"
if bLOAD_FROM_FILE
    D = FileIO.load(D_file, "D")
    firstIndices, dvrhs = FileIO.load(D_file, "firstIndices", "dvrhs")
else
    file = matopen(MatDataFile)
    inD = read(file, "Dij")
    V = read(file, "V")
    close(file)
    println(
        "initial D size: ",
        size(inD),
        "  OAR num (length(V)-1): ",
        length(V) - 1,
        " D nnz: ",
        nnz(inD))
    #println(inD.nzval)
    #sumVec = sum(inD,2)
    nonzeroIdxs = unique(inD.rowval) #finall(sumVec->sumVec>0,sumVec)
    println(
        "nonzero rows: ",
        length(nonzeroIdxs),
        " min: ",
        minimum(nonzeroIdxs),
        " max: ",
        maximum(nonzeroIdxs))
    nb = size(inD, 2)
    D = spzeros(0, nb)
    firstIndices = [] # vector of indices of first voxels for each of the stuctures
    dvrhs = zeros(length(V) - 1)

    # skipping the last strucure that includes dead voxels?
    for k = 1:length(V)#-1
        if k > 1 #size(D,1)>0
            println(size(D))
            global firstIndices = [firstIndices; size(D, 1) + 1]
        end
        oarIndices = Array{Int}(vec(V[k]))
        idxs = intersect(nonzeroIdxs, oarIndices) #findall(in(nonzeroIdxs),oarIndices)
        println(
            "struct: ",
            k,
            " before removing rows: ",
            length(V[k]),
            " after removing rows: ",
            length(idxs),
            " ",
            typeof(nonzeroIdxs),
            " ",
            typeof(oarIndices))
        #println(minimum(oarIndices), " ", maximum(oarIndices), " ")
        appendD = inD[idxs, :]
        rowN, colN = size(appendD)
        println(
            "For organ: ",
            k,
            " appending submat of D of size: ",
            size(appendD))

        if rowN > 0
            global D = [D; appendD]
            if k > 1
                dvrhs[k-1] = floor((1 - ρ[k-1]) * length(V[k]))
                println("DVRHS: ", dvrhs[k-1] , " for organ: ", k-1)
            end
        end
    end
    if bSAVE_FILES
        FileIO.save(
            D_file,
            "D",
            D,
            "firstIndices",
            firstIndices,
            "dvrhs",
            dvrhs)
    end
end
if sum(tmax-t) == 0 || sum(dvrhs) == 0
    dvrhs = []
end
inD = []
println(firstIndices)
@assert(size(γ, 1) == firstIndices[1] - 1)#need to make sure these are the same
#m = maxmin_twostage_subprob.initModel(D,firstIndices,t,dvrhs,β)
#m = maxmin_twostage_subprob.solveModel!(m)
#for
phi_under = ϕ.- δ
phi_under[phi_under.<0].= 0
phi_bar = ϕ.+ δ
phi_bar[phi_bar.>1].= 1


println("Now solving with min phi bar = ", minimum(phi_bar), " min phi_under = " , minimum(phi_under), " max phi_bar = ", maximum(phi_bar), " max phi_under = ", maximum(phi_under))
phi_u_n=[]
phi_b_n=[]
dists=[]
#file_name_gamma=@sprintf("./RS_Dists/Gamma_dist_new_%1.3f.jld2",gamma_const)
file_name_gamma=@sprintf("./RS_Dists/Patient4_Gamma_dist_new_%1.3f.jld2",gamma_const)
#file_name_proj=@sprintf("./Projections/Projection_new_%1.3f_%1.3f.jld2",gamma_const,δ)
file_name_proj=@sprintf("./Projections/Patient4_Projection_new_%1.3f_%1.3f.jld2",gamma_const,δ)
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

time_prof=@elapsed model, htCn, homCn =maxmin_twostage_subprob.robustCuttingPlaneAlg!(D,firstIndices,t,tmax,dvrhs,β,μ,phi_u_n, phi_b_n, dists,λ, ϕ, MAX_HOMOGEN_CONSTR)

@show time_prof

maxmin_twostage_subprob.printDoseVolume(model, t, tmax, !isempty(dvrhs), true) # print out verbose output

xx = value.(model[:x])
g = value.(model[:g])
reg1 = value.(model[:reg1])
reg2 = value.(model[:reg2])
PhysHom=maximum(D[1:firstIndices[1]-1,:]*xx)/minimum(D[1:firstIndices[1]-1,:]*xx)
minPhysDose=minimum(D[1:firstIndices[1]-1,:]*xx)
d = D[1:firstIndices[1]-1,:]*xx
μ_nominal = evaluate_solution.CalculateMinHom(d,firstIndices,dists,ϕ,ϕ)
g_nominal = minimum(ϕ.*d[1:firstIndices[1]-1])
#file_name=@sprintf("results_%1.3f_%.2f_%.3f_%.2f.jld2",β,μ,δ,gamma_const)
#FileIO.save(file_name,"δ",δ,"μ",μ,"β",β,"t",t,"gamma_const",gamma_const,"time_prof",time_prof,"xx",xx,"g",g,"PhysHom",PhysHom,"phi_under",phi_under,"phi_bar",phi_bar)
#summary_file_name="./NewResultsFiles/no_dose_vol_nom_reg.txt" # not recommended that 2 processes print to same file -- check with shimrit
summary_file_name="./ResultsFiles/no_dose_vol_20220811.txt"
open(summary_file_name,"a") do io
    println(io,μ,",",δ,",",gamma_const,",",β,",",λ[1],",",λ[2],",",t,",",g,",",reg1,",",reg2,",",PhysHom,",",minPhysDose,",",time_prof,",",xx,",",μ_nominal,",",g_nominal)
end
