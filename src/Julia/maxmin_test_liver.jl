
module maxmin_test_liver
include("maxmin_twostage_subprob.jl")
#using maxmin_twostage_subprob
#using maxmin_twostage_subprob.initModel, maxmin_twostage_subprob.solveModel
#import maxmin_twostage_subprob: initModel, solveModel
using MAT
using JuMP
using SparseArrays
using FileIO, JLD2
using Printf

const bLOAD_FROM_FILE=false
const bLOAD_FROM_FILE_gamma=false
const bLOAD_FROM_FILE_projection = false
const bSAVE_FILES = false
const bSAVE_DISTORPROJ_FILES = false

#file = matopen("liverEx2.mat")
#ρ = [0.99; 1; 1]
const ρ = [1; 1]# 1]
const t = [40.0 ; 40.0]
const tmax = [40.0 ; 40.0]

#t = [62; 54; 100]
#t=[60; 54; 100]
#tmax = [62; 54; 100]

β = 1e-8
μ = 1.3 #1.25 #1.1
δ=0.1  #0.1 #0.01:0.01:0.1
gamma_const=0.051
if Sys.islinux()
    μ = parse(Float64,ARGS[1])#1.2 #1.25 #1.1
    gamma_const = parse(Float64,ARGS[2])
    δ = parse(Float64,ARGS[3])#0.05
end
α=[0.03145+gamma_const,0.00228,-7.885e-5]
max_γ=0.05+gamma_const
max_dist=10
gamma_func(x) = (x<=max_dist)*(x>0)*(α'*[1;x;x^2])+(x>max_dist)*max_γ #0.04
file = matopen("liverEx_2.mat")
#file = matopen("Patient4_Visit1_16beams_withdeadvoxels.mat") #changed from 13 since it did not cover the PTV
γ = read(file, "neighbors_Mat")
ϕ = read(file, "omf_Vec")

D = []
firstIndices = []
dvrhs = []

if bLOAD_FROM_FILE
    D = FileIO.load("D_Liv_formatted.jld2", "D")
    firstIndices, dvrhs = FileIO.load("D_Liv_formatted.jld2", "firstIndices", "dvrhs")
else
    #file = matopen("Patient4_Visit1_16beams_withdeadvoxels.mat") #changed from 13 since it did not cover the PTV
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
    for k = 1:length(V)-1
        if k > 1 #size(D,1)>0
            println(size(D))
            global firstIndices = [firstIndices; size(D, 1) + 1]
        end
        oarIndices = Array{Int}(vec(V[k]))
        idxs = intersect(nonzeroIdxs, oarIndices) #findall(in(nonzeroIdxs),oarIndices)
        println("struct: ",k," before removing rows: ", length(V[k]), " after removing rows: ", length(idxs), " ", typeof(nonzeroIdxs), " ", typeof(oarIndices))
        #println(minimum(oarIndices), " ", maximum(oarIndices), " ")
        appendD = inD[idxs, :]
        rowN, colN = size(appendD)
        println("For struct: ", k, " appending submat of D of size: ", size(appendD))

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
            "D_Liv_formatted.jld2",
            "D",
            D,
            "firstIndices",
            firstIndices,
            "dvrhs",
            dvrhs)
    end
end
close(file)

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
file_name_gamma=@sprintf("./RS_Dists/Gamma_dist_new_%1.3f.jld2",gamma_const)
file_name_proj=@sprintf("./Projections/Projection_new_%1.3f_%1.3f.jld2",gamma_const,δ)
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

time_prof=@elapsed model, oarCt, homCt =maxmin_twostage_subprob.robustCuttingPlaneAlg(D,firstIndices,t,tmax,dvrhs,β,μ,phi_u_n, phi_b_n, dists,200)
#time_prof=@elapsed model = maxmin_twostage_subprob.parametricSolveIncreasing(D,firstIndices,t,tmax,dvrhs,β,μ,phi_u_n, phi_b_n, dists,200)

@show time_prof

maxmin_twostage_subprob.printDoseVolume(model, t, tmax, !isempty(dvrhs), true) # print out verbose output

xx = value.(model[:x])
g = value.(model[:g])
PhysHom=maximum(D[1:firstIndices[1]-1,:]*xx)/minimum(D[1:firstIndices[1]-1,:]*xx)
minPhysDose=minimum(D[1:firstIndices[1]-1,:]*xx)
#file_name=@sprintf("results_%1.3f_%.2f_%.3f_%.2f.jld2",β,μ,δ,gamma_const)
#FileIO.save(file_name,"δ",δ,"μ",μ,"β",β,"t",t,"gamma_const",gamma_const,"time_prof",time_prof,"xx",xx,"g",g,"PhysHom",PhysHom,"phi_under",phi_under,"phi_bar",phi_bar)
mkpath("ResultsFiles")
summary_file_name="./ResultsFiles/liver_new.txt" # not recommended that 2 processes print to same file -- check with shimrit
open(summary_file_name,"a+") do io
    println(io,μ,",",δ,",",gamma_const,",",β,",",t,",",g,",",PhysHom,",",minPhysDose,",",time_prof,",",xx)
end

end
