module maxmin_test_liver
include("maxmin_twostage_subprob.jl")

using MAT
using JuMP
using SparseArrays
using FileIO, JLD2
using Printf

const bLOAD_FROM_FILE = false
const bLOAD_FROM_FILE_gamma = false
const bLOAD_FROM_FILE_projection = false
const bSAVE_FILES = true
const bSAVE_DISTORPROJ_FILES = false

const MAX_HOMOGEN_CONSTR = 250


file = matopen("liverEx_2.mat")
#const ρ = [0.7; 1; 1] #0.5]# 1]
#const ρ = [1; 0.5; 1] #0.5]# 1]
const ρ = [1; 1; 1]

# In the Liver data the 1st OAR is liver, 2nd OAR is heart
#const t = [60.0 ; 40.0; 60]
const tmax = [60.0 ; 50.0; 60]
const t = tmax
#Heart
#Max Dose (Gy,0.03 cc) ≤ 50Gy ≤ 52 Gy > 52 Gy
#Mean Dose (Gy) ≤ 30 Gy ≤ 31 Gy > 31 Gy
#V40 ≤ 50% ≤ 55% > 55%
#Liver
#Mean Dose(Gy) ≤ 21 Gy ≤ 25 Gy > 25 Gy
#V30 ≤ 30% ≤ 40% > 40%

#λ=0 #unused reg param
β = 0#0.01
μ = 1.45 #1.45 #1.45 #1.25 #1.1  # 1.1 is for physical dose, should be higher for biological
δ = 0 #0.1  #0.1 #0.01:0.01:0.1
gamma_const=1
max_γ=0.05+gamma_const
max_dist=10
α = [0.0138749; 0.0218827; -0.000604924]
gamma_func(x) = (x<=max_dist)*(x>0)*(α'*[1;x;x^2])+(x>max_dist)*max_γ #0.04
γ = read(file, "neighbors_Mat")
ϕ = read(file, "omf_Vec")

n,nn = size(γ)    # to remove later, currently \gamma matrix is incorrect for liver
γ[1:n,1:n].=0

D = []
firstIndices = []
dvrhs = []

D_file="D_Liv_formatted.jld2"
#D_file="Patient4_Visit1_D_formatted.jld2"
if bLOAD_FROM_FILE
    D = FileIO.load(D_file, "D")
    firstIndices, dvrhs = FileIO.load(D_file, "firstIndices", "dvrhs")
else
    #file = matopen("Patient4_Visit1_16beams_withdeadvoxels.mat") #changed from 13 since it did not cover the PTV
    inD = read(file, "Dij")
    V = read(file, "V")
    close(file)
    println("initial D size: ",size(inD),"  OAR num (length(V)-1): ",length(V) - 1," D nnz: ",nnz(inD))
    #println(inD.nzval)
    #sumVec = sum(inD,2)
    nonzeroIdxs = unique(inD.rowval) #finall(sumVec->sumVec>0,sumVec)
    println("nonzero rows: ",length(nonzeroIdxs)," min: ",minimum(nonzeroIdxs)," max: ",maximum(nonzeroIdxs))
    nb = size(inD, 2)
    D = spzeros(0, nb)
    firstIndices = [] # vector of indices of first voxels for each of the stuctures
    dvrhs = zeros(length(t)) #zeros(length(V) - 1)

    # skipping the last strucure that includes dead voxels?
    for k = 1:length(t)+1 #length(V)-1
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
        FileIO.save(D_file,"D",D,"firstIndices",firstIndices,"dvrhs",dvrhs)
    end
end
close(file)

if sum(tmax-t) == 0 || sum(dvrhs) == 0
    dvrhs = []
end
inD = []
println(firstIndices)
@assert(size(γ, 1) == firstIndices[1] - 1) #need to make sure these are the same
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
file_name_gamma=@sprintf("./RS_Dists/Patient4_Visit1_Gamma_dist_new_%1.3f.jld2",gamma_const)
file_name_proj=@sprintf("./Projections/PPatient4_Visit1_rojection_new_%1.3f_%1.3f.jld2",gamma_const,δ)
if bLOAD_FROM_FILE_gamma
    dists = FileIO.load(file_name_gamma,"dists")
else
    dists = []
end
if bLOAD_FROM_FILE_projection*bLOAD_FROM_FILE_gamma
    phi_u_n, phi_b_n = FileIO.load(file_name_proj,"phi_u_n","phi_b_n")
elseif false
    phi_u_n, phi_b_n, dists = maxmin_twostage_subprob.computeProjections(γ, gamma_func, phi_under, phi_bar, dists)
    if bSAVE_DISTORPROJ_FILES
        if !bLOAD_FROM_FILE_gamma
            FileIO.save(file_name_gamma,"dists",dists)
        end
        FileIO.save(file_name_proj,"phi_u_n",phi_u_n,"phi_b_n",phi_b_n)
    end
else
    phi_u_n = ϕ
    phi_b_n = ϕ
end
dvrhs = []
@show t, dvrhs
#@assert(length(t)==length(dvrhs))

#time_prof = @elapsed m, htCn, homCn = maxmin_twostage_subprob.robustCuttingPlaneAlg!(D,firstIndices,t,tmax,[],β,μ,phi_u_n,phi_b_n,dists,[0;0],0,200)

time_prof=@elapsed model, oarCt, homCt =maxmin_twostage_subprob.robustCuttingPlaneAlg!(D,firstIndices,t,tmax,dvrhs,β,μ,phi_u_n, phi_b_n, dists,[0;0],0,MAX_HOMOGEN_CONSTR)
#time_prof=@elapsed model = maxmin_twostage_subprob.parametricSolveIncreasing(D,firstIndices,t,tmax,dvrhs,β,μ,phi_u_n, phi_b_n, dists,200)
#time_prof=@elapsed model,βvec,objValVec,dvCntMat = maxmin_twostage_subprob.parametricSolveDecreasing(D,firstIndices,t,tmax,dvrhs,μ,phi_u_n, phi_b_n, dists,200)
βvec = [β]
objValVec = [JuMP.objective_value(model)]

@show time_prof

maxmin_twostage_subprob.printDoseVolume(model, t, tmax, !isempty(dvrhs), true) # print out verbose output

mkpath("ResultsFiles")
summary_file_name="./ResultsFiles/parametric.jld2" # not recommended that 2 processes print to same file -- check with shimrit
FileIO.save(summary_file_name,"βvec",βvec,"objValVec",objValVec) #,"dvCntMat",dvCntMat)


end
