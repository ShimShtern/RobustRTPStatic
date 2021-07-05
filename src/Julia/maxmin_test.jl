include("maxmin_twostage_subprob.jl")
#using maxmin_twostage_subprob
#using maxmin_twostage_subprob.initModel, maxmin_twostage_subprob.solveModel
#import maxmin_twostage_subprob: initModel, solveModel
using MAT
using JuMP
using SparseArrays
file = matopen("liverEx2.mat")
inD = read(file,"Dij")
V = read(file,"V")
γ = read(file,"neighbors_Mat")
ϕ = read(file,"omf_Vec")
close(file)

println("initial D size: ", size(inD), "  OAR num (length(V)-1): ", length(V)-1, " D nnz: ", nnz(inD))
#println(inD.nzval)

ρ = [0.998 ; 0.998]
#t = [40.0 ; 40.0]
t = [70; 70]
#β = 0.01
β=0
μ = 1.5

#sumVec = sum(inD,2)
nonzeroIdxs = unique(inD.rowval[inD.nzval .!= 0]) #finall(sumVec->sumVec>0,sumVec)
nb = size(inD,2)

D = spzeros(0,nb)
firstIndices = [] # vector of indices of first voxels for each of the stuctures
dvrhs = zeros(length(V)-3)
for k in 1:length(V)-2
    if size(D,1)>0
        println(size(D))
        global firstIndices = [firstIndices; size(D,1)+1]
    end
    oarIndices = Array{Int}(V[k])
    idxs = intersect(nonzeroIdxs,oarIndices) #findall(in(nonzeroIdxs),oarIndices)
    appendD = inD[idxs,:]
    println("appending submat of D of size: ", size(appendD))
    global D = [D; appendD]
    if k > 1
        dvrhs[k-1] = floor((1-ρ[k-1])*length(V[k]))
    end
end
println(firstIndices)
#m = maxmin_twostage_subprob.initModel(D,firstIndices,t,dvrhs,β)
#m = maxmin_twostage_subprob.solveModel!(m)
phi_under = ϕ.-0.25
phi_under[phi_under.<0].= 0
phi_bar = ϕ.+0.3
phi_bar[phi_bar.>1].= 1

println("Now solving with min phi bar = ", minimum(phi_bar), " min phi_under = " , minimum(phi_under), " max phi_under = ", maximum(phi_under))
@time maxmin_twostage_subprob.robustCuttingPlaneAlg(D,firstIndices,t,dvrhs,β,μ,γ,phi_under,phi_bar,30)
