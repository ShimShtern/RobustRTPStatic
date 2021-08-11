include("maxmin_twostage_subprob.jl")
#using maxmin_twostage_subprob
#using maxmin_twostage_subprob.initModel, maxmin_twostage_subprob.solveModel
#import maxmin_twostage_subprob: initModel, solveModel
using MAT
using JuMP
using SparseArrays
using FileIO, JLD2
bLOAD_FROM_FILE=true

#file = matopen("liverEx2.mat")
ρ = [0.998 ; 0.998; 0.998]
#t = [40.0 ; 40.0]
t = [62; 54; 100]
#β = 0.01
β=0
μ = 1.25
file = matopen("Patient4_Visit1_16beams_withdeadvoxels.mat") #changed from 13 since it did not cover the PTV
γ = read(file,"neighbors_Mat")
ϕ = read(file,"omf_Vec")
close(file)
D=[]
firstIndices=[]
if bLOAD_FROM_FILE
    D = FileIO.load("D_formatted.jld2","D")
    firstIndices, dvrhs = FileIO.load("D_formatted.jld2","firstIndices","dvrhs")
else
    file = matopen("Patient4_Visit1_16beams_withdeadvoxels.mat") #changed from 13 since it did not cover the PTV
    inD = read(file,"Dij")
    V = read(file,"V")
    close(file)

    println("initial D size: ", size(inD), "  OAR num (length(V)-1): ", length(V)-1, " D nnz: ", nnz(inD))
    #println(inD.nzval)

    #sumVec = sum(inD,2)
    nonzeroIdxs = unique(inD.rowval) #finall(sumVec->sumVec>0,sumVec)
    println("nonzero rows: ", length(nonzeroIdxs), " min: ", minimum(nonzeroIdxs), " max: ", maximum(nonzeroIdxs))
    nb = size(inD,2)

    D = spzeros(0,nb)
    firstIndices = [] # vector of indices of first voxels for each of the stuctures
    dvrhs = zeros(length(V)-1)
    for k in 1:length(V)-1
        if k> 1 #size(D,1)>0
            println(size(D))
            global firstIndices = [firstIndices; size(D,1)+1]
        end
        oarIndices = Array{Int}(vec(V[k]))
        idxs = intersect(nonzeroIdxs,oarIndices) #findall(in(nonzeroIdxs),oarIndices)
        println("organ: ", k, " before removing rows: ", length(V[k]), " after removing rows: ", length(idxs), " ", typeof(nonzeroIdxs) , " ", typeof(oarIndices))
        #println(minimum(oarIndices), " ", maximum(oarIndices), " ")
        appendD = inD[idxs,:]
        rowN, colN = size(appendD)
        println("For organ: ", k, " appending submat of D of size: ", size(appendD))
        if rowN > 0
            global D = [D; appendD]
            if k > 1
                dvrhs[k-1] = floor((1-ρ[k-1])*length(V[k]))
            end
        end
    end
    FileIO.save("D_formatted.jld2","D",D,"firstIndices",firstIndices,"dvrhs",dvrhs)
end
inD = []
println(firstIndices)
@assert(size(γ,1)==firstIndices[1]-1)#need to make sure these are the same
#m = maxmin_twostage_subprob.initModel(D,firstIndices,t,dvrhs,β)
#m = maxmin_twostage_subprob.solveModel!(m)
phi_under = ϕ.-0.05
phi_under[phi_under.<0].= 0
phi_bar = ϕ.+0.05
phi_bar[phi_bar.>1].= 1

println("Now solving with min phi bar = ", minimum(phi_bar), " min phi_under = " , minimum(phi_under), " max phi_bar = ", maximum(phi_bar), " max phi_under = ", maximum(phi_under))
@time maxmin_twostage_subprob.robustCuttingPlaneAlg(D,firstIndices,t,dvrhs,β,μ,γ,phi_under,phi_bar,200,bLOAD_FROM_FILE)
