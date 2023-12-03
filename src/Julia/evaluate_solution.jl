module evaluate_solution
include("maxmin_twostage_subprob.jl")
using Printf
using FileIO, JLD2


BIG_NUM = 1e6

export EvaluateSolution, CalculateNomPerformance, CalculatePhysPerformance, EvaluateEUD

function CalculateMinHom(d,firstIndices,dists,phi_u_n=[],phi_b_n=[])
    n = length(d)
    #_rowLoc = spzeros(n,1)
    V = 1:firstIndices[1]-1   # PTV
    #lastRow = firstIndices[1]-1
    ptvn = length(V)
    μ=1
    i_max=0
    j_min=0
    max_val=0
    min_val=0
    edge=0
    for i=1:ptvn-1
        for k=1:2
            if k==1
                is=i
                js=i+1:ptvn
            else
                is=i+1:ptvn
                js=i
            end
            phibar = minimum(hcat(vec(phi_u_n[js]*ones(length(is),1)+dists[is,js]),vec(phi_b_n[is]*ones(length(js),1))),dims=2)
            phiun = maximum(hcat(vec(phi_b_n[is]*ones(length(js),1)-dists[js,is]),vec(phi_u_n[js]*ones(length(is),1))),dims=2)
            μ1 = maximum(phibar.*d[is]./(phi_u_n[js].*d[js]))
            μ2 = maximum(phi_b_n[is].*d[is]./(phiun.*d[js]))
            μ=max(μ,μ1,μ2)
        end
    end
    return μ
end

function CalculateNomPerformance(d,firstIndices,phi)
    V = 1:firstIndices[1]-1   # PTV
    g = minimum(d[V].*phi)
    μ =  maximum(d[V].*phi)/g
    return g,μ
end


function CalculatePhysPerformance(d,firstIndices)
    V = 1:firstIndices[1]-1   # PTV
    g = minimum(d[V])
    μ =  maximum(d[V])/g
    return g,μ
end

function EvaluateSolution(d,firstIndices, dists,phi_u_n,phi_b_n)
    if isempty(phi_u_n)  # if not given then initialize phi to unity to solve nominal problem
        phi_u_n = ones(ptvn,1)
        phi_b_n = ones(ptvn,1)
    end
    ptvN, nn = size(dists)
    μ=CalculateMinHom(d,firstIndices,dists,phi_u_n,phi_b_n)
    g=minimum(phi_u_n.*d[1:firstIndices[1]-1]) #CalculateMinBioDose
    return g, μ
end

function EvaluateEUD(d,firstIndices)
    #calculating EUD for physical dose
    α=-10
    PTV_voxel_num=firstIndices[1]-1
    EUD=(sum(d[1:PTV_voxel_num].^α)/PTV_voxel_num).^(1/α)
    return EUD
end
end
