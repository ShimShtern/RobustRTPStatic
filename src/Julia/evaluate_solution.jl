module evaluate_solution
include("maxmin_twostage_subprob.jl")
using Printf
using FileIO, JLD2


BIG_NUM = 1e6

export EvaluateSolution

function CalculateMinHom(d,firstIndices,dists,phi_u_n=[],phi_b_n=[])
    n = length(d)
    #_rowLoc = spzeros(n,1)
    V = 1:firstIndices[1]-1   # PTV
    #lastRow = firstIndices[1]-1
    ptvn = length(V)
    if isempty(phi_u_n)  # if not given then initialize phi to unity to solve nominal problem
        phi_u_n = ones(ptvn,1)
        phi_b_n = ones(ptvn,1)
    end


    μ=1
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
            μ1 = maximum(phibar.*d[is])/minimum(phi_u_n[js].*d[js])
            μ2 = maximum(phi_b_n[is].*d[is])./minimum(phiun.*d[js])
            μ=max(μ,μ1,μ2)
        end
    end
    return μ
end

function EvaluateSolution(d,firstIndices, γ, gamma_const, δ,phi_under,phi_bar, bLOAD_FROM_FILE=false)
    ptvN, nn = size(γ)
    phi_u_n=[]
    phi_b_n=[]
    dists=[]
    file_name=@sprintf("Projection_Gamma_%1.3f_%1.3f.jld2",gamma_const,δ)
    if bLOAD_FROM_FILE
        phi_u_n, phi_b_n, dists = FileIO.load(file_name,"phi_u_n","phi_b_n","dists")
    else
        phi_u_n, phi_b_n, dists = maxmin_twostage_subprob.computeProjections(γ, gamma_const, phi_under, phi_bar)
        FileIO.save(file_name,"phi_u_n",phi_u_n,"phi_b_n",phi_b_n,"dists",dists)
    end
    μ=CalculateMinHom(d,firstIndices,dists,phi_u_n,phi_b_n)
    g=minimum(phi_u_n.*d[1:firstIndices[1]-1]) #CalculateMinBioDose
    return g, μ
end
end
