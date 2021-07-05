module maxmin_twostage_subprob
using MAT
using JuMP
using Gurobi
using SCS
using SparseArrays
using LinearAlgebra
using LightGraphs, SimpleWeightedGraphs
using DataStructures
using SortingAlgorithms
import LightGraphs.Parallel


#using Plots
VIOL_EPS = 1e-4
DBARNZTH = 1e-3
ZNZTH = 1e-3
UNIFORM_GAMMA = 0.04

export initModel, solveModel!, robustCuttingPlaneAlg
#optimizer_constructor = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 10, "verbose" => 0)
#set_optimizer(problem, optimizer_constructor)
# Load Optimizer and model

#
function initModel(Din,firstIndices,t,dvrhs,β=0,phi_u_n=[])
    optimizer = Gurobi.Optimizer #SCS.Optimizer
    m = Model(optimizer)
    set_optimizer_attribute(m, "OutputFlag", 1)
    set_optimizer_attribute(m, "OptimalityTol", 1e-4)

    #println("************ debug **************")
    #println(Din.nzval)
    #println(findnz(Din))
    #println("*********************")

    @variable(m,g)
    n, nb = size(Din)
    if isempty(phi_u_n)  # if not given then initialize phi to unity to solve nominal problem
        phi_u_n = ones(n,1)
        phi_b_n = ones(n,1)
    end
    htn = n - firstIndices[1] + 1 # number of healthy tissue voxels
    @variable(m,x[1:nb]>=0)
    @variable(m,dbar[1:htn]>=0)
    @variable(m,u[1:htn]>=0)
    @variable(m,0<=z[1:htn]<=1)
    ptvn = firstIndices[1]-1
    @assert (n == ptvn+htn)
    # convert sparse matrix to JuMP constraint -- go over nonzero values using .nzval and .nzind
    @constraint(m,cons_ptv[i=1:ptvn], g <= phi_u_n[i]*Din[i,:].nzval'*x[Din[i,:].nzind])

    #@show(cons_ptv[1])

    firstIndices = [firstIndices; n+1]
    @constraint(m,cons_oar[k=1:length(firstIndices)-1,i=firstIndices[k]-ptvn:firstIndices[k+1]-1-ptvn], Din[i+ptvn,:].nzval'*x[Din[i+ptvn,:].nzind]-dbar[i] <= t[k])
    if β > 0 # penalty coefficient of DV related term in objective
        @constraint(m,cons_ind[i=1:htn], dbar[i]^2 <= u[i]*z[i])
    else # if beta coefficient is 0
        for i=1:htn
            JuMP.fix(dbar[i],0,force=true)
            JuMP.fix(z[i],0,force=true)
        end
    end
    if β > 0
        @constraint(m,dos_vol[k=1:length(firstIndices)-1],sum(z[i] for i=firstIndices[k]-ptvn:firstIndices[k+1]-1-ptvn) <= dvrhs[k] )
        @objective(m, Max, g-β*sum(u[i] for i=1:htn))
    else
        @objective(m, Max, g)
    end

    write_to_file(m,"DBGModel.lp")
    return m
end


function computeProjections(γ, phi_under, phi_bar, uniform_gamma=0)
    n, nn = size(γ)
    g = SimpleWeightedGraph(γ) #n,1:n,1:n,
    println("Started all-pairs-shortest path computations......................")
    dists = zeros(n,n)
    if uniform_gamma==0
        fws = Parallel.floyd_warshall_shortest_paths(g)
        dists = fws.dists
    else
        for i=1:n
            dists[i,:]=uniform_gamma*gdistances(g,i)#; sort_alg=RadixSort)
        end
    end
    ###################################### save dists for debug
    file = matopen("dists.mat", "w")
    write(file, "dists", dists)
    close(file)
    ######################################
    println("Finished all-pairs-shortest path computations")
    phi_under_n = phi_under
    phi_bar_n = phi_bar
    for i in 1:n
        phi_under_n[i] = maximum(phi_under[j]-dists[i,j] for j in 1:n)
        phi_bar_n[i] = minimum(phi_bar[j]+dists[i,j] for j in 1:n)
        # if false #i==3131
        #     println("phi_under_n[3131]=", phi_under_n[i], " phi_bar_n[3131]=",phi_bar_n[i])
        #     println(dists[i,i])
        #     println(phi_under[i])
        #     println(phi_under[3140])
        #     println(dists[i,3140])
        #     println(phi_under[3140]-dists[i,3140])
        #     @assert false
        # end
    end
    @assert all(phi_under_n .<= phi_bar_n)
    return phi_under,phi_bar, dists
end


function solveModel!(m,firstIndices,g=[],x=[],dbar=[],u=[],z=[])
    if ~isempty(g)
        JuMP.set_start_value(m[:g],g)
        JuMP.set_start_value(m[:x],x)
        JuMP.set_start_value(m[:dbar],dbar)
        JuMP.set_start_value(m[:u],u)
        JuMP.set_start_value(m[:z],z)
    end

    optimize!(m)
    if termination_status(m) == MOI.OPTIMAL
        println("****************** Optimal Objective Function value: ", JuMP.objective_value(m))
    elseif termination_status(m) == MOI.INFEASIBLE
        println("Infeasible model")
        return
    elseif termination_status(m) == MOI.INFEASIBLE_OR_UNBOUNDED
        println("Infeasible or unbounded model")
        return
    else
        println("Failed to obtain an optimal soln to subprob, status: ", JuMP.termination_status(m))
        return
    end
        #cons = constraint_by_name(m,"cons_oar")
    #n = length(cons)
    z = m[:z]
    zz = value.(z)
    dbar = m[:dbar]
    dbarr = value.(dbar)
    ptvN = length(dbarr)
    oarIdx = firstIndices[1]-1
    begIdx = 1
    firstIndices = firstIndices[2:length(firstIndices)]
    firstIndices = firstIndices.-firstIndices[1].+1
    firstIndices = [firstIndices; ptvN]
    for k=1:length(firstIndices)
        endIdx = firstIndices[k]
        oarDbar = dbarr[begIdx:endIdx]
        oarZ = zz[begIdx:endIdx]
#        nz = length(zz[zz.>10^-3])
        println("OAR k=", k, " Number of nonzero dbar: ", length(oarDbar[oarDbar.> DBARNZTH]), " Number of nonzero z: ",  length(oarZ[oarZ.> ZNZTH]) )
        begIdx = endIdx+1
    end
    return m
end


function robustCuttingPlaneAlg(Din,firstIndices,t,dvrhs,β,μ,γ,phi_under,phi_bar, L=1)
    ptvN, nn = size(γ)
    phi_u_n, phi_b_n, dists = computeProjections(γ, phi_under, phi_bar, UNIFORM_GAMMA)
    m = initModel(Din,firstIndices,t,dvrhs,β,phi_u_n)
    iter=0
    while(true)
        iter=iter+1
        solveModel!(m,firstIndices)
        println("Iteration: ", iter, " Objective value: ", objective_value(m), " *********")
        x = m[:x]
        xx = value.(x)
        d = Din[1:ptvN,:]*xx
        consCollect = SortedMultiDict{Float64,ScalarConstraint}()
        lthviol = VIOL_EPS

        @time for i1 = 1:ptvN
            for i2 = 1:ptvN
                if i1==i2
                    continue
                end
                phibar = min(phi_u_n[i2]+dists[i2,i1],phi_b_n[i1])
                v1 = phibar*d[i1]-μ*phi_u_n[i2]*d[i2]
                phiun = max(phi_u_n[i1]-dists[i1,i2],phi_u_n[i2])
                v2 = phi_b_n[i1]*d[i1]-μ*phiun*d[i2]
                v = max(v1,v2)
                if v > lthviol + VIOL_EPS
                    cons = []
                    if v1 > v2
                        cons = @build_constraint(phibar*sum(Din[i1,j]*x[j] for j in Din[i1,:].nzind)-μ*phi_u_n[i2]*sum(Din[i2,j]*x[j] for j in Din[i2,:].nzind) <= 0)
                    else
                        cons = @build_constraint(phi_b_n[i1]*sum(Din[i1,j]*x[j] for j in Din[i1,:].nzind) - μ*phiun*sum(Din[i2,j]*x[j] for j in Din[i2,:].nzind) <= 0)
                    end
                    insert!(consCollect,v,cons)
                    if length(consCollect) > L
#                        lthviol = min(lthviol,v)
                        delete!((consCollect,startof(consCollect)))
                        lthviol, tmp = first(consCollect)
                    end
                end
            end
        end
        println("Max Violation: ", lthviol)
        if isempty(consCollect) # no violated inequalities found
            println("Terminating cut algorithm..")
            break
        else
            println("Adding ", length(consCollect), " cuts.....................")
        end
        for (v,cons) in exclusive(consCollect,startof(consCollect),pastendsemitoken(consCollect))
            add_constraint(m,cons)
            #println("added constraint ", cons)
        end
    end
    return m
end
end
