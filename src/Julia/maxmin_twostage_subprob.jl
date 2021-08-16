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
using FileIO, JLD2


#using Plots
VIOL_EPS = 1e-4
DBARNZTH = 1e-3
ZNZTH = 1e-3
UNIFORM_GAMMA = 0.04
BIG_OBJ = 1e8


LAZY_CONS_PERORGAN_TH = 5e4
LAZY_CONS_PERORGAN_INIT_NUM = 400
LAZY_CONS_PERORGAN_NUM_PER_ITER = 200
MAX_LAZY_CON_IT = 200

LAZYCONS_TIGHT_TH = 0.2
MAX_VIOL_RATIO_TH = 0.2

MAX_V_CONS = 20 #can be set to infinity
MAX_VIOL_EPS = 1e-4
MAX_VIOL_EPS_INIT = 10

global _D   # save D in order to load rows as needed
#global _rowLoc # 0 - not loaded into optimization problem, rowLoc[i] > 0 indicates row in model constraint matrix
global _V
global _N # voxels by organ not loaded into optimization

export initModel, solveModel!, robustCuttingPlaneAlg, addMostViolated!
#optimizer_constructor = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 10, "verbose" => 0)
#set_optimizer(problem, optimizer_constructor)
# Load Optimizer and model

#
function initModel(Din,firstIndices,t,dvrhs,β=0,phi_u_n=[],xinit=[])

    n, nb = size(Din)
    #_rowLoc = spzeros(n,1)
    global _V = fill(Int[],length(firstIndices)+1)
    global _N = fill(Int[],length(firstIndices)+1)
    _V[1] = 1:firstIndices[1]-1   # PTV
    #lastRow = firstIndices[1]-1
    println("Init Model Start......................")

    optimizer = Gurobi.Optimizer #SCS.Optimizer
    m = Model(optimizer)
    set_optimizer_attribute(m, "OutputFlag", 0)
    set_optimizer_attribute(m, "OptimalityTol", 1e-2)
    ptvn = length(_V[1])
    @variable(m,g)
    if isempty(phi_u_n)  # if not given then initialize phi to unity to solve nominal problem
        phi_u_n = ones(ptvn,1)
        phi_b_n = ones(ptvn,1)
    end
    #htn = n - firstIndices[1] + 1 # number of healthy tissue voxels
    @variable(m,x[1:nb]>=0)
    @constraint(m,cons_ptv[i in _V[1]], g <= phi_u_n[i]*Din[i,:].nzval'*x[Din[i,:].nzind])

        global _D = Din
        # select initial rows according to max violation with given initial soln vector
        if isempty(xinit)
            xinit = 100/nb*ones(nb,1)
        end

        firstIndices = [firstIndices; n+1]
        for k = 1:length(firstIndices)-1
            if firstIndices[k+1]-firstIndices[k] <= LAZY_CONS_PERORGAN_TH
                _V[k+1] =  firstIndices[k]:firstIndices[k+1]-1
            else
                _N[k+1] = firstIndices[k]:firstIndices[k+1]-1
            end
        end
        addMostViolated!(m,LAZY_CONS_PERORGAN_INIT_NUM,xinit,t)

    #println("************ debug **************")
    #println(Din.nzval)
    #println(findnz(Din))
    #println("*********************")

    #ptvn = firstIndices[1]-1
    #@assert (n == ptvn+htn)
    # convert sparse matrix to JuMP constraint -- go over nonzero values using .nzval and .nzind
    #@show(cons_ptv[1])
    #firstIndices = [firstIndices; n+1]
    if β > 0 # penalty coefficient of DV related term in objective
        @variable(m,dbar[k=2:length(_V),i in _V[k]]>=0)
        @variable(m,u[k=2:length(_V),i in _V[k]]>=0)
        @variable(m,0<=z[k=2:length(_V),i in _V[k]]<=1)
        @constraint(m,cons_oar[k=2:length(_V),i in _V[k]], Din[i,:].nzval'*x[Din[i,:].nzind]-dbar[i] <= t[k])
        @constraint(m,cons_ind[k=2:length(_V),i in _V[k]], dbar[i]^2 <= u[i]*z[i])
    else # if beta coefficient is 0
        @show length(_V)

        ##@constraint(m,cons_oar[k=2:length(_V),i in _V[k]], Din[i,:].nzval'*x[Din[i,:].nzind] <= t[k])
        #@constraint(m,cons_oar[k=2:length(_V),i in _V[k]], sum( Din[i,j]*x[j] for j in Din[i,:].nzind) <= t[k-1])

    end
    if β > 0
        @constraint(m,dos_vol[k=2:length(_V)],sum(z[i] for i in _V[k]) <= dvrhs[k] )
        @objective(m, Max, g-β*sum(u[i] for k=2:length(_V),i in _V[k]))
    else
        @objective(m, Max, g)
    end
    println("Init Model End......................")
    #write_to_file(m,"DBGModel.lp")

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
#    file = matopen("dists.mat", "w")
#    write(file, "dists", dists)
#    close(file)
    ######################################
    println("Finished all-pairs-shortest path computations")
    phi_under_n = phi_under
    phi_bar_n = phi_bar
    # for i in 1:n
    #    phi_under_n[i] = maximum(phi_under[j]-dists[i,j] for j in 1:n)
    #    phi_bar_n[i] = minimum(phi_bar[j]+dists[i,j] for j in 1:n)
    #    @assert(phi_under_n[i] <= phi_bar_n[i])
    #end
    phi_under_n = maximum(phi_under*ones(1,n)-dists,dims=1)'
    phi_bar_n = minimum(phi_bar*ones(1,n)+dists,dims=1)'
    @assert(all(phi_under_n .<= phi_bar_n))
    #@assert all(phi_under_n .<= phi_bar_n)
    println("Finished computeProjections")
    return phi_under_n,phi_bar_n, dists
end


function solveModel!(m,firstIndices,g=[],x=[],dbar=[],u=[],z=[])
    doseVol = false
    zVar = variable_by_name(m, "z")
    if zVar != nothing
        doseVol = true
    end

    if ~isempty(g)
        JuMP.set_start_value(m[:g],g)
        JuMP.set_start_value(m[:x],x)
        if doseVol
            JuMP.set_start_value(m[:dbar],dbar)
            JuMP.set_start_value(m[:u],u)
            JuMP.set_start_value(m[:z],z)
        end
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


    if doseVol
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
            println("OAR k=", k, " Number of nonzero dbar: ", length(oarDbar[oarDbar.> DBARNZTH]), " Number of nonzero z: ",  length(oarZ[oarZ.> ZNZTH]) )
            begIdx = endIdx+1
        end
    end
    return m
end

function addMostViolated!(m, n, x, t)
    #adds n most violated constraints
    # m - model
    # n - num constraint per organ to add
    # x - current solution value
    # t - bound vector on organs
    global _N
    global _D
    global _V
    num_const_added_aor = 0
    max_viol_aor = 0.0
    for k in 1:length(_N)-1
        #if length(_V[2])>0
        #    @show _V[2][end], _V[4][end]
        #end
        if length(_N[k+1]) > 0
            indicesToAdd = []
            viol = vec(_D[_N[k+1],:]*x .- t[k])
            n_min=min(n,length(viol))
            violIdxs = partialsortperm(viol,1:n_min,rev=true)
            max_viol_org=viol[violIdxs[1]]
            max_viol_aor=max(max_viol_org,max_viol_aor)
            @show viol[violIdxs[n_min]] viol[violIdxs[1]]
            #@show size(largestViol)
            #largestViol = largestViol[LAZY_CONS_PERORGAN_INIT_NUM]
            ## nextreme(DataStructures.FasterReverse(),LAZY_CONS_PERORGAN_INIT_NUM,viol)
            #violIdxs = findall(viol >= largestViol[length(largestViol)])
            first_not_viol_ind = findfirst((viol[violIdxs] .<= 0.0))
            #initialization of _V[k+1] is neccesary to prevent use of the same pointer
            if length(_V[k+1]) == 0 #if the first add all even if not violated
                indicesToAdd = _N[k+1][violIdxs]
                _V[k+1] = indicesToAdd
                deleteat!(_N[k+1], sort!(violIdxs))
            else
                if  first_not_viol_ind == nothing
                    first_not_viol_ind = n_min+1
                end
                if first_not_viol_ind>1
                    violIdxs = sort!(violIdxs[1:first_not_viol_ind-1])
                    indicesToAdd = _N[k+1][violIdxs]
                    append!(_V[k+1], indicesToAdd) #union!(_V[k+1],indicesToAdd)
                    deleteat!(_N[k+1], violIdxs)
                end
            end
            num_const_added_aor += length(indicesToAdd)
            xx = m[:x]
            @constraint(m,[i in indicesToAdd], sum( _D[i,j]*xx[j] for j in _D[i,:].nzind) <= t[k])
        end
        @show num_const_added_aor
        println("Number of voxels in _V ", length(_V[k+1]), " voxels in _N: ",  length(_N[k+1])  , " OAR: ", k+1 )
        #_rowLoc[indicesToAdd] = lastRow+1:lastRow+length(indicesToAdd)
        #lastRow = lastRow + length(indicesToAdd)
    end
    return max_viol_aor, num_const_added_aor
end

function addHomogenityConstr!(m,consCollect,ptvN,μ,L,phi_u_n,phi_b_n,dists,d)
    global _D
    lthviol = VIOL_EPS
    count_v=spzeros(Int64,ptvN)
    for i1 = 1:(ptvN-1)
        for k=1:2
            if k==1
                is=i1
                js=i1+1:ptvN
            else
                is=i1+1:ptvN
                js=i1
            end
            phibar = minimum(hcat(vec(phi_u_n[js]*ones(length(is),1)+dists[is,js]),vec(phi_b_n[is]*ones(length(js),1))),dims=2)
            phiun = maximum(hcat(vec(phi_b_n[is]*ones(length(js),1)-dists[js,is]),vec(phi_u_n[js]*ones(length(is),1))),dims=2)
            v1 = phibar.*d[is]-μ*phi_u_n[js].*d[js]*ones(length(is),1)
            v2 = phi_b_n[is].*d[is]*ones(length(js),1)-μ*phiun.*d[js].*ones(length(is),1)
            v = maximum(hcat(v1,v2),dims=2)
            inds=partialsortperm(vec(v),1:min(L,ptvN-i1,MAX_V_CONS),rev=true)
            #@show v[inds[1]],v[inds[end]]
            for l in inds
                if v[l] > lthviol + VIOL_EPS
                    if k==1
                        insert!(consCollect,v[l],[is; js[l]])
                        count_v[is] += 1
                        count_v[js[l]] += 1
                    else
                        insert!(consCollect,v[l],[is[l]; js])
                        count_v[is[l]] += 1
                        count_v[js] += 1
                    end
                    if length(consCollect) > L
                        tmp, pair = first(consCollect)
                        count_v[pair[1]] -= 1
                        count_v[pair[2]] -= 1
                        delete!((consCollect,startof(consCollect)))
                        lthviol, tmp = first(consCollect)
                    end
                else
                    break
                end
            end
        end
    end
    maxviol = lthviol
    if !isempty(consCollect)
        maxviol, tmp = last(consCollect)
    end
    return lthviol, maxviol
end

function robustCuttingPlaneAlg(Din,firstIndices,t,dvrhs,β,μ,γ,phi_under,phi_bar, L=1,bLOAD_FROM_FILE=false)
    ptvN, nn = size(γ)
    phi_u_n=[]
    phi_b_n=[]
    dists=[]
    if bLOAD_FROM_FILE
        phi_u_n, phi_b_n, dists = FileIO.load("Projection.jld2","phi_u_n","phi_b_n","dists")
    else
        phi_u_n, phi_b_n, dists = computeProjections(γ, phi_under, phi_bar, UNIFORM_GAMMA)
        FileIO.save("Projection.jld2","phi_u_n",phi_u_n,"phi_b_n",phi_b_n,"dists",dists)
    end
    m = initModel(Din,firstIndices,t,dvrhs,β,phi_u_n)
    iter=0
    stage=1
    sum_num_const_added_aor = LAZY_CONS_PERORGAN_INIT_NUM
    while(true)
        iter=iter+1
        prevObj = BIG_OBJ
        newObj = BIG_OBJ
        prev_viol_aor = 1e6
        max_viol_aor = 1e6
        @timev for it = 1:MAX_LAZY_CON_IT
            solveModel!(m,firstIndices)
            newObj=JuMP.objective_value(m)
            max_viol_aor, num_const_added_aor = addMostViolated!(m, LAZY_CONS_PERORGAN_NUM_PER_ITER, JuMP.value.(m[:x]), t)
            sum_num_const_added_aor+=num_const_added_aor
            @show max_viol_aor
            if ((prevObj-newObj)/prevObj < LAZYCONS_TIGHT_TH && abs(prev_viol_aor-max_viol_aor)/prev_viol_aor < MAX_VIOL_RATIO_TH && max_viol_aor<=MAX_VIOL_EPS_INIT && stage==1) || max_viol_aor<=MAX_VIOL_EPS
                println("Total of ", sum_num_const_added_aor, " constraints")
                println("Reduced the objective by ", (prevObj-newObj)/prevObj)
                println("Constraints violation by ", (prev_viol_aor-max_viol_aor)/prev_viol_aor)
                println("Iter= ", it)
                break
            end
            prev_viol_aor = max_viol_aor
            prevObj = newObj
        end

        println("Iteration: ", iter, " Objective value: ", newObj, " *********")
        consCollect = SortedMultiDict{Float64,Array{Int64,1}}()
        x = m[:x]
        xx = value.(x)
        d = _D[1:ptvN,:]*xx
        @show minimum(d), maximum(d)
        @timev min_hom_viol, max_hom_viol=addHomogenityConstr!(m,consCollect,ptvN,μ,L,phi_u_n,phi_b_n,dists,d)
        if max_hom_viol<=MAX_VIOL_EPS && max_viol_aor<=MAX_VIOL_EPS# no violated inequalities found
            #println("Terminating cut algorithm..")
            break
        elseif !isempty(consCollect)
            stage=1
            println("Max Violation: ", max_hom_viol, " " ,min_hom_viol)
            println("Adding ", length(consCollect), " cuts.....................")
        else
            stage=2
        end

        cons_array = Any[]
        @timev for (v,pair) in exclusive(consCollect,startof(consCollect),pastendsemitoken(consCollect))
            #add_constraint(m,cons)
            i1=pair[1]
            i2=pair[2]
            #shimrit:added both constraints
            phibar = min(phi_u_n[i2]+dists[i2,i1],phi_b_n[i1])
            v1 = phibar*d[i1]-μ*phi_u_n[i2]*d[i2]
            #@show i1, i2, phibar*d[i1], phi_u_n[i2]*d[i2]
            phiun = max(phi_b_n[i1]-dists[i1,i2],phi_u_n[i2])
            v2 = phi_b_n[i1]*d[i1]-μ*phiun*d[i2]
            #@show i1, i2, phi_b_n[i1]*d[i1], phiun*d[i2]
            #if v1 > v2
            if v1 > VIOL_EPS
                push!(cons_array,@build_constraint(phibar*sum(Din[i1,j]*x[j] for j in Din[i1,:].nzind)-μ*phi_u_n[i2]*sum(Din[i2,j]*x[j] for j in Din[i2,:].nzind) <= 0))
            end
            #else
            if v2 > VIOL_EPS
                push!(cons_array,@build_constraint(phi_b_n[i1]*sum(Din[i1,j]*x[j] for j in Din[i1,:].nzind) - μ*phiun*sum(Din[i2,j]*x[j] for j in Din[i2,:].nzind) <= 0))
            end
            #end
            #println("added constraint ", cons)
        end
        lenCons = length(cons_array)
        @timev for i=1:lenCons
            add_constraint(m, cons_array[i])
        end
    end
    return m
end
end
