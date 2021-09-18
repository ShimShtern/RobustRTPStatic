module maxmin_twostage_subprob
using MAT
using JuMP
using Gurobi
#using SCS
using UnicodePlots
using SparseArrays
using LinearAlgebra
using LightGraphs, SimpleWeightedGraphs
using DataStructures
using SortingAlgorithms
using Statistics
import LightGraphs.Parallel
using FileIO, JLD2
using Printf

#using Plots
INITXNORM = 100

VIOL_EPS = 1e-3
DBARNZTH = 1e-4
ZNZTH = 1e-4
BIG_OBJ = 1e8
UNIFORM_GAMMA = true

LAZY_CONS_PERORGAN_TH = 5e4
LAZY_CONS_PERORGAN_INIT_NUM = 400
LAZY_CONS_PERORGAN_NUM_PER_ITER = 200
MAX_LAZY_CON_IT = 500

LAZYCONS_TIGHT_TH = 0.2
MAX_VIOL_RATIO_TH = 0.2

MAX_V_CONS = 10 #can be set to infinity
MAX_VIOL_EPS = 1e-3
MAX_VIOL_EPS_INIT = 10

SURPLUS_VAR_OBJ_NORM = 1  # penalty norm either 1 or 2

BIG_NUM = 1e6

global _D   # save D in order to load rows as needed
#global _rowLoc # 0 - not loaded into optimization problem, rowLoc[i] > 0 indicates row in model constraint matrix
global _V
global _N # voxels by organ not loaded into optimization


export initModel, solveModel!,robustCuttingPlaneAlg, printDoseVolume
#optimizer_constructor = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 10, "verbose" => 0)
#set_optimizer(problem, optimizer_constructor)
# Load Optimizer and model


function initDoseVolume(m, t, tmax, dvrhs)
    println("initDoseVolume....")
    for k =1:length(t)
        if tmax[k] > t[k]
            @variable(m, z[union(_V[k+1],_N[k+1])]==0)
            @constraint(m,dos_vol[k],sum(z[i] for i in _V[k+1] ) + sum(z[i] for i in _N[k+1] ) <= dvrhs[k] )
        end
    end
    return m
end

#
function initModel(Din,firstIndices,t,dvrhs=[],β=0,phi_u_n=[],xinit=[])

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
        xinit = INITXNORM/nb*ones(nb,1)
    end
    dbar = m[:dbar] = Dict()

    firstIndices = [firstIndices; n+1]
    for k = 1:length(firstIndices)-1
        if firstIndices[k+1]-firstIndices[k] <= LAZY_CONS_PERORGAN_TH
            _V[k+1] =  firstIndices[k]:firstIndices[k+1]-1
            #if (!isempty(dvrhs) &&  t[k] < tmax[k] ) ||
            if β>0
                for i in _V[k+1]
                    dbar[i] = @variable(m,lower_bound=0)
                end
                @constraint(m,[i in _V[k+1] ], sum( _D[i,j]*x[j] for j in _D[i,:].nzind) - dbar[i] <= t[k])
            else
                #    if t[k]<tmax[k]
                #        unfix(z[_V[k+1]])
                #        set_lower_bound(z[_V[k+1]],0)
                #        set_upper_bound(z[_V[k+1]],1)
                #        @constraint(m, [i in _V[k+1]], dbar[i] <= (tmax[k]-t[k])*z[k])
                #    end
                # else
                @constraint(m,[i in _V[k+1] ], sum( _D[i,j]*x[j] for j in _D[i,:].nzind) <= t[k])
            end
        else
            _N[k+1] = firstIndices[k]:firstIndices[k+1]-1
        end
    end

    addMostViolated!(m,LAZY_CONS_PERORGAN_INIT_NUM,xinit,t,t,β)
    @show length(_V)

    if β > 0 # penalty coefficient of DV related term in objective
        dbar = m[:dbar]
        @objective(m, Max, g-β*sum(dbar[i]^SURPLUS_VAR_OBJ_NORM for k=2:length(_V),i in _V[k]))
    else
        @objective(m, Max, g)
    end
    println("Init Model End......................")
    #write_to_file(m,"DBGModel.lp")
    return m
end


function computeProjections(γ, phi_under, phi_bar, uniform_gamma = 0)
    n, nn = size(γ)
    g = SimpleWeightedGraph(γ) #n,1:n,1:n,
    println("Started all-pairs-shortest path computations......................")
    dists = zeros(n,n)
    if uniform_gamma==0
        fws = Parallel.floyd_warshall_shortest_paths(g)
        dists = fws.dists
    else
        for i=1:n
            dists[i,:]=gdistances(g,i)#; sort_alg=RadixSort)
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


function solveModel!(m,firstIndices)
    println("In solveModel!")
    @time optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("********** Optimal Objective Function value: ", JuMP.objective_value(m))
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
    #printDoseVolume(m,doseVol)
    return m
end

function addMostViolated!(m, n, x, t, tmax, β, forDVOnly = false)
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
        indicesToAdd = []
        if (length(_N[k+1]) > 0 && !forDVOnly) || tmax[k] > t[k]
            z = []
            voxIdxs = []
            if tmax[k] == t[k]
                voxIdxs = _N[k+1]
            else
            #if tmax[k] > t[k]
                z = m[:z]
                for i in _V[k+1]
                    if is_fixed(z[i])
                        push!(voxIdxs,i)
                    end
                end
            end
            viol = vec(_D[voxIdxs,:]*x .- t[k])
            n_min=min(n,length(viol))
            violIdxs = partialsortperm(viol,1:n_min,rev=true)
            max_viol_org=viol[violIdxs[1]]
            max_viol_aor=max(max_viol_org,max_viol_aor)
            @show viol[violIdxs[n_min]] viol[violIdxs[1]]

            first_not_viol_ind = findfirst((viol[violIdxs] .<= 0.0))
            #initialization of _V[k+1] is neccesary to prevent use of the same pointer
            if length(_V[k+1]) == 0 && t[k]==tmax[k] #if the first add all even if not violated
                indicesToAdd = _N[k+1][violIdxs]
                _V[k+1] = indicesToAdd
                deleteat!(_N[k+1], sort!(violIdxs))
            else
                if  first_not_viol_ind == nothing
                    first_not_viol_ind = n_min+1
                end
                if first_not_viol_ind>1
                    violIdxs = sort!(violIdxs[1:first_not_viol_ind-1])
                    if tmax[k]==t[k]
                        indicesToAdd = _N[k+1][violIdxs]
                        append!(_V[k+1], indicesToAdd) #union!(_V[k+1],indicesToAdd)
                        deleteat!(_N[k+1], violIdxs)
                    else
                        indicesToAdd = voxIdxs[violIdxs]
                    end
                end
            end
            num_const_added_aor += length(indicesToAdd)
            xx = m[:x]
            dbar = []

            if t[k]<tmax[k] || β > 0
                dbar = m[:dbar]
            end
            for l=1:length(indicesToAdd)
                if β>0 || t[k]<tmax[k]
                    dbar[indicesToAdd[l]] = @variable(m,lower_bound=0, start = viol[violIdxs[l]])
                end
                if t[k]<tmax[k]
                    @assert(tmax[k]-t[k] >= 1)
                    unfix(z[indicesToAdd[l]])
                    set_start_value(z[indicesToAdd[l]],1/(tmax[k]-t[k]))
                    set_lower_bound(z[indicesToAdd[l]],0)
                    set_upper_bound(z[indicesToAdd[l]],1)
                end
            end
            if β > 0
                @constraint(m,[i in indicesToAdd], sum( _D[i,j]*xx[j] for j in _D[i,:].nzind) - dbar[i] <= t[k])
                obj = objective_function(m, QuadExpr)
                @objective(m, Max, obj-β*sum(dbar[i]^SURPLUS_VAR_OBJ_NORM for i in indicesToAdd))
            elseif t[k]<tmax[k]
                @constraint(m,[i in indicesToAdd], sum( _D[i,j]*xx[j] for j in _D[i,:].nzind) - dbar[i] <= t[k])
                @constraint(m, [i in indicesToAdd], dbar[i] <= (tmax[k]-t[k])*z[i])
            else
                @constraint(m,[i in indicesToAdd], sum( _D[i,j]*xx[j] for j in _D[i,:].nzind) <= t[k])
            end
        end
        println("Number of voxels in _V ", length(_V[k+1]), " voxels in _N: ",  length(_N[k+1])  , " OAR: ", k+1, " cons added: ", length(indicesToAdd))
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

function robustCuttingPlaneAlg(Din,firstIndices,t,tmax,dvrhs,β,μ, γ, gamma_const, phi_under,phi_bar, L=1, bLOAD_FROM_FILE=false)
    @assert(isempty(dvrhs) || sum(tmax-t)>0)
    ptvN, nn = size(γ)
    phi_u_n=[]
    phi_b_n=[]
    dists=[]
    file_name=@sprintf("Projection_Gamma_%1.3f.jld2",UNIFORM_GAMMA)
    if bLOAD_FROM_FILE
        phi_u_n, phi_b_n, dists = FileIO.load(file_name,"phi_u_n","phi_b_n","dists")
    else
        phi_u_n, phi_b_n, dists = computeProjections(γ, phi_under, phi_bar, UNIFORM_GAMMA)
        FileIO.save(file_name,"phi_u_n",phi_u_n,"phi_b_n",phi_b_n,"dists",dists)
    end
    m = initModel(Din,firstIndices,t,dvrhs,β,phi_u_n)
    iter=0
    addedDVCons = false
    stage = 1
    sum_num_const_added_aor = LAZY_CONS_PERORGAN_INIT_NUM  # this single counter is initialized to the initial number per organ??
    while(true)
        iter=iter+1
        prevObj = BIG_OBJ
        newObj = BIG_OBJ
        prev_viol_aor = BIG_NUM
        max_viol_aor = BIG_NUM
        num_const_added_wdv = BIG_NUM
        @time for it = 1:MAX_LAZY_CON_IT
            solveModel!(m,firstIndices)
            newObj=JuMP.objective_value(m)
            if ((prevObj-newObj)/prevObj < LAZYCONS_TIGHT_TH && abs(prev_viol_aor-max_viol_aor)/prev_viol_aor < MAX_VIOL_RATIO_TH && max_viol_aor<=MAX_VIOL_EPS_INIT && stage==1) || (max_viol_aor<=MAX_VIOL_EPS && (isempty(dvrhs) || stage ==1 || num_const_added_wdv == 0))
                println("Terminating lazy cons loop at It= ", it)
                flush(stdout)
                break
            end
            if stage == 2 && !isempty(dvrhs)
                # add constraints with dbar variables and t RHS only for organs k with tmax[k] > t[k]
                max_viol_dev , num_const_added_wdv = addMostViolated!(m, LAZY_CONS_PERORGAN_NUM_PER_ITER, JuMP.value.(m[:x]), t, tmax,β,true)
                println("Max violation wrt to t in _V: ", max_viol_dev, " Num of cons with dev added: ", num_const_added_wdv)
            end
            max_viol_aor, num_const_added_aor = addMostViolated!(m, LAZY_CONS_PERORGAN_NUM_PER_ITER, JuMP.value.(m[:x]), tmax, tmax,β)
            sum_num_const_added_aor+=num_const_added_aor
            @show max_viol_aor, num_const_added_aor
            prev_viol_aor = max_viol_aor
            prevObj = newObj
        end
        println("Total of ", sum_num_const_added_aor, " lazy constraints added so far")
        println("Reduced the objective by ", (prevObj-newObj)/prevObj)
        println("Constraints violation by ", (prev_viol_aor-max_viol_aor)/prev_viol_aor)
        println("Outer loop Iter=", iter, " Objective value: ", newObj, " *********")
        consCollect = SortedMultiDict{Float64,Array{Int64,1}}()
        x = m[:x]
        xx = value.(x)
        d = _D[1:ptvN,:]*xx
        @show minimum(d), maximum(d)
        @time min_hom_viol, max_hom_viol=addHomogenityConstr!(m,consCollect,ptvN,μ,L,phi_u_n,phi_b_n,dists,d)
        if stage==2 && max_hom_viol<=MAX_VIOL_EPS && max_viol_aor<=MAX_VIOL_EPS && (isempty(dvrhs) || num_const_added_wdv == 0 )# no violated inequalities found
            println("Terminating cut algorithm.. iter=", iter)
            break
        elseif !isempty(consCollect)
            stage=1
            println("Max Violation: ", max_hom_viol, " " ,min_hom_viol)
            println("Adding ", length(consCollect), " cuts.....................")
        else
            if stage == 1
                println("Switching to stage 2 ******")
                #zVar = variable_by_name(m, "z")
                if !isempty(dvrhs) && !addedDVCons
                    initDoseVolume(m,t,tmax,dvrhs)
                    addedDVCons = true
                end
            end
            stage=2
            #set_optimizer_attribute(m, "OptimalityTol", 1e-3)  # tighten tolerance in 2nd phase of the algorithm
        end

        cons_array = Any[]
        @time for (v,pair) in exclusive(consCollect,startof(consCollect),pastendsemitoken(consCollect))
            i1=pair[1]
            i2=pair[2]
            #shimrit:added both constraints
            phibar = min(phi_u_n[i2]+dists[i2,i1],phi_b_n[i1])
            v1 = phibar*d[i1]-μ*phi_u_n[i2]*d[i2]
            #@show i1, i2, phibar*d[i1], phi_u_n[i2]*d[i2]
            phiun = max(phi_b_n[i1]-dists[i1,i2],phi_u_n[i2])
            v2 = phi_b_n[i1]*d[i1]-μ*phiun*d[i2]
            if v1 > VIOL_EPS
                push!(cons_array,@build_constraint(phibar*sum(Din[i1,j]*x[j] for j in Din[i1,:].nzind)-μ*phi_u_n[i2]*sum(Din[i2,j]*x[j] for j in Din[i2,:].nzind) <= 0))
            end
            if v2 > VIOL_EPS
                push!(cons_array,@build_constraint(phi_b_n[i1]*sum(Din[i1,j]*x[j] for j in Din[i1,:].nzind) - μ*phiun*sum(Din[i2,j]*x[j] for j in Din[i2,:].nzind) <= 0))
            end
            #end
            #println("added constraint ", cons)
        end
        lenCons = length(cons_array)
        @time for i=1:lenCons
            add_constraint(m, cons_array[i])
        end
    end
    #printDoseVolume(m, t, tmax, !isempty(dvrhs), true) # print out verbose output
    return m
end

function printDoseVolume(m, t = [], tmax = [], doseVol = false, verbose = false)
    conicForm = false

    zVar = variable_by_name(m, "z")
    gVar = m[:g] #variable_by_name(m,"g")
    xVar = m[:x]
    #if zVar != nothing
    #    conicForm = true
    #end
    x = value.(xVar)
    println("Min PTV bio dose, g=",value(gVar))

    if verbose
        for k=1:length(_V)
            dose = _D[[_V[k];_N[k]],:]*x
            dev = []
            zVar = []
            if !isempty(t) && !isempty(tmax) && k>=2 && tmax[k-1] > t[k-1]
                dev = max.(dose .- t[k-1],0)
                zVar = dev./(tmax[k-1]-t[k-1])
            end
            numDev = 0
            if !isempty(dev)
                numDev = count(i->(i>DBARNZTH),dev)
            end
            sumZVar = 0
            if !isempty(zVar)
                sumZVar = sum(zVar)
                histogram(zVar[zVar .> ZNZTH])
                #savefig("zlikehist.png")
            end
            println("Structure: ", k, " Max (phys) dose: ", maximum(dose), " Min (phys) dose: ", minimum(dose), " 99%-tile: ", quantile!(dose,0.99), " vox num exceeding t: " , numDev, " sum(zVar):", sumZVar)
        end
    end

    if doseVol
        dbarVar = m[:dbar] #variable_by_name(m,"dbar")
        zVar = m[:z]
        zVarVal = []
        if zVar != nothing
            zVarVal = value.(zVar)
            #I = axes(zVarVal,Axis{1})
            nzZVal = zVarVal.data[zVarVal.data .> ZNZTH]
            histogram(nzZVal)
            #savefig("zhist.png")
        else
            println("z var not created")
        end

        for k=2:length(_V)
            if t[k-1] < tmax[k-1]
                indices = _V[k]
                @show length(dbarVar)
                dev = max.(_D[indices,:]*x.-t[k-1],0)
                nzNum = 0
                nzDvar = 0
                nzZvar = 0
                for i = 1:length(indices)
                    if dev[i] > DBARNZTH
                        nzNum+=1
                        if value(dbarVar[indices[i]]) > DBARNZTH
                            nzDvar+=1
                        end
                    end
                    if zVar!=nothing && zVarVal[indices[i]] > ZNZTH
                        nzZvar+=1
                    end
                end
                println("OAR k=", k, " Number of deviations from t: ", nzNum)
                println("OAR k=", k, " Number of nonzero dbar: ", nzDvar)
                if zVar != nothing
                    println("OAR k=", k, " Number of nonzero z: ", nzZvar, " sum z: ", sum(zVarVal.data))
                end
            end
        end
    end
    return
end

end
