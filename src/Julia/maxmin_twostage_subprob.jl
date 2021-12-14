module maxmin_twostage_subprob
using MAT
using JuMP
using Gurobi
#using SCS
#using UnicodePlots
using SparseArrays
using LinearAlgebra
using LightGraphs, SimpleWeightedGraphs
using SortingAlgorithms
using Statistics
using DataStructures
import LightGraphs.Parallel
using FileIO, JLD2
using Printf

const __DEBUG = false
#using Plots
const INITXNORM = 100

const VIOL_EPS = 1e-2
const DBARNZTH = 1e-4
const ZNZTH = 1e-4
const BIG_OBJ = 1e8
const UNIFORM_GAMMA = true

const LAZY_CONS_PERORGAN_TH = 1e5 #5e4 # 5e3 #1e5 #5e4
const LAZY_CONS_PERORGAN_INIT_NUM = 400
const LAZY_CONS_PERORGAN_NUM_PER_ITER = 200
const MAX_LAZY_CON_IT = 500

const LAZYCONS_TIGHT_TH = 0.2
const MAX_VIOL_RATIO_TH = 0.2

const MAX_V_CONS = 10 #can be set to infinity
const MAX_VIOL_EPS = 1e-2
const MAX_VIOL_EPS_INIT = 10

const SURPLUS_VAR_OBJ_NORM = 1  # penalty norm either 1 or 2

const BIG_NUM = 1e6

const OPTIMALITY_TOL = 1e-5

global _D   # save D in order to load rows as needed
#global _rowLoc # 0 - not loaded into optimization problem, rowLoc[i] > 0 indicates row in model constraint matrix
global _V
global _N # voxels by organ not loaded into optimization


export initModel, solveModel!,robustCuttingPlaneAlg, printDoseVolume, computeProjections, parametricSolveIncreasing
#optimizer_constructor = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 10, "verbose" => 0)
#set_optimizer(problem, optimizer_constructor)
# Load Optimizer and model


function initDoseVolume(m, t, tmax, dvrhs)
    println("initDoseVolume....")
    for k =1:length(t)
        if tmax[k] > t[k]
            @variable(m, z[[_V[k+1];_N[k+1]]]==0)
            @constraint(m,dos_vol[k],sum(z[i] for i in _V[k+1] ) + sum(z[i] for i in _N[k+1] ) <= dvrhs[k] )
            println("Adding dose vol cons for organ k=", k, " dvrhs[k]=", dvrhs[k])
        end
    end
    return m
end

#
function initModel(Din,firstIndices,t,tmax,dvrhs,β,phi_u_n,λ=[0;0],ϕ_nom=0,xinit=[])

    n, nb = size(Din)
    #_rowLoc = spzeros(n,1)
    global _V = fill(Int[],length(firstIndices)+1)
    global _N = fill(Int[],length(firstIndices)+1)
    _V[1] = 1:firstIndices[1]-1   # PTV
    #lastRow = firstIndices[1]-1
    if __DEBUG
        println("Init Model Start......................")
    end
    optimizer = Gurobi.Optimizer #SCS.Optimizer
    m = Model(optimizer)
    set_optimizer_attribute(m, "OutputFlag", 0)
    set_optimizer_attribute(m, "OptimalityTol", OPTIMALITY_TOL)
    ptvn = length(_V[1])
    @variable(m,g)
    if isempty(phi_u_n)  # if not given then initialize phi to unity to solve problem with physical dose
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
        oarIdxs = firstIndices[k]:firstIndices[k+1]-1
        if tmax[k] > t[k] && !isempty(dvrhs)
            @variable(m, z[oarIdxs]==0)
            @constraint(m,dos_vol[k],sum(z[i] for i in oarIdxs) <= dvrhs[k] )
        end

        if firstIndices[k+1]-firstIndices[k] <= LAZY_CONS_PERORGAN_TH
            _V[k+1] =  oarIdxs#firstIndices[k]:firstIndices[k+1]-1
            #if (!isempty(dvrhs) &&  t[k] < tmax[k] ) ||
            if β>0
                for i in _V[k+1]
                    dbar[i] = @variable(m,lower_bound=0,upper_bound=tmax[k]-t[k])
                end
                @constraint(m,[i in _V[k+1] ], sum( _D[i,j]*x[j] for j in _D[i,:].nzind) - dbar[i] <= t[k])
            else
                if t[k] < tmax[k] && !isempty(dvrhs)
                    unfix(z[_V[k+1]])
                    set_lower_bound(z[_V[k+1]],0)
                    set_upper_bound(z[_V[k+1]],1)
                    @constraint(m, [i in _V[k+1]], dbar[i] <= (tmax[k]-t[k])*z[k])
                    @constraint(m,[i in _V[k+1] ], sum( _D[i,j]*x[j] for j in _D[i,:].nzind) - dbar[i] <= t[k])
                else
                     @constraint(m,[i in _V[k+1] ], sum( _D[i,j]*x[j] for j in _D[i,:].nzind) <= t[k])
                end
            end
        else
            _N[k+1] = oarIdxs #firstIndices[k]:firstIndices[k+1]-1
        end
    end

    addMostViolated!(m,LAZY_CONS_PERORGAN_INIT_NUM,xinit,t,tmax,β,false,isempty(dvrhs))
    if __DEBUG
        @show length(_V)
    end
	if λ[1] > 0
		@expression(m, reg1, λ[1]*sum(x[i] for i=1nb))
	else
		@expression(m, reg1, 0)
	end
	if λ[2] > 0
		@variable(m,g_nom)
		print(size(ϕ_nom))
		@constraint(m,cons_ptv_reg2[i in _V[1]], g_nom <= ϕ_nom[i]*Din[i,:].nzval'*x[Din[i,:].nzind])
		@expression(m, reg2, λ[2]*g_nom)
	else
		@expression(m, reg2, 0)
	end

    if β > 0 # penalty coefficient of DV related term in objective
        dbar = m[:dbar]
#        @objective(m, Max, g-β*sum(dbar[i]^SURPLUS_VAR_OBJ_NORM for k=2:length(_V),i in _V[k]))
        @objective(m, Max, g + reg1 + reg2 -β*sum(dbar[i]^SURPLUS_VAR_OBJ_NORM for k=2:length(_V),i in _V[k]))
    else
        @objective(m, Max, g + reg1 + reg2)
    end
    println("Init Model End......................")
    #write_to_file(m,"DBGModel.lp")
    return m
end


function computeProjections(γ, gamma_func, phi_under, phi_bar,dists=[])
    n, nn = size(γ)
    g = SimpleWeightedGraph(γ) #n,1:n,1:n,
#    println("Started all-pairs-shortest pth computations")
	if dists==[]
    	dists = zeros(n,n)
    	if UNIFORM_GAMMA==0
        	fws = Parallel.floyd_warshall_shortest_aths(g)
        	dists = fws.dist
    	else
        	for i=1:n
            	dists[i,:]=gdistances(g,i)#; sort_alg=RadixSort)
				for j=1:n
					dists[i,j]=gamma_func(dists[i,j])
				end
        	end
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
    phi_under_n = maximum(phi_under*ones(1,n)-dists,dims=1)'
    phi_bar_n = minimum(phi_bar*ones(1,n)+dists,dims=1)'
	@show minimum(phi_bar_n-phi_under_n)
    @assert(all(phi_under_n .<= phi_bar_n))
    println("Finished computeProjections")
    return phi_under_n,phi_bar_n, dists
end


function solveModel!(m,firstIndices)
    if __DEBUG
        println("In solveModel!")
    end
    @time optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        if __DEBUG
            println("********** Optimal Objective Function value: ", JuMP.objective_value(m))
        end
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

function getValidBetaInterval(m)
    dbar = m[:dbar]
    report = lp_sensitivity_report(m)
    minMaxBeta = Inf
    maxMinBeta = 0
    for i=1:length(dbar)
        dbarlo, dbarhigh = report[dbar[i]]
        if dbarhigh < minMaxBeta
            minMaxBeta = dbarhigh
        end
        if dbarlo > maxMinBeta
            maxMinBeta = dbarlo
        end
    end
    return maxMinBeta,minMaxBeta
end

function addMostViolated!(m, n, x, t, tmax, β, fromVOnly = false, noDVCons = false)
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
        if (length(_N[k+1]) > 0 && !fromVOnly) || tmax[k] > t[k]
            z = []
            #if tmax[k] == t[k]
            #if !fromVOnly
            voxIdxs = []
            if !fromVOnly
                voxIdxs = _N[k+1]
            end
            #else
            #else #if tmax[k] > t[k]
            if tmax[k] > t[k] && !noDVCons
                z = m[:z]
                if fromVOnly
                    for i in _V[k+1]
                        if is_fixed(z[i])
                            push!(voxIdxs,i)
                        end
                    end
                end
            end

            viol = vec(_D[voxIdxs,:]*x .- t[k])
            n_min=min(n,length(viol))
            violIdxs = partialsortperm(viol,1:n_min,rev=true)
            max_viol_org=viol[violIdxs[1]]
            max_viol_aor=max(max_viol_org,max_viol_aor)
            if __DEBUG
                @show viol[violIdxs[n_min]] viol[violIdxs[1]]
            end
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
                    if !fromVOnly # tmax[k]==t[k]
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
                    dbar[indicesToAdd[l]] = @variable(m,lower_bound=0, upper_bound=tmax[k]-t[k],start = viol[violIdxs[l]])
                    if t[k]<tmax[k] && !noDVCons
                        @assert(tmax[k]-t[k] >= 1)
                        unfix(z[indicesToAdd[l]])
                        set_start_value(z[indicesToAdd[l]],1/(tmax[k]-t[k]))
                        set_lower_bound(z[indicesToAdd[l]],0)
                        set_upper_bound(z[indicesToAdd[l]],1)
                    end
                end
            end
            if β > 0
                obj = objective_function(m, QuadExpr)
                #@objec#tive(m, M#ax, ob#j-β*sum(dbar[i]^SURPLUS_VAR_OBJ_NORM for i in indicesToAdd))
                @objective(m, Max, obj-β*sum(dbar[i] for i in indicesToAdd))
            end

            if t[k]<tmax[k] && (!noDVCons || β > 0)
                @constraint(m,[i in indicesToAdd], sum( _D[i,j]*xx[j] for j in _D[i,:].nzind) - dbar[i] <= t[k])
            else
                @constraint(m,[i in indicesToAdd], sum( _D[i,j]*xx[j] for j in _D[i,:].nzind) <= t[k])
            end
            if t[k]<tmax[k] && !noDVCons
                @constraint(m, [i in indicesToAdd], dbar[i] <= (tmax[k]-t[k])*z[i])
            end
        end
        if __DEBUG
            println("Number of voxels in _V ", length(_V[k+1]), " voxels in _N: ",  length(_N[k+1])  , " OAR: ", k+1, " cons added: ", length(indicesToAdd))
        end
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


function evaluateDevNum(dbar)
    cntVec = zeros(length(_V)-1,1)
    for k=2:length(_V)
        for i in _V[k]
            if (dbar[i] > DBARNZTH)
                cntVec[k-1] += 1
            end
        end
    end
    return cntVec
end

# parametricSolveIncreasing -
function parametricSolveIncreasing(Din,firstIndices,t,tmax,dvrhs,β,μ, phi_u_n, phi_b_n, dists, ϕ, L=1, βstop=BIG_NUM)
    m = initModel(Din,firstIndices,t,tmax,[],β,phi_u_n)
    βPrev = β
    while (β<βstop)
        m, htCn, homCn = robustCuttingPlaneAlg(Din,firstIndices,t,tmax,[],β,μ, phi_u_n, phi_b_n, dists, [0;0], 0, 200,m)
        dbar = m[:dbar]
        devVec = evaluateDevNum(dbar) # function that returns deviation vector with dimension of OAR num
        if (all(devVec.<dvrhs))
            break
        end
        βLb,βUb = getValidBetaInterval(m)
        if (htCn+homCn > 0 && βLb > βPrev) # if generated cuts and lower bound for validity of basis exceeds the previous beta
            β = βLb - eps
        else
            β = βUb + eps
        end
    end
    return m
end

function robustCuttingPlaneAlg(Din,firstIndices,t,tmax,dvrhs,β,μ, phi_u_n, phi_b_n, dists, λ, ϕ, L=1, m=[])

    @assert(isempty(dvrhs) || sum(tmax-t)>0)

    ptvN = firstIndices[1]-1


    if isempty(m)
        m = initModel(Din,firstIndices,t,tmax,dvrhs,β,phi_u_n,λ,ϕ)
    end
    iter=0
    addedDVCons = false
    stage = 1
    sum_num_const_added_aor = 0 #LAZY_CONS_PERORGAN_INIT_NUM  # this single counter is initialized to the initial number per organ??
    sum_num_const_added_hom = 0

    while(true)
        iter=iter+1
        prevObj = BIG_OBJ
        newObj = BIG_OBJ
        prev_viol_aor = BIG_NUM
        max_viol_aor = BIG_NUM
        num_const_added_wdv = BIG_NUM
        #@time
        for it = 1:MAX_LAZY_CON_IT
            @time solveModel!(m,firstIndices)
            newObj=JuMP.objective_value(m)
              ########### exit inner loop if relative decrease in obj fun less than threshold or absolute diff and in the first stage or otherwise if violation less than threshold without dv cons or in dose vol
            if ((prevObj-newObj)/prevObj < LAZYCONS_TIGHT_TH && abs(prev_viol_aor-max_viol_aor)/prev_viol_aor < MAX_VIOL_RATIO_TH && max_viol_aor<=MAX_VIOL_EPS_INIT && stage==1)   #&& (isempty(dvrhs) || stage ==1 || num_const_added_wdv == 0))
                println("Terminating lazy cons loop at It= ", it, "  Infeas reduction: ", (prev_viol_aor-max_viol_aor)/prev_viol_aor,
                " Obj red %: ", (prevObj-newObj)/prevObj)
                flush(stdout)
                break
            end
            #if stage == 2 && !isempty(dvrhs)
                # add constraints with dbar variables and t RHS only for organs k with tmax[k] > t[k]
            #    max_viol_dev , num_const_added_wdv = addMostViolated!(m, LAZY_CONS_PERORGAN_NUM_PER_ITER, JuMP.value.(m[:x]), t, tmax,β,true)
            #    println("Max violation wrt to t in _V: ", max_viol_dev, " Num of cons with dev added: ", num_const_added_wdv)
            max_viol_aor, num_const_added_aor = addMostViolated!(m, LAZY_CONS_PERORGAN_NUM_PER_ITER, JuMP.value.(m[:x]), t, tmax,β,false,isempty(dvrhs))
            if max_viol_aor<=MAX_VIOL_EPS
                println("Terminating lazy cons loop at It= ", it)
                break
            end
            #else
            #    max_viol_aor, num_const_added_aor = addMostViolated!(m, LAZY_CONS_PERORGAN_NUM_PER_ITER, JuMP.value.(m[:x]), tmax, tmax,β)
            #end
            sum_num_const_added_aor += num_const_added_aor
            if __DEBUG
                @show max_viol_aor, num_const_added_aor
            end
            prev_viol_aor = max_viol_aor
            prevObj = newObj
        end
#        println("Total of ", sum_num_const_added_aor, " lazy constraints added so far")
#        println("Reduced the objective by ", (prevObj-newObj)/prevObj)
#        println("Constraints violation by ", (prev_viol_aor-max_viol_aor)/prev_viol_aor)
        println("Outer loop Iter=", iter, " Objective value: ", newObj, " OAR lazy cons: ", sum_num_const_added_aor, " Hom lazy cons: ",sum_num_const_added_hom)
        consCollect = SortedMultiDict{Float64,Array{Int64,1}}()
        x = m[:x]
        xx = value.(x)
        d = _D[1:ptvN,:]*xx
        if __DEBUG
            @show minimum(d), maximum(d)
        end
        #@time
        min_hom_viol, max_hom_viol=addHomogenityConstr!(m,consCollect,ptvN,μ,L,phi_u_n,phi_b_n,dists,d)
        if stage==2 && max_hom_viol<=MAX_VIOL_EPS && max_viol_aor<=MAX_VIOL_EPS #&& (isempty(dvrhs) || num_const_added_wdv == 0 )# no violated inequalities found
            println("Terminating cut algorithm.. iter=", iter)
            break
        elseif !isempty(consCollect)
            stage=1
            if __DEBUG
                println("Max Violation: ", max_hom_viol, " " ,min_hom_viol, " Adding ", length(consCollect), " cuts..")
            end
        else
            if stage == 1
                println("Switching to stage 2 ******")
                #if !isempty(dvrhs) && !addedDVCons
                #    initDoseVolume(m,t,tmax,dvrhs)
                #    addedDVCons = true
                #end
            end
            stage=2
            #set_optimizer_attribute(m, "OptimalityTol", 1e-3)  # tighten tolerance in 2nd phase of the algorithm
        end

        cons_array = Any[]
        for (v,pair) in exclusive(consCollect,startof(consCollect),pastendsemitoken(consCollect))
        #@time
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
        #@time
        for i=1:lenCons
            add_constraint(m, cons_array[i])
        end
        sum_num_const_added_hom += lenCons
    end # robust cuts loop
    #printDoseVolume(m, t, tmax, !isempty(dvrhs), true) # print out verbose output
    return m, sum_num_const_added_aor, sum_num_const_added_hom
end

function printDoseVolume(m, t = [], tmax = [], doseVol = false, verbose = false)
    conicForm = false
    #zVar = m[:z] #variable_by_name(m, "z")
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
                nzZ = zVar[zVar .> ZNZTH]
                #if !isempty(nzZ)
                #    histogram(nzZ)
                #end
                #savefig("zlikehist.png")
            end
            println("Structure: ", k, " Max (phys) dose: ", maximum(dose), " Min (phys) dose: ", minimum(dose), " 99%-tile: ", quantile!(dose,0.99), " vox num exceeding t: " , numDev, " sum dbar/(tmax-t) :", sumZVar)
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
            if !isempty(nzZVal)
            #    histogram(nzZVal)
            end
            #savefig("zhist.png")
        else
            println("z var not created")
        end

        for k=2:length(_V)
            if t[k-1] < tmax[k-1]
                indices = _V[k]
                if __DEBUG
                    @show length(dbarVar)
                end
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
