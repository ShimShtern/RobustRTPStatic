module maxmin_twostage_subprob
#using MAT
using JuMP
#using CPLEX
using Gurobi
#using SCS
#using UnicodePlots
using DataFrames
using SparseArrays
using LinearAlgebra
using LightGraphs, SimpleWeightedGraphs, GraphBLASInterface, SuiteSparseGraphBLAS,LightGraphsGraphBLAS
using SortingAlgorithms
using Statistics
using DataStructures
import LightGraphs.Parallel
using FileIO, JLD2
using Printf

const __ALLCONS_NO_GENERATION = true #adds all homogen constraint, false generate constraints

const NO_DEBUG = 0
const DEBUG_LOW = 1
const DEBUG_MED = 2
const DEBUG_HIGH = 3
const __DEBUG = DEBUG_LOW #NO_DEBUG #NO_DEBUG  #true   #DEBUG_LOW#
const __SOLVER_DEBUG = 0 #1

#using Plots
const INITXNORM = 100
const VIOL_EPS = 1e-4 #allowed violation for homogeneity constraints
const INFEAS_TOL = 1e-5
const DBARNZTH = 2e-6
const ZNZTH = 1e-4
const BIG_OBJ = 1e8
const UNIFORM_GAMMA = true

const LAZY_CONS_PERORGAN_TH = 5e3 #5e4 # 5e3 #1e5 #5e4
global LAZY_CONS_PERORGAN_INIT_NUM = 0 # 2000 #2000 #4000000 #4000 #4000000		#maximum number of constraints for OAR added in model initialization

const LAZY_CONS_PERORGAN_NUM_PER_ITER = 400 #maximum number of constraints for OAR added in each iteration
const MAX_LAZY_CON_IT = 1e6		#maximum number of iterations done for adding OAR constraints

const LAZYCONS_TIGHT_TH = 0.01
const MAX_VIOL_RATIO_TH = 0.01


MAX_V_CONS = 5 #1000000000 #10 #100000000 #Inf #10 #can be set to infinity
const MAX_VIOL_EPS = 1e-2 #oar constraint violation allowed NOT USED
const MAX_VIOL_EPS_INIT = 10  #initial oar constraint violation allowed in phase I / to generate homogeneity constraints

const BETA_EPS = 5e-9 # used for parametricSolve routines
const BETA_DELTA_NZ_TH = BETA_EPS-eps()
const BISECTION_GAP = 50

const SMALL_EPS = 1e-16 # used for sensitivity / basis validity interval calculations


const DEV_VAR_OBJ_NORM = 1  # penalty norm either 1 or 2

const BIG_NUM = 1e6

const OPTIMALITY_TOL = 1e-6
const FEASIBILITY_TOL = 1e-6

const W_th = 0.5 # bisection weight

const BASIS_EPS = 1e-10 #determines minimal nonzero value to include in basis

global const GRB_ENV = Gurobi.Env()

global _D   # save D in order to load rows as needed
#global _rowLoc # 0 - not loaded into optimization problem, rowLoc[i] > 0 indicates row in model constraint matrix
global _V
global _N # voxels by organ not loaded into optimization
#global _dbarIdxToGrbIdx
#global _GrbIdxToDbarIdx

export initModel, solveModel!,robustCuttingPlaneAlg!, printDoseVolume, computeProjections, parametricSolveIncreasing, CVAR_solve
#optimizer_constructor = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 10, "verbose" => 0)
#set_optimizer(problem, optimizer_constructor)
# Load Optimizer and model

struct GlobalIndexes
	 _V::Vector{Vector{Int}}
      _N::Vector{Vector{Int}}
      #_dbarIdxToGrbIdx::Dict{Int64,Cint}
	  #_GrbIdxToDbarIdx::Dict{Cint,Int64}
end

struct Solution
	 	gvalue::Float64
       basicXIdxs::Vector{Int64}
       basicDeltaIdxs::Vector{Int64}
	   XValue::Vector{Float64}
	   DbarDict::Dict{Int64,Float64}
end

struct Solution2
       basicXIdxs::Vector{Int64}
       basicDeltaIdxs::Vector{Int64}
	   XValue::Vector{Float64}
	   DbarDict::Dict{Int64,Vector{AffExpr}}
end

function clearGlobalRunData()
	_V = nothing
	_N = nothing
#	_dbarIdxToGrbIdx = nothing
#	_GrbIdxToDbarIdx = nothing
	GC.gc()
end

function getGlobals()
	global _V
	global _N
	#global _dbarIdxToGrbIdx
	#global _GrbIdxToDbarIdx
	IndexStr=GlobalIndexes(deepcopy(_V),deepcopy(_N))#,deepcopy(_dbarIdxToGrbIdx),deepcopy(_GrbIdxToDbarIdx))
	return IndexStr
end

function initGlobals(Din,firstIndices)
	global _V = fill(Int[],length(firstIndices)+1)
    global _N = fill(Int[],length(firstIndices)+1)
	global _D = Din
	n, nb = size(Din)
	_V[1] = 1:firstIndices[1]-1   # PTV
	for k = 1:length(firstIndices)-1
        _V[k+1] = firstIndices[k]:firstIndices[k+1]-1
    end
end


function GetSol(m)
	basicXIdxs, basicDeltaIdxs = getBasisXandDelta(m)
	dbarVarDict=m[:dbar]
	XValue=value.(m[:x])
	dbarValues=[value(dbarVarDict[i]) for i in keys(dbarVarDict)]
	dbarDict=Dict(keys(dbarVarDict) .=> dbarValues)
	gvalue=value(m[:g])
	Sol = Solution(gvalue,basicXIdxs,collect(basicDeltaIdxs),XValue,dbarDict)
	return Sol
end

function initModel(Din, firstIndices, t, tmax, dvrhs, β, phi_u_n, λ=[0;0], ϕ_nom=0, xinit=[], alwaysCreateDeltaVars=false, idxOfGreaterThanCons=1)
    n, nb = size(Din)
    #_rowLoc = spzeros(n,1)
    global _V = fill(Int[],length(firstIndices)+1)
    global _N = fill(Int[],length(firstIndices)+1)
    _V[1] = 1:firstIndices[1]-1   # PTV
	@assert(β>=0 && all(t<=tmax))
    if __DEBUG >= DEBUG_LOW
        println("Init Model Start......................")
		@show t,tmax,β
    end
#	global _dbarIdxToGrbIdx = Dict{Int64,Cint}()
#	global _GrbIdxToDbarIdx = Dict{Cint,Int64}()
    #m = Model(() ->SCS.Optimizer())
	#m = Model(()-> CPLEX.Optimizer())
	#m = Model(() -> Gurobi.Optimizer(GRB_ENV))
	m = direct_model(Gurobi.Optimizer(GRB_ENV))
	#grb = unsafe_backend(m)

	MOI.set(m, MOI.Silent(), true)
	MOI.set(m, MOI.NumberOfThreads(), 1)
	#set_optimizer_attribute(m,"CPX_PARAM_LPMETHOD",4)
	#set_optimizer_attribute(m,"CPX_PARAM_BAREPCOMP",1e-4)
	if __SOLVER_DEBUG > 0
    	set_optimizer_attribute(m, "OutputFlag", __SOLVER_DEBUG)
	end
    set_optimizer_attribute(m, "OptimalityTol", OPTIMALITY_TOL)
    set_optimizer_attribute(m, "FeasibilityTol", FEASIBILITY_TOL)
	set_optimizer_attribute(m, "NumericFocus", 3)
#	set_optimizer_attribute(m, "Quad", 1)
#	set_optimizer_attribute(m, "MarkowitzTol", 0.05)

    ptvn = length(_V[1])
    @variable(m,g)

    if isempty(phi_u_n)  # if not given then initialize phi to unity to solve problem with physical dose
        phi_u_n = ones(ptvn,1)
        phi_b_n = ones(ptvn,1)
    end
    #htn = n - firstIndices[1] + 1 # number of healthy tissue voxels
    @variable(m,x[1:nb]>=0)
    @constraint(m,cons_ptv[i in _V[1]], phi_u_n[i]*Din[i,:].nzval'*x[Din[i,:].nzind]-g >= 0)
    global _D = Din
    # select initial rows according to max violation with given initial soln vector
    if isempty(xinit)
        xinit = INITXNORM/nb*ones(nb,1)
    end
    dbar = m[:dbar] = Dict()

    firstIndices = [firstIndices; n+1]
    for k = 1:length(firstIndices)-1
        oarIdxs = firstIndices[k]:firstIndices[k+1]-1
        #if tmax[k] > t[k] && !isempty(dvrhs)
		#	@error("this mode with z vars is not supported")
        #    @variable(m, z[oarIdxs]==0)
        #    @constraint(m,dos_vol[k],dvrhs[k] - sum(z[i] for i in oarIdxs) >= 0)
        #end
        if firstIndices[k+1]-firstIndices[k] <= LAZY_CONS_PERORGAN_TH
            _V[k+1] =  oarIdxs   #firstIndices[k]:firstIndices[k+1]-1
            if t[k] < tmax[k]
                @assert(β>0 || !isempty(dvrhs) || alwaysCreateDeltaVars)
                for i in _V[k+1]
					if alwaysCreateDeltaVars
						error("the mode alwaysCreateDeltaVars is not currently supported")
						dbar[i] = @variable(m) #,lower_bound=0,upper_bound=tmax[k]-t[k])
					    fix(dbar[i],0)
					else
                    	dbar[i] = @variable(m,lower_bound=0,upper_bound=tmax[k]-t[k])
						#grbIdx = Gurobi.column(grb,index(dbar[i]))
						#_dbarIdxToGrbIdx[i] = grbIdx
						#_GrbIdxToDbarIdx[grbIdx] = i
					end
                end
                @constraint(m,[i in _V[k+1] ], t[k]-sum( _D[i,j]*x[j] for j in _D[i,:].nzind) + dbar[i] >= 0)
                if !isempty(dvrhs)
					error("the mode z vars is not currently supported")
                    unfix(z[_V[k+1]])
                    set_lower_bound(z[_V[k+1]],0)
                    set_upper_bound(z[_V[k+1]],1)
                    @constraint(m, [i in _V[k+1]], (tmax[k]-t[k])*z[k]-dbar[i] >= 0)
                end
            else # t[k] = tmax[k]
				if k == idxOfGreaterThanCons
            		@constraint(m,[i in _V[k+1] ], sum( _D[i,j]*x[j] for j in _D[i,:].nzind) <= t[k])
				else
					@constraint(m,[i in _V[k+1] ], t[k]-sum( _D[i,j]*x[j] for j in _D[i,:].nzind) >= 0)
				end
            end
        else
            _N[k+1] = oarIdxs #firstIndices[k]:firstIndices[k+1]-1
        end
    end

    infeas, numAdded = addMostViolated!(m, LAZY_CONS_PERORGAN_INIT_NUM, xinit, t, tmax, β, alwaysCreateDeltaVars)
    if __DEBUG >= DEBUG_MED
        @show length(_V), infeas, numAdded
    end

	if λ[1] > 0   # ???
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
	@show  β, alwaysCreateDeltaVars
    if β > 0 || alwaysCreateDeltaVars # penalty coefficient of DV related term in objective
        dbar = m[:dbar]
#        @objective(m, Max, g-β*sum(dbar[i]^SURPLUS_VAR_OBJ_NORM for k=2:length(_V),i in _V[k]))
		println("setting up problem for β=",β)
		scaling = ComputeScaling(Din, firstIndices, t,tmax)
        #@objective(m, Max, g*scaling -β*sum(dbar[i]^DEV_VAR_OBJ_NORM for k=2:length(_V) if tmax[k-1]>t[k-1] for i in _V[k] ) + reg1 + reg2)
		#in case the quadratic slows it down
		@objective(m, Max, g*scaling -β*sum(dbar[i] for k=2:length(_V) if tmax[k-1]>t[k-1] for i in _V[k] ) + reg1 + reg2)
    else
		scaling = ComputeScaling(Din, firstIndices, t,tmax)
		if scaling==0
        	@objective(m, Max, g + reg1 + reg2)
		else
			@objective(m, Max, g*scaling + reg1 + reg2)
		end
    end
	if __DEBUG >= DEBUG_LOW
    	println("Init Model End......................")
	end
	@show length(_V[2]) length(m[:dbar])
    return m
end



function computeProjections(γ, gamma_func, phi_under, phi_bar,dists=[])
    n, nn = size(γ)
#    println("Started all-pairs-shortest pth computations")
    @assert(tr(γ)==0)
	if dists==[]
    	dists = zeros(n,n)
    	if UNIFORM_GAMMA==0
            g = SimpleWeightedGraph(γ) #n,1:n,1:n,
        	fws = Parallel.floyd_warshall_shortest_aths(g)
        	dists = fws.dist
    	else
            #GrB_init(GrB_NONBLOCKING)
            ##broadcast!(.!=,γ,0)
            #bg = BLASGraph(γ)#{Int64}(n)
			#println(γ[1:12,1:12])
            g = SimpleGraph(γ)
        	for i=1:n
				comp_dist = gdistances(g,i)
				dists[i,:]=gamma_func.(comp_dist)#; sort_alg=RadixSort)
        	end
			#println(dists[1:12,1:12])
		end
    end
    ###################################### save dists for debug
#    file = matopen("dists.mat", "w")
#    write(file, "dists", dists)
#    close(file)
    ######################################
    println("Finished all-pairs-shortest path computations")
    #phi_under_n = phi_under
    #phi_bar_n = phi_bar
    phi_under_n = maximum(phi_under*ones(1,n)-dists,dims=1)'
    phi_bar_n = minimum(phi_bar*ones(1,n)+dists,dims=1)'
	#@show minimum(phi_bar_n-phi_under_n)
    if __DEBUG >= DEBUG_LOW
		for i = 1:length(phi_under_n)
			if phi_under_n[i] > phi_bar_n[i] || phi_under_n[i]<phi_under[i] || phi_bar_n[i]>phi_bar[i]
				val, ind=findmax(phi_under-dists[i,:]',1)
				org_dist=gdistances(g,i)[ind]
				phi_dist=gamma_func(org_dist)
				phi_dist2=dists[i,ind]
				@printf("i=%d,ind=%d,phi_i=%f,phi_ind=%f,dist=%f,phi_dist1=%f,phi_dist2=%f,bphi=%f\n",i, ind,phi_under[i],phi_under[ind],org_dist,phi_dist,phi_dist2,phi_bar[ind])
				error("")
				#
				#@printf("i=%d,ind=%d,dist=%d,phi_dist1=%f,phi_dist2=%f,uphi,=%f\n",i, ind,org_dist,phi_dist,phi_dist2,phi_under[ind])
				#val, ind=findmin(phi_bar+dists[i,:],2)
				#org_dist=gdistances(g,i)[ind]
				#phi_dist=gamma_func(org_dist)
				#phi_dist2=dists[i,ind]
				println("phi_under_n[i]=", phi_under_n[i], " phi_bar_n[i]=", phi_bar_n[i])
				#error("")
			end
		end
	end
    @assert(all(phi_under_n .<= phi_bar_n))
    println("Finished computeProjections")
    return phi_under_n,phi_bar_n, dists
end


function solveModel!(m,firstIndices)
    if __DEBUG >= DEBUG_LOW
        println("In solveModel!")
    end
    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        if __DEBUG >= DEBUG_LOW
            println("********** Solved model with cons num: ", num_constraints(m, AffExpr,MOI.LessThan{Float64}) + num_constraints(m, AffExpr,MOI.GreaterThan{Float64}) , " Optimal Objective Function value: ", JuMP.objective_value(m))
        end
    elseif termination_status(m) == MOI.INFEASIBLE
        error("Infeasible model")
    elseif termination_status(m) == MOI.INFEASIBLE_OR_UNBOUNDED
        error("Infeasible or unbounded model")
    else
        println("Failed to obtain an optimal soln to subprob, status: ", JuMP.termination_status(m))
        return m
    end
    #printDoseVolume(m,doseVol)
    return m
end

function getValidBetaInterval(m, t, tmax)
	println("In getValidBetaInterval..., dual_status(m): ", dual_status(m))
	global _V
	dbarDic = m[:dbar]
	lub = Inf64
	glb = -Inf64
	for (key,dbarVar) in dbarDic
		lb = MOI.get(m, Gurobi.VariableAttribute("SAObjLow"), dbarVar)
		if lb > glb
			glb = lb
		end
		ub = MOI.get(m, Gurobi.VariableAttribute("SAObjUp"), dbarVar)
		if ub < lub
			lub = ub
		end
	end
	return glb,lub
	#deltaDec,deltaInc
end




function addMostViolated!(m, n, x, t, tmax, β, alwaysAddDbarVars = false, idxOfLessThanCons=1)
	#adds n most violated constraints
	# m - model
	# n - num constraint per organ to add
	# x - current solution value
	# t - bound vector on organs
	global _N
	global _D
	global _V
	#    x = value.(m[:x])
	# grb = unsafe_backend(m)

	num_const_added_oar = 0
	max_viol_oar = 0.0
	for k in 1:length(_N)-1
		indicesToAdd = []
		if length(_N[k+1]) > 0 #|| tmax[k] > t[k]
			z = []
			voxIdxs = _N[k+1]
			viol = vec(_D[voxIdxs,:]*x .- t[k])
			n_min=min(n,length(viol))
			violIdxs = partialsortperm(viol,1:n_min,rev=true)
			max_viol_org = viol[violIdxs[1]]
			max_viol_oar = max(max_viol_org,max_viol_oar)
			if __DEBUG >= DEBUG_LOW
				@show viol[violIdxs[n_min]] viol[violIdxs[1]]
			end
			first_not_viol_ind = findfirst((viol[violIdxs] .<= 0.0))
			if  first_not_viol_ind == nothing
				first_not_viol_ind = n_min+1
			end
			#initialization of _V[k+1] is neccesary to prevent use of the same pointer
			if length(_V[k+1]) == 0 #&& t[k]==tmax[k] #if the first add all even if not violated
				indicesToAdd = _N[k+1][violIdxs]
				_V[k+1] = indicesToAdd
				deleteat!(_N[k+1], sort!(violIdxs))
			else
				if first_not_viol_ind > 1  #
					#    if first_not_viol_ind>1
					violIdxs = sort!(violIdxs[1:first_not_viol_ind-1])
					indicesToAdd = _N[k+1][violIdxs]
					append!(_V[k+1], indicesToAdd) #union!(_V[k+1],indicesToAdd)
					deleteat!(_N[k+1], violIdxs)
					#    end
				end
			end
			num_const_added_oar += length(indicesToAdd)
			xx = m[:x]
			dbar = []
			if t[k]<tmax[k] || β > 0 || alwaysAddDbarVars
				dbar = m[:dbar]
			end

			if (β>0 || alwaysAddDbarVars) && t[k]<tmax[k]
				for l=1:length(indicesToAdd)
					dbarIdx = indicesToAdd[l]
					dbar[dbarIdx] = @variable(m,lower_bound=0, upper_bound=tmax[k]-t[k],start = viol[violIdxs[l]])
					#grbIdx = Gurobi.column(grb,index(dbar[dbarIdx]))
					#_dbarIdxToGrbIdx[dbarIdx] = grbIdx
					#_GrbIdxToDbarIdx[grbIdx] = dbarIdx
				end
				@assert(length(_V[k+1])==length(dbar))
			end

			if t[k] < tmax[k] && (β > 0 || alwaysAddDbarVars)
				if __DEBUG >= DEBUG_MED
					println("addMostViolated! - adding constraints with dbar vars, k=", k)
				end
				@constraint(m,[i in indicesToAdd], t[k]-sum( _D[i,j]*xx[j] for j in _D[i,:].nzind) + dbar[i] >= 0)
				obj = objective_function(m, AffExpr) # NG - previously AffExpr - error
				#@show obj #, indicesToAdd
				#@objective(m, Max, obj-β*sum(dbar[i]^DEV_VAR_OBJ_NORM for i in indicesToAdd))
				@objective(m, Max, obj-β*sum(dbar[i] for i in indicesToAdd))
			elseif k==idxOfLessThanCons
				if __DEBUG >= DEBUG_MED
					println("addMostViolated! - adding <= constraint for oar: ", k)
				end
				@constraint(m,[i in indicesToAdd], sum( _D[i,j]*xx[j] for j in _D[i,:].nzind)-t[k] <= 0)
			else
				@constraint(m,[i in indicesToAdd], t[k]-sum( _D[i,j]*xx[j] for j in _D[i,:].nzind) >= 0)
			end
			#end
			if __DEBUG >= DEBUG_LOW
				println("Number of voxels in _V ", length(_V[k+1]), " voxels in _N: ",  length(_N[k+1])  , " OAR: ", k, " cons added: ", length(indicesToAdd))
			end
		end
	end
	if __DEBUG >= DEBUG_LOW
		@show length(_V[2]) length(m[:dbar])
	end
	return max_viol_oar, num_const_added_oar
end


# addNominalHomogenConstr!
# In the nominal case compare with g instead of going through all pairs
function addNominalHomogenConstr!(m,d,ϕ,μ,L=1)
    global _D
    ptvN = length(d)
	g = m[:g]
	v = ones(length(d),1)
	if !__ALLCONS_NO_GENERATION || any(d .> 0)
    	v  = ϕ.*d .- μ*value(g)
	end
    v = max.(v,0)
    idxs=partialsortperm(vec(v),1:min(L,ptvN),rev=true)
    vv = v[idxs]
    vIdxs = findall(>(VIOL_EPS),vv)
	@assert(!__ALLCONS_NO_GENERATION || length(idxs) == ptvN)
	#@assert(!__ALLCONS_NO_GENERATION || !isempty(vIdxs))
	maxViol = 0
	numConstr = 0
	x = m[:x]
	if !isempty(vIdxs)
		if L < length(vIdxs)
			vIdxs = vIdxs[1:L;]
		end
    	idxs = idxs[vIdxs]
		if __DEBUG >= DEBUG_HIGH
			@show idxs, vIdxs, size(_D), length(d)
		end
    	vv = vv[vIdxs]
    	@constraint(m, [i in idxs] ,μ*g-ϕ[i]*sum(_D[i,j]*x[j] for j in _D[i,:].nzind) >= 0)
    	numConstr = length(idxs)
		maxViol = first(vv)
	end
	if __DEBUG >= DEBUG_LOW
		println("In addNominalHomogenConstr!.. numConstr: ", numConstr)
	end
#    if !isempty(vv)
#    end
    return maxViol, numConstr #v[last(idxs)]
end

function findViolatingPairs!(m,consCollect,μ,L,phi_u_n,phi_b_n,dists,d)
    global _D
    ptvN = length(d)
    lthviol = VIOL_EPS

    #count_v=spzeros(Int64,ptvN)
    for i1 = 1:(ptvN-1)
        for k=1:2
            if k==1
                is=i1
                js=i1+1:ptvN
            else
                is=i1+1:ptvN
                js=i1
            end
			v = ones(max(length(js),length(is)),1)
			if !__ALLCONS_NO_GENERATION || any(d .> 0)
				phi_u_n_js_ones = vec(phi_u_n[js]*ones(length(is),1))
				phi_b_n_is_ones = vec(phi_b_n[is]*ones(length(js),1))
            	phibar = minimum(hcat(phi_u_n_js_ones + dists[is,js],phi_b_n_is_ones) ,dims=2)
            	phiun = maximum(hcat(phi_b_n_is_ones - dists[js,is],phi_u_n_js_ones),dims=2)
            	v1 = phibar.*d[is]-μ*phi_u_n[js].*d[js]*ones(length(is),1)
            	v2 = phi_b_n[is].*d[is]*ones(length(js),1)-μ*phiun.*d[js].*ones(length(is),1)
            	v = maximum(hcat(v1,v2),dims=2)
			end
			indexEnd = min(L,ptvN-i1,MAX_V_CONS)
			@assert(typeof(indexEnd)==Int64)
			inds = 1:indexEnd
			if !__ALLCONS_NO_GENERATION
            	inds=partialsortperm(vec(v),inds,rev=true)
			end
            #@show v[inds[1]],v[inds[end]]
            for l in inds
                if v[l] >= lthviol #+ VIOL_EPS
                    if k==1
                        insert!(consCollect,v[l],[is; js[l]])
                        #count_v[is] += 1
                        #count_v[js[l]] += 1
                    else
                        insert!(consCollect,v[l],[is[l]; js])
                        #count_v[is[l]] += 1
                        #count_v[js] += 1
                    end
                    if length(consCollect) > L
                        tmp, pair = first(consCollect)
                        #count_v[pair[1]] -= 1
                        #count_v[pair[2]] -= 1
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
    return maxviol
end


function evaluateDevNumNoDbar(x, t, tmax)
	#oarIdxs = find(x->x > 0,tmax-t)
	#x = value.(m[:x])
	violNum = zeros(length(t))
	violSum = zeros(length(t))
	for k = 1:length(t)
	  	if tmax[k] > t[k]
			dose = _D[[_V[k+1];_N[k+1]],:]*x
			violNum[k] += count(x->x > t[k] + DBARNZTH, dose)
			if violNum[k]>0
				violSum[k] += sum(dose[i]-t[k] for i in 1:length(dose) if dose[i]>t[k]+DBARNZTH)
			end
		end
	end
	return violNum, violSum
end

function evaluateDevNum(m, t, tmax)
	dbar = m[:dbar]
	@show dbar
    cntVec = zeros(length(_V)-1,1)
	sumVec = zeros(length(_V)-1,1)
	objFunc = objective_function(m)
    for k=2:length(_V)
        if t[k-1]<tmax[k-1]
            for i in _V[k]
				if __DEBUG >= DEBUG_LOW
					if i == first(_V[k])
						println("***************** β from objective function: ", coefficient(objFunc,dbar[i]))
					end
					@assert(!is_fixed(dbar[i]))
				end
				#@show k, i, t, tmax
                if (value(dbar[i]) > DBARNZTH)
                    cntVec[k-1] += 1
					sumVec[k-1] += value.dbar[i]
                end
            end
        end
    end
    return cntVec, sumVec
end


function betaBisection!(m, betaLb, betaUb, dvrhs, Din, firstIndices, t, tmax, μ, phi_u_n, phi_b_n, dists, L, dvlhsLB=[])
	βlb = betaLb
	βub = betaUb
	dvlhs = zeros(length(t))
	idx = findall(x->x>0,tmax-t)
	idx = idx[1]
	dvlhsUB = 0
	ResultsArray=[]
	while βub > BETA_EPS #any(dvrhs-dvlhs > BISECTION_GAP)

		if !isempty(dvlhsLB)
			W = (dvrhs[idx]-dvlhsUB)/(dvlhsLB-dvlhsUB)
			@assert(W>0)
			if W>=1-W_th
				W=1-W_th
			elseif W<=W_th
				W=W_th
			end
		else
			W=0.5
		end
			println("betaLb=", βlb, " betaUb=", βub, " W=", W, " violNumLb=", dvlhsLB, " violNumUB=",dvlhsUB)
		βmid = βub*(1-W)+βlb*W
		println("********* In betaBisection, current βmid=",βmid)
		time1 = @elapsed m, htCn, homCn = robustCuttingPlaneAlg!(Din,firstIndices,t,tmax,[],βmid,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m)
		#dbar = m[:dbar]
		x = value.(m[:x])
		time2 = @elapsed dvlhs, devSum = evaluateDevNumNoDbar(x,t,tmax) # function that returns deviation vector with dimension of OAR num
		time = time1 + time2
		push!(ResultsArray,[1,βmid,value(m[:g]),dvlhs[idx],devSum[idx],value.(m[:x]),time])
		if dvlhs[idx] - BISECTION_GAP < dvrhs[idx] && dvlhs[idx] >= dvrhs[idx]
			βlb = βmid
			break
		elseif devSum[idx]>dvrhs[idx]*(tmax[idx]-t[idx])
			βlb = βmid
			if !isempty(dvlhsLB)
				dvlhsLB=dvlhs[idx]
			end
		elseif dvrhs[idx] > dvlhs[idx]
			βub = βmid
			if !isempty(dvlhsLB)
				dvlhsUB=dvlhs[idx]
			end
		else # dvlhs[idx] > dvrhs[idx]
			βlb = βmid
			if !isempty(dvlhsLB)
				dvlhsLB=dvlhs[idx]
			end
		end
	end
	return  βlb, βub, ResultsArray
end

# getBasisXandDelta
function getBasisXandDelta(model)
	println(termination_status(model))
	println(MOI.get.(model, MOI.PrimalStatus()))
	x = model[:x]
	IsBasic=(value.(x).>BASIS_EPS)
	basicXIdxs=[]
	basicDeltaIdxs=[]
	if sum(IsBasic)>0
		try
			baseX = MOI.get.(model, MOI.VariableBasisStatus(), x)
			basicXIdxs = findall(x->x == MOI.BASIC,baseX)
		catch
			println("catch 1")
			basicXIdxs = findall(IsBasic)
		end
		#baseX = MOI.get(model, MOI.ConstraintBasisStatus(), LowerBoundRef(x))
		#basicXIdxs = findall(x->x == MOI.NONBASIC,baseX) # nonbasic LB constraint means basic variable
    	#	dbar = m[:dbar] = Dict()
  		basicDeltaIdxs = SortedSet{Int64}()
		dbarVarDict = model[:dbar]
		try
			for (key,val) in dbarVarDict
				baseDelta = MOI.get(model, MOI.VariableBasisStatus(), val)
				if baseDelta == MOI.BASIC
					insert!(basicDeltaIdxs,key)
				end
			end
		catch
			println("catch 2")
			for (key,val) in dbarVarDict
				if value(val)>BASIS_EPS
					insert!(basicDeltaIdxs,key)
				end
			end
		end
	end
	return basicXIdxs, basicDeltaIdxs
end

# basesAdjacent - assume idx vectors are sorted
function basesAdjacent(Sol1,Sol2,m)
	bNeighbors=true
	diffX1=setdiff(Sol1.basicXIdxs,Sol2.basicXIdxs)
	diffX2=setdiff(Sol2.basicXIdxs,Sol1.basicXIdxs)
	diffDelta1=setdiff(Sol1.basicDeltaIdxs,Sol2.basicDeltaIdxs)
	diffDelta2=setdiff(Sol2.basicDeltaIdxs,Sol1.basicDeltaIdxs)
	total_basis_diff1=(length(diffX1)+length(diffDelta1))
	total_basis_diff2=(length(diffDelta2)+length(diffX2))
	if (total_basis_diff1>1 || total_basis_diff2>1)
		bNeighbors=false
		return bNeighbors
	else #compute slacks
		#slacks for dbar upper bound
		SlackDbarIndexes1=[]
		SlackDbarIndexes2=[]
		for k=2:length(_V)
			for i in intersect(Sol1.basicDeltaIdxs,_V[k])
				if tbar[k-1]-t[k-1]-Sol1.DbarDict[i]>1e-10
					append!(SlackDbarIndexes1,i)
				end
			end
			for i in intersect(Sol2.basicDeltaIdxs,_V[k])
				if tbar[k-1]-t[k-1]-Sol2.DbarDict[i]>1e-10
					append!(SlackDbarIndexes2,i)
				end
			end
		end
		diffSlack1=setdiff(SlackDbarIndexes1,SlackDbarIndexes2)
		diffSlack2=setdiff(SlackDbarIndexes2,SlackDbarIndexes1)
		total_basis_diff1 += length(diffSlack1)
		total_basis_diff2 += length(diffSlack2)
		if (total_basis_diff1>1 || total_basis_diff2>1)
			bNeighbors=false
			return bNeighbors
		end
		#slacks for homogeneity
		vars = all_variables(m)
		cons = all_constraints(m, AffExpr, MOI.LessThan{Float64})
		xx=m[:x]
		dbar=m[:dbar]
		dind=sort(collect(keys(dbar)))
		inds1=indexin(Sol1.basicDeltaIdxs,dind)
		inds2=indexin(Sol2.basicDeltaIdxs,dind)
		gg=m[:g]
		EPS_BASIS=1e-10
		k=0
		coef_x=zeros(length(xx))
		for con in cons
			k += 1
			@show(k)
			for i=1:length(xx)
				coef_x[i]=normalized_coefficient(con,xx[i])
			end
			rhs = normalized_rhs(con)
			value1=rhs-sum(coef_x[k]*Sol1.XValue[k] for k in Sol1.basicXIdxs)#-sum(b[inds1[k]]*Sol1.DbarDict[Sol1.basicDeltaIdxs[k]][2] for k in 1:length(inds1) )
			value2=rhs-sum(coef_x[k]*Sol2.XValue[k] for k in Sol2.basicXIdxs)#-sum(b[inds2[k]]*Sol2.DbarDict[Sol2.basicDeltaIdxs[k]][2] for k in 1:length(inds2) )
			if (value1>EPS_BASIS && value2<=EPS_BASIS)
				total_basis_diff1 += 1
			elseif (value1<=EPS_BASIS && value2>EPS_BASIS)
				total_basis_diff2 += 1
			end
			if (total_basis_diff1>1.0) || (total_basis_diff2>1.0)
				bNeighbors=false
				break
			end
		end
	end
	cons = all_constraints(m, AffExpr, MOI.GreaterThan{Float64})
	k=0
	for con in cons
		coef_g=normalized_coefficient(con,gg)
		if coef_g!=0
			k += 1
			@show(k)
			for i=1:length(xx)
				coef_x[i]=normalized_coefficient(con,xx[i])
			end
			rhs = normalized_rhs(con)
			value1=rhs-coef_g*Sol1.gvalue-sum(coef_x[k]*Sol1.XValue[k] for k in Sol1.basicXIdxs)#-sum(b[inds1[k]]*Sol1.DbarDict[Sol1.basicDeltaIdxs[k]][2] for k in 1:length(inds1) )
			value2=rhs-coef_g*Sol2.gvalue-sum(coef_x[k]*Sol2.XValue[k] for k in Sol1.basicXIdxs)#-sum(b[inds1[k]]*Sol1.DbarDict[Sol1.basicDeltaIdxs[k]][2] for k in 1:length(inds1) )
			if (-value1>EPS_BASIS && -value2<=EPS_BASIS)
				total_basis_diff1 += 1
			elseif (-value1<=EPS_BASIS && -value2>EPS_BASIS)
				total_basis_diff2 += 1
			end
			if total_basis_diff1>1.0 || total_basis_diff2>1.0
				bNeighbors=false
				return bNeighbors
			end
		end
	end
	return bNeighbors
end

function getMaxLessThanConsDual(m,t=nothing)
	lamUb = 0
	cons = all_constraints(m, AffExpr, MOI.LessThan{Float64})
	for con in cons
		lam = dual(con)
		#if t!=nothing && value(con) > t + INFEAS_TOL
		#	violNum +=1
		#end
		if (lam > INFEAS_TOL)
			error("positive dual var value for <= cons, lam = ", lam)
		elseif -lam > lamUb
			lamUb = -lam
		end
	end
	return lamUb #,violNum
end

function ComputeScaling(Din, firstIndices, t,tmax)
	#compute number of voxels in organ in which t_max>t
	if firstIndices[end]<size(Din)[1]+1
		firstIndicesnew = [firstIndices;size(Din)[1]+1]
	else
		firstIndicesnew = firstIndices
	end
	scaling = 0
	K_indeces=2:length(firstIndicesnew)
	K_indeces=K_indeces[tmax.>t]
	@show K_indeces
	if length(K_indeces)>0
		scaling += sum(firstIndicesnew[k]-firstIndicesnew[k-1] for k=K_indeces)
	end
	@show scaling
	return scaling
end

function AddBudgetConstraint!(m, Din, firstIndices, dvrhs, t ,tmax, μ, phi_u_n, phi_b_n, dists, budget_limit, L)
	global _V
	solveModel!(m,firstIndices)
	obj_cur = objective_value(m)
	obj_prev = obj_cur+1
	dbarVarDict = m[:dbar]
	@show length(dbarVarDict) length(_V[2])
	K_indices=1:length(t)
	K_indices=K_indices[tmax.>t]
	temp_V=deepcopy(_V)
	@constraint(m,Budget[k in K_indices],sum(dbarVarDict[i] for i in _V[k+1]) <= budget_limit[k])
	solveModel!(m,firstIndices)
	while obj_cur-obj_prev<0
		m, htCn, homCn = robustCuttingPlaneAlg!(Din, firstIndices, t, tmax, [], 0.0, μ, phi_u_n, phi_b_n, dists, [0;0], 0, L, m, true)
		dbarVarDict = m[:dbar]
		for k in K_indices
			for i in setdiff(_V[k+1],temp_V[k+1])
				set_normalized_coefficient(Budget[k],dbarVarDict[i],1.0)
			end
		end
		temp_V=deepcopy(_V)
		solveModel!(m,firstIndices)
		obj_prev = obj_cur
		obj_cur = objective_value(m)
	end
	return m, obj_cur
end

function addMissingDoseVolume!(m,t,tmax)
	global _V
	global _D
	dbar=m[:dbar]
	x=m[:x]
	xx=value.(x)
	for k=1:length(t)
		if tmax[k]>t[k]
			dbarkeys=keys(dbar)
			newkeys=setdiff(_V[k+1],dbarkeys)
			for i in newkeys
				viol = max(sum( _D[i,j]*xx[j] for j in _D[i,:].nzind)-t[k],0)
				dbar[i] = @variable(m,lower_bound=0, upper_bound=tmax[k]-t[k],start = viol)
			end
			dbar=m[:dbar]
			@show length(dbar) length(_V[k+1])
			@constraint(m,[i in newkeys], t[k]-sum( _D[i,j]*x[j] for j in _D[i,:].nzind) + dbar[i] >= 0)
		end
	end
end

function parametricSolveIncreasing(Din,firstIndices,t,tmax,dvrhs,μ, phi_u_n, phi_b_n, dists, L=1)
    global _V
    global _N

	ResultsArray=[]
	println("##################")
	println("Compute upper bound for beta")
	println("##################")

	scaling = Float64(ComputeScaling(Din, firstIndices, t,tmax))
	time1=@elapsed m = initModel(Din, firstIndices, t, t, [], 0, phi_u_n, [0;0], 0, [], false)
	g=m[:g]
	set_objective_coefficient(m, g, scaling)
	time2=@elapsed m, htCn, homCn = robustCuttingPlaneAlg!(Din,firstIndices,t,t,[],0,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m) #,true)
	time=time1+time2
	UBg = value(m[:g])
	UBx = value.(m[:x])
	betaUb = getMaxLessThanConsDual(m)
	push!(ResultsArray,[0,betaUb,UBg,0,0,UBx,time])
	m = nothing
	clearGlobalRunData()

	idx = findall(x->x>0,tmax-t)
	idx = idx[1] #assumes only one organ has dose volume constraint
	budget_limit = zeros(length(t))

	println("##################")
	println("Compute lower bound for beta")
	println("##################")
	budget_limit[idx] = (tmax[idx]-t[idx])*scaling
	time1=@elapsed m = initModel(Din, firstIndices, tmax, tmax, [], 0, phi_u_n, [0;0], 0, [], false)
	g=m[:g]
	set_objective_coefficient(m, g, scaling)
	time2=@elapsed m, htCn, homCn = robustCuttingPlaneAlg!(Din, firstIndices, tmax, tmax, [], 0, μ, phi_u_n, phi_b_n, dists, [0;0], 0, L, m, false)
	time=time1+time2
	LBg = value(g)
	LBx = value.(m[:x])
	devNum, devSum = evaluateDevNumNoDbar(LBx, t, tmax)
	violNumLb = devNum[idx]
	obj_lb=objective_value(m)
	betaLb = 0
	addMissingDoseVolume!(m,t,tmax)
	#computer tighter lower bound using budget constraint
	if tmax[idx]<Inf && false #skipping this cause it takes too long
		# dual multipliers do not give a LB in this case
		budget_limit[idx]=dvrhs[idx]*(tmax[idx]-t[idx])
		if devSum[idx]>budget_limit[idx]
			m, obj_lb_new = AddBudgetConstraint!(m, Din, firstIndices, dvrhs, t ,tmax, μ, phi_u_n, phi_b_n, dists, budget_limit, L)
			LBGNew = value(m[:g])
			LBx = value.(m[:x])
			devNum, devSum = evaluateDevNumNoDbar(LBx, t, tmax)
			violNumLb = devNum[idx]
			@assert(budget_limit[idx]-devSum[idx]<1e-5)
			betaLb = (LBg-LBGNew)*scaling/(devSum[idx]-budget_limit[idx])
			LBg=LBGNew
			betaLb = max(-dual(m[:Budget][idx]),betaLb)
			delete(m,m[:Budget][idx])
		end
	end
	println("betaLb=", betaLb, " betaUb=", betaUb," violNumLb=", violNumLb)
	push!(ResultsArray,[0,betaLb,LBg,devNum[idx],devSum[idx],LBx,time])

	##
	#set_optimizer_attribute(m,"Method",4)
	βLb, βUb, ResultsArrayTemp = betaBisection!(m, betaLb, betaUb, dvrhs, Din, firstIndices, t, tmax, μ, phi_u_n, phi_b_n, dists, L, violNumLb)
	append!(ResultsArray,ResultsArrayTemp)
	println("************** after betaBisection, β = ", βLb)
	dbar=m[:dbar]
	for i in keys(dbar)
		set_objective_coefficient(m, dbar[i], -βLb)
	end
	set_optimizer_attribute(m,"Method",5)
	set_optimizer_attribute(m,"OutputFlag",1)
	time1=@elapsed m, htCn, homCn = @time maxmin_twostage_subprob.robustCuttingPlaneAlg!(Din,firstIndices,t,tmax,[],βLb,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m)
	time2=@elapsed solveModel!(m,firstIndices)
	time=time1+time2
	x= value.(m[:x])
	SolNew=GetSol(m)
	devVecNew, devSumNew = maxmin_twostage_subprob.evaluateDevNumNoDbar(x,t,tmax) # function that returns deviation vector with dimension of OAR num
	dbar = m[:dbar]
	new_num_dbar=length(dbar)
	NumConstOld = num_constraints(m,AffExpr, MOI.LessThan{Float64})+num_constraints(m,AffExpr, MOI.GreaterThan{Float64})
	devVecOld = devVecNew
	old_num_dbar = new_num_dbar
	βnew = βLb
	βprev= βLb
    iter = 0
	while devVecNew[idx] >= dvrhs[idx] || βnew==βprev
		iter+=1
        BETA_EPS = 0
		dbar = m[:dbar]
		for i in keys(dbar)
			set_objective_coefficient(m, dbar[i], -βnew)
		end
		time1=@elapsed m, htCn, homCn =  maxmin_twostage_subprob.robustCuttingPlaneAlg!(Din,firstIndices,t,tmax,[],βnew,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m)
		dbar = m[:dbar]
		new_num_dbar=length(dbar)
		x=value.(m[:x])
		devVecNew, devSumNew = maxmin_twostage_subprob.evaluateDevNumNoDbar(x,t,tmax) # function that returns deviation vector with dimension of OAR num
		NumConstNew=num_constraints(m,AffExpr, MOI.LessThan{Float64})+num_constraints(m,AffExpr, MOI.GreaterThan{Float64})
		println(termination_status(model))
		time2=@elapsed SolNew=GetSol(m)
		time=time+time1+time2
		if (new_num_dbar>old_num_dbar || NumConstNew>NumConstOld)
			time1=@elapsed is_adj=basesAdjacent(Sol1,Sol2,m)
			time=time+time1
			if (βnew>βprev && maximum(devVecOld-devVecNew)>1) || !is_adj
				βnew=βprev
				continue
			end
		end
		if βnew!=βprev
			push!(ResultsArray,[2,βnew,value(m[:g]),devVecNew[idx], devSumNew[idx],value.(m[:x]),time])
		end

		time=@elapsed β1,β2 = maxmin_twostage_subprob.getValidBetaInterval(m,t,tmax)
		βLbB = -β2
		βUbB = -β1
		βprev = βnew
		if devVecNew[idx] <= dvrhs[idx]
			βnew=max(βLbB-BETA_EPS,βLb)
			βUB=βprev
		else
			βnew=min(βUbB+BETA_EPS,βUb)
		end
		devVecOld=devVecNew
		SolOld=SolNew;
		NumConstOld=NumConstNew
		old_num_dbar=new_num_dbar
	end
    return m,βnew,ResultsArray
end

function CVAR_solve(Din,firstIndices,t,tmax,dvrhs,μ, phi_u_n, phi_b_n, dists, L=1,method="all")
	#=based on Romeijn, H. Edwin, Ravindra K. Ahuja, James F. Dempsey, Arvind Kumar, and Jonathan G. Li.
	 "A novel linear programming approach to fluence map optimization for intensity modulated radiation therapy
	 treatment planning." Physics in Medicine & Biology 48, no. 21 (2003): 3521.=#
	m = initModel(Din, firstIndices, tmax, tmax, [], 0, phi_u_n, [0;0], 0, [], false)
	#add CVAR constraints
	MOI.set(m, MOI.NumberOfThreads(), 10)
	htCn=0
	homCn=0
	xVar = m[:x]
	q=length(t)
	@variable(m,dvvar[1:q])
	@constraint(m,con[i=1:q],dvvar[i]<=t[i])

	isinit=fill(true,q)
	w = m[:w] = Dict()
	𐤃 = m[:𐤃] = Dict()
	indices = fill(Int[],length(t))
	if method=="iter"
		threshold=t
		m, htCn, homCn = @time robustCuttingPlaneAlg!(Din,firstIndices,tmax,tmax,[],0,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m) #,true)
	elseif method=="all"
		threshold=-1*t #just a negative number, adds all CVAR variables in the first iteration
		solveModel!(m,firstIndices)
	end
	iter = 0
	xx_new = value.(xVar)
	xx_old = xx_new.-1
	while norm(xx_new-xx_old)>1e-10
		d = Din*xx_new
		xx_old = xx_new
		iter = iter + 1
		println("start iteration ",iter)
		for i=1:q
			if t[i]<tmax[i]
				OARsize=firstIndices[i+1]-firstIndices[i]
				if isinit[i]
					indices[i] = firstIndices[i]:firstIndices[i+1]-1
				else
					indices[i] = setdiff(indices[i],keys(w))
				end
				println("adding varibales and constraints")
				added_indices=findall(d[indices[i]].>threshold[i])
				for j=added_indices
					𐤃[j] = @expression(m, sum(Din[j,k]*xVar[k] for k in Din[j,:].nzind))
					w[j] = @variable(m,lower_bound=0)
					@constraint(m,w[j]>=𐤃[j]-dvvar[i])
					set_normalized_coefficient(con[i],w[j],1/dvrhs[i])
				end
				println("added ",length(added_indices)," constraints")
			end
		end
		#solveModel!(m,firstIndices)
		#xx=value.(m[:x])
		#dvvarvalue=value.(m[:dvvar])
		println("Start resolve")
		m, htCn, homCn = @time robustCuttingPlaneAlg!(Din,firstIndices,tmax,tmax,[],0,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m) #,true)
		xx_new = value.(xVar)
		println("End resolve")
	end
	return m, htCn, homCn
end

function unfixDeltaVariables!(m,t,tmax)
	dbarVarDict = m[:dbar]
	for k=2:length(_V)
		if tmax[k-1]>t[k-1]
			for i in _V[k]
				if is_fixed(dbarVarDict[i])
					unfix( dbarVarDict[i])
					set_lower_bound( dbarVarDict[i],0)
					set_upper_bound( dbarVarDict[i],tmax[k-1]-t[k-1])
				else
					error("delta var not fixed in unfixDeltaVariables!")
				end
			end
		end
	end
end

function robustCuttingPlaneAlg!(Din,firstIndices,t,tmax,dvrhs,β,μ, phi_u_n, phi_b_n, dists, λ, ϕ, numHomogenCuts=1, m=nothing, alwaysCreateDeltaVars=false, numHomogenCutsPerVox = MAX_V_CONS, numOARConsInit = LAZY_CONS_PERORGAN_INIT_NUM)


    MAX_V_CONS = numHomogenCutsPerVox
    global LAZY_CONS_PERORGAN_INIT_NUM = numOARConsInit

    δ=maximum(phi_b_n-phi_u_n)./2

    println("In robustCuttingPlaneAlg!, δ=", δ, " numHomogenCuts=", numHomogenCuts, " numHomogenCutsPerVox=", numHomogenCutsPerVox, " numOARConsInit=", numOARConsInit)

    @assert(isempty(dvrhs) || sum(tmax-t)>0)
    ptvN = firstIndices[1]-1
    global model = m
    if m == nothing
        global model = initModel(Din,firstIndices,t,tmax,dvrhs,β,phi_u_n,λ,ϕ)
    else
        dbarVarDict = model[:dbar]
		#@show dbarVarDict
        #for (key,val) in dbarVarDict
		for k=2:length(_V)
			#@show k
			if tmax[k-1]>t[k-1]
				for i in _V[k]
	        		set_objective_coefficient(model,dbarVarDict[i], -β)
        		end
    		end
		end
	end
	iter=0
    stage = 1
    sum_num_const_added_oar = 0 #LAZY_CONS_PERORGAN_INIT_NUM  # this single counter is initialized to the initial number per organ??
    sum_num_const_added_hom = 0

    while(true)
        iter=iter+1
        prevObj = BIG_OBJ
        newObj = BIG_OBJ
        num_const_added_oar = 0 #BIG_NUM
        prev_viol_oar = BIG_NUM
        max_viol_oar = BIG_NUM
        #@time
        xVar = model[:x]
        x = []

		if (sum_num_const_added_hom > 0 && newObj == BIG_OBJ) ||  !__ALLCONS_NO_GENERATION
			for it = 1:MAX_LAZY_CON_IT
				println("outer iteration=",iter," inner iteration=",it, " max_viol_oar=",max_viol_oar)
				global model
				solveModel!(model,firstIndices)
				newObj = JuMP.objective_value(model)
				x = value.(xVar)
				#set_optimizer_attribute(model,"CPX_PARAM_LPMETHOD",1)
				num_const_added_oar = 0
				prev_viol_oar = max_viol_oar
				println(solution_summary(model))
				println(termination_status(model))
				simplex_iter=simplex_iterations(model)
				if !__ALLCONS_NO_GENERATION
					max_viol_oar, num_const_added_oar = addMostViolated!(model, LAZY_CONS_PERORGAN_NUM_PER_ITER, x, t, tmax,β,alwaysCreateDeltaVars)
				end
				if ( (prevObj-newObj)/prevObj < LAZYCONS_TIGHT_TH && max_viol_oar<=MAX_VIOL_EPS_INIT && stage == 1 ) || num_const_added_oar == 0 #|| __ALLCONS_NO_GENERATION  #abs(prev_viol_oar-max_viol_oar)/prev_viol_oar < MAX_VIOL_RATIO_TH && max_viol_oar<=MAX_VIOL_EPS_INIT && stage==1)   #&& (isempty(dvrhs) || stage ==1 || num_const_added_wdv == 0))
					println("Terminating lazy cons loop at It= ", it, "  Infeas reduction: ", (prev_viol_oar-max_viol_oar)/prev_viol_oar," Obj red %: ", (prevObj-newObj)/prevObj, " phase: ", stage, #" last solve simplex iter: ", simplex_iterations(model),
					" Simplex iter: ", simplex_iter)
					flush(stdout)
					break
				end
				if __DEBUG >= DEBUG_LOW
					println("max_viol_oar= ", max_viol_oar, "  num_const_added_oar =", num_const_added_oar)
				end
				#if num_const_added_oar == 0#max_viol_oar <= MAX_VIOL_EPS
				#    println("Terminating lazy cons loop at It= ", it)
				#    break
				#end
				sum_num_const_added_oar += num_const_added_oar
				if __DEBUG >= DEBUG_LOW
					@show max_viol_oar, num_const_added_oar
				end
				prevObj = newObj
			end
		end
#        println("Total of ", sum_num_const_added_oar, " lazy constraints added so far")
#        println("Reduced the objective by ", (prevObj-newObj)/prevObj)
#        println("Constraints violation by ", (prev_viol_oar-max_viol_oar)/prev_viol_oar)
        consCollect = SortedMultiDict{Float64,Array{Int64,1}}()
		d = zeros(ptvN)
		if !__ALLCONS_NO_GENERATION #newObj != BIG_OBJ
        	d = _D[1:ptvN,:]*x
			if __DEBUG >= DEBUG_LOW
	            @show minimum(d), maximum(d)
	        end
		end
        #@time
        #min_hom_viol=0
        constr_num = 0
        max_viol_hom=0
        if δ==0
			if sum_num_const_added_hom == 0 || !__ALLCONS_NO_GENERATION
				solveModel!(model,firstIndices)
            	max_viol_hom, constr_num = addNominalHomogenConstr!(model,d,phi_b_n,μ,numHomogenCuts)
			end
        else
			if sum_num_const_added_hom == 0 || !__ALLCONS_NO_GENERATION
            	max_viol_hom = findViolatingPairs!(model,consCollect,μ,numHomogenCuts,phi_u_n,phi_b_n,dists,d)
				constr_num = length(consCollect)
			end
        end
		sum_num_const_added_hom += constr_num
        if stage==2 && constr_num == 0#max_viol_hom<=MAX_VIOL_EPS && max_viol_oar<=MAX_VIOL_EPS #&& (isempty(dvrhs) || num_const_added_wdv == 0 )# no violated inequalities found
            println("Terminating cut algorithm.. outer iter=", iter, " Objective value: ", newObj)
            break
        elseif constr_num > 0 #|| max_viol_hom > MAX_VIOL_EPS
            #stage=1
            if __DEBUG >= DEBUG_LOW
                println("Max Violation: ", max_viol_hom, "  Adding ", length(consCollect), " cuts..")
            end
        else
            if stage == 1 && __DEBUG >= DEBUG_LOW
                println("Switching to stage 2 ******")
            end
            stage = 2
            #set_optimizer_attribute(m, "OptimalityTol", 1e-3)  # tighten tolerance in 2nd phase of the algorithm
        end
		if __ALLCONS_NO_GENERATION
			stage = 2
		end

        cons_array = Array{AbstractConstraint}(undef,length(consCollect)*2) #Any[]
		#affFunc =  Array{MOI.ScalarAffineFunction}(undef,length(consCollect)*2)
		#rhs = Array{MOI.MathOptInterface.LessThan{Float64}}(0,length(consCollect)*2)
		vIdx = 0
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
            if v1 > VIOL_EPS || __ALLCONS_NO_GENERATION
                #push!(cons_array,@build_constraint(phibar*sum(Din[i1,j]*xVar[j] for j in Din[i1,:].nzind) - μ*phi_u_n[i2]*sum(Din[i2,j]*xVar[j] for j in Din[i2,:].nzind) <= 0))
				vIdx += 1
				cons_array[vIdx] = @build_constraint(phibar*sum(Din[i1,j]*xVar[j] for j in Din[i1,:].nzind) - μ*phi_u_n[i2]*sum(Din[i2,j]*xVar[j] for j in Din[i2,:].nzind) <= 0)
				#affFunc[vIdx] = @expression(phibar*sum(Din[i1,j]*xVar[j] for j in Din[i1,:].nzind) - μ*phi_u_n[i2]*sum(Din[i2,j]*xVar[j] for j in Din[i2,:].nzind))
            end
            if v2 > VIOL_EPS || __ALLCONS_NO_GENERATION
				vIdx +=1
				cons_array[vIdx] = @build_constraint(phi_b_n[i1]*sum(Din[i1,j]*xVar[j] for j in Din[i1,:].nzind) - μ*phiun*sum(Din[i2,j]*xVar[j] for j in Din[i2,:].nzind) <= 0)
                #push!(cons_array,@build_constraint(phi_b_n[i1]*sum(Din[i1,j]*xVar[j] for j in Din[i1,:].nzind) - μ*phiun*sum(Din[i2,j]*xVar[j] for j in Din[i2,:].nzind) <= 0))
				#affFunc[vIdx] = @expression(phi_b_n[i1]*sum(Din[i1,j]*xVar[j] for j in Din[i1,:].nzind) - μ*phiun*sum(Din[i2,j]*xVar[j] for j in Din[i2,:].nzind))
            end
            if max(v1,v2) > max_viol_hom
                max_viol_hom = max(v1,v2)
            end
            #end
            #println("added constraint ", cons)
        end
        lenCons = vIdx #length(cons_array)
		cons_array = cons_array[1:lenCons]
		consCollect = []

		if __DEBUG >= DEBUG_LOW
			println("Constructed constraints, now adding to model ", lenCons, " constraints")
		end
		#@constraint(model, [i=1:lenCons],cons_array[i])
        for i=1:lenCons
           add_constraint(model, cons_array[i])
	    end

        sum_num_const_added_hom += lenCons
        println("Outer loop Iter=", iter, " Objective value: ", newObj, " OAR lazy cons: ", sum_num_const_added_oar, " Hom lazy cons: ",sum_num_const_added_hom, " max_viol_oar= ",max_viol_oar, " max_viol_hom= ", max_viol_hom)
    end # main loop
    #printDoseVolume(m, t, tmax, !isempty(dvrhs), true) # print out verbose output
    return model, sum_num_const_added_oar, sum_num_const_added_hom
end

function printDoseVolume(m, t = [], tmax = [], doseVol = false, verbose = false)
    conicForm = false
    #zVar = m[:z] #variable_by_name(m, "z")
    gVar = m[:g] #variable_by_name(m,"g")
    xVar = m[:x]
	println("Num of constraints, less than: ", JuMP.num_constraints(m, AffExpr,MOI.LessThan{Float64}), " Num of constraints, greater than: ", JuMP.num_constraints(m, AffExpr,MOI.GreaterThan{Float64}))
    println("Min PTV bio dose, g=",value(gVar))
    x = value.(xVar)

	println("Physical dose summary (including for TV):")
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
            println("Structure: ", k, " Max dose: ", maximum(dose), " Mean dose: ", mean(dose),  " Min dose: ", minimum(dose), " 99%-tile: ", quantile!(dose,0.99), " 90%-tile: ", quantile!(dose,0.90), " 70%-tile: ", quantile!(dose,0.70)  , " vox num exceeding t: " , numDev, " sum dbar/(tmax-t) :", sumZVar)
        end
    end

    if !doseVol
        dbarVar = m[:dbar] #variable_by_name(m,"dbar")
        zVar = variable_by_name(m,"z") #m[:z]
        zVarVal = []
        if zVar != nothing
            zVarVal = value.(zVar)
            #I = axes(zVarVal,Axis{1})
            nzZVal = zVarVal.data[zVarVal.data .> ZNZTH]
            if !isempty(nzZVal)
                histogram(nzZVal)
            end
            savefig("zhist.png")
        else
            println("z var not created")
        end

        for k=2:length(_V)
            if t[k-1] < tmax[k-1]
                indices = _V[k]
                if __DEBUG >= DEBUG_LOW
                    @show length(dbarVar)
                end
                dev = max.(_D[indices,:]*x.-t[k-1],0)
                nzNum = 0
                nzDvar = 0
                nzZvar = 0
                for i = 1:length(indices)
                    if dev[i] > DBARNZTH
                        nzNum+=1
                        if value.(dbarVar[indices[i]]) > DBARNZTH
                            nzDvar+=1
                        end
                    end
                    if zVar!=nothing && zVarVal[indices[i]] > ZNZTH
                        nzZvar+=1
                    end
                end
                println("Struct k=", k, " Number of deviations from t: ", nzNum)
                println("Struct k=", k, " Number of nonzero dbar: ", nzDvar)
                if zVar != nothing
                    println("Struct k=", k, " Number of nonzero z: ", nzZvar, " sum z: ", sum(zVarVal.data))
                end
            end
        end
    end
    return
end

end
