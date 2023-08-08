include("maxmin_twostage_subprob.jl")

using MAT
using DelimitedFiles
using JuMP
using SparseArrays
using FileIO, JLD2
using Printf

const bLOAD_FROM_FILE = true
const bLOAD_FROM_FILE_gamma = true
const bLOAD_FROM_FILE_projection = true
const bSAVE_FILES = false
const bSAVE_DISTORPROJ_FILES = false

#file = matopen("liverEx_2.mat")
#const ρ = [0.7; 1; 1] #0.5]# 1]
#const ρ = [1; 0.5; 1] #0.5]# 1]
# In the Liver data the 1st OAR is liver, 2nd OAR is heart
#const t = [60.0 ; 40.0; 60]
#const tmax = [60.0 ; 50.0; 60]


# brain

#λ=0 #unused reg param
β = 0
μ = 1.125 #1.45 #1.45 #1.25 #1.1  # 1.1 is for physical dose, should be higher for biological
δ = 0 #0.1  #0.1 #0.01:0.01:0.1
gamma_const=1

#n,nn = size(γ)    # to remove later, currently \gamma matrix is incorrect for liver
#γ[1:n,1:n].=0


include("BrainDScript.jl")
include("BrainPhiScript.jl")
global _V
global _N
L=200
scaling = Float64(maxmin_twostage_subprob.ComputeScaling(Din, firstIndices, t,tmax))
m = maxmin_twostage_subprob.initModel(Din, firstIndices, t, t, [], 0, phi_u_n, [0;0], 0, [], false)
g=m[:g]
set_objective_coefficient(m, g, scaling)
m, htCn, homCn = @time maxmin_twostage_subprob.robustCuttingPlaneAlg!(Din,firstIndices,t,t,[],0,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m) #,true)
initialG = value(m[:g])
initialX = value.(m[:x])
betaUb = maxmin_twostage_subprob.getMaxLessThanConsDual(m)
m = nothing
_V = nothing
_N = nothing
GC.gc()

idx = findall(x->x>0,tmax-t)
idx = idx[1]
budget_limit = zeros(length(t))


budget_limit[idx] = (tmax[idx]-t[idx])*scaling
m = maxmin_twostage_subprob.initModel(Din, firstIndices, tmax, tmax, [], 0, phi_u_n, [0;0], 0, [], false)
g=m[:g]
set_objective_coefficient(m, g, scaling)
m, htCn, homCn = maxmin_twostage_subprob.robustCuttingPlaneAlg!(Din, firstIndices, tmax, tmax, [], 0, μ, phi_u_n, phi_b_n, dists, [0;0], 0, L, m, false)
LBG = value(m[:g])
initialX = value.(m[:x])
devVec, devSum = maxmin_twostage_subprob.evaluateDevNumNoDbar(m, t, tmax)
violNumLb = devVec[idx]
obj_lb=objective_value(m)
betaLb = 0


budget_limit[idx]=dvrhs[idx]*(tmax[idx]-t[idx])
maxmin_twostage_subprob.addMissingDoseVolume!(m,t,tmax)
#=
#this takes too long!
budget_limit[idx]=dvrhs[idx]*(tmax[idx]-t[idx])
maxmin_twostage_subprob.addMissingDoseVolume!(m,t,tmax)
m, obj_lb_new = maxmin_twostage_subprob.AddBudgetConstraint!(m, Din, firstIndices, dvrhs, t ,tmax, μ, phi_u_n, phi_b_n, dists, budget_limit, L)
IndexStr=maxmin_twostage_subprob.getGlobals()
LBGNew = value(m[:g])
devVecNew, devSumNew = maxmin_twostage_subprob.evaluateDevNumNoDbar(m, t, tmax)
violNumLb = devVecNew[idx]
@assert(budget_limit[idx]-devSumNew[idx]<1e-5)
betaLb = max(-dual(m[:Budget][idx]),betaLb)=#


W = dvrhs[idx]/violNumLb
println("betaLb=", betaLb, " betaUb=", betaUb," violNumLb=", violNumLb)
dbar=m[:dbar]
for i in keys(dbar)
	set_objective_coefficient(m, dbar[i], -betaLb)
end
optimize!(m)
#delete(m,m[:Budget][idx])
#set_optimizer_attribute(m,"Method",4)
#optimize!(m)

global βLB
global βUB
βLB, βUB, ResultsArray = maxmin_twostage_subprob.betaBisection!(m, betaLb, betaUb, dvrhs, Din, firstIndices, t, tmax, μ, phi_u_n, phi_b_n, dists, L, violNumLb)

#βLB = 72.22957647838726
#βUB = 74.53033859716632
#β = 72.57469079620412
set_optimizer_attribute(m, "NumericFocus", 3)
basicXIdxsPrev = []
basicDeltaIdxsPrev = []
for i in keys(dbar)
	set_objective_coefficient(m, dbar[i], -βLB)
end
m, htCn, homCn = @time maxmin_twostage_subprob.robustCuttingPlaneAlg!(Din,firstIndices,t,tmax,[],βLB,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m)
devVecNew, devSumNew = maxmin_twostage_subprob.evaluateDevNumNoDbar(m,t,tmax) # function that returns deviation vector with dimension of OAR num
devVecOld = devVecNew
global mm
mm=m;
global βnew = βLB
global βprev= βLB
while devVecNew[idx] >= dvrhs[idx] || βnew==βprev
	global mm
	global βnew
	global βprev
	global βLB
	global βUB
	BETA_EPS = 5e-9
	dbar = mm[:dbar]
	for i in keys(dbar)
		set_objective_coefficient(mm, dbar[i], -βnew)
	end
	mm, htCn, homCn = @time maxmin_twostage_subprob.robustCuttingPlaneAlg!(Din,firstIndices,t,tmax,[],βnew,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,mm)
	devVecNew, devSumNew = maxmin_twostage_subprob.evaluateDevNumNoDbar(mm,t,tmax) # function that returns deviation vector with dimension of OAR num
	if βnew>βprev && devVecOld-devVecNew>1
		βnew=βprev
		continue
	end
	SolNew=maxmin_twostage_subprob.GetSol(mm)

	deltaDec,deltaInc = maxmin_twostage_subprob.getValidBetaInterval(mm,t,tmax)
	βLbB = βnew + deltaDec
	βUbB = βnew + deltaInc #+ deltaInc
	βprev = βnew
	if devVecNew[idx] <= dvrhs[idx]
		βnew=max(βLbB-BETA_EPS,βLB)
		βUB=βprev
	else
		βnew=min(βUbB+BETA_EPS,βUB)
	end
	basicXIdxsPrev =basicXIdxs
	basicDeltaIdxsPrev = basicDeltaIdxs
	devVecOld=devVecNew
end

struct Solution2
       basicXIdxs::Vector{Int64}
       basicDeltaIdxs::Vector{Int64}
	   XValue::Vector{Float64}
	   DbarDict::Dict{Int64,Vector{AffExpr}}
end

function bNeighbors=CheckNeighbors(Sol1,Sol2,m)
	global _N
	global _D
	global _V
	#Sol.x
	#Sol.baseX
	#Sol.dbarVarDict
	#Sol.basicDeltaIdxs
	diffX1=setdiff(Sol1.basicXIdxs,Sol2.basicXIdxs)
	diffX2=setdiff(Sol2.basicXIdxs,Sol1.basicXIdxs)
	diffDelta1=setdiff(Sol1.basicDeltaIdxs,Sol2.basicDeltaIdxs)
	diffDelta2=setdiff(Sol2.basicDeltaIdxs,Sol1.basicDeltaIdxs)

	if (length(diffX1)+length(diffDelta1))>1) || (length(diffDelta2)+length(diffX2))>1)
		bNeighbors=false
	else #compute slacks
		#slacks for dbar upper bound
		SlackDbarIndexes1=[]
		SlackDbarIndexes2=[]
		for k=1:length(_V)
			for i in intersect(Sol1.basicDeltaIdxs,_V[k])
				if tbar[k]-t[k]-Sol1.dbarVarDict[i]>1e-10
					append!(SlackDbarIndexes1,i)
				end
			end
			for i in intersect(Sol2.basicDeltaIdxs,_V[k])
				if tbar[k]-t[k]-Sol2.dbarVarDict[i]>1e-10
					append!(SlackDbarIndexes2,i)
				end
			end
		end
		#slacks for homogeneity
		all_variables(m)
		cons = all_constraints(m, AffExpr, MOI.LessThan{Float64})
		for con in cons
		end
	end
end

function GetSol(m)
	basicXIdxs, basicDeltaIdxs = maxmin_twostage_subprob.getBasisXandDelta(mm)
	dbarVarDict=m[:dbar]
	XValue=value.(m[:x])
	dbarValues=[[dbarVarDict[i],value(dbarVarDict[i])] for i in keys(dbarVarDict)]
	dbarDict=Dict(keys(dbarVarDict) .=> dbarValues)
	Sol = Solution2(basicXIdxs,collect(basicDeltaIdxs),XValue,dbarDict)
	return Sol
end
