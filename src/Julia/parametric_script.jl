include("maxmin_twostage_subprob.jl")

using MAT
using JuMP
using SparseArrays
using FileIO, JLD2
using Printf

const bLOAD_FROM_FILE = true
const bLOAD_FROM_FILE_gamma = false
const bLOAD_FROM_FILE_projection = false
const bSAVE_FILES = false
const bSAVE_DISTORPROJ_FILES = false

#file = matopen("liverEx_2.mat")
#const ρ = [0.7; 1; 1] #0.5]# 1]
#const ρ = [1; 0.5; 1] #0.5]# 1]
# In the Liver data the 1st OAR is liver, 2nd OAR is heart
#const t = [60.0 ; 40.0; 60]
#const tmax = [60.0 ; 50.0; 60]


# brain
file = matopen("Patient4_Visit1_16beams_withdeadvoxels.mat") #changed from 13 since it did not cover the PTV
 const ρ = [0.99; 1; 1]
 const t= [60; 54; 100]
 const tmax = [62; 54; 100]

#λ=0 #unused reg param
β = 0
μ = 1.15 #1.45 #1.45 #1.25 #1.1  # 1.1 is for physical dose, should be higher for biological
δ = 0 #0.1  #0.1 #0.01:0.01:0.1
gamma_const=0.051

gamma_func(x) = 1*(x>0) #(x<=max_dist)*(x>0)*(α'*[1;x;x^2])+(x>max_dist)*max_γ #0.04
γ = read(file, "neighbors_Mat")
ϕ = read(file, "omf_Vec")

n,nn = size(γ)    # to remove later, currently \gamma matrix is incorrect for liver
γ[1:n,1:n].=0

D = []
firstIndices = []
dvrhs = []

#D_file="D_Liv_formatted.jld2"
D_file="Patient4_Visit1_D_formatted.jld2"
if bLOAD_FROM_FILE
    D = FileIO.load(D_file, "D")
    firstIndices, dvrhs = FileIO.load(D_file, "firstIndices", "dvrhs")
else
    #file = matopen("Patient4_Visit1_16beams_withdeadvoxels.mat") #changed from 13 since it did not cover the PTV
    inD = read(file, "Dij")
    V = read(file, "V")
    close(file)
    println("initial D size: ",size(inD),"  OAR num (length(V)-1): ",length(V) - 1," D nnz: ",nnz(inD))
    #println(inD.nzval)
    #sumVec = sum(inD,2)
    nonzeroIdxs = unique(inD.rowval) #finall(sumVec->sumVec>0,sumVec)
    println("nonzero rows: ",length(nonzeroIdxs)," min: ",minimum(nonzeroIdxs)," max: ",maximum(nonzeroIdxs))
    nb = size(inD, 2)
    D = spzeros(0, nb)
    firstIndices = [] # vector of indices of first voxels for each of the stuctures
    dvrhs = zeros(length(t)) #zeros(length(V) - 1)

    # skipping the last strucure that includes dead voxels?
    for k = 1:length(t)+1 #length(V)-1
        if k > 1 #size(D,1)>0
            println(size(D))
            global firstIndices = [firstIndices; size(D, 1) + 1]
        end
        oarIndices = Array{Int}(vec(V[k]))
        idxs = intersect(nonzeroIdxs, oarIndices) #findall(in(nonzeroIdxs),oarIndices)
        println("struct: ",k," before removing rows: ", length(V[k]), " after removing rows: ", length(idxs), " ", typeof(nonzeroIdxs), " ", typeof(oarIndices))
        #println(minimum(oarIndices), " ", maximum(oarIndices), " ")
        appendD = inD[idxs, :]
        rowN, colN = size(appendD)
        println("For struct: ", k, " appending submat of D of size: ", size(appendD))
        if rowN > 0
            global D = [D; appendD]
            if k > 1
                dvrhs[k-1] = floor((1 - ρ[k-1]) * length(V[k]))
                println("DVRHS: ", dvrhs[k-1] , " for organ: ", k-1)
            end
        end
    end
    if bSAVE_FILES
        FileIO.save(D_file,"D",D,"firstIndices",firstIndices,"dvrhs",dvrhs)
    end
end
close(file)

if sum(tmax-t) == 0 || sum(dvrhs) == 0
    dvrhs = []
end
inD = []
println(firstIndices)
@assert(size(γ, 1) == firstIndices[1] - 1)#need to make sure these are the same
#m = maxmin_twostage_subprob.initModel(D,firstIndices,t,dvrhs,β)
#m = maxmin_twostage_subprob.solveModel!(m)
#for
phi_under = ϕ.- δ
phi_under[phi_under.<0].= 0
phi_bar = ϕ.+ δ
phi_bar[phi_bar.>1].= 1


println("Now solving with min phi bar = ", minimum(phi_bar), " min phi_under = " , minimum(phi_under), " max phi_bar = ", maximum(phi_bar), " max phi_under = ", maximum(phi_under))
phi_u_n=[]
phi_b_n=[]
dists=[]
file_name_gamma=@sprintf("./RS_Dists/Patient4_Visit1_Gamma_dist_new_%1.3f.jld2",gamma_const)
file_name_proj=@sprintf("./Projections/PPatient4_Visit1_rojection_new_%1.3f_%1.3f.jld2",gamma_const,δ)

phi_u_n = ϕ
phi_b_n = ϕ

Din=D
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
	m, obj_lb_new = maxmin_twostage_subprob.AddBudgetConstraint!(m, Din, firstIndices, dvrhs, t ,tmax, μ, phi_u_n, phi_b_n, dists, budget_limit, L)
	LBGNew = value(m[:g])
	devVecNew, devSumNew = maxmin_twostage_subprob.evaluateDevNumNoDbar(m, t, tmax)
	violNumLb = devVecNew[idx]
	@assert(budget_limit[idx]-devSumNew[idx]<1e-5)
	betaLb = max(-dual(m[:Budget][idx]),betaLb)


	W = dvrhs[idx]/violNumLb
	println("betaLb=", betaLb, " betaUb=", betaUb," violNumLb=", violNumLb)
	dbar=m[:dbar]
	for i in keys(dbar)
		set_objective_coefficient(m, dbar[i], -betaLb)
	end
	optimize!(m)
	delete(m,m[:Budget][idx])
	set_optimizer_attribute(m,"Method",4)
	optimize!(m)

	global βLB
	global βUB
	β, βLB, βUB = maxmin_twostage_subprob.betaBisection!(m, betaLb, betaUb, dvrhs, Din, firstIndices, t, tmax, μ, phi_u_n, phi_b_n, dists, L, violNumLb)


	set_optimizer_attribute(m, "NumericFocus", 3)
	basicXIdxsPrev = []
	basicDeltaIdxsPrev = []
	for i in keys(dbar)
		set_objective_coefficient(m, dbar[i], -β)
	end
	m, htCn, homCn = @time maxmin_twostage_subprob.robustCuttingPlaneAlg!(Din,firstIndices,t,tmax,[],β,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m)
	devVecNew, devSumNew = maxmin_twostage_subprob.evaluateDevNumNoDbar(m,t,tmax) # function that returns deviation vector with dimension of OAR num
	global mm
	mm=m;
	global βnew
	βnew=β
	global βprev=β
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
		if βnew>βprev && devVecNew-devVecOld>1
			βnew=βprev
			continue
		end
		basicXIdxs, basicDeltaIdxs = maxmin_twostage_subprob.getBasisXandDelta(mm)
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
