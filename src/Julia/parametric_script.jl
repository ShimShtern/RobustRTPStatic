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

	β = maxmin_twostage_subprob.betaBisection!(m, betaLb, betaUb, dvrhs, Din, firstIndices, t, tmax, μ, phi_u_n, phi_b_n, dists, L, violNumLb)
	set_optimizer_attribute(m, "NumericFocus", 3)
	basicXIdxsPrev = []
	basicDeltaIdxsPrev = []
	#β=72.229
	m, htCn, homCn = @time maxmin_twostage_subprob.robustCuttingPlaneAlg!(Din,firstIndices,t,tmax,[],β,μ,phi_u_n,phi_b_n,dists,[0;0],0,L,m)
	devVecNew, devSumNew = maxmin_twostage_subprob.evaluateDevNumNoDbar(m,t,tmax) # function that returns deviation vector with dimension of OAR num
	basicXIdxs, basicDeltaIdxs = maxmin_twostage_subprob.getBasisXandDelta(m)
	deltaDec,deltaInc = maxmin_twostage_subprob.getValidBetaInterval(m,t,tmax)
