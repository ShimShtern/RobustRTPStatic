include("maxmin_twostage_subprob.jl")
include("evaluate_solution.jl")
include("evaluateOAR_solution.jl")


using MAT
using DelimitedFiles
using JuMP
using SparseArrays
using FileIO, JLD2
using Printf


const bLOAD_FROM_FILE = false
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
ρ = [0.99; 1; 1]
t = [60; 54; 100]
tmax = [62; 54; 100]

β = 0
μ = 1.125 #1.45 #1.45 #1.25 #1.1  # 1.1 is for physical dose, should be higher for biological
δ = 0 #0.1  #0.1 #0.01:0.01:0.1
gamma_const=1.0
include("BrainDScript.jl")
include("BrainPhiScript.jl")
println("Now solving with min phi bar = ", minimum(phi_bar), " min phi_under = " , minimum(phi_under), " max phi_bar = ", maximum(phi_bar), " max phi_under = ", maximum(phi_under))

L=200
time = @elapsed m,βnew,ResultsArray=maxmin_twostage_subprob.parametricSolveIncreasing(Din,firstIndices,t,tmax,dvrhs,μ, phi_u_n, phi_b_n, dists,L)
outfile=@sprintf("./ResultsFiles/ChangingBeta_delta%.2f_gamma%.2f_mu%.3f.csv",δ,gamma_const,μ)
f = open(outfile, "w")
writedlm(f, ResultsArray,",")
close(f)

xx = value.(m[:x])
g = value.(m[:g])
d = D*xx
mu_calc=evaluate_solution.CalculateMinHom(d,firstIndices,dists,phi_u_n,phi_b_n)
minPhysDose, PhysHom = evaluate_solution.CalculatePhysPerformance(d,firstIndices)
g_nominal, μ_nominal = evaluate_solution.CalculateNomPerformance(d,firstIndices,ϕ)
violProp = evaluateOAR_solution.EvaluateOARSolution(d,firstIndices, t,tmax)
EUD = evaluate_solution.EvaluateEUD(d,firstIndices)
outfile="./ResultsFiles/ChangingBeta.csv"
open(outfile,"a") do io
    println(io,μ,",",δ,",",gamma_const,",",t,",",tmax,",",1 .- ρ,",",mu_calc,",",g,",",PhysHom,",",minPhysDose,",",time,",",xx,",",μ_nominal,",",g_nominal,",", EUD,",",violProp')
end
