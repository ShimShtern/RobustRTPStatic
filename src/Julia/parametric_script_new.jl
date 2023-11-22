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
m,βnew,ResultsArray=maxmin_twostage_subprob.parametricSolveIncreasing(Din,firstIndices,t,tmax,dvrhs,μ, phi_u_n, phi_b_n, dists,L)
outfile=sprintf("./ResultsFiles/ChangingBeta_delta%.2f_gamma%.2f_mu%.3f.csv",δ,gamma_const,μ)
f = open(outfile, "w")
writedlm(f, ResultsArray,",")
close(f)
