include("maxmin_twostage_subprob_nogen2.jl")
include("evaluate_solution.jl")
#using maxmin_twostage_subprob
#using maxmin_twostage_subprob.initModel, maxmin_twostage_subprob.solveModel
#import maxmin_twostage_subprob: initModel, solveModel
using MAT
using JuMP
using SparseArrays
using FileIO, JLD2
using Printf
using LinearAlgebra
bLOAD_FROM_FILE= false 
bLOAD_FROM_FILE_gamma= false
bLOAD_FROM_FILE_projection = false 
bSAVE_FILES = false 
bSAVE_DISTORPROJ_FILES = false

MAX_HOMOGEN_CONSTR = 0
MAX_V_CONS = 0
LAZY_CONS_PERORGAN_INIT_NUM = 0

#file = matopen("liverEx2.mat")
#ρ = [0.99; 1; 1]
#ρ = [1; 1; 1]
#t = [40.0 ; 40.0]
#t = [62; 54; 100]
#t=[60; 54; 100]
ρ = [1; 1; 1]
t = [62; 54; 100]
tmax = [62; 54; 100]
#β = 0.01
β = 0 #1e-8 # coeficient do deal with dose volume constraint
λ = [0; 0] # coeficient for regularizer preventing multiple solutions
#the first parameter coincides with the sum(x_i) regularizer
#the first parameter coincides with the g_nom - minimum nominal dose regularizer
δ=0 #0.08 #0.01:0.01:0.1
μ = 1.125 #1.25 #1.1
gamma_const=1
if Sys.islinux()
    #μ = parse(Float64,ARGS[1])#1.2 #1.25 #1.1
    #gamma_const = parse(Float64,ARGS[2])
    #δ = parse(Float64,ARGS[3])#0.05
    MAX_HOMOGEN_CONSTR = parse(Int64,ARGS[1])
    MAX_V_CONS = parse(Int64,ARGS[2])
    LAZY_CONS_PERORGAN_INIT_NUM = parse(Int64,ARGS[3])
end
#for patient 4

include("BrainDScript.jl")
include("BrainPhiScript.jl")
println(dvrhs)
time_prof=@elapsed model, htCn, homCn =maxmin_twostage_subprob.robustCuttingPlaneAlg!(D,firstIndices,t,tmax,dvrhs,β,μ,phi_u_n, phi_b_n, dists,λ, ϕ, MAX_HOMOGEN_CONSTR,nothing, false, MAX_V_CONS, LAZY_CONS_PERORGAN_INIT_NUM)

@show time_prof
maxmin_twostage_subprob.printDoseVolume(model, t, tmax, !isempty(dvrhs), true) # print out verbose output

xx = value.(model[:x])
g = value.(model[:g])
mu_calc=evaluate_solution.CalculateMinHom(d,firstIndices,dists,phi_u_n,phi_b_n)
reg1 = value.(model[:reg1])
reg2 = value.(model[:reg2])
d = D[1:firstIndices[1]-1,:]*xx
minPhysDose, PhysHom = evaluate_solution.CalculatePhysPerformance(d,firstIndices)
g_nominal, μ_nominal = evaluate_solution.CalculateNomPerformance(d,firstIndices,ϕ)
EUD = evaluate_solution.EvaluateEUD(d,firstIndices)
#file_name=@sprintf("results_%1.3f_%.2f_%.3f_%.2f.jld2",β,μ,δ,gamma_const)
#FileIO.save(file_name,"δ",δ,"μ",μ,"β",β,"t",t,"gamma_const",gamma_const,"time_prof",time_prof,"xx",xx,"g",g,"PhysHom",PhysHom,"phi_under",phi_under,"phi_bar",phi_bar)
#summary_file_name="./NewResultsFiles/no_dose_vol_nom_reg.txt" # not recommended that 2 processes print to same file -- check with shimrit
summary_file_name="./ResultsFiles/no_dose_vol_20230811.txt"
open(summary_file_name,"a") do io
    println(io,μ,",",δ,",",gamma_const,",",β,",",λ[1],",",λ[2],",",t,",",g,",",mu_calc,",",reg1,",",reg2,",",PhysHom,",",minPhysDose,",",time_prof,",",xx,",",μ_nominal,",",g_nominal,",",EUD)
end
