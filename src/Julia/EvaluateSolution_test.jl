include("evaluate_solution.jl")
include("evaluateOAR_solution.jl")
include("maxmin_twostage_subprob.jl")
using MAT
using JuMP
using XLSX
using SparseArrays
using FileIO, JLD2
using Printf
bLOAD_FROM_FILE=true
bSAVE_FILES = false
bLOAD_FROM_FILE_gamma=true
bLOAD_FROM_FILE_projection = true

#solutions imput and output path
InputXLS="./ResultsFiles/20230906_no_dose_vol_no_zeros.xlsx" #no_dose_vol_b_no_zeros.xlsx"#"./ResultsFiles/nominal_solutions.xlsx"#"./ResultsFiles/no_dose_vol_box.xlsx"#"./ResultsFiles/no_dose_vol_b_no_zeros.xlsx"#"no_dose_vol.xlsx" #
InputSheet="no_dose_vol"#"nominal_solutions"#"no_dose_vol_box"#"no_dose_vol"#"nominal_solutions"#
summary_file_name="./ResultsFiles/20230919_solution_evaluation.txt"#"./ResultsFiles/solution_evaluation_nominal.txt" #"./ResultsFiles/solution_evaluation_box.txt"#"solution_evaluation_nominal.txt" #

ρ = [0.99;1;1]#[0.99 ; 1; 1]
t = [60; 54; 100]
tmax = [62; 54; 100]
β = 0 #1e-8
include("BrainDScript.jl")
maxmin_twostage_subprob.initGlobals(Din,firstIndices)
println("after read MAT")

println("before read XLS")
Ɣ_start = 1
Ɣ_end = 1
δ_end = 0.1
δ_start = Dict(0 => 0.02, 0.01 => 0.02, 0.02 => 0.01, 0.03 => 0.01, 0.04 => 0.00, 0.05 => 0.00, 1=>0.00)

XLSX.openxlsx(InputXLS,mode="r") do xf
    sh = xf[InputSheet]
    #sh2= xf["no_dose_vol"]
	new_r = 1
    for r in XLSX.eachtablerow(sh)
        number_of_organs=length(ρ)
        number_of_beamlets=size(D,2)
        t_max=zeros(number_of_organs,1)
        for i=1:number_of_organs
            t_max[i] = r[6+i]
        end
        x=zeros(number_of_beamlets,1)
        for i=1:number_of_beamlets
			x[i]=r[(6+number_of_organs+6+i)]
        end

		μ_origin=r[1]
		δ_origin=r[2]
		Ɣ_origin=r[3]
		@show μ_origin, δ_origin, Ɣ_origin
		if sum(x)!=0
			d = D*x
			global dists=[]
			for Ɣ in Ɣ_start:0.01:Ɣ_end
				global dists
				global gamma_const=Ɣ
				file_name_gamma=@sprintf("./RS_Dists/Gamma_dist_new_%1.3f.jld2",gamma_const)
				if δ_origin<δ_start[gamma_const]
					continue
				end
				for delta in δ_start[gamma_const]:0.01:δ_end#δ_origin#
					#if δ_origin==δ
					#	continue
					#end
					global dists
					global δ=delta
					include("BrainPhiScript.jl")
					new_r += 1
					#@show new_r
					time_prof=@elapsed g, μ=evaluate_solution.EvaluateSolution(d,firstIndices,dists,phi_u_n,phi_b_n)
					#sh2[new_r,:]=vec(hcat(μ_origin,t',δ,gamma_const,μ,g))
					open(summary_file_name,"a") do io
						println(io,μ_origin,",",δ_origin,",",Ɣ_origin,",",t_max',",",δ,",",gamma_const,",",μ,",",g)
					end
				end
			end#
			#time_prof=@elapsed violProp=evaluateOAR_solution.EvaluateOARSolution(d,firstIndices,t,t_max)
			#=time_prof=@elapsed violNum, violSum=maxmin_twostage_subprob.evaluateDevNumNoDbar(x, t, tmax)
			println(violNum)
			println(OARSize)
			violProp=violNum./OARSize
			open(summary_file_name,"a") do io
				println(io,μ_origin,",",δ_origin,",",Ɣ_origin,",",t',",",t_max',",",violProp)
			end=#
		else
			#=for gamma_const in Ɣ_start:0.01:Ɣ_end
				for δ in δ_start[gamma_const]:0.01:δ_end
					open(summary_file_name,"a") do io
						println(io,μ_origin,",",δ_origin,",",Ɣ_origin,",",t',",",δ,",",gamma_const,",","NaN",",",0)
					end
				end
			end=#
			violProp=zeros(length(t),1)
			open(summary_file_name,"a") do io
				println(io,μ_origin,",",δ_origin,",",Ɣ_origin,",",t',",",t_max',",",violProp)
			end
		end
    end
end
