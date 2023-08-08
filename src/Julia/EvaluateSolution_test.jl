include("evaluate_solution.jl")
include("evaluateOAR_solution.jl")
include("maxmin_twostage_subprob.jl")
using MAT
using JuMP
using XLSX
using SparseArrays
using FileIO, JLD2
using Printf
bLOAD_FROM_FILE=false
bSAVE_FILES = true
bLOAD_FROM_FILE_gamma=true
bLOAD_FROM_FILE_projection = true

#solutions imput and output path
InputXLS="C:/Users/shimrit/Dropbox (MIT)/technion/research/Robust Radio Therapy (All)/Robust Radiotherapy/Results/NewGammaModel/20220818_no_dose_vol_no_zeros.xlsx" #no_dose_vol_b_no_zeros.xlsx"#"./ResultsFiles/nominal_solutions.xlsx"#"./ResultsFiles/no_dose_vol_box.xlsx"#"./ResultsFiles/no_dose_vol_b_no_zeros.xlsx"#"no_dose_vol.xlsx" #
InputSheet="no_dose_vol"#"nominal_solutions"#"no_dose_vol_box"#"no_dose_vol"#"nominal_solutions"#
summary_file_name="./ResultsFiles/20220818_solutionOAR_evaluation.txt"#"./ResultsFiles/solution_evaluation_nominal.txt" #"./ResultsFiles/solution_evaluation_box.txt"#"solution_evaluation_nominal.txt" #

ρ = [0.99;1;1]#[0.99 ; 1; 1]
t = [60; 54; 100]
tmax = [62; 54; 100]
β = 0 #1e-8
include("BrainDScript.jl")
maxmin_twostage_subprob.initGlobals(Din,firstIndices)
println("after read MAT")

println("before read XLS")
Ɣ_start = 0.05
Ɣ_end = 0.05
δ_end = 0.1
δ_start = Dict(0 => 0.05, 0.01 => 0.03, 0.02 => 0.02, 0.03 => 0.02, 0.04 => 0.01, 0.05 => 0.01, 0.06 => 0, 0.07 => 0, 0.08 => 0, 0.09 => 0, 0.1 => 0)

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
			#=for gamma_const in Ɣ_start:0.01:Ɣ_end
				global dists
				file_name_gamma=@sprintf("./RS_Dists/Gamma_dist_new_%1.3f.jld2",gamma_const)
				if bLOAD_FROM_FILE_gamma
					global dists = FileIO.load(file_name_gamma,"dists")
				else
					global dists=[]
				end
				#if δ_origin<δ_start[gamma_const]
				#	continue
				#end
				for δ in δ_start[gamma_const]:0.01:δ_end#δ_origin#
					if δ_origin==δ
						continue
					end
					global dists
					file_name_proj=@sprintf("./Projections/Projection_new_%1.3f_%1.3f.jld2",gamma_const,δ)
					new_r += 1
					#@show new_r
					phi_under = ϕ.- δ
					phi_under[phi_under.<0].= 0
					phi_bar = ϕ.+ δ
					phi_bar[phi_bar.>1].= 1
					println("Now solving with min phi bar = ", minimum(phi_bar), " min phi_under = " , minimum(phi_under), " max phi_bar = ", maximum(phi_bar), " max phi_under = ", maximum(phi_under))
					phi_u_n=[]
					phi_b_n = []
					if bLOAD_FROM_FILE_projection*bLOAD_FROM_FILE_gamma
						phi_u_n, phi_b_n = FileIO.load(file_name_proj,"phi_u_n","phi_b_n")
					else
						#for patient 4
						α=[0.03145+gamma_const,0.00228,-7.885e-5]
						max_γ=0.05+gamma_const
						max_dist=10
						gamma_func(x) = (x<=max_dist)*(x>0)*(α'*[1;x;x^2])+(x>max_dist)*max_γ #0.04

						global dists
				        phi_u_n, phi_b_n, dists = maxmin_twostage_subprob.computeProjections(γ, gamma_func, phi_under, phi_bar)
						if !bLOAD_FROM_FILE_gamma
				            FileIO.save(file_name_gamma,"dists",dists)
				        end
				        FileIO.save(file_name_proj,"phi_u_n",phi_u_n,"phi_b_n",phi_b_n)
				    end
					time_prof=@elapsed g, μ=evaluate_solution.EvaluateSolution(d,firstIndices,dists,phi_u_n,phi_b_n)
					#sh2[new_r,:]=vec(hcat(μ_origin,t',δ,gamma_const,μ,g))
					open(summary_file_name,"a") do io
						println(io,μ_origin,",",δ_origin,",",Ɣ_origin,",",t_max',",",δ,",",gamma_const,",",μ,",",g)
					end
				end
			end=#
			#time_prof=@elapsed violProp=evaluateOAR_solution.EvaluateOARSolution(d,firstIndices,t,t_max)
			time_prof=@elapsed violNum, violSum=maxmin_twostage_subprob.evaluateDevNumNoDbar(x, t, tmax)
			println(violNum)
			println(OARSize)
			violProp=violNum./OARSize
			open(summary_file_name,"a") do io
				println(io,μ_origin,",",δ_origin,",",Ɣ_origin,",",t',",",t_max',",",violProp)
			end
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
