include("evaluate_solution.jl")
using MAT
using JuMP
using XLSX
using SparseArrays
using FileIO, JLD2
using Printf
bLOAD_FROM_FILE=true
bLOAD_FROM_FILE2=true

ρ = [1;1;1]#[0.99 ; 1; 1]
β = 0 #1e-8
file = matopen("Patient4_Visit1_16beams_withdeadvoxels.mat") #changed from 13 since it did not cover the PTV
γ = read(file, "neighbors_Mat")
ϕ = read(file, "omf_Vec")
close(file)
println("after read MAT")

D = []
firstIndices = []

if bLOAD_FROM_FILE
    D = FileIO.load("D_formatted.jld2", "D")
    firstIndices, dvrhs =
    FileIO.load("D_formatted.jld2", "firstIndices", "dvrhs")
    println("after read jld2")
else
    file = matopen("Patient4_Visit1_16beams_withdeadvoxels.mat") #changed from 13 since it did not cover the PTV
    inD = read(file, "Dij")
    V = read(file, "V")
    close(file)

    println(
        "initial D size: ",
        size(inD),
        "  OAR num (length(V)-1): ",
        length(V) - 1,
        " D nnz: ",
        nnz(inD),
    )
    #println(inD.nzval)

    #sumVec = sum(inD,2)
    nonzeroIdxs = unique(inD.rowval) #finall(sumVec->sumVec>0,sumVec)
    println(
        "nonzero rows: ",
        length(nonzeroIdxs),
        " min: ",
        minimum(nonzeroIdxs),
        " max: ",
        maximum(nonzeroIdxs),
    )
    nb = size(inD, 2)

    D = spzeros(0, nb)
    firstIndices = [] # vector of indices of first voxels for each of the stuctures
    dvrhs = zeros(length(V) - 1)

    for k = 1:length(V)-1
        if k > 1 #size(D,1)>0
            println(size(D))
            global firstIndices = [firstIndices; size(D, 1) + 1]
        end
        oarIndices = Array{Int}(vec(V[k]))
        idxs = intersect(nonzeroIdxs, oarIndices) #findall(in(nonzeroIdxs),oarIndices)
        println(
            "struct: ",
            k,
            " before removing rows: ",
            length(V[k]),
            " after removing rows: ",
            length(idxs),
            " ",
            typeof(nonzeroIdxs),
            " ",
            typeof(oarIndices),
        )
        #println(minimum(oarIndices), " ", maximum(oarIndices), " ")
        appendD = inD[idxs, :]
        rowN, colN = size(appendD)
        println(
            "For organ: ",
            k,
            " appending submat of D of size: ",
            size(appendD),
        )
        if rowN > 0
            global D = [D; appendD]
            if k > 1
                dvrhs[k-1] = floor((1 - ρ[k-1]) * length(V[k]))
                println("DVRHS: ", dvrhs[k-1] , " for organ: ", k-1, )
            end
        end
    end
    FileIO.save(
        "D_formatted.jld2",
        "D",
        D,
        "firstIndices",
        firstIndices,
        "dvrhs",
        dvrhs,
    )
end
println("before read XLS")
if sum(dvrhs) == 0
    dvrhs = []
end
inD = []
println(firstIndices)
Ɣ_start = 0.06
Ɣ_end = 0.06
δ_end = 0
δ_start = Dict(0 => 0.05, 0.01 => 0.03, 0.02 => 0.02, 0.03 => 0.02, 0.04 => 0.01, 0.05 => 0.01, 0.06 => 0, 0.07 => 0, 0.08 => 0, 0.09 => 0, 0.1 => 0)

summary_file_name="solution_evaluation.txt"
@assert(size(γ, 1) == firstIndices[1] - 1)#need to make sure these are the same
XLSX.openxlsx("no_dose_vol.xlsx",mode="rw") do xf
    sh = xf["no_dose_vol"]
    #sh2= xf["no_dose_vol"]
	new_r = 1
    for r in XLSX.eachtablerow(sh)
        number_of_organs=length(ρ)
        number_of_beamlets=size(D,2)
        t=zeros(number_of_organs,1)
        for i=1:number_of_organs
            t[i] = r[4+i]
        end
        t_max=t
        x=zeros(number_of_beamlets,1)
        for i=1:number_of_beamlets
			x[i]=r[(4+number_of_organs+4+i)]
        end
		
		μ_origin=r[1]
		δ_origin=r[2]
		Ɣ_origin=r[3]
			
		if sum(x)!=0
			d = D*x
			for gamma_const in Ɣ_start:0.01:Ɣ_end
				for δ in δ_start[gamma_const]:0.01:δ_end
					new_r += 1
					@show new_r
					phi_under = ϕ.- δ
					phi_under[phi_under.<0].= 0
					phi_bar = ϕ.+ δ
					phi_bar[phi_bar.>1].= 1
					println("Now solving with min phi bar = ", minimum(phi_bar), " min phi_under = " , minimum(phi_under), " max phi_bar = ", maximum(phi_bar), " max phi_under = ", maximum(phi_under))
					time_prof=@elapsed g, μ=evaluate_solution.EvaluateSolution(d,firstIndices,γ,gamma_const,δ,phi_under,phi_bar,bLOAD_FROM_FILE2)
					#sh2[new_r,:]=vec(hcat(μ_origin,t',δ,gamma_const,μ,g))
					open(summary_file_name,"a") do io
						println(io,μ_origin,",",δ_origin,",",Ɣ_origin,",",t',",",δ,",",gamma_const,",",μ,",",g)
					end
				end
			end
		else
			for gamma_const in Ɣ_start:0.01:Ɣ_end
				for δ in δ_start[gamma_const]:0.01:δ_end
					open(summary_file_name,"a") do io
						println(io,μ_origin,",",δ_origin,",",Ɣ_origin,",",t',",",δ,",",gamma_const,",","NaN",",",0)
					end
				end
			end
		end
    end
end


