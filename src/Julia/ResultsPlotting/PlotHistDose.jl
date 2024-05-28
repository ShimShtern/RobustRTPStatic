
using MAT
using JuMP
using SparseArrays
using FileIO, JLD2, XLSX
using Printf
using LinearAlgebra
using Plots
using Measures


bLOAD_FROM_FILE=false
bSAVE_FILES=false

gamma=0.04
delta=0.14

global x_nominal=[]
global mu_nominal=2
global g_nominal=0
global x_box=[]
global mu_box=2
global g_box=0
global x_spatial=[]
global mu_spatial=2
global g_spatial=0
t = [60; 54; 100]
tmax = [62; 54; 100]
ρ = [1; 1; 1]

#path relative to running from ./src/julia/ and not current location
include("../BrainDScript.jl")
PlotFolder = "./ResultsFiles/Plots/Patient4/"
InputXLS="./ResultsFiles/Patient4_Visit1_20230906_no_dose_vol_no_zeros.xlsx"
InputSheet="no_dose_vol"
XLSX.openxlsx(InputXLS,mode="r") do xf
    sh = xf[InputSheet]
    for r in XLSX.eachtablerow(sh)
        number_of_organs=length(ρ)
        number_of_beamlets=size(D,2)
        t_max=zeros(number_of_organs,1)
        for i=1:number_of_organs
            t_max[i] = r[6+i]
        end
        phys_hom=r[(10+number_of_organs)]
        phys_dose=r[(11+number_of_organs)]
        x=zeros(number_of_beamlets,1)
        for i=1:number_of_beamlets
            x[i]=r[(6+number_of_organs+6+i)]
        end

        μ_origin=r[1]
        δ_origin=r[2]
        Ɣ_origin=r[3]
        if δ_origin==0 && Ɣ_origin==1
            if phys_hom<=1.1 && phys_dose>g_nominal#mu_nominal>μ_origin
                println(δ_origin," ",Ɣ_origin," ",μ_origin," ",mu_nominal)
                global x_nominal = x
                global mu_nominal = μ_origin
                global g_nominal = phys_dose
            end
        elseif Ɣ_origin==1 && δ_origin==delta
            if phys_hom<=1.1 && phys_dose>g_box #mu_box>μ_origin
                println(δ_origin," ",Ɣ_origin," ",μ_origin," ",mu_nominal)
                global x_box = x
                global mu_box = μ_origin
                global g_box = phys_dose
            end
        elseif Ɣ_origin==gamma && δ_origin==delta
            if phys_hom<=1.1 && phys_dose>g_spatial #mu_spatial>μ_origin
                println(δ_origin," ",Ɣ_origin," ",μ_origin," ",mu_nominal)
                global x_spatial = x
                global mu_spatial = μ_origin
                global g_spatial=phys_dose
            end
        end
    end
end

IndsT=1:(firstIndices[1]-1)
dose_nominal_T=[46;sort(D[IndsT,:]*x_nominal,dims=1)]
dose_box_T=[46;sort(D[IndsT,:]*x_box,dims=1)]
dose_spatial_T=[46;sort(D[IndsT,:]*x_spatial,dims=1)]
nominal_b_range = range(60, 67, length=29)
box_b_range = range(45, 52, length=29)
freqT=reverse(IndsT)/(firstIndices[1]-1)
freqT=[1;freqT]*100

Plots.plot(dose_nominal_T,freqT,label="Nominal",ylabel="Relative Volume (%)", xlabel="Dose (Gy)",color="blue",linestyle=:solid,# fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_spatial_T,freqT,label="Spatial", color="green",linestyle=:dash,# fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_box_T,freqT,label="Box", color="red",linestyle=:dot,legend=:bottom, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=box_b_range
Plots.savefig(@sprintf("./%s/Patient4_Visit1_PTVDoseHist.png",PlotFolder))


IndsB=firstIndices[1]:(firstIndices[2]-1)
dose_nominal_B=[0;sort(D[IndsB,:]*x_nominal,dims=1);62]
dose_box_B=[0;sort(D[IndsB,:]*x_box,dims=1);62]
dose_spatial_B=[0;sort(D[IndsB,:]*x_spatial,dims=1);62]
freq=(reverse(IndsB).-(firstIndices[1]+1))/(firstIndices[2]-firstIndices[1])
freqB=[1;freq;0]*100

Plots.plot(dose_nominal_B,freqB,label="Nominal",ylabel="Relative Volume (%)", xlabel="Dose (Gy)",color="blue",linestyle=:solid, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_spatial_B,freqB,label="Spatial", color="green", linestyle=:dash,#fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_box_B,freqB,label="Box", color="red",linestyle=:dot,legend=:topright, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=box_b_range
Plots.savefig(@sprintf("./%s/Patient4_Visit1_BrainDoseHist.png",PlotFolder))


#IndsE=[firstIndices[1]:(firstIndices[2]-1);firstIndices[3]:size(D,1)]
IndsE=firstIndices[3]:size(D,1)
dose_nominal_E=[0;sort(D[IndsE,:]*x_nominal,dims=1)]
dose_box_E=[0;sort(D[IndsE,:]*x_box,dims=1)]
dose_spatial_E=[0;sort(D[IndsE,:]*x_spatial,dims=1)]
freq=reverse(1:length(IndsE))/length(IndsE)
freqE=[1;freq]*100

Plots.plot(dose_nominal_E,freqE,label="Nominal",ylabel="Relative Volume (%)", xlabel="Dose (Gy)",color="blue",linestyle=:solid, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_spatial_E,freqE,label="Spatial", color="green", linestyle=:dash,#fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_box_E,freqE,label="Box", color="red",linestyle=:dot,legend=:topright, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=box_b_range
Plots.savefig(@sprintf("./%s/Patient4_Visit1_ExternalDoseHist.png",PlotFolder))


IndsC=firstIndices[2]:firstIndices[3]-1
dose_nominal_C=[0;sort(D[IndsC,:]*x_nominal,dims=1)]
dose_box_C=[0;sort(D[IndsC,:]*x_box,dims=1)]
dose_spatial_C=[0;sort(D[IndsC,:]*x_spatial,dims=1)]
freq=reverse(1:length(IndsC))/length(IndsC)
freqC=[1;freq]*100

Plots.plot(dose_nominal_C,freqC,label="Nominal",ylabel="Relative Volume (%)", xlabel="Dose (Gy)",color="blue",linestyle=:solid, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_spatial_C,freqC,label="Spatial", color="green", linestyle=:dash,#fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_box_C,freqC,label="Box", color="red",linestyle=:dot,legend=:topright, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=box_b_range
Plots.savefig(@sprintf("./%s/Patient4_Visit1_ChaismDoseHist.png",PlotFolder))


Plots.plot(dose_nominal_E,freqE,label="Exterior - Nominal",color="red",linestyle=:solid, #fill=true, fillalpha=0.2,
        xtickfont=font(16),
        ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_spatial_E,freqE,label="Exterior - Spatial", color="red",linestyle=:dash, #fill=true, fillalpha=0.2,
        xtickfont=font(16),
        ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_box_E,freqE,label="Exterior - Box", color="red",linestyle=:dot,legend=:topright, #fill=true, fillalpha=0.2,
        xtickfont=font(16),
        ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=box_b_range
Plots.plot!(dose_nominal_B,freqB,label="Brain - Nominal",ylabel="Relative Volume (%)", xlabel="Dose (Gy)",color="blue",linestyle=:solid, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_spatial_B,freqB,label="Brain - Spatial", color="blue",linestyle=:dash, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_box_B,freqB,label="Brain - Box", color="blue",linestyle=:dot,legend=:topright, #fill=true, fillalpha=0.2,
    xtickfont=font(16),
    ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=box_b_range
Plots.plot!(dose_nominal_C,freqC,label="Chaism - Nominal",ylabel="Relative Volume (%)", xlabel="Dose (Gy)",color="green",linestyle=:solid, #fill=true, fillalpha=0.2,
        xtickfont=font(16),
        ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_spatial_C,freqC,label="Chaism - Spatial", color="green",linestyle=:dash, #fill=true, fillalpha=0.2,
        xtickfont=font(16),
        ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_box_C,freqC,label="Chaism - Box", color="green",linestyle=:dot,#legend=:outertopright, #fill=true, fillalpha=0.2,
        xtickfont=font(16),
        ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=box_b_range
Plots.plot!(dose_nominal_T,freqT,label="PTV - Nominal",ylabel="Relative Volume (%)", xlabel="Dose (Gy)",color="black",linestyle=:solid,# fill=true, fillalpha=0.2,
            xtickfont=font(16),
            ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_spatial_T,freqT,label="PTV - Spatial", color="black",linestyle=:dash,# fill=true, fillalpha=0.2,
            xtickfont=font(16),
            ytickfont=font(16),xguidefontsize=16,yguidefontsize=16,legendfontsize=12)#,bins=nominal_b_range
Plots.plot!(dose_box_T,freqT,label="PTV - Box", color="black",linestyle=:dot, #fill=true, fillalpha=0.2,
            xlims = (0,105),
            xtickfont=font(16),
            ytickfont=font(16),
            xguidefontsize=16,yguidefontsize=16,legendfontsize=12,size=(800,500),
            ylabel="Relative Volume (%)", xlabel="Dose (Gy)",
            left_margin=4mm,bottom_margin=4mm)#,bins=box_b_range
Plots.savefig(@sprintf("./%s/Patient4_Visit1_OARHist.png",PlotFolder))
