using MAT
using FileIO, JLD2
using Printf
using Plots
using ColorSchemes
using DataFrames
using Statistics
using QuantileRegressions
using FileIO, JLD2
#import PyPlot
using PyPlot
using PyCall
using StatsBase
using GLM
using LinearAlgebra
using LaTeXStrings
using JuMP, Gurobi

CodeDir = ".\\"
DataDir ="..\\..\\data\\"
PlotDir = ".\\ResultsFiles\\Plots\\"
ResultDir = ".\\data"
#Patient = "Patient1" #for figure 1
Patient = "Patient4" #for figure 4
if Patient=="Patient1"
    data_file_name="Patient1_Visit1_16beams_refpointpercent50_value26_notincludingdeadvoxels_20230905"
    Visit = "Visit1"
    extra = "refpointpercent50_value26"
elseif Patient=="Patient4"
    data_file_name="Patient4_Visit1_16beams_refpointpercent50_value26_notincludingdeadvoxels_20230719"
    Visit = "Visit1"
    extra = "refpointpercent50_value26"
else
    error("No such patient")
end
include(string(CodeDir, "maxmin_twostage_subprob.jl"))

bLOAD_FROM_FILE = false
bLOAD_DATAFRAME_FROM_FILE = true
CREATE_PLOTS = true

file_name_data_frame =
    @sprintf("%s/%s/%s_%s_%s_Dist_data.jld2",ResultDir,Patient,Patient, Visit, extra)
if bLOAD_DATAFRAME_FROM_FILE
    df = FileIO.load(file_name_data_frame, "df")
    df2 =FileIO.load(file_name_data_frame, "dfSummary")
else
    file_name_gamma = @sprintf(
        "%s/%s/RS_Dists/%s_%s_Gamma_dist_%s_forPloting.jld2",
        ResultDir,
        Patient,
        Patient,
        Visit,
        extra
    )
    if bLOAD_FROM_FILE
        dists = FileIO.load(file_name_gamma, "dists")
        ϕ = FileIO.load(file_name_gamma, "ϕ")
    else
        if Patient=="Patient4" || Patient=="Patient1"
            file = matopen(string(DataDir,Patient,"/",data_file_name,".mat")) #changed from 13 since it did not cover the PTV
        elseif Patient=="liver"
            file = matopen(string(DataDir,"liverEx_2.mat"))
        else
            file = matopen(string(DataDir,Patient,"/",Patient,"_",Visit,"_16beams_",extra,".mat"))
        end
        ϕ = read(file, "omf_Vec")
        γ = read(file, "neighbors_Mat")
        phi_u_n = []
        phi_b_n = []
        dists = []
        gamma_func(x) = x
        γ[diagind(γ)].=0
        phi_u_n, phi_b_n, dists = maxmin_twostage_subprob.computeProjections(
            γ,
            gamma_func,
            ϕ,
            ϕ,
            dists,
        )
        FileIO.save(file_name_gamma, "dists", dists, "ϕ", ϕ)
    end
    phi_u_n = []
    phi_b_n = []
    m, n = size(dists)
    println(m, ' ', n, ' ', n * (n - 1) / 2)
    distVec = zeros(Int32(n * (n - 1) / 2), 1)
    phiDist = zeros(Int32(n * (n - 1) / 2), 1)
    k = 1
    for i = 1:n
        for j = i+1:n
            distVec[k] = dists[i, j]
            phiDist[k] = abs(ϕ[j] - ϕ[i])
            global k = k + 1
        end
    end
    dists = []
    ϕ = []


    df = DataFrame(Dist = vec(distVec), PhiDist = vec(phiDist))
    df2 = combine(
    groupby(df, :Dist),
    :PhiDist => mean => :MeanPhiDist,
    :PhiDist => (t -> quantile(t, 0.80)) => :PhiDistPercent80,
    :PhiDist => (t -> quantile(t, 0.85)) => :PhiDistPercent85,
    :PhiDist => (t -> quantile(t, 0.90)) => :PhiDistPercent90,
    :PhiDist => (t -> quantile(t, 0.95)) => :PhiDistPercent95,
    :PhiDist => (t -> quantile(t, 0.96)) => :PhiDistPercent96,
    :PhiDist => (t -> quantile(t, 0.97)) => :PhiDistPercent97,
    :PhiDist => (t -> quantile(t, 0.98)) => :PhiDistPercent98,
    :PhiDist => (t -> quantile(t, 0.99)) => :PhiDistPercent99,
    :PhiDist => std => :StdPhiDist,
    :PhiDist => maximum => :MaxPhiDist,
    nrow => :count,
    )
    sort!(df2,[:Dist])
    #rename!(df2,:PhiDist_function => :PhiDistPercent98)
    FileIO.save(file_name_data_frame,
    "df",
    df,
    "dfSummary",
    df2,
    )
end

println("Calculate upper bound function")
gamma_coef_file_name=@sprintf("%s/%s/%s_%s_%s_Gamma_dist_func.jld2",ResultDir,Patient,Patient, Visit, extra)
rounding=10^7
max_dist=10
df3 = df2[df2.Dist .<= max_dist,:]
X=hcat(ones(max_dist),log.(df3.Dist),df3.Dist)
y=df3.PhiDistPercent98
max_y_val=maximum(df2.PhiDistPercent98)
m_ub = Model(Gurobi.Optimizer)
@variable(m_ub,α[1:3])
@constraint(m_ub, X*α .>= y)
@constraint(m_ub, X*α .<= X[max_dist,:]'*α)
@constraint(m_ub,X[max_dist,:]'*α>=max_y_val)
@objective(m_ub, Min, sum((X*α-y).^2))
optimize!(m_ub)
D2=df2[:, "Dist"]
coeff=ceil.(value.(α)*rounding)/rounding
calcPhiDistNew=X*coeff
println(calcPhiDistNew)
rounding_new=10^(ceil(-log10(minimum(calcPhiDistNew[2:length(calcPhiDistNew)]-calcPhiDistNew[1:length(calcPhiDistNew)-1]))))
max_ϕ_dist=ceil(max_y_val*rounding_new)/rounding_new
max_ϕ_dist_init=ceil(maximum(calcPhiDistNew)*rounding_new)/rounding_new
max_ϕ_dist=max(max_ϕ_dist,max_ϕ_dist_init)
add_intercept=(max_ϕ_dist-max_ϕ_dist_init)
println(add_intercept)
if add_intercept>0
    coeff[1] += add_intercept
    calcPhiDistNew=X*coeff
    calcPhiDistNew=vcat(calcPhiDistNew,max_ϕ_dist*ones(length(D2)-max_dist))
end
gamma_func_str=@sprintf("gamma_func(x)=ceil((dot(hcat(1,log(x),x)*(x<=%d)*(x>0),[%.10f;%.7f;%.7f])+(x>%d)*%.10f)*%g)/%g",max_dist,coeff[1],coeff[2],coeff[3],max_dist,max_ϕ_dist,rounding_new,rounding_new)
eval(Meta.parse(gamma_func_str))
calcPhiDistNew=gamma_func.(D2)
println(calcPhiDistNew[1:11])
FileIO.save(gamma_coef_file_name,"gamma_func",gamma_func_str,"reg_coeff",coeff,"max_ϕ_dist",max_ϕ_dist,"max_dist",max_dist)

if CREATE_PLOTS
    println("before plotting")

    namess = [
    "Mean",
    "PhiDistPercent80",
    "PhiDistPercent85",
    "PhiDistPercent90",
    "PhiDistPercent95",
    "PhiDistPercent98",
    "MaxPhiDist",
    ]
    labels = ["Mean", "80%", "85%", "90%", "95%", "98%", "Max"]
    linestyles = [:solid, :solid, :dash, :dashdot, :dot, :solid, :solid]


    Plots.plot(
    df2[:, "Dist"],
    df2[:, "MeanPhiDist"],
    label = "Mean",
    xlabel = "Distance between voxels",
    ylabel = "Difference in radiosensitivity",
    legend = :best,
    )


    for i = 2:length(namess)
        c = get(ColorSchemes.rainbow, i ./ Base.length(namess))
        Plots.plot!(
        df2[:, "Dist"],
        df2[:, namess[i]],
        linestyle = linestyles[i],
        color = c,
        label = labels[i],
        xlim = (1, 60),
        legendfont = font(10),
        ytickfont = font(10),
        xtickfont = font(10),
        )
    end
    Plots.savefig(@sprintf("%s/%s/%s_%s_%s_Dist_to_deltaphi.png", PlotDir,Patient,Patient, Visit, extra))


    Plots.plot!(
    df2[:, "Dist"],
    calcPhiDistNew,
    linestyle = :dot,
    color = :black,
    label = "\$\\Gamma\$ w/ \$\\gamma=0\$",
    xlim = (1, 60),
    legendfont = font(10),
    ytickfont = font(10),
    xtickfont = font(10),
    )
    if Patient=="Patient4"
        shift=0.04
    elseif Patient=="Patient1"
        shift=0.025
    end
    Plots.plot!(
        df2[:, "Dist"],
        calcPhiDistNew.+shift,
        linestyle = :dash,
        color = :black,
        label = @sprintf("\$\\Gamma\$ w/ \$\\gamma=%1.2f\$",shift),
        xlim = (1, 60),
        legendfont = font(10),
        ytickfont = font(10),
        xtickfont = font(10),
        )
    Plots.savefig(@sprintf("%s/%s/%s_%s_%s_Dist_to_deltaphi_ub_ref%s.png", PlotDir,Patient,Patient, Visit, extra,max_dist))

    Plots.plot(
    df2[:, "Dist"],
    df2[:, "MeanPhiDist"],
    label = "Mean",
    xlabel = "Distance between voxels",
    ylabel = "Difference in radiosensitivity",
    legend = :best,
    )
    for i = 2:length(namess)
        currc = get(ColorSchemes.rainbow, i ./ Base.length(namess))
        Plots.plot!(
        df2[:, "Dist"],
        df2[:, namess[i]],
        linestyle = linestyles[i],
        color = currc,
        label = labels[i],
        xlim = (1, max_dist),
        legendfont = font(10),
        ytickfont = font(10),
        xtickfont = font(10),
        )
    end
    Plots.savefig(
    @sprintf(
    "%s/%s/%s_%s_%s_Dist_to_deltaphi_%sneighbors.png",
    PlotDir,
    Patient,
    Patient,
    Visit,
    extra,
    max_dist
    )
    )

    stepsize = 0.002
    phiDist=df[:,:PhiDist]
    distVec=df[:,:Dist]
    bins = 0:stepsize:maximum(phiDist)+stepsize
    nbins = Int64(floor((maximum(phiDist) + stepsize) / stepsize))
    h1 = fit(StatsBase.Histogram, vec(phiDist)[vec(distVec).==1], bins)
    h2 = fit(StatsBase.Histogram, vec(phiDist)[vec(distVec).==2], bins)
    h3 = fit(StatsBase.Histogram, vec(phiDist)[vec(distVec).==3], bins)
    hall = fit(StatsBase.Histogram, vec(phiDist), bins)
    currc = get(ColorSchemes.rainbow, 0.75)
    Plots.plot(
    collect(hall.edges[1])[1:nbins],
    hall.weights ./ sum(hall.weights),
    seriestype = :steps,
    color = currc,
    fill = true,
    label = "all",
    xlabel = "Difference in radiosensitivity",
    ylabel = "Frequency",
    legendfont = font(10),
    linewidth = 2,
    thickness_scaling = 1,
    opacity = 0.3,
    )
    currc = get(ColorSchemes.rainbow, 1.0)
    Plots.plot!(
    collect(h1.edges[1])[1:nbins],
    h1.weights ./ sum(h1.weights),
    seriestype = :steps,
    color = currc,
    label = "dist=1",
    linewidth = 2,
    thickness_scaling = 1,
    legendfont = font(10),
    )
    currc = get(ColorSchemes.rainbow, 0.5)
    Plots.plot!(
    collect(h2.edges[1])[1:nbins],
    h2.weights ./ sum(h2.weights),
    seriestype = :steps,
    color = currc,
    label = "dist=2",
    linestyle = :dash,
    linewidth = 2,
    thickness_scaling = 1,
    legendfont = font(10),
    )
    currc = get(ColorSchemes.rainbow, 0.25)
    Plots.plot!(
    collect(h3.edges[1])[1:nbins],
    h3.weights ./ sum(h3.weights),
    seriestype = :steps,
    color = currc,
    label = "dist=3",
    linestyle = :dot,
    linewidth = 2,
    thickness_scaling = 1,
    legendfont = font(10),
    ytickfont = font(10),
    xtickfont = font(10),
    )

    Plots.savefig(@sprintf("%s/%s/%s_%s_%s_PhiDistHist.png", PlotDir,Patient,Patient, Visit, extra))

    println("finished plotting")
end
