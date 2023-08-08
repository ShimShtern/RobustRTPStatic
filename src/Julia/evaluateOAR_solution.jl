module evaluateOAR_solution
include("maxmin_twostage_subprob.jl")
using Printf
using FileIO, JLD2


BIG_NUM = 1e6

export EvaluateOARSolution


function EvaluateOARSolution(d,firstIndices, t,tmax)
    R=length(tmax)
    violProp=zeros(R,1)
    for r in 1:R
        if tmax[r]>t[r]
            violProp[r]=sum(d[firstIndices[r]:firstIndices[r+1]-1].>(t[r]+1e-5))/(firstIndices[r+1]-firstIndices[r])
        end
    end
    return violProp
end
end
