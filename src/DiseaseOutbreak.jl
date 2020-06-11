module DiseaseOutbreak

include("distancing.jl")
export makeCitizens, infectCitizens!, followCitizens!
export uninfected, newInfections, cumulativeInfections, infectionSpread

include("dynamics.jl")
export Dynamics, Population

include("evolution.jl")
export evolve

include("sirx.jl")
include("sei3r.jl")
include("sir.jl")
export SEI3R
export SIRX
export SIR
export getParams, initialize
export caseModel

include("casemodel-sirx.jl")
export CaseModelFitResult
export fitCaseModel, summary, fitted, estimatedStates

include("casemodel-sei3r.jl")
export CaseModelFitResult
export fitCaseModel, estimatedStates, predictCases

include("casemodel-sir.jl")
export CaseModelFitResult
export fitCaseModel, summary, fitted, estimatedStates

include("sirplot.jl")
export plotEvolution, plotFit
export pyplotEvolution, pyplotFit

include("helper.jl")
export logit, invlogit
end # module
