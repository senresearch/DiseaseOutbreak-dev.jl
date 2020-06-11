using Optim
using Calculus

include("helper.jl")
include("sir.jl")

struct CaseModelFitSIR <: CaseModelFitResult
    info::Matrix{Float64}
    inputs::NamedTuple
    fit::Optim.MultivariateOptimizationResults
end

function poissonLogLik( x::Vector{Int64}, λ::Vector{Float64})
    return sum(poissonLogLik.(x,λ))
end

function poissonLogLik(x::Int64,λ::Float64)
    return logpdf(Poisson(λ),x)
end

"""
function fitCaseModel(cases::Vector{Int64},N::Int64,
    ρ::Float64,d0::SIR)
"""
function fitCaseModel(cases::Vector{Int64},N::Int64,
    ρ::Float64,ϕ0::Float64,d0::SIR)

    ntime = length(cases)

    function f(param::Vector{Float64})
        d = SIR(invlogit(param[2] + logit(d0.β)), d0.α)
        s = estimatedStates(ntime,N,(cases[1]/N)/invlogit(param[1]),d)
        caseIntensity = s.I*invlogit(param[1])
        # println(caseIntensity)
        return -poissonLogLik(diff(cases),caseIntensity[2:end])
        # return caseIntensity
    end

    fit = optimize(f,[logit(ρ),logit(ϕ0)],NelderMead())
    d = SIR(invlogit(fit.minimizer[2] + logit(d0.β)),d0.α )
    info = Calculus.hessian(f,fit.minimizer)
    return (cases=cases,N=N,ρ=invlogit(fit.minimizer[1]),
            d=d,info=info,fit=fit)
end

###########
"""
function fitCaseModel(cases::Vector{Int64},x::Matrix{Float64},
    N::Int64,ρ::Float64,ϕ0::Float64,E::Float64,d0::SEI3R)

Fit a model to confirmed cases where the β1 transmission parameter is allowed to depend on covariates x.  The dependence is through a logit link function.
"""
function fitCaseModel(cases::Vector{Int64},x::Matrix{Float64},
    N::Int64,ρ::Float64,E::Float64,d0::SIR)

    d = deepcopy(d0)
    nc = size(x,2)

    ntime = length(cases)
    # dvec = fill(d0,ntime)

    function f(param::Vector{Float64})
        ϕ = hcat(param[3:end])
        betachange = (x*ϕ)[:,1]
        s = estimatedStates(N,cases[1],exp(param[2]),
                            invlogit(param[1]),betachange,d)
        caseIntensity = s.I1*invlogit(param[1])
        # println(caseIntensity)
        return -poissonLogLik(diff(cases),caseIntensity[2:end])
        # return caseIntensity
    end

    fit = optimize(f,[logit(ρ);log(E);zeros(nc)],NelderMead())
    # d.β[1] = invlogit(fit.minimizer[2] + logit(d0.β[1]))
    info = Calculus.hessian(f,fit.minimizer)
    return (cases=cases,N=N,ρ=invlogit(fit.minimizer[1]),
            E0=exp(fit.minimizer[3]),d=d0,info=info,fit=fit)
end

"""
predictCases(fit,ntime)

Predict number of cases using a fit object.
"""
function predictCases(fit,ntime::Int64)
    s = estimatedStates(ntime,fit.N,(fit.cases[1]/fit.N)/fit.ρ,fit.d)
    return cumsum(s.I*fit.ρ)
end

function predictCases(fit,x::Matrix{Float64})
    ϕ = hcat(fit.fit.minimizer[3:end])
    betachange = (x*ϕ)[:,1]
    s = estimatedStates(fit.N,(fit.cases[1]/fit.N)/fit.ρ,
        betachange,fit.d)
    return cumsum(s.I*fit.ρ)
end
