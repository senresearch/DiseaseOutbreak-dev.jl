## implement Markov Chain with SEIR
include("dynamics.jl")

struct SIR <: Dynamics
    α::Float64
    β::Float64
end

# default constructor
SIR() = SIR(0.0,0.0)

# change in a day
"""
change(s::Vector{Float64},d::SIR)

Retutn the change in the state of the population in a day

s = state vector (S,I,R)
d = SIR dynamics parameters
"""
function change(s::Vector{Float64},d::SIR)
    # this is done purely for readability of the formula
    st = NamedTuple{(:S,:I,:R)}(s)
    # N = sum(s) # population size
    S = -( d.β*st.I )*st.S
    I =  ( d.β*st.I )*st.S - d.α*st.R
    R = d.α*st.R
    return [S,I,R]
end

function initialize(I::Float64,d::SIR)
    state = zeros(nstates(d))
    state[2] = I
    state[1] = 1-I
    return state
end

function nstates(d::SIR)
    return 3
end

function stateNames(d::SIR)
    return ["S" "I" "R"]
end

function getParams(α::Float64,β::Float64,d::SIR)
    return SIR(α,β)
end

"""
estimatedStates( nt::Int64, N::Int64, I0::Float64, d::SIR)

nt: time points
N: population size
I0: proportion infected at start
d: SIR parameters
"""
function estimatedStates(nt::Int64,N::Int64,I0::Float64,d::SIR)
    s0 = initialize(I0,d) # initial state
    (d,ds) = evolve(N*1.0,s0,d,nt) # evolve over time
    s1 = DataFrame(N.*d',[:S,:I,:R]) # data frame with states over time
    return s1
end


## betachange is how the β1 parameter is changed on the logit scale
function estimatedStates( N::Int64, E0::Float64,  betachange::Vector{Float64}, d::SEI3R)
    s0 = initialize(E0,d)
    (d,ds) = evolve(N*1.0,s0,betachange,d)
    s1 = DataFrame(N.*d',[:S,:I,:R])
    return s1
end
