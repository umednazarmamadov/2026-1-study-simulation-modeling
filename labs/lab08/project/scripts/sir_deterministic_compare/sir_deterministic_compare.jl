using DrWatson
@quickactivate "project"
using ConcurrentSim
using ResumableFunctions
using Distributions
using DataFrames
using Random, StatsPlots

mutable struct SIRPerson2
    id::Int64
    status::Symbol
end

mutable struct SIRModel2
    sim::ConcurrentSim.Simulation
    beta::Float64
    c::Float64
    gamma::Float64
    ta::Array{Float64}
    Sa::Array{Int64}
    Ia::Array{Int64}
    Ra::Array{Int64}
    allIndividuals::Array{SIRPerson2}
    deterministic::Bool
end

function infection_update2!(sim, m)
    push!(m.ta, ConcurrentSim.now(sim))
    push!(m.Sa, m.Sa[end] - 1)
    push!(m.Ia, m.Ia[end] + 1)
    push!(m.Ra, m.Ra[end])
end

function recovery_update2!(sim, m)
    push!(m.ta, ConcurrentSim.now(sim))
    push!(m.Sa, m.Sa[end])
    push!(m.Ia, m.Ia[end] - 1)
    push!(m.Ra, m.Ra[end] + 1)
end

@resumable function live2(env::ConcurrentSim.Simulation, individual::SIRPerson2, m::SIRModel2)
    while individual.status == :S
        @yield timeout(env, rand(Exponential(1 / m.c)))
        alter = individual
        while alter == individual
            N = length(m.allIndividuals)
            index = rand(DiscreteUniform(1, N))
            alter = m.allIndividuals[index]
        end
        if alter.status == :I
            if rand(Uniform(0, 1)) < m.beta
                individual.status = :I
                infection_update2!(env, m)
            end
        end
    end
    if individual.status == :I
        if m.deterministic
            @yield timeout(env, 1 / m.gamma)
        else
            @yield timeout(env, rand(Exponential(1 / m.gamma)))
        end
        individual.status = :R
        recovery_update2!(env, m)
    end
end

function run_model(u0, p, tmax, deterministic)
    (S, I, R) = u0
    N = S + I + R
    (beta, c, gamma) = p
    sim = ConcurrentSim.Simulation()
    allIndividuals = SIRPerson2[]
    for i = 1:S
        push!(allIndividuals, SIRPerson2(i, :S))
    end
    for i = (S+1):(S+I)
        push!(allIndividuals, SIRPerson2(i, :I))
    end
    for i = (S+I+1):N
        push!(allIndividuals, SIRPerson2(i, :R))
    end
    ta = Float64[0.0]
    Sa = Int64[S]
    Ia = Int64[I]
    Ra = Int64[R]
    m = SIRModel2(sim, beta, c, gamma, ta, Sa, Ia, Ra, allIndividuals, deterministic)
    [@process live2(m.sim, individual, m) for individual in m.allIndividuals]
    ConcurrentSim.run(m.sim, tmax)
    result = DataFrame()
    result[!, :t] = m.ta
    result[!, :S] = m.Sa
    result[!, :I] = m.Ia
    result[!, :R] = m.Ra
    return result
end

tmax = 40.0
u0 = [990, 10, 0]
p = [0.05, 10.0, 0.25]

Random.seed!(1234)
data_stoch = run_model(u0, p, tmax, false)

Random.seed!(1234)
data_det = run_model(u0, p, tmax, true)

mkpath(plotsdir())

plot(data_stoch.t, data_stoch.I, label = "Стохастическое выздоровление", xlab = "Время", ylab = "Инфицированные I")
plot!(data_det.t, data_det.I, label = "Детерминированное выздоровление")
title!("Сравнение стохастической и детерминированной длительности болезни")

savefig(plotsdir("sir_deterministic_compare.png"))

println("Готово! Сравнение сохранено.")
