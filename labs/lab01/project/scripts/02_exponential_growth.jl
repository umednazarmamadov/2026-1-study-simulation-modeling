using DrWatson
@quickactivate "project"
using DifferentialEquations
using DataFrames
using Plots
using JLD2
using BenchmarkTools

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

function exponential_growth!(du, u, p, t)
    a = p.a
    du[1] = a * u[1]
end

param_grid = Dict(
    :u0 => [[1.0]],
    :a => [0.1, 0.3, 0.5, 0.8, 1.0],
    :tspan => [(0.0, 10.0)],
    :solver => [Tsit5()],
    :saveat => [0.1],
)

all_params = dict_list(param_grid)

p2 = plot(size=(800, 500))
for params in all_params
    prob = ODEProblem(exponential_growth!, params[:u0], params[:tspan], (a=params[:a],))
    sol = solve(prob, params[:solver]; saveat=params[:saveat])
    plot!(p2, sol.t, first.(sol.u), label="a = $(params[:a])", lw=2)
end

plot!(p2,
    xlabel="Время, t",
    ylabel="Популяция, u(t)",
    title="Параметрическое исследование: влияние a",
    legend=:topleft
)
savefig(plotsdir(script_name, "parametric_scan.png"))
println("Готово! Графики сохранены.")
