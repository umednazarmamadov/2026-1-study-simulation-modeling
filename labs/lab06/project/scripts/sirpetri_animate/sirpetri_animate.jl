using DrWatson
@quickactivate "project"
include(srcdir("SIRPetri.jl"))
using .SIRPetri
using DataFrames, Plots

beta = 0.3
gamma = 0.1
tmax = 100.0

net, u0, states = build_sir_network(beta, gamma)
df = simulate_deterministic(net, u0, (0.0, tmax), saveat = 0.2, rates = [beta, gamma])

anim = @animate for i in 1:nrow(df)
    row = df[i, :]
    bar(
        ["S", "I", "R"],
        [row.S, row.I, row.R],
        legend = false,
        ylims = (0, 1000),
        xlabel = "Группа",
        ylabel = "Численность",
        title = "Время = $(round(row.time, digits=1))",
    )
end

mkpath(plotsdir())
gif(anim, plotsdir("sir_animation.gif"), fps = 10)

println("Анимация сохранена в plots/sir_animation.gif")
