# # Модель Лотки-Вольтерры
# Цель: Исследовать динамику системы хищник-жертва.
#
# ## Инициализация проекта и загрузка пакетов
using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots
using DataFrames
using JLD2

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# ## Определение модели
function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, γ, δ = p
    du[1] = α * x - β * x * y
    du[2] = δ * x * y - γ * y
end

# ## Параметры модели
α = 0.1
β = 0.02
γ = 0.4
δ = 0.01
u0 = [40.0, 9.0]
tspan = (0.0, 200.0)
p = [α, β, γ, δ]

# ## Решение системы
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
sol = solve(prob, Tsit5(), saveat=0.5)

# ## Визуализация
plot(sol.t, [first.(sol.u), last.(sol.u)],
    label=["Жертвы" "Хищники"],
    xlabel="Время", ylabel="Популяция",
    title="Модель Лотки-Вольтерры",
    lw=2, legend=:topright)
savefig(plotsdir(script_name, "lotka_volterra.png"))

df = DataFrame(t=sol.t, prey=first.(sol.u), predator=last.(sol.u))
println("Первые 5 строк:")
println(first(df, 5))
@save datadir(script_name, "lotka_volterra_results.jld2") df
println("Готово!")
