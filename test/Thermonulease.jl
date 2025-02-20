using MutationEntropy:gamma
using MutationEntropy
using Test
using Optim

# Compute the Gamma matrix
eta1, gamma1 = 1, 1
eta2, gamma2 = 2, 2
coordinates = read_coordinates(pkgdir(MutationEntropy, "data", "1EY0.pdb"))
Γ1 = gamma(coordinates, eta1, gamma1)
Γ2 = gamma(coordinates, eta2, gamma2)
Γ = Γ1 + Γ2

# Read ddG data of Rosetta
ddG_path = "/data1/xsqr/entropy/data"
single_ddG = MutationEntropy.read_ddGs(ddG_path)

# Read the PAE data
task_file_path = "/data1/xsqr/entropy/MutationEntropy/data/task"
WT_pae = MutationEntropy.read_single_pae("thermonulease")
paes = MutationEntropy.read_paes(task_file_path)

# Read the experimental ddG data
using CSV
using DataFrames

ddG_exp = CSV.read("/data1/xsqr/entropy/data/Thermonuclease.csv", DataFrame)
ddG_exp = filter(row -> !ismissing(row[:ddG]), ddG_exp)
# ddG_exp[(ddG_exp.position .== 172) .& (ddG_exp.mutation .== "S"), :ddG]

## Process the data
rho = 1
initial_A = [1.0] 
lower_bound = [-20000.0] 
upper_bound = [20000.0] 

# 使用 Optim.jl 进行优化
result = optimize(
    (A_vec) -> MutationEntropy.objective_function(A_vec, task_file_path, rho, single_ddG, ddG_exp, Γ, WT_pae, paes),
    lower_bound,
    upper_bound,
    initial_A,
    Fminbox(LBFGS())
)

# 获取最优 A 值和最大 PCC
optimal_A::Float64 = Optim.minimizer(result)[1]
max_pcc::Float64 = -Optim.minimum(result)
filtered_ddG_exp, ΔΔGs, r_ddGs = MutationEntropy.calculate_ddgs(task_file_path, single_ddG, Γ, WT_pae, paes, ddG_exp, rho, A)
println(cor(r_ddGs, filtered_ddG_exp), " ", cor(ΔΔGs, r_ddGs))

# 打印结果
println("Optimal A: ", optimal_A)
println("Max PCC: ", max_pcc)

## Plot the results
using CairoMakie

A = optimal_A


using Statistics
cor(r_ddGs, filtered_ddG_exp)
cor(ΔΔGs, r_ddGs)

fig = Figure()
ax = Axis(fig[1, 1])
# scatter!(ax, filtered_ddG_exp, ΔΔGs, marker=:circle, markersize=8, label="ΔΔG_exp")
scatter!(ax, filtered_ddG_exp, r_ddGs, color=:red, marker=:circle, markersize=5, label="ΔΔG")
fig
