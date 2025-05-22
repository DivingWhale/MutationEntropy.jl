using MutationEntropy
using Test
using Optim

# Compute the Gamma matrix
# function compute_gamma()
#     eta1, gamma1 = 1, 1
#     eta2, gamma2 = 2, 2
#     coordinates = read_coordinates(pkgdir(MutationEntropy, "data", "1EY0.pdb"))
#     Γ1 = MutationEntropy.gamma(coordinates, eta1, gamma1)
#     Γ2 = MutationEntropy.gamma(coordinates, eta2, gamma2)
#     return Γ = Γ1 + Γ2
# end
# Γ = compute_gamma()

# Read ddG data of Rosetta
ddG_path = "../data"
single_ddG = MutationEntropy.read_ddGs(ddG_path)

# Read the PAE data
task_file_path = "data/task"
WT_pae = MutationEntropy.read_single_pae("data", "thermonulease")
paes = MutationEntropy.read_paes(task_file_path)

# Read the experimental ddG data
using CSV
using DataFrames

ddG_exp = CSV.read("data/Thermonuclease.csv", DataFrame)
ddG_exp = filter(row -> !ismissing(row[:ddG]), ddG_exp)
# ddG_exp[(ddG_exp.position .== 172) .& (ddG_exp.mutation .== "S"), :ddG]

## Process the data
using Statistics

# Define grid search parameters
A_range = range(-100.0, 0.0, step=1)
rho_range = [x for x in 0.0:0.1:3.0 if abs(x - 2.0) > 1e-10]  # Exclude rho = 2

# Initialize variables to store best results
best_pcc = -Inf
best_A = nothing
best_rho = nothing
best_ddgs = nothing
best_r_ddgs = nothing
best_filtered_exp = nothing

# Perform grid search
for A_test in A_range, rho_test in rho_range
    # println("A: ", A_test, " rho: ", rho_test)
    filtered_ddG_exp, ΔΔGs, r_ddGs = calculate_ddgs(
        task_file_path, single_ddG, Γ, WT_pae, paes, ddG_exp, 
        rho_test, A_test
    )
    
    current_pcc = cor(filtered_ddG_exp, ΔΔGs)
    
    if current_pcc > best_pcc
        global best_pcc = current_pcc
        global best_A = A_test
        global best_rho = rho_test
        global best_ddgs = ΔΔGs
        global best_r_ddgs = r_ddGs
        global best_filtered_exp = filtered_ddG_exp
    end
end

# Print results
println("Best parameters found:")
println("A: ", best_A)
println("rho: ", best_rho)
println("Max PCC: ", best_pcc)
println("Correlation between ΔΔG_exp and r_ddGs: ", cor(best_filtered_exp, best_r_ddgs))

## Plot the results
using CairoMakie

function plot_correlation(x, y, title)
    fig = Figure()
    ax = Axis(fig[1, 1], title=title)
    scatter!(ax, x, y, color=:blue, marker=:circle, markersize=5)
    ax.xlabel = "ΔΔG_Rosetta"
    ax.ylabel = "ΔΔS"
    fig
end

plot_correlation(best_r_ddgs, best_ddgs-best_r_ddgs, "Correlation between ΔΔG_Rosetta and ΔΔS")
plot_correlation(best_ddgs, best_filtered_exp, "Correlation between ΔΔG_pred and ΔΔG_exp")