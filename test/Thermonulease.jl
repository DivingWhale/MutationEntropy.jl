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
    filtered_ddG_exp, ΔΔGs, r_ddGs = MutationEntropy.calculate_ddgs(
        task_file_path, single_ddG, Γ, WT_pae, paes, ddG_exp, 
        rho_test, A_test
    )
    
    current_pcc = cor(filtered_ddG_exp, ΔΔGs)
    
    if current_pcc > best_pcc
        best_pcc = current_pcc
        best_A = A_test
        best_rho = rho_test
        best_ddgs = ΔΔGs
        best_r_ddgs = r_ddGs
        best_filtered_exp = filtered_ddG_exp
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

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, best_filtered_exp, best_r_ddgs, color=:red, marker=:circle, markersize=5, label="ΔΔG")
# scatter!(ax, best_filtered_exp, best_ddgs, color=:blue, marker=:circle, markersize=5, label="ΔΔG_pred")
fig
