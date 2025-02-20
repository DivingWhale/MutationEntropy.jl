phi(r, eta, kappa) = exp(-(r/eta)^kappa)

function gamma(coordinates::AbstractVector{Vector{Float64}}, eta::Int, gamma::Int)
    natoms = length(coordinates)
    Γ = zeros(natoms, natoms)
    for i = 1:natoms, j = i+1:natoms
        Γ[i, j] = -phi(norm(coordinates[i] - coordinates[j]), eta, gamma)
    end
    Γ += Γ'
    for i = 1:natoms
        Γ[i, i] = -sum(view(Γ, :, i))
    end
    return Γ
end

function read_coordinates(filename::String)
    WT_structure = read(filename, PDBFormat)
    CA_atoms = collectatoms(WT_structure, sel"name CA")
    coordinates = getfield.(CA_atoms, :coords)
    return coordinates
end

function read_xvg(filename::String)
    result = Tuple{Int, Float64}[]
    for line in eachline(filename)
        if startswith(line, "#") || startswith(line, "@") || isempty(line)
            continue
        end
        parts = split(line)
        push!(result, (parse(Int, parts[1]), parse(Float64, parts[2])))
    end
    return result
end

function msf(Γ::AbstractMatrix)
    return diag(pinv(Γ))
end

function read_single_pae(mutation::String)
    data_path = "/data1/af3/af_output/"
    pae_files = glob("$(mutation)_*", data_path)
    pae_data = nothing
    num_files = 0
    for pae_file in pae_files
        top1_files = glob("seed-*_sample-0", pae_file)
        if !isempty(top1_files)
            top1 = top1_files[1]
            raw_data = JSON.parsefile(joinpath(top1, "confidences.json"))
            pae_any = raw_data["pae"]
            pae = map(x -> map(Float64, x), pae_any)
            pae =  hcat(pae...)'
          
            if pae_data === nothing
                pae_data = zeros(Float64, size(pae))
            end
           
            pae_data .+= pae
            num_files += 1
        end
    end
    if num_files > 0
        return pae_data ./ num_files
    else
        return nothing
    end
end

function read_paes(task_file_path::String)
    paes = Dict{String, Matrix{Float64}}()
    open(task_file_path, "r") do file
        for line in eachline(file)
            mutation = string(strip(line))
            println(mutation)
            pae = read_single_pae(mutation)
            if pae !== nothing
                paes[mutation] = pae
            end
        end
    end
    return paes
end

function ΔΔS(position::Int64, rho::Int64, Γ::Matrix{Float64}, PAE_mut::Matrix{Float64}, PAE_wt::Matrix{Float64})
    indices = findall(x -> isapprox(x, 0.0, atol=1e-6), Γ[position-87, :])
    ΔΔS = 0.0
    for i in indices
        ΔΔS += PAE_mut[position, i]^(2-rho) - PAE_wt[position, i]^(2-rho)
    end
    ΔΔS = ΔΔS ./ length(indices)
    return ΔΔS
end

function ΔΔG_prime(A::Float64, ΔΔS::Float64, ΔΔG::Float64)
    return ΔΔG + A * ΔΔS
end

function calculate_ddgs(task_file_path::String, single_ddG::Dict{String, Float64},  Γ::Matrix{Float64}, WT_pae::Matrix{Float64}, paes::Dict{String, Matrix{Float64}}, ddG_exp::DataFrame, rho::Int64, A::Float64)
    num_muts = countlines(task_file_path)
    ΔΔGs = zeros(num_muts)
    filtered_ddG_exp = zeros(num_muts)
    r_ddGs = zeros(num_muts)
    count = 0
    open(task_file_path, "r") do file
        for line in eachline(file)
            mutation = uppercase(strip(line))
            position = parse(Int, match(r"\d+", mutation).match)
            if haskey(single_ddG, mutation)
                count += 1
                ddG = single_ddG[mutation]
                pae = paes[string(strip(line))]
                ΔΔS = MutationEntropy.ΔΔS(position, rho, Γ, pae, WT_pae)
                ΔΔGs[count] = ΔΔG_prime(A, ΔΔS, ddG)
                ddG_e = mean(ddG_exp[(ddG_exp.position .== position) .& (ddG_exp.mutation .== string(mutation[end])), :ddG])
                filtered_ddG_exp[count] = ddG_e
                r_ddGs[count] = ddG
            else
                continue
            end
        end
    end
    filtered_ddG_exp = filtered_ddG_exp[1:count]
    ΔΔGs = ΔΔGs[1:count]
    r_ddGs = r_ddGs[1:count]
    return filtered_ddG_exp, ΔΔGs, r_ddGs
end

function objective_function(A_vec::Vector{Float64}, task_file_path::String, rho::Int64, single_ddG::Dict{String, Float64}, ddG_exp::DataFrame, Γ::Matrix{Float64}, WT_pae::Matrix{Float64}, paes::Dict{String, Matrix{Float64}})::Float64
    A::Float64 = A_vec[1] # 解包向量中的第一个元素
    filtered_ddG_exp, ΔΔGs, _ = calculate_ddgs(task_file_path, single_ddG, Γ, WT_pae, paes, ddG_exp, rho, A)
    if all(x -> x == 0, filtered_ddG_exp) || all(x -> x == 0, ΔΔGs)
        return Inf # 如果有一个向量全0，返回无穷大来避免计算
    else
        return -cor(filtered_ddG_exp, ΔΔGs)
    end
end