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

function read_pae(mutation::String)
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
        
            # 如果你需要确保它是矩阵，并且之后使用convert，则可以使用 hcat 和 转置来保证它是一个Matrix{Float64}
            pae =  hcat(pae...)'
          
            if pae_data === nothing
               
                println("Initializing pae_data with the shape:")
                println(size(pae))
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