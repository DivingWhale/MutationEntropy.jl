const RES_DICT = Dict("ALA" => "A", "ARG" => "R", "ASN" => "N", "ASP" => "D", "CYS" => "C", "GLN" => "Q", "GLU" => "E", "GLY" => "G", "HIS" => "H", "ILE" => "I", "LEU" => "L", "LYS" => "K", "MET" => "M", "PHE" => "F", "PRO" => "P", "SER" => "S", "THR" => "T", "TRP" => "W", "TYR" => "Y", "VAL" => "V")
const AA = collect(values(RES_DICT))

# 简单的缓存机制
const FILE_CACHE = Dict{String, Tuple{Float64, Float64}}()  # path -> (ddG, mtime)
const DIR_CACHE = Dict{String, Tuple{Dict{String, Float64}, Float64}}()  # path -> (ddGs, mtime)
const VARIANT_CACHE = Dict{Tuple{String, String}, Tuple{Float64, Float64}}() # (base_dir, variant_name) -> (ddG, mtime)

# 清理缓存的辅助函数
function clear_cache!()
    empty!(FILE_CACHE)
    empty!(DIR_CACHE)
    empty!(VARIANT_CACHE)
end

function read_ddG_file(ddGPath::String)
    # 检查缓存
    if haskey(FILE_CACHE, ddGPath)
        cached_ddG, cached_mtime = FILE_CACHE[ddGPath]
        current_mtime = stat(ddGPath).mtime
        if current_mtime == cached_mtime
            return cached_ddG
        end
    end
    
    # 读取文件
    ddG = 0.0
    open(ddGPath, "r") do file
        lines = readlines(file)
        
        if isempty(lines)
            @warn "Empty ddG file: $ddGPath"
            return 0.0
        end
        
        wt_scores = Float64[]
        mut_scores = Float64[]

        for line in lines
            parts = split(line)
            if length(parts) < 4
                continue # 跳过格式不正确的行
            end

            # 第三列包含 WT_ 或 MUT_ 标识符
            if occursin("WT_", parts[3])
                push!(wt_scores, parse(Float64, parts[4]))
            elseif occursin("MUT_", parts[3])
                push!(mut_scores, parse(Float64, parts[4]))
            end
        end

        if isempty(wt_scores)
            @warn "No WT scores found in $ddGPath"
            return 0.0
        end

        if isempty(mut_scores)
            @warn "No mutant scores found in $ddGPath"
            return 0.0
        end

        native_score = mean(wt_scores)
        mut_ene = mean(mut_scores)
        
        ddG = mut_ene - native_score
    end
    
    # 缓存结果
    FILE_CACHE[ddGPath] = (ddG, stat(ddGPath).mtime)
    return ddG
end

function read_ddg_rosetta(file_path::String)
    if !isdir(file_path)
        throw(ErrorException("The directory $file_path does not exist."))
    end

    # 检查目录缓存
    if haskey(DIR_CACHE, file_path)
        cached_ddGs, cached_mtime = DIR_CACHE[file_path]
        current_mtime = stat(file_path).mtime
        if current_mtime == cached_mtime
            return cached_ddGs
        end
    end

    # 支持单点突变（如 A95G）和多点突变（如 A95G_D206A）
    subDirs = filter(x -> isdir(x) && occursin(r"^[A-Z]\d+[A-Z](?:_[A-Z]\d+[A-Z])*$", basename(x)), readdir(file_path, join=true))
    subDirSorted = sort(subDirs, by=x -> begin
        # 提取第一个突变位点的位置进行排序
        m = match(r"[A-Z](\d+)[A-Z]", basename(x))
        m !== nothing ? parse(Int, m.captures[1]) : 0
    end)
    ddGs = Dict{String, Float64}()

    for subDir in subDirSorted
        ddGFiles = glob("*.ddg", subDir)   
        ddG = 0.0   
        if !isempty(ddGFiles)
            ddG = read_ddG_file(ddGFiles[1])
            ddGs[basename(subDir)] = ddG
        end
    end
    
    # 缓存结果
    DIR_CACHE[file_path] = (ddGs, stat(file_path).mtime)
    return ddGs
end

"""
    read_ddg_foldx(path::String)

Read FoldX ddG data from a directory of FoldX output files.
"""
function read_ddg_foldx(path::String)
    foldx_ddgs = Dict{String, Float64}()
    
    # Iterate through folders in the current directory
    for folder_name in readdir(path)
        folder_path = joinpath(path, folder_name)

        # Check if it's a directory and not "slogs"
        if isdir(folder_path) && folder_name != "slogs"
            fxout_file = joinpath(folder_path, "Average_1EY0.fxout")
            
            # Check if the fxout file exists
            if isfile(fxout_file)
                # Read the last line and extract the third column
                lines = readlines(fxout_file)
                if !isempty(lines)
                    last_line = lines[end]
                    columns = split(last_line)
                    if length(columns) >= 3
                        ddg_value = parse(Float64, columns[2])
                        
                        # Add the data to our dictionary
                        # Assuming folder_name is the mutation name, e.g., "A12G"
                        foldx_ddgs[uppercase(folder_name)] = ddg_value
                    end
                end
            end
        end
    end
    return foldx_ddgs
end

"""
    read_ddg_pythia(path::String)

Read Pythia ddG data from a CSV file.
"""
function read_ddg_pythia(path::String)
    pythia_ddgs = Dict{String, Float64}()
    
    # Read the CSV file
    if isfile(path)
        df = CSV.read(path, DataFrame, header=false)
        
        # Iterate through each row
        for row in eachrow(df)
            mutation_str = string(row[1])  # First column is mutation
            ddg_value = row[2]             # Second column is ddG
            
            # Convert mutation format from "K_1_A" to "k88a" (lowercase)
            # Pythia uses 1-based PDB numbering, need to convert to full sequence numbering
            # PDB position 1 corresponds to full sequence position 88
            parts = split(mutation_str, "_")
            if length(parts) == 3
                wild_type = lowercase(parts[1])
                pdb_position = parse(Int, parts[2])
                mutant = lowercase(parts[3])
                
                # Convert PDB position to full sequence position
                # PDB position 1 -> sequence position 88
                full_seq_position = pdb_position + 87
                
                mutation_clean = wild_type * string(full_seq_position) * mutant
                pythia_ddgs[mutation_clean] = ddg_value
            end
        end
    else
        error("Pythia CSV file not found at: $path")
    end
    
    return pythia_ddgs
end

"""
    read_ddg_stabilityoracle(path::String)

Read StabilityOracle prediction data from a CSV file.
"""
function read_ddg_stabilityoracle(path::String)
    """
    Read StabilityOracle prediction data from a CSV file.
    Expected format: pdb_code,chain_id,mutation,exp_ddG,pred_ddG
    Note: PDB numbering + 82 = biological numbering
    """
    stabilityoracle_ddgs = Dict{String, Float64}()
    
    if !isfile(path)
        error("StabilityOracle file not found: $path")
    end
    
    # Read the CSV file
    df = CSV.read(path, DataFrame)
    
    # Check if required columns exist
    required_cols = ["pdb_code", "chain_id", "mutation", "pred_ddG"]
    for col in required_cols
        if !(col in names(df))
            error("Required column '$col' not found in StabilityOracle data")
        end
    end
    
    # Process each row
    for row in eachrow(df)
        # Skip rows with missing prediction values
        if ismissing(row.pred_ddG)
            continue
        end
        
        # Parse mutation string (e.g., K6A -> wild_type=K, position=6, mutation=A)
        mutation_str = string(row.mutation)
        if length(mutation_str) < 3
            println("Warning: Invalid mutation format: $mutation_str. Skipping...")
            continue
        end
        
        wild_type = mutation_str[1]
        mutation_aa = mutation_str[end]
        pdb_position = parse(Int, mutation_str[2:end-1])
        
        # Convert PDB position to biological position (PDB + 82 = biological)
        bio_position = pdb_position + 82
        
        # Create mutant name in the same format as experimental data (wild_type + position + mutation)
        mutant_name = lowercase(string(wild_type)) * string(bio_position) * lowercase(string(mutation_aa))
        
        stabilityoracle_ddgs[mutant_name] = Float64(row.pred_ddG)
    end
    
    return stabilityoracle_ddgs
end

"""
    read_single_ddG(base_dir::String, variant_name::String)

Reads the ddG value for a single, specified variant from its subdirectory, using a dedicated cache for single variants.

# Arguments
- `base_dir::String`: The base directory containing all variant subdirectories.
- `variant_name::String`: The name of the variant subdirectory (e.g., "A123G" for single mutations or "A123G_D206A" for multiple mutations).

# Returns
- `Union{Float64, Nothing}`: The calculated ddG value if found, otherwise `nothing`.
"""
function read_single_ddG(base_dir::String, variant_name::String)::Union{Float64, Nothing}
    # 支持单点突变和多点突变格式验证
    if !occursin(r"^[A-Z]\d+[A-Z](?:_[A-Z]\d+[A-Z])*$", variant_name)
        @warn "Invalid variant name format: $variant_name"
        return nothing
    end

    variant_dir = joinpath(base_dir, variant_name)
    cache_key = (base_dir, variant_name)

    # 检查目录是否存在
    if !isdir(variant_dir)
        return nothing # 如果目录不存在，则该变体不存在，静默返回
    end

    # 查找 .ddg 文件以检查修改时间
    ddg_files = glob("*.ddg", variant_dir)
    if isempty(ddg_files)
        @warn "No .ddg file found for variant $variant_name in $variant_dir"
        return nothing
    end
    ddg_path = first(ddg_files)
    current_mtime = stat(ddg_path).mtime

    # 检查缓存
    if haskey(VARIANT_CACHE, cache_key)
        cached_ddG, cached_mtime = VARIANT_CACHE[cache_key]
        if current_mtime == cached_mtime
            return cached_ddG
        end
    end

    # 如果缓存未命中或文件已更新，则读取文件
    ddG = read_ddG_file(ddg_path)

    # 存储到缓存
    VARIANT_CACHE[cache_key] = (ddG, current_mtime)
    
    return ddG
end