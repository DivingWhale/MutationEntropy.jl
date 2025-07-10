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

# 定义 extract_parts 函数
"""
    extract_parts(genetic_string::String)

Extracts the numeric and alphabetic parts from a genetic string.

# Arguments
- `genetic_string::String`: The input genetic string.

# Returns
- `Tuple`: A tuple containing the numeric part and the alphabetic part extracted from the genetic string.
           If no match is found, returns `(nothing, nothing)`.
"""
function extract_parts(genetic_string::String)
    m= match(r"(\d+)([A-Za-z]*)", genetic_string)
    return m !== nothing ? (parse(Int, m.captures[1]), m.captures[2]) : (nothing, nothing)
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
        WT_flag = true
        lines = readlines(file)
        native_score = 0.0
        for i in 1:3:length(lines)
            cur_lines = lines[i:min(i+2, length(lines))]
            if isempty(cur_lines)
                break
            end
            if WT_flag
                native_score = mean([parse(Float64, split(line)[4]) for line in cur_lines])
                WT_flag = false
            else
                mut_ene = mean([parse(Float64, split(line)[4]) for line in cur_lines])
                ddG = mut_ene - native_score
            end
        end
    end
    
    # 缓存结果
    FILE_CACHE[ddGPath] = (ddG, stat(ddGPath).mtime)
    return ddG
end

function read_ddGs(file_path::String)
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

    subDirs = filter(x -> isdir(x) && occursin(r"^[A-Z]\d+[A-Z]$", basename(x)), readdir(file_path, join=true))
    subDirSorted = sort(subDirs, by=x -> begin
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
    read_single_ddG(base_dir::String, variant_name::String)

Reads the ddG value for a single, specified variant from its subdirectory, using a dedicated cache for single variants.

# Arguments
- `base_dir::String`: The base directory containing all variant subdirectories.
- `variant_name::String`: The name of the variant subdirectory (e.g., "A123G").

# Returns
- `Union{Float64, Nothing}`: The calculated ddG value if found, otherwise `nothing`.
"""
function read_single_ddG(base_dir::String, variant_name::String)::Union{Float64, Nothing}
    # 确保变体名称格式与 read_ddGs 函数的逻辑一致
    if !occursin(r"^[A-Z]\d+[A-Z]$", variant_name)
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