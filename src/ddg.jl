const RES_DICT = Dict("ALA" => "A", "ARG" => "R", "ASN" => "N", "ASP" => "D", "CYS" => "C", "GLN" => "Q", "GLU" => "E", "GLY" => "G", "HIS" => "H", "ILE" => "I", "LEU" => "L", "LYS" => "K", "MET" => "M", "PHE" => "F", "PRO" => "P", "SER" => "S", "THR" => "T", "TRP" => "W", "TYR" => "Y", "VAL" => "V")
const AA = collect(values(RES_DICT))

# 简单的缓存机制
const FILE_CACHE = Dict{String, Tuple{Float64, Float64}}()  # path -> (ddG, mtime)
const DIR_CACHE = Dict{String, Tuple{Dict{String, Float64}, Float64}}()  # path -> (ddGs, mtime)

# 清理缓存的辅助函数
function clear_cache!()
    empty!(FILE_CACHE)
    empty!(DIR_CACHE)
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
                native_score = mean([parse(Float64, split(line)[4]) - parse(Float64, split(line)[8]) for line in cur_lines])
                WT_flag = false
            else
                mut_ene = mean([parse(Float64, split(line)[4]) - parse(Float64, split(line)[8]) for line in cur_lines])
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