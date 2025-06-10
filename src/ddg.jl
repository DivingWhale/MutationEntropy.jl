const RES_DICT = Dict("ALA" => "A", "ARG" => "R", "ASN" => "N", "ASP" => "D", "CYS" => "C", "GLN" => "Q", "GLU" => "E", "GLY" => "G", "HIS" => "H", "ILE" => "I", "LEU" => "L", "LYS" => "K", "MET" => "M", "PHE" => "F", "PRO" => "P", "SER" => "S", "THR" => "T", "TRP" => "W", "TYR" => "Y", "VAL" => "V")
const AA = collect(values(RES_DICT))

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

# 定义 read_single_ddG 函数
"""
    read_single_ddG(file_path::String, column_indices::Vector{Int})

Reads the ddG files in the specified directory and returns a dictionary of ddG values.

# Arguments
- `file_path::String`: The path to the directory containing the ddG files.
- `column_indices::Vector{Int}`: List of column indices to be used for ddG calculations.

# Returns
- `Dict`: A dictionary containing the ddG values. The keys are in the format 'index_residue', and the values are dictionaries containing the ddG values for each residue.

# Raises
- `ErrorException`: If the specified directory does not exist.
"""
# function read_single_ddG(file_path::String, column_indices::Vector{Int})
#     if !isdir(file_path)
#         throw(ErrorException("The directory $file_path does not exist."))
#     end

#     subDirs = filter(isdir, readdir(file_path, join=true))
#     subDirSorted = sort(subDirs, by=x -> parse(Int, split(basename(x), '_')[1]))
#     single_ddG = Dict{String, Dict{String, Float64}}()

#     for subDir in subDirSorted
#         ddGFiles = glob("*.ddg", subDir)
#         ddGs = Dict{String, Float64}()
#         if !isempty(ddGFiles)
#             ddGPath = ddGFiles[1]
#             open(ddGPath, "r") do file
#                 native_residue = split(basename(dirname(ddGPath)), '_')[2]
#                 native_index = parse(Int, split(basename(dirname(ddGPath)), '_')[1]) + 34
#                 WT_flag = true
#                 lines = readlines(file)
#                 native_score = 0.0
#                 for i in 1:3:length(lines)
#                     cur_lines = lines[i:min(i+2, length(lines))]
#                     if isempty(cur_lines)
#                         break
#                     end
#                     if WT_flag
#                         native_score = mean([parse(Float64, split(line)[4]) for line in cur_lines])
#                         adjustment = mean([parse(Float64, split(line)[index]) for index in column_indices for line in cur_lines])
#                         native_score -= adjustment
#                         ddGs[native_residue] = 0.0
#                         WT_flag = false
#                     else
#                         current_residue = split(split(cur_lines[1])[3], '_')[2]
#                         _, res = extract_parts(string(current_residue))
#                         mut_ene = mean([parse(Float64, split(line)[4]) for line in cur_lines])
#                         adjustment = mean([parse(Float64, split(line)[index]) for index in column_indices for line in cur_lines])
#                         ddG = mut_ene - native_score - adjustment
#                         ddGs[RES_DICT[res]] = ddG
#                     end
#                 end
#                 native_res = "$native_index$native_residue"
#                 single_ddG[native_res] = ddGs
#             end
#         end
#     end
#     return single_ddG
# end

function read_ddG_file(ddGPath::String)
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
    return ddG
end

function read_ddGs(file_path::String)
    if !isdir(file_path)
        throw(ErrorException("The directory $file_path does not exist."))
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
    return ddGs
end