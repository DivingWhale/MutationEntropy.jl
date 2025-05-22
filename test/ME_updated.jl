using Distributed

# Add worker processes if none are available yet
if nprocs() == 1
    # You can adjust the number of processes based on your system resources
    # Using 'Sys.CPU_THREADS' will use all available logical cores
    println("Adding 3 worker processes...")
    addprocs(3)
    println("Worker setup complete. Total processes: $(nprocs())")
end

# Load the required modules on all processes
@everywhere using MutationEntropy

data_path = "data"

task_file = joinpath(data_path, "task")
task_lines = readlines(task_file)

# Read the available mutations from task file or directory structure
mutations = ["a140e", "a140g", "a140v", "a142c", "a142f", "a142g", "a142v", "a151g", "a151t"]

# Example: Using the new process_multiple_mutations function with parallel processing
println("Starting parallel processing of multiple mutations...")

# Process all mutations for a specific round in parallel
results = MutationEntropy.process_multiple_mutations(
    data_path, 
    mutations, 
    1,  # Process only round 1 
    parallel=true   # Enable parallel processing
)

println("Processing complete. Processed $(length(results)) mutations.")

# Get distance map for a specific round
d_matrix = MutationEntropy.get_dist_map(data_path, "a140e", 1)

# Find residues within distance using the distance matrix
nearby_residues = MutationEntropy.find_residues_within_distance(140, d_matrix)

# Read PAE matrix for a specific round
pae = MutationEntropy.read_pae(data_path, "a140e", 1)

# Example of using the low pLDDT function
low_confidence_residues = MutationEntropy.get_low_plddt_residues("a140e", 1, data_path, threshold=85.0)

# Calculate strain for a mutation
mutation = "A140E"
alpha = 0.5
strain = MutationEntropy.calculate_strain(alpha, pae, d_matrix, mutation)

println("Strain value for mutation $mutation with alpha=$alpha: $strain")

# Example of processing for a range of residues
actual_range = 83:231

for k in actual_range
    # This is just a placeholder to demonstrate usage - modify as needed
    if k == actual_range.start
        println("Processing residues in range $actual_range...")
    end
    
    # Example of accessing data for a specific residue
    if d_matrix[k, k] == 0
        println("Processing residue $k")
    end
end
