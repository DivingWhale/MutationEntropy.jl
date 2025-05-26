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
mutations = [mutant, "a140g", "a140v", "a142c", "a142f", "a142g", "a142v", "a151g", "a151t"]

# Example: Using the new process_multiple_mutations function with parallel processing
println("Starting parallel processing of multiple mutations...")

# Precompute dimensions to optimize processing
MutationEntropy.precompute_dimensions(data_path, mutations, 1)

# Process all mutations for a specific round in parallel
results = MutationEntropy.process_multiple_mutations(
    data_path, 
    mutations, 
    1,  # Process only round 1
    precache=true,  # Use the precomputed dimensions
    parallel=true   # Enable parallel processing
)

println("Processing complete. Processed $(length(results)) mutations.")

# Get residues within distance for a specific round
nearby_residues = MutationEntropy.find_residues_within_distance(140, "data", 1, mutant)

# Read PAE matrix for a specific round
wt_pae = MutationEntropy.read_paes("data", "thermonulease", 1)

using MutationEntropy
## Test the different alpha values
data_path = "data"
mutant = "a140e"

actual_range = 83:231
alpha = 2.0

# Use the new collect_strains function to calculate average strain values
# for the range of residues (actual_range) over 20 rounds
wt_avg_S_values = MutationEntropy.collect_strains(
    data_path,                # Base directory for data
    "thermonulease",          # Protein name
    alpha=alpha,                # Alpha parameter for strain calculation
    residue_range=actual_range, # Range of residues to analyze
    num_rounds=20,            # Number of rounds to average
    verbose=true              # Print progress information
)

mutant_avg_S_values = MutationEntropy.collect_strains(
    data_path,                # Base directory for data
    mutant,                  # Mutation name
    alpha=alpha,                # Alpha parameter for strain calculation
    residue_range=actual_range, # Range of residues to analyze
    num_rounds=20,            # Number of rounds to average
    verbose=true              # Print progress information
)

mutant_ME = MutationEntropy.calculate_ME(mutant_avg_S_values, wt_avg_S_values)

# Plot ME vs Distance using the CairoMakie implementation
fig = MutationEntropy.plot_MEvsDist(mutant_ME, data_path, mutant)