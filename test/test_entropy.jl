using Test
using MutationEntropy
using LinearAlgebra
using JSON
using DataFrames
using BioStructures

@testset "MutationEntropy.jl Tests" begin

    @testset "Φ function tests" begin
        # Test basic functionality of the Φ function
        @test MutationEntropy.Φ(0.0, 20, 7) ≈ 1.0
        @test MutationEntropy.Φ(20.0, 20, 7) ≈ exp(-1.0)
        @test 0.0 < MutationEntropy.Φ(10.0, 20, 7) < 1.0
        
        # Test with different parameters
        @test MutationEntropy.Φ(0.0, 13, 10) ≈ 1.0
        @test MutationEntropy.Φ(13.0, 13, 10) ≈ exp(-1.0)
    end

    @testset "Γ matrix function tests" begin
        # Create simple test coordinates
        test_coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ]
        
        # Test Γ matrix computation
        Γ_matrix = MutationEntropy.Γ(test_coords, 20, 7)
        
        # Check dimensions
        @test size(Γ_matrix) == (3, 3)
        
        # Check symmetry
        @test Γ_matrix ≈ Γ_matrix'
        
        # Check that diagonal elements are negative sum of off-diagonal elements
        for i in 1:3
            @test Γ_matrix[i, i] ≈ -sum(Γ_matrix[i, 1:end .!= i])
        end
        
        # Check that row sums are approximately zero (due to construction)
        for i in 1:3
            @test sum(Γ_matrix[i, :]) ≈ 0.0 atol=1e-10
        end
    end

    @testset "compute_Γ function tests" begin
        # Create test coordinates
        test_coords = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0]
        ]
        
        # Test compute_Γ function
        Γ_combined = MutationEntropy.compute_Γ(test_coords)
        
        # Check dimensions
        @test size(Γ_combined) == (4, 4)
        
        # Check symmetry
        @test Γ_combined ≈ Γ_combined'
        
        # Verify it's a sum of two Γ matrices
        Γ1 = MutationEntropy.Γ(test_coords, 20, 7)
        Γ2 = MutationEntropy.Γ(test_coords, 13, 10)
        @test Γ_combined ≈ Γ1 + Γ2
    end

    @testset "msf function tests" begin
        # Create a simple test matrix
        test_matrix = [2.0 -1.0; -1.0 2.0]
        
        # Test msf function
        msf_result = MutationEntropy.msf(test_matrix)
        
        # Check that result is a vector with correct length
        @test length(msf_result) == 2
        @test all(msf_result .> 0)  # MSF values should be positive
        
        # Test with singular matrix (should still work with pinv)
        singular_matrix = [1.0 1.0; 1.0 1.0]
        msf_singular = MutationEntropy.msf(singular_matrix)
        @test length(msf_singular) == 2
    end

    @testset "read_xvg function tests" begin
        # Create a temporary xvg file for testing
        temp_xvg_path = tempname() * ".xvg"
        
        test_content = """
        # This is a comment
        @ title "Test XVG"
        @ xaxis label "Time (ps)"
        @ yaxis label "Value"
        
        1    1.5
        2    2.3
        3    1.8
        4    2.1
        """
        
        write(temp_xvg_path, test_content)
        
        try
            # Test reading the xvg file
            result = MutationEntropy.read_xvg(temp_xvg_path)
            
            # Check that we get the expected number of data points
            @test length(result) == 4
            
            # Check the data types
            @test all(x -> isa(x, Tuple{Int, Float64}), result)
            
            # Check specific values
            @test result[1] == (1, 1.5)
            @test result[2] == (2, 2.3)
            @test result[3] == (3, 1.8)
            @test result[4] == (4, 2.1)
            
        finally
            # Clean up temporary file
            isfile(temp_xvg_path) && rm(temp_xvg_path)
        end
    end

    @testset "ΔΔS function tests" begin
        # Create test matrices
        test_Γ = [
            -2.0  1.0  1.0  0.0;
             1.0 -2.0  1.0  0.0;
             1.0  1.0 -2.0  0.0;
             0.0  0.0  0.0  0.0
        ]
        
        test_PAE_mut = rand(4, 4) * 5  # Random PAE values between 0-5
        test_PAE_wt = rand(4, 4) * 5
        
        # Test ΔΔS calculation
        position = 4  # Position with zero connections
        rho = 1.5
        
        # Create test data for the new parameters
        datadir = tempdir()  # Use temp directory for test
        mutation = "TestMutation123"  # Valid non-empty mutation string
        
        # For testing, we'll create mock distance matrices
        test_dist_matrix = rand(10, 10) .+ 1.0  # Ensure positive distances
        for i in 1:10
            test_dist_matrix[i, i] = 0.0  # Self-distance is 0
        end
        
        # Create entropy configuration for new interface - now uses vector of matrices
        test_wt_dist_matrices = [test_dist_matrix]  # Single round for testing
        config = MutationEntropy.EntropyConfig(datadir; wt_identifier="test", wt_dist_matrices=test_wt_dist_matrices)
        
        # We'll pass the WT distance matrices directly to avoid file I/O in tests
        # For testing, use single-element vectors to match new multi-round structure
        test_PAE_mut_vector = [test_PAE_mut]
        test_PAE_wt_vector = [test_PAE_wt]
        α = 2.5  # Test α parameter
        result = MutationEntropy.ΔΔS(position, rho, α, test_Γ, test_PAE_mut_vector, test_PAE_wt_vector, 
                                     mutation, config)
        
        # Result should be a real number
        @test isa(result, Real)
        
        # Test with offset
        config_offset = MutationEntropy.EntropyConfig(datadir; wt_identifier="test", wt_dist_matrices=test_wt_dist_matrices, offset=1)
        result_offset = MutationEntropy.ΔΔS(position + 1, rho, α, test_Γ, test_PAE_mut_vector, test_PAE_wt_vector, 
                                           mutation, config_offset)
        # Handle NaN case properly - if both are NaN, they should be considered equal for this test
        @test (isnan(result) && isnan(result_offset)) || result ≈ result_offset
    end

    @testset "ΔΔG_prime function tests" begin
        A = 0.5
        ΔΔS = 2.0
        ΔΔG = 1.0
        
        result = MutationEntropy.ΔΔG_prime(A, ΔΔS, ΔΔG)
        
        # Check expected calculation
        @test result ≈ ΔΔG + A * ΔΔS
        @test result ≈ 2.0
        
        # Test with different values
        result2 = MutationEntropy.ΔΔG_prime(0.0, ΔΔS, ΔΔG)
        @test result2 ≈ ΔΔG
    end

    @testset "parse_mutation_position function tests" begin
        # Test various mutation string formats
        @test MutationEntropy.parse_mutation_position("A123B") == 123
        @test MutationEntropy.parse_mutation_position("VAL42ALA") == 42
        @test MutationEntropy.parse_mutation_position("G1A") == 1
        @test MutationEntropy.parse_mutation_position("PRO256GLY") == 256
    end

    @testset "get_experimental_ddg function tests" begin
        # Create test DataFrame
        test_df = DataFrame(
            position = [10, 10, 20, 20, 30],
            mutation = ["A", "A", "V", "V", "G"],
            ddG = [1.5, 1.7, -0.8, -0.6, 2.3]
        )
        
        # Test getting experimental ddG
        result1 = MutationEntropy.get_experimental_ddg(test_df, 10, "A")
        @test result1 ≈ 1.6  # Mean of 1.5 and 1.7
        
        result2 = MutationEntropy.get_experimental_ddg(test_df, 20, "V")
        @test result2 ≈ -0.7  # Mean of -0.8 and -0.6
        
        result3 = MutationEntropy.get_experimental_ddg(test_df, 30, "G")
        @test result3 ≈ 2.3  # Single value
    end

    @testset "read_mutations_from_file function tests" begin
        # Create a temporary file with mutation list
        temp_file_path = tempname()
        
        test_mutations = """
        A123B
        VAL42ALA
        G1A
        PRO256GLY
        """
        
        write(temp_file_path, test_mutations)
        
        try
            result = MutationEntropy.read_mutations_from_file(temp_file_path)
            
            @test length(result) == 4
            @test result[1] == "A123B"
            @test result[2] == "VAL42ALA"
            @test result[3] == "G1A"
            @test result[4] == "PRO256GLY"
            
        finally
            isfile(temp_file_path) && rm(temp_file_path)
        end
    end

    @testset "Integration tests" begin
        # Test that functions work together
        test_coords = [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],  # Typical CA-CA distance
            [0.0, 3.8, 0.0]
        ]
        
        # Test full workflow from coordinates to MSF
        Γ_matrix = MutationEntropy.compute_Γ(test_coords)
        msf_values = MutationEntropy.msf(Γ_matrix)
        
        @test length(msf_values) == 3
        @test all(msf_values .> 0)
        
        # Test PAE processing workflow
        test_PAE = [
            1.0 2.0 3.0;
            2.0 1.0 2.5;
            3.0 2.5 1.0
        ]
        
        # Calculate ΔΔS for all positions
        rho = 1.5
        # Create test parameters for ΔΔS
        test_datadir = tempdir()
        test_mutation = "TestMutation"
        test_dist_matrix = rand(3, 3) .+ 1.0
        for i in 1:3
            test_dist_matrix[i, i] = 0.0
        end
        
        # Create entropy configuration for new interface
        test_wt_dist_matrices = [test_dist_matrix]  # Single round for testing
        config = MutationEntropy.EntropyConfig(test_datadir; wt_identifier="test", wt_dist_matrices=test_wt_dist_matrices)
        
        for pos in 1:3
            # Use vector format for new multi-round structure
            test_PAE_vector = [test_PAE]
            test_PAE_wt_vector = [test_PAE * 0.9]
            α = 2.5  # Test α parameter
            result = MutationEntropy.ΔΔS(pos, rho, α, Γ_matrix, test_PAE_vector, test_PAE_wt_vector, 
                                         test_mutation, config)
            @test isa(result, Real)
        end
    end

end