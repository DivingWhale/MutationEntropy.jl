using Test
using MutationEntropy
using LinearAlgebra
using JSON
using DataFrames
using BioStructures

@testset "MutationEntropy.jl Tests" begin

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
        test_PAE_mut = rand(4, 4) * 5  # Random PAE values between 0-5
        test_PAE_wt = rand(4, 4) * 5
        test_dist_mut = rand(4, 4) * 20
        test_dist_wt = rand(4, 4) * 20
        
        # Test ΔΔS calculation
        position = 4
        rho = 1.5
        α = 2.0
        offset = 0
        datadir = tempdir()
        mutation = "A4G"

        data = MutationData([test_PAE_wt], [test_PAE_mut], [test_dist_wt], [test_dist_mut], mutation)
        params = EntropyParams(position, rho, α, offset, false, 90.0, datadir)

        result = ΔΔS(params, data)
        
        # Result should be a real number
        @test isa(result, Real)
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

end
