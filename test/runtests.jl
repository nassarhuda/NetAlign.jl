using NetworkAlignment
using MatrixNetworks
using Test
using Random
using Statistics
using SparseArrays
using LinearAlgebra

all_tests = [
             "networkalignment"]

for ti = 1:length(all_tests)
    t = all_tests[ti]
    test_name = join(["$(t)", "_test",".jl"])
    @show test_name
    test_path = joinpath(dirname(@__FILE__), test_name)
    include(test_path)
end
