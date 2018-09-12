#!/usr/bin/env julia

using FourierFlows, Requires, Test

# Run tests

testtime = @elapsed begin

@testset "TracerAdvDiff" begin
  include("test_traceradvdiff.jl")
end

end
println("Total test time: ", testtime)
