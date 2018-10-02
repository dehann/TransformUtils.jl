using TransformUtils
using Base: Test


@testset "test in-place and efficient inverse SE3..." begin

dest = zeros(3)
src = [1.0;0;0]
wTb = SE3([1.0;0;0],SO3(Matrix{Float64}(LinearAlgebra.I,3,3)))

doinversetransform!(dest, wTb, src)

@test norm(dest[1:3]) < 1e-10


# using BenchmarkTools
# @btime for i in 1:1000 doinversetransform!(dest, wTb, src); end
# Profile.clear()
# @profiler for i in 1:1000000 doinversetransform!(dest, wTb, src); end



end
