# profile/test conversion utils

using TransformUtils
using Test
using BenchmarkTools, ProfileView

R = SO3(0)
E = Euler(0)

RR = convert(SO3, so3(randn(3)))
E = convert(Euler, RR)

# pvt: 371.464 ms (5000000 allocations: 534.06 MiB)
@btime for i in 1:10^6 convert!(R, E); end
@profiler for i in 1:10^6 convert!(R, E); end






T1 = SE3(zeros(3),R)
T2 = SE3(0)

# pvt: 520.225 ns (20 allocations: 1.39 KiB - norm()
@btime T3 = T1 \ T2

# pvt: 519.044 μs (20000 allocations: 1.36 MiB)
@btime for i in 1:10^3 T3 = T1 \ T2; end
using Profile
Profile.clear()
@profile for i in 1:10^3 T3 = T1 \ T2; end

Juno.profiler()
Juno.profiletree()




@testset "From: quaternion" begin
@test
end



function fast_norm_sonar(u::::Vector{Float64})
  s = 0.0
  @fastmath @inbounds @simd for i in 1:n
      s += u[i]^2
  end
  @fastmath @inbounds return sqrt(s/n)
end






import Base: normalize!

function normalize!(q::Quaternion, tol=1e-6)

  s = q.s^2
  @fastmath @inbounds @simd for i in 1:3
      s += q.v[i]^2
  end
  @fastmath mag = sqrt(s)
  @fastmath q.v ./= mag
  @fastmath q.s /= mag
  nothing
end



@code_warntype normalize!(q)

@code_llvm normalize!(q)

# code_native(normalize!, (Quaternion,))
@code_native normalize!(q)


@btime normalize!(q) # 17ns


@btime for i in 1:1000 normalize!(q); end

@time normalize!(q)






q = convert(Quaternion, so3(randn(3)))
E = convert(Euler, convert(SO3, so3(randn(3))))



@btime convert!(q, E)

#
Profile.clear()
@profiler convert!(q, E)


Profile.clear()
@profiler for i in 1:10^6 convert!(q, E); end




@code_warntype convert!(q, E)

# julia --track-allocation=user testInPlaceConvert.jl
