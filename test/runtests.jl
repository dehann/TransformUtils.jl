using TransformUtils
using Base.Test

print("[TEST] constructors SO3, Quaternion, AngleAxis... ")
q = Quaternion(0)
ss = so3(0)
R = SO3(0)
A = AngleAxis(0)
T = SE3(0)
SE3(T.t, q)
SE3(T.t, R)
SE3(T.t, A)
println("[SUCCESS]")


print("[TEST] comparison functions for SO3, Quaternion, AngleAxis... ")
@test compare(Quaternion(0),Quaternion(0))
@test !compare(Quaternion(0),Quaternion(0,[1;0;0]))
@test compare(SO3(0),SO3(0))
IrotSO3 = SO3(0)
for i in 1:1000
  r = 0.1*randn(3)
  if norm(r) > 1e-10
    @test !compare(IrotSO3, convert(SO3,so3( r )) )
  end
end
@test compare(AngleAxis(0), AngleAxis(0))
println("[SUCCESS]")

print("[TEST] trivial case quaterion -> SO3 -> Euler -> quaternion... ")
q = Quaternion(0)
R = convert(SO3, q)
E = convert(Euler, R)
@test compare(q,convert(Quaternion, E))
# @show SO3(0) * Quaternion(0) * so3(0) * AngleAxis(0)
@test compare(SO3(0), SO3(0) * Quaternion(0) * so3(0) * AngleAxis(0))
println("[SUCCESS]")

println("[TEST] convert functions ")
warn("Need better coverage on convert function tests")
# print("[SUCSSESS]")

print("[TEST] compare SO3 and quaternion rotations... ")
q = Quaternion(0)
R = convert(SO3,q)
A = convert(AngleAxis, q)

for dAA in
    [AngleAxis(pi/4, [1;0;0]);
    so3(0.1*randn(3));
    so3(0.1*randn(3));
    so3(0.1*randn(3));
    convert(Quaternion, Euler(0.1*randn(3)));
    convert(SO3,Euler(0.1*randn(3)))  ]
  @show dAA
  q1 = q*dAA
  R1 = R*dAA
  AA1 = A*dAA
  @test compare( q1, convert(Quaternion, R1) , tol=1e-6)
  @test compare( R1, convert(SO3, q1) , tol=1e-6)
  @test compare( q1, convert(Quaternion,AA1) , tol=1e-6)
  @test compare( R1, convert(SO3,AA1) , tol=1e-6)
end
println("[SUCCESS]")


print("[TEST] basic SE3 mechanics... ")
# SE3 tests
a = SE3([0;0;0],SO3(0))
b = SE3([1;0;0],SO3(0))
@test compare(a.R,b.R)
c = a*b
@test compare(b,c)

ap = SE3(a.t, a.R*so3(0.1*randn(3)))
@test !compare(ap,a)

back = ap*b*SE3(zeros(3),transpose(ap.R))
@test compare(SO3(0),back.R)
println("[SUCCESS]")

print("[TEST] SE3 ⊕ and ⊖ mechanics... ")
include("se3DevTesting.jl")
println("[SUCCESS]")


println("[TEST] previous discovered issues")
va = SE3(zeros(3),Euler(0,0,0.0))*SE3(zeros(3),Euler(pi/4,0,0))
@test abs(TransformUtils.convert(Euler, va.R).R - pi/4) < 1e-4

va = SE3(zeros(3),Euler(0,0,pi/2))*SE3(zeros(3),Euler(pi/4,0,0))
q = TransformUtils.convert(Quaternion, va.R)
@show TransformUtils.convert(Euler, q)
ce = TransformUtils.convert(Euler, va.R)
@test abs(ce.Y - pi/2) < 1e-8
@test abs(ce.R - pi/4) < 1e-8







# end
