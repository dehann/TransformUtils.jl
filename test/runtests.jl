using TransformUtils
using Base.Test

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
  @test compare( q1, convert(Quaternion, R1) )
  @test compare( R1, convert(SO3, q1) )
  @test compare( q1, convert(Quaternion,AA1) )
  @test compare( R1, convert(SO3,AA1) )
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
