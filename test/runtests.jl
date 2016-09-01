using TransformUtils
using Base.Test

print("[TEST] comparison functions... ")
@test compare(Quaternion(0),Quaternion(0))
@test !compare(Quaternion(0),Quaternion(0,[1;0;0]))
@test compare(SO3(0),SO3(0))
println("[SUCCESS]")

print("[TEST] trivial case quaterion -> SO3 -> Euler -> quaternion... ")
q = Quaternion(0)
R = convert(SO3, q)
E = convert(Euler, R)
qtest = convert(Quaternion, E)

@test compare(q,qtest)
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

back = ap*b*SE3(zeros(3),T(ap.R))
@test compare(SO3(0),back.R)
println("[SUCCESS]")
