using TransformUtils
using Base.Test

print("[TEST] trivial case quaterion -> SO3 -> Euler -> quaternion... ")
q = Quaternion(0)
R = convert(SO3, q)
E = convert(Euler, R)
qtest = convert(Quaternion, E)
qiq = q*q_conj(qtest)

@test 0.99 <= qiq.s <= 1.0 && norm(qiq.v) < 1e-14
println("[SUCCESS]")

print("[TEST] basic SE3 mechanics and comparisons... ")
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
