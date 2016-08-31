using TransformUtils
using Base.Test


q = Quaternion(0)
R = convert(SO3, q)
E = convert(Euler, R)
qtest = convert(Quaternion, E)

qiq = q*q_conj(qtest)

@test 0.99 <= qiq.s <= 1.0
