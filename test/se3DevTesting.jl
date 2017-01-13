# dev test oplus and ominus for SE3

using TransformUtils
using Base.Test


xi = SE3(0)

Dx = SE3([1.0;0.0;0.0], Euler(0.0,0.0,0.0))

xj = xi ⊕ Dx
@test compare(xj, Dx)
@test compare(xi, xj ⊖ Dx)

println("Compare operations on two different SE3 transforms")

Dx1 = SE3(randn(3), convert(SO3, so3(0.25*randn(3))))
Dx2 = SE3(randn(3), convert(SO3, so3(0.25*randn(3))))

xi = SE3(0) ⊕ Dx1
xj = SE3(0) ⊕ Dx2
Dx = SE3(0) ⊖ xi ⊕ xj # inverse(xi) ⊕ xj

@show xi
@show xj
@show Dx
@test compare( xi, xj ⊖ Dx )
@test compare( xj, xi ⊕ Dx )
@test compare( Dx, inverse(xi) ⊕ xj)




















#
