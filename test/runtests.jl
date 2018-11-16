using TransformUtils
using LinearAlgebra
using Test

@testset "Constructors SO3, Quaternion, AngleAxis... " begin
  q = Quaternion(0)
  ss = so3(0)
  R = SO3(0)
  A = AngleAxis(0)
  T = SE3(0)
  SE3(T.t, q)
  SE3(T.t, R)
  SE3(T.t, A)
  E = Euler(0)

  RfromE = SO3(E)
  RfromQ = SO3(q)
  RfromA = SO3(A)
  Rfromss = SO3(ss)

  QfromE = Quaternion(E)
  QfromR = Quaternion(R)
  QfromA = Quaternion(A)
  Qfromss = Quaternion(ss)

  # AfromE = AngleAxis(E) # Convert not implemented yet
  AfromR = AngleAxis(R)
  AfromQ = AngleAxis(q)
  # Afromss = AngleAxis(ss) # Convert not implemented yet

  EfromR = Euler(R)
  EfromQ = Euler(q)
  # EfromA = Euler(A) # Convert not implemented yet
  # Efromss = Euler(ss) # Convert not implemented yet

end

@testset "Ensure basic quaternion operations hold" begin
  q1 = Quaternion(0)
  q2 = convert(Quaternion, TransformUtils.AngleAxis(pi/2, [0,0,1] ))
  q3 = Quaternion(0,[0,0,1])
  qD = deepcopy(q2)
  qmY = Quaternion(1/sqrt(2), [0; -1/sqrt(2); 0])
  qmYt = convert(Quaternion, TransformUtils.AngleAxis(-pi/2, [0;1.0;0]) )

  @test !compare(q1, qD)
  @test compare(q2, qD)
  @test !compare(q3, qD)
  @test compare(qmY, qmYt)
  @test qmYt.s >= 0

  q1D = q1 * qD
  @test compare(q2, q1D)
  q2D = q2 * qD
  @test compare(q3, q2D)
  # special case q2 === qD
  qD2 = qD * q2
  @test compare(q3, qD2)
  # ensure local frame rotations are multiplied on the right
  q1DmY = q1 * qD * qmY
  @test norm(rotate(q1DmY, [1,0,0.0]) - [0,0,1.0]) < 1e-12
  q1mYD = q1 * qmY * qD
  @test norm(rotate(q1mYD, [1,0,0.0]) - [0,1.0,0]) < 1e-12
end

@testset "Comparison functions for SO3, Quaternion, AngleAxis" begin
  @test compare(Quaternion(0),Quaternion(0))
  @test !compare(Quaternion(0),Quaternion(0,[1;0;0]))
  @test compare(SO3(0),SO3(0))
  IrotSO3 = SO3(0)
  for i in 1:1 #000
    r = 0.1*randn(3)
    if norm(r) > 1e-10
      @test !compare(IrotSO3, convert(SO3,so3( r )) )
    end
  end
  @test compare(AngleAxis(0), AngleAxis(0))
end

@testset "Trivial case quaterion -> SO3 -> Euler -> quaternion" begin
  q = Quaternion(0)
  R = convert(SO3, q)
  E = convert(Euler, R)
  @test compare(q,convert(Quaternion, E))
  # @show SO3(0) * Quaternion(0) * so3(0) * AngleAxis(0)
  @test compare(SO3(0), SO3(0) * Quaternion(0) * so3(0) * AngleAxis(0))
end

@warn "Need better coverage on convert function tests"

@testset "Compare SO3 and quaternion rotations" begin
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
end

 @testset "Basic SE3 mechanics" begin
  # SE3 tests
  a = SE3([0;0;0],SO3(0))
  b = SE3([1;0;0],SO3(0))
  @test compare(a.R,b.R)
  c = a*b
  @test compare(b,c)

  ap = SE3(a.t, a.R*so3(0.1*randn(3)))
  @test !compare(ap,a)

  back = ap*b*SE3(zeros(3), transpose(ap.R) )
  @test compare(SO3(0),back.R)
end

@testset "SE3 ⊕ and ⊖ mechanics" begin
  include("se3DevTesting.jl")
end

@testset "Previous discovered issues" begin
  va = SE3(zeros(3),Euler(0,0,0.0))*SE3(zeros(3),Euler(pi/4,0,0))
  @test abs(TransformUtils.convert(Euler, va.R).R - pi/4) < 1e-4

  va = SE3(zeros(3),Euler(0,0,pi/2))*SE3(zeros(3),Euler(pi/4,0,0))
  q = TransformUtils.convert(Quaternion, va.R)
  @show TransformUtils.convert(Euler, q)
  ce = TransformUtils.convert(Euler, va.R)
  @test abs(ce.Y - pi/2) < 1e-8
  @test abs(ce.R - pi/4) < 1e-8
end

@testset "test SO3 logmap..." begin

global R0 = SO3(0)
@test norm(vee(convert(so3, R0))) < 1e-8

for i in 1:100
  val = logmap(convert(SO3,so3(begin global a = 1e-10*randn(3) end))) - a
  @test norm(val) < 1e-8
end

for i in 1:100
  val = logmap(convert(SO3,so3(begin global a = 0.0001*randn(3) end))) - a
  @test norm(val) < 1e-8
end

for i in 1:100
  val = logmap(convert(SO3,so3(begin global a = 0.5*randn(3) end))) - a
  @test norm(val) < 1e-8
end

end

include("testEfficientSE3.jl")





# end
