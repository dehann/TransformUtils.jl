# weighted product of gaussian terms on rotational manifold.

using TransformUtils
# using JuMP
using NLsolve

import Base: convert

# means = zeros(4,2)
# means[1,:] = 1.0


using DrakeVisualizer, GeometryTypes, Colors, CoordinateTransformations, Rotations, StaticArrays


DrakeVisualizer.any_open_windows() || DrakeVisualizer.new_window();

vis = Visualizer()


function drawquat!(vis::DrakeVisualizer.Visualizer,
            q::Rotations.Quat,
            qsym::Symbol;
            length=1.,
            radius=0.01,
            color=RGBA(1.,0,0,0.5) )
    #
    v = SVector(length,0,0)
    line = PolyLine([SVector(0.,0,0), LinearMap(q)(v)], radius=radius)
    linedata = GeometryData(line, color)
    setgeometry!(vis[qsym], linedata)
    nothing
end



function convert(::Type{Rotations.Quat}, q::TransformUtils.Quaternion)
    Rotations.Quat(q.s, q.v...)
end


function distRiemann(q1::Quaternion, q2::Quaternion)
  dR = convert(SO3, q_conj(q1)*q2 )
  norm(logmap(dR))/sqrt(2)
end

# distRiemann(q1,q2)



function weightedmean(
            qArrMu::Vector{Quaternion};
            covariances::Vector{Array{Float64,2}}=Vector{Array{Float64,2}}(),
            tol::Float64=1e-7,
            allowediters::Int=30,
            docovariance::Bool=false  )
  #
  covariances = length(covariances) == length(qArrMu) ? covariances : [eye(3) for i in 1:length(qArrMu)]
  muq = deepcopy(qArrMu[1])
  S = zeros(3)
  invC = zeros(3,3)
  epsi = 1.0
  i = 0
  while epsi > tol && i < allowediters
    @show i+=1
    j=0
    S[:] = 0.0
    invC[:] = 0.0
    for q in qArrMu
      j += 1
      invCov = covariances[j] \ eye(3)
      invC += invCov
      S[:] +=  deltaso3vee(muq, q) #invC *
    end
    dq = 0.5*S #invC \ S
    newmuq = muq * so3(dq)
    epsi = abs(muq.s - newmuq.s) + norm(muq.v - newmuq.v)
    muq = newmuq
  end

  if docovariance
      totCov = invC \ eye(3)
      return muq, totCov
  else
      return muq
  end
end
# a = q1
# b = q2
# c = deepcopy(a)
# Sctoa = deltaso3vee(c, a)
# Sctob = deltaso3vee(c, b)
# S = 0.5*(Sctoa + Sctob)
# d = c * so3(S)
# epsi = abs(c.s - d.s) + norm(c.v-d.v)
# c = d


function quatGuassianProd(q1::Quaternion, cov1, q2::Quaternion, cov2)
  # aa = convert(TransformUtils.AngleAxis, q_conj(q2)*q3)
  inform1 = inv(cov1)
  inform2 = inv(cov2)
  inform = inform1 + inform2

  qD = q_conj(q1) * q2
  s = convert(so3, qD)
  sx = inform \ inform2 * vee(s)
  q1*convert(Quaternion, so3(sx)), inform
end



# test quaternion rotation and commutation conventions are clear
using Base: Test

q1 = Quaternion(0)
q2 = convert(Quaternion, TransformUtils.AngleAxis(pi/2, [0,0,1] ))
q3 = Quaternion(0,[0,0,1])
qD = deepcopy(q2)
qmY = Quaternion(1/sqrt(2), [0; -1/sqrt(2); 0])
qmYt = convert(Quaternion, TransformUtils.AngleAxis(-pi/2, [0;1.0;0]) )

@testset "Ensure basic quaternion operations hold" begin

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

# test the Guassian product

cov1 = eye(3)
cov2 = eye(3)


q12, cov12 = quatGuassianProd(q1, cov1, q2, cov2)

drawquat!(vis, convert(Quat, q1), :q1)
drawquat!(vis, convert(Quat, q2), :q2)
drawquat!(vis, convert(Quat, q12), :q12, color=RGBA(0,1.0,0,0.5))



quatGuassianProd(q1*qmY, cov1, q2, cov2)






sigmas = rand(1,2)
sigmas = ones(1,2)

q1 = Quaternion(0)
q2 = convert(Quaternion, so3(randn(3)))



# analytic example from literature, M Moakher, Means and averaging in the group of rotations

R1 = convert(SO3, q1)
R2 = convert(SO3, q2)

R = R1.R*sqrtm(R1.R'*R2.R)

q12 = convert(Quaternion, SO3(real.(R)))
q12s = convert(Rotations.Quat, q12)



using Optim

# looking at direct optimization
function getdists(X, q1, q2)
    qq = Quaternion(X[1],X[2:4])
    @show unitconst = sum(X[1:4].^2) - 1.0

    # sum([distRiemann(q1,qq); distRiemann(q2,qq)].^2) + X[5]*unitconst
    distRiemann(q2,qq) + X[5]*unitconst
end

X = [1;0;0;0]
distRiemann(q2,Quaternion(X[1],X[2:4]))

result = optimize((x) -> getdists(x, q1, q2), [1.,0,0,0,1000])

q12A = result.minimizer
sum(q12A[1:4].^2)
q12 = Quaternion(q12A[1],q12A[2:4])

q12s = convert(Quat, q12)

# m = Model()
#
# @variable(m, 0.0 <= w <=1.0 )
# @variable(m, -1.0 <= x <=1.0 )
# @variable(m, -1.0 <= y <=1.0 )
# @variable(m, -1.0 <= z <=1.0 )
#
# @constraint(m, w^2 + x^2 + y^2 + z^2 == 1)
# @objective(m, Min, getdists(getvalue(w),getvalue(x),getvalue(y),getvalue(z), q1, q2) )
# print(m)



# Dehann working

q1s = convert(Rotations.Quat, q1)
q2s = convert(Rotations.Quat, q2)


q12 = weightedmean([q1,q2], docovariance=true)
q12s = convert(Rotations.Quat, q12[1])



## Visualize for all

drawquat!(vis, q1s, :q1)
drawquat!(vis, q2s, :q2, color=RGBA(0,1.,0,0.5))

drawquat!(vis, q12s, :q12, color=RGBA(1.,1.,0,0.5))


# q12t = convert(Rotations.Quat, d)
# drawquat!(vis, q12t, :q12t, color=RGBA(0.,0.,1,0.5))








#
