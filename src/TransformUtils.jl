__precompile__(true)

module TransformUtils

using LinearAlgebra

import Base: convert, promote_rule, *, \, vec
import LinearAlgebra: transpose, normalize, normalize!

export
  Quaternion,
  AngleAxis,
  AxisAngle,
  SO3,
  so3,
  SE3,
  Euler,
  skew,
  vee!,
  vee,
  vec,
  setso3!,
  *,
  transpose,
  matrix,
  matrix!,
  inverse,
  compare,
  normalize!,
  normalize,
  q_conj,
  q_conj!,
  convert,
  convert!,
  rotate!,
  rotate,
  wrapRad,
  logmap,
  rightJacExmap,
  rightJacExmapinv,
  deltaso3vee,
  # Should be good to go
  veeAngleAxis,
  veeQuaternion,
  veeEuler,
  A_invB,
  ominus,
  oplus,
  ⊖,
  ⊕,
  \,
  doinversetransform!,

  # type aliases
  FloatInt,
  VectorFloatInt,

  # basic Taylor exported to try against Pade version in expm
  expmOwn,
  expmOwn1,
  expmOwnT,
  expmOwn2,
  expmOwn3,
  expmOwn4,


  # TODO -- needs refactoring
  R,
  SE2,
  se2vee,
  se2vee!



  function skew(v::Array{Float64,1})
      S = zeros(3,3)
      S[1,:] = [0., -v[3], v[2]]
      S[2,:] = [v[3],0.,-v[1]]
      S[3,:] = [-v[2],v[1],0.]
      return S
  end

const FloatInt = Union{Float64,Int}
const VectorFloatInt = Union{Vector{Float64},Vector{Int}}


mutable struct Quaternion
    s::Float64
    v::Array{Float64,1}
    Quaternion() = new()
    Quaternion(s::FloatInt) = new(1.0,zeros(3))
    Quaternion(s::FloatInt,v::VectorFloatInt) = new(s,v)
end

mutable struct AngleAxis
    theta::Float64
    ax::Array{Float64,1}
    AngleAxis() = new()
    AngleAxis(s::FloatInt) = new(0,[1;0;0])
    AngleAxis(s::FloatInt,v::VectorFloatInt) = new(s,v)
end

const AxisAngle = AngleAxis

mutable struct SO3
    R::Array{Float64,2}
    SO3() = new()
    SO3(dummy::FloatInt) = new(Matrix{Float64}(LinearAlgebra.I, 3,3))
    SO3(r::Array{Float64,2}) = new(r)
end

mutable struct so3
    S::Array{Float64,2}
    so3() = new()
    so3(s::FloatInt) = new(zeros(3,3))
    so3(v::VectorFloatInt) = new(skew(v))
    so3(S::Array{Float64,2}) = new(S)
end

function setso3!(s::so3, v::Union{Array, SubArray})
  s.S[1,2] = -v[3]
  s.S[1,3] = v[2]
  s.S[2,3] = -v[1]
  s.S[2,1] = v[3]
  s.S[3,1] = -v[2]
  s.S[3,2] = v[1]
  nothing
end

mutable struct Euler
    R::Float64
    P::Float64
    Y::Float64
    # these will be inconsistent with R, P, Y -- not for public consumption...
    fastconvert::Quaternion
    Euler() = new()
    Euler(s::FloatInt) = new(0.0,0.0,0.0, Quaternion(0))
    Euler(r::FloatInt,p::FloatInt,y::FloatInt) = new(r,p,y, Quaternion(0))
    Euler(r::FloatInt,p::FloatInt,y::FloatInt,q::Quaternion) = new(r,p,y, q)
    Euler(v::VectorFloatInt) = new(v[1],v[2],v[3], Quaternion(0))
    Euler(v::VectorFloatInt, q::Quaternion) = new(v[1],v[2],v[3], q)
end


mutable struct SE3
  R::SO3
  t::Vector{Float64}
  SE3() = new()
  SE3(dummy::FloatInt) = new(SO3(0.0), zeros(3))
  SE3(t::VectorFloatInt, r::SO3) = new(r,t)
  SE3(v::VectorFloatInt, E::Euler) = new(convert(SO3,E), v[1:3] )
  SE3(v::VectorFloatInt, aa::AngleAxis) = new(convert(SO3,aa), v[1:3] ) # maybe replace with promote via dispatch
  SE3(v::VectorFloatInt, q::Quaternion) = new(convert(SO3,q), v[1:3] ) # maybe replace with promote via dispatch
  SE3(M::Array{Float64,2}) = new(SO3(M[1:3,1:3]), vec(M[1:3,4]) )
end

function fast_norm_TU(u)
  # dest[1] = ...
  n = length(u)
  T = eltype(u)
  s = zero(T)
  @fastmath @inbounds @simd for i in 1:n
      s += u[i]^2
  end
  @fastmath @inbounds return sqrt(s)
end

function normalize!(q::Quaternion, tol=1e-6)
  # Quaternion(s, [x,y,z])
  s = q.s^2
  @fastmath @inbounds @simd for i in 1:3
      s += q.v[i]^2
  end
  @fastmath mag1 = 1.0/sqrt(s)
  # @fastmath @inbounds for i in 1:3
  #   q.v[i] *= mag1
  # end
  @fastmath @inbounds q.v .*= mag1
  @fastmath q.s *= mag1
  nothing
end

function normalize(q::Quaternion, tol=0.00001)
  qq = deepcopy(q)
  normalize!(qq,tol)
  return qq
end



function normalize(v::Array{Float64,1})
  return v / fast_norm_TU(v)
end


function matrix!(M::Array{Float64,2}, a::SE3)
  M[1:3,1:3] = a.R.R
  M[1:3,4] = a.t
  nothing
end
# Return 4x4 matrix version of SE3 transform, also see matrix! for in place version
function matrix(a::SE3)
  T = Matrix{Float64}(LinearAlgebra.I,4,4)
  matrix!(T, a)
  return T
end


function q_conj!(q::Quaternion)
    normalize!(q)
    q.v = -q.v
    nothing
end
function q_conj(q::Quaternion)
    qq = deepcopy(q)
    q_conj!(qq)
    return qq
end


transpose(a::SO3) = SO3(collect(a.R')) # TODO use adjoint instead of collect
inverse(a::SO3) = transpose(a)




# put result back in a
# function A_invB!(a::SE3, b::SE3)
#  = SE3( ( matrix(b)' \ (matrix(a)') )' )
# end

function *(a::SO3, b::SO3)
  return SO3(a.R*b.R)
end



function *(a::SE3, b::SE3)
  return SE3(vec(a.R.R*b.t + a.t), a.R*b.R)
end

function *(q1::Quaternion, q2::Quaternion)
  ee  = [q1.s; q1.v] * [q2.s; q2.v]'
  w = ee[1,1] - ee[2,2] - ee[3,3] - ee[4,4]
  x = ee[1,2] + ee[2,1] + ee[3,4] - ee[4,3]
  y = ee[1,3] - ee[2,4] + ee[3,1] + ee[4,2]
  z = ee[1,4] + ee[2,3] - ee[3,2] + ee[4,1]
  if (w < 0.0)
    return Quaternion(-w, -[x; y; z])
  end
  return Quaternion(w, [x; y; z])
end

#mangled type products, return first type or nearest Group type (doesn't return an Algebra)
*(a::SO3, bq::Quaternion) = a*convert(SO3,bq)
*(a::SO3, b::so3) = a*convert(SO3,b)
*(a::SO3, baa::AngleAxis) = a*convert(SO3,baa)

*(aq::Quaternion, b::SO3) = aq*convert(Quaternion,b)
*(aq::Quaternion, b::so3) = aq*convert(Quaternion, convert(SO3,b) )
*(aq::Quaternion, baa::AngleAxis) = aq*convert(Quaternion, baa)

*(a::AngleAxis, b::AngleAxis) = convert(AngleAxis,convert(Quaternion,a)*convert(Quaternion,b))
*(a::AngleAxis, bq::Quaternion) = convert(AngleAxis, convert(Quaternion,a)*bq)
*(a::AngleAxis, b::SO3) = convert(AngleAxis, convert(Quaternion,a)*b )
*(aa::AngleAxis, b::so3) = aa*convert(SO3,b)

*(a::so3,b::so3) = convert(SO3,a)*convert(SO3,b)
*(a::so3, b::SO3) = convert(SO3,a)*b
*(a::so3, bq::Quaternion) = convert(Quaternion, convert(SO3,a) )*bq
*(a::so3, b::AngleAxis) = convert(AngleAxis, convert(SO3,a))*b

inverse(a::SE3) = SE3( matrix(a) \ Matrix{Float64}(LinearAlgebra.I,4,4) )
# TODO -- optimize this, and abstraction is wrong here
# Xj = Xi ⊕ ΔX
# Xj ⊖ ΔX = Xi
# Xi \ Xj = ΔX   # ⊖ Xi ⊕ Xj = ΔX
A_invB(a::SE3, b::SE3) = SE3( collect( ( matrix(b)' \ (matrix(a)') )' )  )
ominus(xi::SE3, xj::SE3) = A_invB(xi,xj)
oplus(xi::SE3, xj::SE3) = xi*xj
⊖(xi::SE3, xj::SE3) = A_invB(xi,xj)
⊕(xi::SE3, Δx::SE3) = xi*Δx
# SE3(0) ⊖ xi ⊕ xj = Δx
\(xi::SE3, xj::SE3) = SE3( matrix(xi) \ matrix(xj) )

\(qi::Quaternion,qj::Quaternion) = q_conj(qi) * qj

function doinversetransform!(dst::Vector{Float64}, aTb::SE3, src::Vector{Float64})::Nothing
  mul!(dst, transpose(aTb.R.R), src)  # At_mul_B!(dst, aTb.R.R, src)
  @fastmath @inbounds for i in 1:3, j in 1:3
    dst[i] -= aTb.R.R[j,i]*aTb.t[i]
  end
  nothing
end
# bR

# comparison functions

compare(a::SO3, b::SO3; tol::Float64=1e-14) = norm((a.R*transpose(b.R))-Matrix{Float64}(LinearAlgebra.I,3,3)) < tol

# function compare(a::SE3, b::SE3; tol::Float64=1e-14)
#   norm(a.t-b.t) < tol ? nothing : return false
#   return compare(a.R,b.R, tol=tol)
# end
function compare(a::SE3,b::SE3; tol::Float64=1e-10)
  norm(a.t-b.t)<tol && TransformUtils.compare(a.R, b.R, tol=tol)
end
function compare(a::Quaternion, b::Quaternion; tol::Float64=1e-14)
  qiq = a*q_conj(b)
  return tol <= qiq.s <= 1.0+tol && norm(qiq.v) < tol
end
function compare(a::AngleAxis,b::AngleAxis; tol::Float64=1e-14)
  aTb = q_conj(convert(Quaternion,a))*b
  return compare(Quaternion(0),aTb, tol=tol)
end


# convert functions

function convert(::Type{Quaternion}, v::VectorFloatInt)
    return Quaternion(v[1],v[2:4])
end

function convert(::Type{Quaternion}, aa::AngleAxis)
    v = normalize(aa.ax)
    theta = aa.theta/2
    w = cos(theta)
    x = v[1] * sin(theta)
    y = v[2] * sin(theta)
    z = v[3] * sin(theta)
    return Quaternion(w, [x, y, z])
end

function convert(::Type{AngleAxis}, q::Quaternion)
    theta = acos(q.s) * 2.0
    if norm(q.v)==0
      return AngleAxis(theta, [1.,0,0])
    end
    return AngleAxis(theta, normalize(q.v))
end

function convert!(R::SO3, q::Quaternion)
  # q.s = q.s;
  # x = q.v[1];
  # y = q.v[2];
  # z = q.v[3];
  @inbounds begin
    if q.s < 0.0
      q.s = -q.s
      q.v[1:3] = -q.v[1:3]
    end
    # nrm = sqrt(q.s^2+sum(q.v[1:3].^2))
    # if (nrm < 0.999)
    #   println("q2C -- not a unit quaternion nrm = $(nrm)")
    #   R = Matrix{Float64}(LinearAlgebra.I, 3,3)
    # else
      normalize!(q)
      # nrm = 1.0/nrm
      # q.s *= nrm
      # q.v .*= nrm
      # # x = x*nrm
      # # y = y*nrm
      # # z = z*nrm
      w2, x2, y2, z2 = q.s*q.s, q.v[1]*q.v[1], q.v[2]*q.v[2], q.v[3]*q.v[3]
      # y2 = y*y
      # z2 = z*z
      # w2 =
      xy = 2.0*q.v[1]*q.v[2]
      xz = 2.0*q.v[1]*q.v[3]
      yz = 2.0*q.v[2]*q.v[3]
      wx = 2.0*q.s*q.v[1]
      wy = 2.0*q.s*q.v[2]
      wz = 2.0*q.s*q.v[3]
      #build matrix
      # R = zeros(3,3)
      R.R[1,1] = w2+x2-y2-z2
      R.R[1,2] = xy-wz
      R.R[1,3] = xz+wy
      R.R[2,1] = xy+wz
      R.R[2,2] = w2-x2+y2-z2
      R.R[2,3] = yz-wx
      R.R[3,1] = xz-wy
      R.R[3,2] = yz+wx
      R.R[3,3] = w2-x2-y2+z2
    # end
  end
  # return SO3(R)
  nothing
end
function convert(::Type{SO3}, q::Quaternion)
  R = SO3(0)
  convert!(R,q)
  return R
end

function convert(::Type{Quaternion}, S::SO3)
  @inbounds begin
    rot = S.R
    tr = rot[1,1] + rot[2,2] + rot[3,3];

    if (tr > 0.0)
      S = sqrt(tr+1.0) * 2.0; # S=4*qw
      qw = 0.25 * S;
      qx = (rot[3,2] - rot[2,3]) / S;
      qy = (rot[1,3] - rot[3,1]) / S;
      qz = (rot[2,1] - rot[1,2]) / S;
    else
      if ((rot[1,1] > rot[2,2]) && (rot[1,1] > rot[3,3]))
        S = sqrt(1.0 + rot[1,1] - rot[2,2] - rot[3,3]) * 2.0; # S=4*qx
        qw = (rot[3,2] - rot[2,3]) / S;
        qx = 0.25 * S;
        qy = (rot[1,2] + rot[2,1]) / S;
        qz = (rot[1,3] + rot[3,1]) / S;
      else
        if (rot[2,2] > rot[3,3])
          S = sqrt(1.0 + rot[2,2] - rot[1,1] - rot[3,3]) * 2.0 # S=4*qy
          qw = (rot[1,3] - rot[3,1]) / S;
          qx = (rot[1,2] + rot[2,1]) / S;
          qy = 0.25 * S;
          qz = (rot[2,3] + rot[3,2]) / S;
        else
          S = sqrt(1.0 + rot[3,3] - rot[1,1] - rot[2,2]) * 2.0 # S=4*qz
          qw = (rot[2,1] - rot[1,2]) / S
          qx = (rot[1,3] + rot[3,1]) / S
          qy = (rot[2,3] + rot[3,2]) / S
          qz = 0.25 * S;
        end
      end
    end
  end

  if (qw >= 0)
    q = Quaternion(qw,[qx;qy;qz])
  else
    q = Quaternion(-qw,-[qx;qy;qz])
  end
  return q
end

function convert(::Type{SO3}, aa::AngleAxis)
  convert(SO3,convert(Quaternion, aa))
end

function convert(::Type{AngleAxis}, s::SO3)
  convert(AngleAxis, convert(Quaternion, s))
end


function convert!(Eu::Euler, q::Quaternion)
  #Function to convert quaternion to Euler angles, from Titterton and Weston
  #04, following quaternion convention of Titterton

  # a = q.s
  # b = q.v[1]
  # c = q.v[2]
  # d = q.v[3]

  Eu.R = atan(2.0*(q.v[2]*q.v[3] + q.v[1]*q.s), q.s*q.s - q.v[1]^2 - q.v[2]^2 + q.v[3]^2)
  # Eu.R = -atan(2.0*(q.v[2]*q.v[3] - q.v[1]*q.s), q.s*q.s - q.v[1]^2 - q.v[2]^2 + q.v[3]^2)
  #-atan(2.0*(q.v[2]*q.v[3] - q.v[1]*q.s),1.0-2.0*(q.v[1]^2+q.v[2]^2)); # -atan(2.0*(q.v[2]*q.v[3] - q.v[1]*q.s), q.s*q.s - q.v[1]^2 - q.v[2]^2 + q.v[3]^2) #numerically more stable
  Eu.P = asin(2.0*(q.v[2]*q.s - q.v[1]*q.v[3]));
  Eu.Y =  atan(2.0*(q.v[1]*q.v[2] + q.v[3]*q.s), q.s^2 + q.v[1]^2 - q.v[2]^2 - q.v[3]^2)
  # Eu.Y = -atan(2.0*(q.v[1]*q.v[2] - q.v[3]*q.s),1.0-2.0*(q.v[2]^2+q.v[3]^2)); # -atan(2.0*(q.v[1]*q.v[2] - q.v[3]*q.s), q.s^2 + q.v[1]^2 - q.v[2]^2 - q.v[3]^2)

  nothing
end
function convert(::Type{Euler}, q::Quaternion)
  Eu = Euler(0)
  convert!(Eu, q)
  return Eu
end

function convert(::Type{so3}, G::SO3)
  so3(logmap(G))
end

function convert(::Type{so3}, q::Quaternion)
  so3(logmap(convert(SO3,q)))
end

function convert(::Type{Quaternion}, alg::so3)
  convert(Quaternion, convert(SO3, alg))
end

function convert(::Type{Euler}, R::SO3)
  convert(Euler, convert(Quaternion, R))
end

function convert!(q::Quaternion, E::Euler)::Nothing

  @fastmath halfroll = 0.5*E.R;
  @fastmath halfpitch = 0.5*E.P;
  @fastmath halfyaw = 0.5*E.Y;
  @fastmath sin_r2 = sin(halfroll);
  @fastmath sin_p2 = sin(halfpitch);
  @fastmath sin_y2 = sin(halfyaw);
  @fastmath cos_r2 = cos(halfroll);
  @fastmath cos_p2 = cos(halfpitch);
  @fastmath cos_y2 = cos(halfyaw);

  @fastmath q.s = cos_r2 * cos_p2 * cos_y2 + sin_r2 * sin_p2 * sin_y2;
  @fastmath @inbounds q.v[1] = sin_r2 * cos_p2 * cos_y2 - cos_r2 * sin_p2 * sin_y2;
  @fastmath @inbounds q.v[2] = cos_r2 * sin_p2 * cos_y2 + sin_r2 * cos_p2 * sin_y2;
  @fastmath @inbounds q.v[3] = cos_r2 * cos_p2 * sin_y2 - sin_r2 * sin_p2 * cos_y2;

  # Enforce positive scalar quaternions following conversion from Euler angles
  @inbounds if (q.s<0.0)
    q.s = -q.s
    q.v[1:3] = -q.v[1:3];
  end
  # q = q./norm(q);
  normalize!(q)
  # return Quaternion(q[1],q[2:4])
  nothing
end
  # mutable struct HalfAngles
  #   halfang::Vector{Float64}
  #   halfsin::Vector{Float64}
  #   halfcos::Vector{Float64}
  #   HalfAngles() = new()
  #   HalfAngles(::Int) = new(zeros(3), zeros(3), zeros(3))
  # end
  # unsafe_reuseE2Q::Vector{HalfAngles}
  # thr_reuse = E.unsafe_reuseE2Q[Threads.threadid()]
  # thr_reuse.halfang[1] = 0.5*E.R;
  # thr_reuse.halfang[2] = 0.5*E.P;
  # thr_reuse.halfang[3] = 0.5*E.Y;
  # thr_reuse.halfang[1] = 0.5*E.R;
  # thr_reuse.halfang[2] = 0.5*E.P;
  # thr_reuse.halfang[3] = 0.5*E.Y;
  # thr_reuse.halfsin[1] = sin(thr_reuse.halfang[1]);
  # thr_reuse.halfsin[2] = sin(thr_reuse.halfang[2]);
  # thr_reuse.halfsin[3] = sin(thr_reuse.halfang[3]);
  # thr_reuse.halfcos[1] = cos(thr_reuse.halfang[1]);
  # thr_reuse.halfcos[2] = cos(thr_reuse.halfang[2]);
  # thr_reuse.halfcos[3] = cos(thr_reuse.halfang[3]);
  # @fastmath q.s = thr_reuse.halfcos[1] * thr_reuse.halfcos[2] * thr_reuse.halfcos[3] + thr_reuse.halfsin[1] * thr_reuse.halfsin[2] * thr_reuse.halfsin[3];
  # @fastmath @inbounds q.v[1] = thr_reuse.halfsin[1] * thr_reuse.halfcos[2] * thr_reuse.halfcos[3] - thr_reuse.halfcos[1] * thr_reuse.halfsin[2] * thr_reuse.halfsin[3];
  # @fastmath @inbounds q.v[2] = thr_reuse.halfcos[1] * thr_reuse.halfsin[2] * thr_reuse.halfcos[3] + thr_reuse.halfsin[1] * thr_reuse.halfcos[2] * thr_reuse.halfsin[3];
  # @fastmath @inbounds q.v[3] = thr_reuse.halfcos[1] * thr_reuse.halfcos[2] * thr_reuse.halfsin[3] - thr_reuse.halfsin[1] * thr_reuse.halfsin[2] * thr_reuse.halfcos[3];
  # @fastmath q.s = thr_reuse.halfcos[1] * thr_reuse.halfcos[2] * thr_reuse.halfcos[3]
  # @fastmath @inbounds q.v[1] = thr_reuse.halfsin[1] * thr_reuse.halfcos[2] * thr_reuse.halfcos[3]
  # @fastmath @inbounds q.v[2] = thr_reuse.halfcos[1] * thr_reuse.halfsin[2] * thr_reuse.halfcos[3]
  # @fastmath @inbounds q.v[3] = thr_reuse.halfcos[1] * thr_reuse.halfcos[2] * thr_reuse.halfsin[3]
  # @fastmath q.s += thr_reuse.halfsin[1] * thr_reuse.halfsin[2] * thr_reuse.halfsin[3];
  # @fastmath @inbounds q.v[1] -= thr_reuse.halfcos[1] * thr_reuse.halfsin[2] * thr_reuse.halfsin[3];
  # @fastmath @inbounds q.v[2] += thr_reuse.halfsin[1] * thr_reuse.halfcos[2] * thr_reuse.halfsin[3];
  # @fastmath @inbounds q.v[3] -= thr_reuse.halfsin[1] * thr_reuse.halfsin[2] * thr_reuse.halfcos[3];


function convert(::Type{Quaternion}, E::Euler)
  q = Quaternion(0)
  convert!(q, E)
  return q
end

function convert(::Type{SO3}, E::Euler)
  return convert(SO3,convert(Quaternion, E))
end

function convert!(R::SO3, E::Euler)
  # TODO -- refactor to inplace operations
  convert!(E.fastconvert, E)
  convert!(R, E.fastconvert)
  nothing
end

# type fastconvert!
#
# end

function rotate!(q1::Quaternion, v1::Array{Float64,1})
    #v = (q1*Quaternion(0.0,v1)*q_conj(q1)).v
    R = convert(SO3, q1);
    v = R.R * v1
    for i = 1:3 v1[i] = v[i] end
    nothing
end
function rotate(q1::Quaternion, v1::Array{Float64,1})
    vv1 = deepcopy(v1)
    rotate!(q1, vv1)
    return vv1
end


function expmOwn(alg::so3)
  v_norm = sqrt(alg.S[1,2]^2 + alg.S[1,3]^2 + alg.S[2,3]^2)
  I = Matrix{Float64}(LinearAlgebra.I, 3,3)
  if (v_norm>1e-6)
    #1E-6 is chosen because double numerical LSB is around 1E-18 for nominal values [-pi..pi] and (1E-7)^2 is at 1E-14, but 1E-7 rad/s is 0.02 deg/h
    R = I + sin(v_norm)/v_norm*alg.S + (1.0-cos(v_norm))/(v_norm^2)*(alg.S^2)
  else
    R = I + alg.S + 0.5*(alg.S^2)
  end
  return R
end

function logmap(grp::SO3)
  R = grp.R
  tra = LinearAlgebra.tr(R)
  if abs(tra+1.0) < 1e-10
    if abs(R[3,3]+1.0) > 1e-10
      ret = (pi / sqrt(2.0+2.0*R[3,3] )) * [R[1,3]; R[2,3]; (1.0+R[3,3])]
    else
      if abs(grp.R[1,1]+1.0) > 1e-10
        ret = (pi / sqrt(2.0+2.0*R[2,2])) * [R[1,2]; (1.0+R[2,2]); R[3,2]]
      else
        ret = (pi / sqrt(2.0+2.0*R[1,1])) * [(1.0+R[1,1]); R[2,1]; R[3,1]]
      end
    end
  else
    tra = round(tra, digits=12)
    tr_3 = tra-3.0
    # TODO: this sign might have been the wrong way round, changed on 10/10/2018 for Taylor series when tr_3 small
    ## also introduced abs(tr_3) rather than tr_3 < 1e-7, which looks to be a bug
    if (abs(tr_3) < -1e-7)
      magnitude = 0.5 - tr_3*tr_3/12.0
    else
      theta = acos((tra-1.0)/2.0)
      magnitude = theta/(2.0*sin(theta)) # 0.5/sinc(theta)
    end
    ret = magnitude*[ (R[3,2]-R[2,3]); (R[1,3]-R[3,1]); (R[2,1]-R[1,2])]
  end
  return ret
end

function logmap(q::Quaternion)
  logmap(convert(SO3, q))
end

# Chikjain's book p??
function rightJacExmap(alg::so3)
  Gam = alg.S
  v_norm = sqrt(alg.S[1,2]^2 + alg.S[1,3]^2 + alg.S[2,3]^2)
  I = Matrix{Float64}(LinearAlgebra.I, 3,3)
  if (v_norm>1e-7)
    Jr = I + (v_norm-sin(v_norm))/(v_norm^3)*(Gam^2) - (1.0-cos(v_norm))/(v_norm^2+0.0)*Gam
  else
    Jr = I
  end
  return Jr
end

function rightJacExmapinv(alg::so3)
  Gam = alg.S
  v_norm = sqrt(alg.S[1,2]^2 + alg.S[1,3]^2 + alg.S[2,3]^2)
  I = Matrix{Float64}(LinearAlgebra.I, 3,3)
  if (v_norm>1e-7)
    brace = (  1.0/(v_norm^2) - (1.0+cos(v_norm))/(2.0*v_norm*sin(v_norm))  )*Gam
    Jr_inv = I + 0.5*Gam + brace*Gam
  else
    Jr_inv = I
  end
  return Jr_inv
end

# function computeRatesFromRPY(Q, prevQ, dt)
#   if (dt<1):
#     # Assume oneside rectangular derivative, we could use a trapezoidal
#     # approach to obtain a more accurate estimate of the rotation rate
#     # TODO: Improve derivative scheme
#     S = (    np.dot(q2R(Q),np.transpose(q2R(prevQ)))  - np.identity(3)    ) / dt;
#     bw = np.matrix([[-S[2,1].item()],[-S[0,2].item()],[-S[1,0].item()]])
#   else
#     bw = zeros(3,1);
#   end
#   return bw
# end

# See Julia's implementation the matrix exponential function -- expm!()
# https://github.com/acroy/julia/blob/0bce8951f19e8f47d361e73b5b5bd6283926c01b/base/linalg/dense.jl
expmOwnT(M::Array{Float64,2}) = (Matrix{Float64}(LinearAlgebra.I,size(M,1),size(M,1)) + 0.5*M)*((Matrix{Float64}(LinearAlgebra.I,size(M,1),size(M,1)) - 0.5*M) \ Matrix{Float64}(LinearAlgebra.I,size(M,1),size(M,1)));
expmOwn1(M::Array{Float64,2}) = Matrix{Float64}(LinearAlgebra.I,size(M,1),size(M,1)) + M;
expmOwn2(M::Array{Float64,2}) = Matrix{Float64}(LinearAlgebra.I,size(M,1),size(M,1)) + M + 0.5*(M^2);
function expmOwn3(M::Array{Float64,2})
  M052 = 0.5*M^2;
  return Matrix{Float64}(LinearAlgebra.I,size(M,1),size(M,1)) + M + M052 + M052*M/3.0;
end
function expmOwn4(M::Array{Float64,2})
  M052 = 0.5*M^2;
  return Matrix{Float64}(LinearAlgebra.I,size(M,1),size(M,1)) + M + M052 + M052*M/3.0 + M052*M052/6.0
end

function convert(::Type{SO3}, alg::so3)
  v = vee(alg)
  nv = norm(v)
  if nv < 1e-3
    return SO3(exp(alg.S))
  else
    invnv = 1.0/nv
    return SO3(Matrix{Float64}(LinearAlgebra.I, 3,3) + invnv*(sin(nv)*alg.S + invnv*(1.0-cos(nv))*(alg.S^2) ) )
  end
end



function wrapRad(th::Float64)
  th = rem(th, 2pi) # returns in range (-2pi,2pi)
  # now move to range [-pi, pi)
  if th >= pi
    th -= 2.0*pi
  elseif th < -pi
    th += 2.0*pi
  end
  return th
end
# rem(pi, 2pi)
# rem(4pi, 2pi)
# rem(4.5pi, 2pi)
# rem(-4.5pi, 2pi)
# rem(7pi, 2pi)
# rem(-7pi, 2pi)
# rem(-pi, 2pi)

R(th::Float64) = [[cos(th);-sin(th)]';[sin(th);cos(th)]'];
R(;x::Float64=0.0,y::Float64=0.0,z::Float64=0.0) = convert(SO3, so3([x,y,z]))

function vee!(rv::Vector{Float64}, alg::so3)
  rv[1] = -alg.S[2,3]
  rv[2] = alg.S[1,3]
  rv[3] = -alg.S[1,2]
  nothing
end

function vee(alg::so3)
   rv = zeros(3)
   @inbounds vee!(rv, alg)#[-alg.S[2,3], alg.S[1,3], -alg.S[1,2]]
   return rv
end

function vee!(rv::Vector{Float64}, q::Quaternion)
  rv[1] = q.s
  rv[2:4] = q.v[1:3]
  nothing
end
function vee(q::Quaternion)
   rv = zeros(4)
   @inbounds vee!(rv, q)
   return rv
end

# vectorize SE3 group to translation and AngleAxis numbers
# dim7 out: [xyz, theta, axyz]
function veeAngleAxis(G::SE3)
  aa = convert(AngleAxis, G.R)
  v = zeros(7)
  v[1:3] = G.t
  v[4] = aa.theta
  v[5:7] = aa.ax[:]
  return v
end

function vec(q::Quaternion)
  v = zeros(4)
  v[1] = q.s
  v[2:4] = q.v[:]
  return v
end

# vectorize SE3 group to translation and Quaternion numbers
# dim7 out: [xyz, cos(th/2), sin(th/2)xyz]
function veeQuaternion(G::SE3)
  q = convert(Quaternion, G.R)
  v = zeros(7)
  v[1:3] = G.t
  v[4:7] = vee(q)
  # v[4] = q.s
  # v[5:7] = q.v[:]
  return v
end

# vectorize SE3 group to translation and Quaternion numbers
# dim7 out: [xyz, cos(th/2), sin(th/2)xyz]
function veeEuler(G::SE3)
  E = convert(Euler, G.R)
  v = zeros(6)
  v[1:3] = G.t
  v[4] = E.R
  v[5] = E.P
  v[6] = E.Y

  return v
end

function deltaso3vee(aQb::Quaternion,aQc::Quaternion)
  dq = aQb\aQc
  s = logmap(dq)
end



# TODO -- Change to type and overload the operators
# TODO -- uncomment
# function vee!(rv::Vector{Float64,1},T::SE2)
#   rv[1] = T.t[1]
#   rv[2] = T.t[2]
#   rv[3] = wrapRad(atan(-T.R.R[1,2], T.R.R[1,1]))
#   nothing
# end
# function *(a::SE2, b::SE2)
#   return SE2(R.R*b.R, vec(a.R.R*b.t + a.t))
# end

function SE2(X::Array{Float64,1})
    T = Matrix{Float64}(LinearAlgebra.I, 3,3)
    T[1:2,1:2] = R(X[3])
    T[1,3] = X[1]
    T[2,3] = X[2]
    return T
end

function se2vee!(retval::Array{Float64,1}, T::Array{Float64,2})
    retval[1] = T[1,3]
    retval[2] = T[2,3]
    retval[3] = wrapRad(atan(-T[1,2], T[1,1]))
    nothing
end

function se2vee(T::Array{Float64,2})
    retval = zeros(3)
    se2vee!(retval, T)
    return retval
end


# TODO Switch to using SE(2) oplus
# DX = [transx, transy, theta]
function addPose2Pose2!(retval::Array{Float64,1}, x::Array{Float64,1}, dx::Array{Float64,1})
  X = SE2(x)
  DX = SE2(dx)
  se2vee!(retval, X*DX)
  nothing
end
function addPose2Pose2(x::Array{Float64,1}, dx::Array{Float64,1})
    retval = zeros(3)
    addPose2Pose2!(retval, x, dx)
    return retval
end



end #module
