module TransformUtils

import Base: convert, promote_rule, *

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
  *,
  normalize!,
  normalize,
  q_conj,
  q_conj!,
  convert,
  rotate!,
  rotate,
  wrapRad,
  logmap,
  rightJacExmap,
  rightJacExmapinv,
  # Should be good to go
  veeAngleAxis,
  veeQuaternion,

  # basic Taylor exported to try against Pade version in expm
  expmOwn,
  expmOwn1,
  expmOwnT,
  expmOwn2,
  expmOwn3,
  expmOwn4,


  # needs refactoring
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

typealias FloatInt Union{Float64,Int}

type Quaternion
    s::Float64
    v::Array{Float64,1}
    Quaternion() = new()
    Quaternion(s::FloatInt) = new(1.0,zeros(3))
    Quaternion(s::FloatInt,v::Array{FloatInt,1}) = new(s,v)
end

type AngleAxis
    theta::Float64
    ax::Array{Float64,1}
    AngleAxis() = new()
    AngleAxis(s::FloatInt) = new(1.0,zeros(3))
    AngleAxis(s::FloatInt,v::Array{FloatInt,1}) = new(s,v)
end

typealias AxisAngle AngleAxis

type SO3
    R::Array{Float64,2}
    SO3() = new()
    SO3(dummy::FloatInt) = new(eye(3))
    SO3(r::Array{Float64,2}) = new(r)
end

type so3
    S::Array{Float64,2}
    so3() = new()
    so3(s::FloatInt) = new(zeros(3,3))
    so3(v::Vector{Float64}) = new(skew(v))
    so3(S::Array{Float64,2}) = new(S)
end

type Euler
    R::Float64
    P::Float64
    Y::Float64
    Euler() = new()
    Euler(s::FloatInt) = new(0.0,0.0,0.0)
    Euler(r::FloatInt,p::FloatInt,y::FloatInt) = new(r,p,y)
end


type SE3
  R::SO3
  t::Vector{Float64}
  SE3() = new()
  SE3(dummy::FloatInt) = new(SO3(0.0), zeros(3))
  SE3(r::SO3, t::Vector{FloatInt}) = new(r,t)
  SE3(v::Vector{FloatInt}, E::Euler) = new(v[1:3],convert(SO3,E))
  SE3(v::Vector{FloatInt}, aa::AngleAxis) = new(v[1:3],convert(SO3,aa))
  SE3(v::Vector{FloatInt}, q::Quaternion) = new(v[1:3],convert(SO3,q))
end

function normalize!(q::Quaternion, tol=0.00001)
    mag2 = sum(q.v.^2) + q.s^2
    if abs(mag2 - 1.0) > tol
        mag = sqrt(mag2)
        q.v = q.v ./ mag
        q.s = q.s / mag
    end
    nothing
end
function normalize(q::Quaternion, tol=0.00001)
  qq = deepcopy(q)
  normalize!(qq,tol)
  return qq
end

function normalize(v::Array{Float64,1})
  return v / norm(v)
end


function *(a::SO3, b::SO3)
  return SO3(a.R*b.R)
end

function *(a::SE3, b::SE3)
  return SE3(R.R*b.R, vec(a.R.R*b.t + a.t))
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


function convert(::Type{Quaternion}, v::Array{Float64,1})
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

function convert(::Type{SO3}, q::Quaternion)
    w = q.s;
    x = q.v[1];
    y = q.v[2];
    z = q.v[3];
    if w < 0.0
        w,x,y,z = -w, -x, -y, -z
    end
    nrm = sqrt(w*w+x*x+y*y+z*z)
    if (nrm < 0.999)
        println("q2C -- not a unit quaternion nrm = $(nrm)")
        R = eye(3)
    else
        nrm = 1.0/nrm
        w = w*nrm
        x = x*nrm
        y = y*nrm
        z = z*nrm
        x2 = x*x
        y2 = y*y
        z2 = z*z
        w2 = w*w
        xy = 2.0*x*y
        xz = 2.0*x*z
        yz = 2.0*y*z
        wx = 2.0*w*x
        wy = 2.0*w*y
        wz = 2.0*w*z
        #build matrix
        R = zeros(3,3)
        R[1,1] = w2+x2-y2-z2
        R[1,2] = xy-wz
        R[1,3] = xz+wy
        R[2,1] = xy+wz
        R[2,2] = w2-x2+y2-z2
        R[2,3] = yz-wx
        R[3,1] = xz-wy
        R[3,2] = yz+wx
        R[3,3] = w2-x2-y2+z2
    end
    return SO3(R)
end

function convert(::Type{Quaternion}, S::SO3)
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


function convert(::Type{Euler}, q::Quaternion)
  #Function to convert quaternion to Euler angles, from Titterton and Weston
  #04, following quaternion convention of Titterton
  #UJ UAV Research 11/05/09

  a = q.s
  b = q.v[1]
  c = q.v[2]
  d = q.v[3]

  phi = -atan2(2.0*(c*d - b*a),1.0-2.0*(b^2+c^2)); # -atan2(2.0*(c*d - b*a), a*a - b*b - c*c + d*d) #numerical more stable
  theta = -asin(-2.0*(c*a - b*d));
  psi = -atan2(2.0*(b*c - d*a),1.0-2.0*(c^2+d^2)); # -atan2(2.0*(b*c - d*a), a*a + b*b - c*c - d*d)

  return  Euler(phi,theta,psi)
end

function convert(::Type{so3}, G::SO3)
  so3(logmap(G))
end

function convert(::Type{Euler}, R::SO3)
  convert(Euler, convert(Quaternion, R))
end

function convert(::Type{Quaternion}, E::Euler)
  # Using fixed frame rotation scheme, as used in MIT libbot
  q = zeros(4)

  halfroll = E.R/2.0;
  halfpitch = E.P/2.0;
  halfyaw = E.Y/2.0;

  sin_r2 = sin (halfroll);
  sin_p2 = sin (halfpitch);
  sin_y2 = sin (halfyaw);

  cos_r2 = cos (halfroll);
  cos_p2 = cos (halfpitch);
  cos_y2 = cos (halfyaw);

  q[1] = cos_r2 * cos_p2 * cos_y2 + sin_r2 * sin_p2 * sin_y2;
  q[2] = sin_r2 * cos_p2 * cos_y2 - cos_r2 * sin_p2 * sin_y2;
  q[3] = cos_r2 * sin_p2 * cos_y2 + sin_r2 * cos_p2 * sin_y2;
  q[4] = cos_r2 * cos_p2 * sin_y2 - sin_r2 * sin_p2 * cos_y2;

  q = q./norm(q);

  # Enforce positive scalar quaternions following conversion from Euler angles
  if (q[1]<0)
    q = -q;
  end
  return Quaternion(q[1],q[2:4])
end


function rotate!(q1::Quaternion, v1::Array{Float64,1})
    #v = (q1*Quaternion(0.0,v1)*q_conj(q1)).v
    R = convert(SO3, q1);
    v = R.R * v1
    [v1[i] = v[i] for i = 1:3]
    nothing
end
function rotate(q1::Quaternion, v1::Array{Float64,1})
    vv1 = deepcopy(v1)
    rotate!(q1, vv1)
    return vv1
end


function expmOwn(alg::so3)
  v_norm = sqrt(alg.S[1,2]^2 + alg.S[1,3]^2 + alg.S[2,3]^2)
  I = eye(3)
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
  tr = trace(R)
  if abs(tr+1.0) < 1e-10
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
    tr_3 = tr-3.0
    if (tr_3<-1e-7)
      theta = acos((tr-1.0)/2.0)
      magnitude = theta/(2.0*sin(theta))
    else
      magnitude = 0.5 - tr_3*tr_3/12.0
    end
    ret = magnitude*[ (R[3,2]-R[2,3]); (R[1,3]-R[3,1]); (R[2,1]-R[1,2])]
  end
  return ret
end

# Chikjain's book p??
function rightJacExmap(alg::so3)
  Gam = alg.S
  v_norm = sqrt(alg.S[1,2]^2 + alg.S[1,3]^2 + alg.S[2,3]^2)
  I = eye(3)
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
  I = eye(3)
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
expmOwnT(M::Array{Float64,2}) = (eye(size(M,1)) + 0.5*M)*((eye(size(M,1)) - 0.5*M) \ eye(size(M,1)));
expmOwn1(M::Array{Float64,2}) = eye(size(M,1)) + M;
expmOwn2(M::Array{Float64,2}) = eye(size(M,1)) + M + 0.5*(M^2);
function expmOwn3(M::Array{Float64,2})
  M052 = 0.5*M^2;
  return eye(size(M,1)) + M + M052 + M052*M/3.0;
end
function expmOwn4(M::Array{Float64,2})
  M052 = 0.5*M^2;
  return eye(size(M,1)) + M + M052 + M052*M/3.0 + M052*M052/6.0
end

function convert(::Type{SO3}, alg::so3)
  return SO3(expm(alg.S))
end


function wrapRad(th::Float64)
  if th >= pi
    th -= 2.0*pi
  end
  if th < -pi
    th += 2.0*pi
  end
  return th
end

R(th::Float64) = [[cos(th);-sin(th)]';[sin(th);cos(th)]'];

function vee!(rv::Vector{Float64}, alg::so3)
  rv[1] = -alg.S[2,3]
  rv[2] = alg.S[1,3]
  rv[3] = -alg.S[1,2]
  nothing
end

function vee(alg::so3)
   rv = zeros(3)
   vee!(rv, alg)#[-alg.S[2,3], alg.S[1,3], -alg.S[1,2]]
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

# vectorize SE3 group to translation and Quaternion numbers
# dim7 out: [xyz, cos(th/2), sin(th/2)xyz]
function veeQuaternion(G::SE3)
  q = convert(Quaternion, G.R)
  v = zeros(7)
  v[1:3] = G.t
  v[4] = q.s
  v[5:7] = q.v[:]
  return v
end



# TODO -- Change to type and overload the operators
# TODO -- uncomment
# function vee!(rv::Vector{Float64,1},T::SE2)
#   rv[1] = T.t[1]
#   rv[2] = T.t[2]
#   rv[3] = wrapRad(atan2(-T.R.R[1,2], T.R.R[1,1]))
#   nothing
# end
# function *(a::SE2, b::SE2)
#   return SE2(R.R*b.R, vec(a.R.R*b.t + a.t))
# end

function SE2(X::Array{Float64,1})
    T = eye(3)
    T[1:2,1:2] = R(X[3])
    T[1,3] = X[1]
    T[2,3] = X[2]
    return T
end

function se2vee!(retval::Array{Float64,1}, T::Array{Float64,2})
    retval[1] = T[1,3]
    retval[2] = T[2,3]
    retval[3] = wrapRad(atan2(-T[1,2], T[1,1]))
    nothing
end

function se2vee(T::Array{Float64,2})
    retval = zeros(3)
    se2vee!(retval, T)
    return retval
end



end #module
