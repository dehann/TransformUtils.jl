__precompile__(true)

module TransformUtils

using LinearAlgebra
using ManifoldsBase
using Manifolds
using StaticArrays
using DocStringExtensions

import Base: convert, promote_rule, *, \, vec
import LinearAlgebra: transpose, adjoint, normalize, normalize!
import ManifoldsBase: vee, vee!, hat, hat!

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
  logmap_SO2,
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
  se2vee!,
  TU,
  veePose3,
  veePose

const TU = TransformUtils
const LA = LinearAlgebra


include("LazyLieManifolds.jl")
include("Legacy.jl")
include("LegacySpecialEuclidean.jl")




end #module
