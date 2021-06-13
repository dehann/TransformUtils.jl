# legacy TU SpecialEuclidean






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

function SE2(X::AbstractArray{P,1}) where P <: Real
    T = Matrix{P}(LinearAlgebra.I, 3,3)
    T[1:2,1:2] = R(X[3])
    T[1,3] = X[1]
    T[2,3] = X[2]
    return T
end

function se2vee!(retval::AbstractArray{<:Real,1}, T::AbstractArray{<:Real,2})
    retval[1] = T[1,3]
    retval[2] = T[2,3]
    retval[3] = wrapRad(atan(-T[1,2], T[1,1]))
    nothing
end

function se2vee(T::AbstractArray{<:Real,2})
    # retval = zeros(3)
    # se2vee!(retval, T)
    # return retval
    return [T[1,3], T[2,3], wrapRad(atan(-T[1,2], T[1,1]))]
end



function veePose3(s::SE3)
  TransformUtils.veeEuler(s)
end
function veePose(s::SE3)
  TransformUtils.veeEuler(s)
end