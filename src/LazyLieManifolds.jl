
export se, SE


# Tr0, R0 = SA[0.0; 0], SA[1.0 0; 0 1]
# p0 = ProductRepr(Tr0, R0)

# t0, w0 = [0.0; 0],  get_vector(M_R, R0, 0, B)
# x0 = ProductRepr(t0, w0)
# p0_ = exp(G, p0, x0)  # or,  exp(M, p0, x0)

# t1a, w1a = [1.0; 0],  get_vector(M_R, R0, 0  , B)
# x1a = ProductRepr(t1a, w1a)   # Lie algebra element
# p1a = exp(G, p0, x1a)         # Lie group element


"""
    $TYPEDEF

Lie Group representation for `LieGroups.SpecialEuclideanGroup(N)`

Related

[`se`](@ref)
"""
struct SE{N, T}
  p::T # Group
end


function SE(R::Matrix{T}, t::Vector{T}) where T<:Real
    N = length(t)
    @assert size(R,1) == size(R,2) == N

    _t = SVector{N}(t)
    _R = SMatrix{N,N}(R)

    p = ArrayPartition(_t, _R)
    
    return SE{N, typeof(P)}(p)
end

manifold(::SE{N,T}) where {N,T} = SpecialEuclideanGroup(N)


function _identity_SE(N::Int, ::Type{T}=Float64) where {T<:Real}
  p1 = zeros(SVector{N, T})
  p2 = SMatrix{N,N,T}(one(T)I)
  o = ArrayPartition(p1, p2)
  SE{N, typeof(o)}(o)
end


"""
    $TYPEDEF

Lie algebra representation for `LieGroups.SpecialEuclideanGroup(N)`

Related

[`SE`](@ref)
"""
struct se{N, T}
  X::T # Tangent Vector (Lie algebra, skew simetric matrix)
end



function se(X::Vector)
    N = length(X)
    
    for D in 2:3
        if manifold_dimension(SpecialEuclideanGroup(D)) == N
            M = SpecialEuclideanGroup(D)
            _X = hat(LieAlgebra(M), X)
            return se{D, typeof(_X)}(_X)
        end
    end
    error("Only se(2) and se(3) suppored")
end

manifold(::se{N,T}) where {N,T} = SpecialEuclideanGroup(N)

# p0 = SE([1.0 0.0; 0 1.0], [0.0, 0.0])

function ManifoldsBase.vee(X::se) 
    M = manifold(X)
    return vee(LieAlgebra(M), X.X)
end 

# we must use
# seVec = Manifolds.hat(M, coords)

function ⊕(p::SE, X)
    M = manifold(p)
    typeof(p)(exp(M, p.p, X))
end

⊕(X, p::SE) =  typeof(p)(compose(manifold(p), exp(manifold(p), X), p.p))

# I_SO2 = ProductRepr(SA[0.0, 0.0], SA[1.0 0.0; 0 1.0]);


# M = SpecialEuclidean(2)
# X1 = hat(M, I_SO2, SA[1., 0, pi/2]);
# X2 = hat(M, I_SO2, SA[1., 0, -pi/2]);

# p1 = p0 ⊕ X1

# p2 = p1 ⊕ X

##

#