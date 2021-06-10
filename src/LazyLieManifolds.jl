
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

Lie Group representation for `Manifolds.SpecialEuclidean(N)`

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

    p = ProductRepr(_t, _R)
    
    return SE{N, typeof(P)}(p)
end

manifold(::SE{N,T}) where {N,T} = SpecialEuclidean(N)


function _identity_SE(N::Int, ::Type{T}=Float64) where {T<:Real}
  p1 = zeros(SVector{N, T})
  p2 = SMatrix{N,N,T}(one(T)I)
  o = ProductRepr(p1, p2)
  SE{N, typeof(o)}(o)
end


"""
    $TYPEDEF

Lie algebra representation for `Manifolds.SpecialEuclidean(N)`

Related

[`SE`](@ref)
"""
struct se{N, T}
  X::T # Tangent Vector (Lie algebra, skew simetric matrix)
end



function se(X::Vector)
    N = length(X)
    
    for D in 2:3
        dim = manifold_dimension(SpecialEuclidean(D))
        if dim == N
            M = SpecialEuclidean(D)
            e = _identity_SE(dim, eltype(X)).p
            _X = hat(M, e, X)
        
            return se{D, typeof(_X)}(_X)
        end
    end
    error("Only se(2) and se(3) suppored")
end

manifold(::se{N,T}) where {N,T} = SpecialEuclidean(N)

# p0 = SE([1.0 0.0; 0 1.0], [0.0, 0.0])

function ManifoldsBase.vee(X::se) 
    M = manifold(se)
    return vee(M, Identity(M, X.X), X.X)
end 

# we must use
# seVec = Manifolds.hat(M, coords)

function ⊕(p::SE, X::ProductRepr)
    M = manifold(p)
    e = identity(M, p.p)
    typeof(p)(compose(M, p.p, exp(M, e, X)))
end

⊕(X::ProductRepr, p::SE) =  typeof(p)(compose(manifold(p), exp(manifold(p), I_SO2, X), p.p))

# I_SO2 = ProductRepr(SA[0.0, 0.0], SA[1.0 0.0; 0 1.0]);


# M = SpecialEuclidean(2)
# X1 = hat(M, I_SO2, SA[1., 0, pi/2]);
# X2 = hat(M, I_SO2, SA[1., 0, -pi/2]);

# p1 = p0 ⊕ X1

# p2 = p1 ⊕ X

##

#