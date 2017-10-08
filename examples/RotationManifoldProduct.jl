# weighted product of gaussian terms on rotational manifold.

using TransformUtils


# means = zeros(4,2)
# means[1,:] = 1.0

sigmas = rand(1,2)
sigmas = ones(1,2)

q1 = Quaternion(0)
q2 = convert(Quaternion, so3(randn(3)))


function weightedmean(
            qArrMu::Vector{Quaternion};
            covariances::Vector{Array{Float64,2}}=Vector{Array{Float64,2}}(0),
            tol::Float64=1e-7,
            allowediters::Int=30  )
  #
  covariances = length(covariances) == length(qArrMu) ? covariances : [eye(3) for i in 1:length(qArrMu)]
  dq = deepcopy(qArr[1])
  S = zeros(3)
  i = 0
  epsi = 1.0
  while epsi > tol && i < allowediters
    i+=1
    S[:] = 0.0
    for q in qArr
      S[:] += covariances[i] \ deltaso3vee(dq, q)
    end
  end

  return 
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

Vector{Array{Float64,2}}(0)

#
