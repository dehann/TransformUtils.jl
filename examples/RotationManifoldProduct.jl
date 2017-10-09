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
            covariances::Vector{Array{Float64,2}}=Vector{Array{Float64,2}}(),
            tol::Float64=1e-7,
            allowediters::Int=30  )
  #
  covariances = length(covariances) == length(qArrMu) ? covariances : [eye(3) for i in 1:length(qArrMu)]
  muq = deepcopy(qArrMu[1])
  S = zeros(3)
  invC = zeros(3,3)
  epsi = 1.0
  i = 0
  while epsi > tol && i < allowediters
    i+=1
    j=0
    S[:] = 0.0
    invC[:] = 0.0
    for q in qArrMu
      j += 1
      invCov = covariances[j] \ eye(3)
      invC += invCov
      S[:] += invC * deltaso3vee(muq, q)
    end
    dq = invC \ S
    newmuq = muq * so3(dq)
    epsi = abs(muq.s - newmuq.s) + norm(muq.v - newmuq.v)
    muq = newmuq
  end

  return muq
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

q2

weightedmean([q1,q2])


#
