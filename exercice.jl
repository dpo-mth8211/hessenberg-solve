using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
  #x = similar(b)
  ### votre code ici ; ne rien modifier d'autre
  # ...
  ###
  n = size(R, 1)
  x = zeros(n)

  x[n] = b[n] / R[n, n]

  for i in n:-1:1
    x[i] = (b[i] - dot(R[i, i+1:end], x[i+1:end])) / R[i, i]
  end

  return x
end


# 2. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système Hessenberg supérieur Hx = b à l'aide de rotations ou
#    réflexions de Givens et d'une remontée triangulaire.
#    Votre fonction peut modifier H et b si nécessaire.
#    Il n'est pas nécessaire de garder les rotations en mémoire et la
#    fonction ne doit pas les renvoyer.
#    Seul le cas réel sera testé ; pas le cas complexe.
function hessenberg_solve(H::UpperHessenberg, b)
  ### votre code ici ; ne rien modifier d'autre
  # ...
  # x = ...
  ###
  #x = similar(b) 

  x = b
  m = size(H, 1)
  n = size(b, 2)
  u = Vector{typeof(zero(eltype(H)))}(undef, m) # for last rotated col of H
  copyto!(u, 1, H, m * (m - 1) + 1, m) # u .= H[:,m]

  cs = Vector{Tuple{real(eltype(u)),eltype(u)}}(undef, length(u)) # store Givens rotations

  for k = m:-1:2
    # Givens Algorithm
    f1 = u[k]
    g1 = H[k, k-1]
    rho = sqrt(f1*f1 + g1*g1)

    c = f1/rho
    s = g1/rho

    cs[k] = (c, s)

    for i = 1:n
      x[k, i] /= rho
      t_1 = s * x[k, i]
      t_2 = c * x[k, i]

      for j = 1:k-2
        x[j, i] -= u[j] * t_2 + H[j, k-1] * t_1
      end
      x[k-1, i] -= u[k-1] * t_2 + (H[k-1, k-1]) * t_1
    end

    for j = 1:k-2
      u[j] = H[j, k-1] * c - u[j] * s'
    end

    u[k-1] = (H[k-1, k-1]) * c - u[k-1] * s'
  end

  for i = 1:n
    tau_1 = x[1, i] / u[1]

    for j = 2:m
      tau_2 = x[j, i]
      c, s = cs[j]
      x[j-1, i] = c * tau_1 + s * tau_2
      tau_1 = c * tau_2 - s'tau_1
    end

    x[m, i] = tau_1
  end

  return x
end


# vérification
using Test
for n ∈ (10, 20, 30)
  A = rand(n, n)
  b = rand(n)
  R = UpperTriangular(A)
  x = backsolve(R, b)
  @test norm(R * x - b) ≤ sqrt(eps()) * norm(b)
  H = UpperHessenberg(A)
  x = hessenberg_solve(copy(H), copy(b))
  @test norm(H * x - b) ≤ sqrt(eps()) * norm(b)
end
