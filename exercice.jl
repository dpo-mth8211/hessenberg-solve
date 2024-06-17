using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
  x = similar(b)
  ### votre code ici ; ne rien modifier d'autre
  # ...
  ###
  n = length(b)
  for i in n:-1:1
    x[i] = (b[i] - sum(R[i, j] * x[j] for j in i+1:n; init=0.0)) / R[i, i]
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
  n = size(H, 1)
  for j in 1:n-1
    for i in j+1:n
      if H[i, j] != 0.0
        r = sqrt(H[j, j]^2 + H[i, j]^2)
        c = H[j, j] / r
        s = -H[i, j] / r

        for k in j:n
          temp = c * H[j, k] - s * H[i, k]
          H[i, k] = s * H[j, k] + c * H[i, k]
          H[j, k] = temp
        end

        temp = c * b[j] - s * b[i]
        b[i] = s * b[j] + c * b[i]
        b[j] = temp
      end
    end
  end

  x = similar(b)
  for i in n:-1:1
    x[i] = (b[i] - sum(H[i, j] * x[j] for j in i+1:n; init=0.0)) / H[i, i]
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
