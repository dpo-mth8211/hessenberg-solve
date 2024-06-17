using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
  ### votre code ici ; ne rien modifier d'autre
  # Pour faire du code efficace on fait une copie de b dans x
  n = length(b)
  x = copy(b)
  for i in n:-1:1
      @assert R[i, i] != 0
      for j=i+1:n
        x[i] -= R[i, j] * x[j]
      end
      x[i] /= R[i, i]
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
  n = size(H, 1)
  for i in 1:n-1
      a = H[i, i]
      b_val = H[i+1, i]
      r = sqrt(a^2 + b_val^2)
      # Definir les coefficients de la rotation de Givens
      c = a / r
      s = b_val / r

      # Appliquer les rotations de Givens à H
      for j in i:n
          t1 = H[i, j]
          t2 = H[i+1, j]
          H[i, j] = c * t1 + s * t2
          H[i+1, j] = -s * t1 + c * t2
      end

      # Appliquer les rotations de Givens à b aussi
      t1 = b[i]
      t2 = b[i+1]
      b[i] = c * t1 + s * t2
      b[i+1] = -s * t1 + c * t2
  end

  # Maintenant on peut résoudre le système triangulaire supérieur
  x = backsolve(UpperTriangular(H), b)
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
