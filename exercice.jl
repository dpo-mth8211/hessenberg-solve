using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
  x = similar(b)
  ### votre code ici ; ne rien modifier d'autre
  # ...
  ###
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
