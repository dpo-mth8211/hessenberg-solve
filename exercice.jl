using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
  x = similar(b)
  ### votre code ici ; ne rien modifier d'autre
  n = length(b)

  # Premier pivot
  x[end] = b[end]/R[n,n]

  # Solution pour le reste. Reviens à faire b - produit scalaire tronqué R*x
  for i = n-1:-1:1
    x[i] = (b[i]-R[i,i+1:n]'*x[i+1:n])/R[i,i]
  end
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

# Fonction pour avoir les valeurs de sinus/cosinus des rotations de Givens
function givens_homemade(a, b)
  r = sqrt(a^2 + b^2)
  c = a / r
  s = b / r
  return c, s
end

function hessenberg_solve(H::UpperHessenberg, b)
  ### votre code ici ; ne rien modifier d'autre

  # Nombre de colonnes
  _,n = size(H)

  # Initialisation de la matrice Q
  Qtot = Matrix{Float64}(I,n,n)

  # Boucle pour la triangularisation de H.
  for i = 1:n-1
      Q_i = Matrix{Float64}(I,n,n)
      c,s = givens_homemade(H[i,i],H[i+1,i])
      Q_i[i,i:i+1] = [c s]
      Q_i[i+1,i:i+1] = [-s c]
      H = Q_i*H
      Qtot = Q_i*Qtot
  end
  # Valeur de R & mettre le bon type
  R = UpperTriangular(H)

  # Résolution Système Rx = c
  c = Qtot*b
  x = backsolve(R, c)
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
