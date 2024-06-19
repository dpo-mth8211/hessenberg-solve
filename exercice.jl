using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
  x = similar(b)
  ### votre code ici ; ne rien modifier d'autre
  n = length(b)
# Pour chacune des colonnes de R en partant de la dernière:
  for i in n:-1:1
    total = 0.0
# Pour chacune des valeurs sur les lignes (Donc, on fait un produit (lignes R) * (x))
    for j in i+1:n  
      total += R[i,j] * x[j]
    end
# Calcul des valeurs de x[i] tel que R[i,i] * x[i] + total = b[i] ---> x[i] = (b[i] - total)/ R[i,j]
    x[i] = (b[i] - total) / R[i,i] 
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
  n = size(H,1)
  # Accès à chacune des valeurs de la matrice A:
    for i in 1:n-1
      for j = i+1:n
  # Pour chacun des coefficients de la matrice A non-nul: 
  # Calcul des coefficients c,s,p du système 2x2 avec c = x/p et s = y/p (Note de cours diapo 32/58)
        if H[j,i] != 0
          p = sqrt( H[i,j]^2 + H[j,i]^2 )
          c = H[i,i]/p
          s = H[j,i]/p
  # Application des rotations de Givens sur H:
          for k in i:n
            val = H[i,k]
            H[i,k] = c * H[i,k] - s * H[j,k]
            H[j,k] = s * val    + c * H[j,k]
          end
  # Même application mais sur le vecteur b:
          val = b[i]
          b[i] = c * b[i] - s * b[j]
          b[j] = s * val  + c * b[j]
        end
      end
    end
  # Remontée triangulaire en suivant le même principe que pour la fonction backsolve:
    x = similar(b)
    for i = n:-1:1
      x[i] = b[i]
      for j = i+1:n
        x[n] -= H[i,j] *x[j]
      end
      x[i]/= H[i,i]
    end
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
