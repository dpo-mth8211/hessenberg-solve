using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
  x = similar(b)  # création d'un vecteur de la meme taille que b
  m=size(R,1)     # nombre de ligne et de colonne de R pour Un R carree
 for i=m:-1:1     # loop permettant de navigué les lignes de R commençant de la derniere ligne à la première 
  x[i] = (b[i] - dot(R[i, i+1:end], x[i+1:end])) / R[i, i]   # Bacskolve implementation 
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
  m=size(H,1)
  for i=1:m-1    
      j=i+1
    ########Calcul des parametres de la matrice de rotation############
      r=norm(H[:,i])
      c=H[i,i]/r
      s=H[j,i]/r
    ########Appliquer les rotations de givens tout en s'assurant de modifier seulement les parametres non nuls 
       for k = i:m 
        Hik = c * H[i, k] + s * H[j, k]  
        Hjk = -s * H[i, k] + c * H[j, k]
        H[i, k] = Hik
        H[j, k] = Hjk
      end
    ##########Modification de la matrice B############################################
    bi = c * b[i] + s * b[j]
    bj = -s * b[i] + c * b[j]
    b[i] = bi
    b[j] = bj
  end 
  ############Compatibilité avec la fonction backsolve##############
    R=UpperTriangular(H)
  ############Résolution en utilisant la fonction backsolve#########
  x = backsolve(R, b)
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
