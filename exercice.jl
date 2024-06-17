using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.

function backsolve(R::UpperTriangular, b)
    x = similar(b)
    n = length(b)
    
    for i in n:-1:1
        x[i] = (b[i] - dot(R[i, i+1:end], x[i+1:end])) / R[i, i]
    end
    return x
end
###
#------------------------------------------------------------------------
#------------------------------------------------------------------------
# 2. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système Hessenberg supérieur Hx = b à l'aide de rotations ou
#    réflexions de Givens et d'une remontée triangulaire.
#    Votre fonction peut modifier H et b si nécessaire.
#    Il n'est pas nécessaire de garder les rotations en mémoire et la
#    fonction ne doit pas les renvoyer.
#    Seul le cas réel sera testé ; pas le cas complexe.

#function hessenberg_solve(H::UpperHessenberg{Float64, Matrix{Float64}}, b::Vector{Float64})
function hessenberg_solve(H::UpperHessenberg, b)
  n = size(H, 1)

    # Apply Givens rotations to transform H into an upper triangular matrix
    for i in 1:n-1
        r = sqrt.(H[i,i]^2 + H[i+1,i]^2)
        c, s = givens_rotation(H[i, i], H[i+1, i])

          for k in i:n
              temp_Hik = c * H[i, k] + s * H[i+1, k]
              H[i+1, k] = -s * H[i, k] + c * H[i+1, k]
              H[i, k] = temp_Hik
          end

          temp_bi = c * b[i] + s * b[i+1]
          b[i+1] = -s * b[i] + c * b[i+1]
          b[i] = temp_bi
    end

    # Now H is upper triangular, solve Hx = b using back substitution
    x = backsolve(UpperTriangular(H), b)
    return x
end
    # n = size(H,1)
    # for i in 1:n-1
    #     r = sqrt.(H[i,i]^2 + H[i+1,i]^2)
    #     c = H[i,i] / r
    #     s = H[i+1,i] / r

    #     H[i:i+1,i:end] .= [c s; -s c] * H[i:i+1,i:end]#appliquer Givens rotation sur H
    #     b[i:i+1] .= [c s; -s c] * b[i:i+1]#appliquer Givens rotation sur b

    # end
    # H est maintenat triangulaire superieure on peut appeler la fonction backsolve
    #x = backsolve(UpperTriangular(H), b)
#     x = similar(b)
#     x = backsolve(H, b)
#     return x
# end
###
#------------------------------------------------------------------------
function givens_rotation(a, b)
  r = hypot(a, b)
  c = a / r
  s = b / r
  return c, s
end
#------------------------------------------------------------------------
# vérification
using Test
for n ∈ (10, 20, 30)
#for n ∈ (10, 20, 30, 40, 50, 60)#enlever apres
  A = rand(n, n)
  b = rand(n)
  R = UpperTriangular(A)
  x = backsolve(R, b)
  @test norm(R * x - b) ≤ sqrt(eps()) * norm(b)
  H = UpperHessenberg(A)
  x = hessenberg_solve(copy(H), copy(b))
  @test norm(H * x - b) ≤ sqrt(eps()) * norm(b)

  if norm(H * x - b) ≤ sqrt(eps()) * norm(b)#enlever a la remise
    println("Test Passed")#enlever a ala remis
  end

end

