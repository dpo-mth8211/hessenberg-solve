using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
    x = copy(b)

    # m == n
    m,n = size(R)
    for j = n:-1:1
        x[j] /= R[j,j]
        for i = j-1:-1:1
            x[i] -= (R[i,j] * x[j])
         end
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
    x = similar(b)
    m,n = size(H)
    for i = 1:n-1
        ρ = H[i,i]^2 + H[i+1,i]^2
        c = H[i,i] / ρ
        s = H[i+1,i] / ρ
       
        givens_rotation_mat = [c s; -s c]
        H[i:i+1,i:end] .= givens_rotation_mat * H[i:i+1,i:end]
        b[i:i+1] .= givens_rotation_mat * b[i:i+1]
    end

    R = UpperTriangular(H)
    return backsolve(R,b)
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

