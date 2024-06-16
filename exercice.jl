using LinearAlgebra

# 1. Modifiez la fonction suivante pour qu'elle renvoie la solution x
#    du système triangulaire supérieur Rx = b.
#    Votre fonction ne doit modifier ni R ni b.
function backsolve(R::UpperTriangular, b)
    x = similar(b)
    n = length(b)
    for i in n:-1:1
        x[i] = (b[i] - sum(R[i, j] * x[j] for j in i+1:n)) / R[i, i]
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
    n = size(H, 1)
    H = copy(H)
    b = copy(b)
    
    for j in 1:n-1
        for i in n:-1:j+1
            a = H[i-1, j]
            b_elem = H[i, j]
            r = hypot(a, b_elem)
            if r != 0
                c = a / r
                s = -b_elem / r

                # Appliquer la rotation de Givens aux lignes i-1 et i de H
                G = [c -s; s c]
                H[i-1:i, j:end] = G * H[i-1:i, j:end]

                # Appliquer la rotation de Givens aux lignes i-1 et i de b
                b[i-1:i] = G * b[i-1:i]
            end
        end
    end
    
    # Maintenant que H est triangulaire supérieure, on utilise la substitution arrière
    x = backsolve(UpperTriangular(H), b)
    
    return x
end

# Fonction auxiliaire pour calculer l'hypoténuse
function hypot(a, b)
    return sqrt(a^2 + b^2)
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
