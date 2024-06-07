using LinearAlgebra

function backsolve(R::UpperTriangular, b)
    x = copy(b)
    m,n = size(R)

    # x contains the right-hand side values of the equation Rx = b that
    # we calculate iteratively.
    # At first it holds only the value of b, but each time we calculate
    # the value of a new variable x[j] (starting with x[n]), we shift
    # the value of that variable to the right-hand side.
    # After that, all that remains on the right-hand side is R[j,j] * x[j]
    # so we need to divide the right-hand side by R[j,j]
    for j = n:-1:1
        x[j] /= R[j,j]
        for i = j-1:-1:1
            x[i] -= (R[i,j] * x[j])
         end
    end    
    return x
end

function hessenberg_solve(H::UpperHessenberg, b)
    x = similar(b)
    m,n = size(H)
    for i = 1:n-1
        ρ = H[i,i]^2 + H[i+1,i]^2
        c = H[i,i] / ρ
        s = H[i+1,i] / ρ
        givens_rotation_mat = [c s; -s c]
        
        # Apply givens rotation in-place on the line of each element on the
        # diagonal together with the element directly below the diagonal.
        H[i:i+1,i:end] .= givens_rotation_mat * H[i:i+1,i:end]
        b[i:i+1] .= givens_rotation_mat * b[i:i+1]
    end

    # After applying Givens rotations on each element below the diagonal
    # H is upper triangular. We can solve with a simple backsolve.
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

