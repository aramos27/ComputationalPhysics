using GLMakie
using LinearAlgebra

#=
We want to solve Δϕ = 0 by considering the matrix formulation
We use finite differences to rewrite the equation as:
ϕ_(i+1,j) + ϕ_(i-1,j) + ϕ_(i,j+1) + ϕ_(i,j-1) - 4ϕ_(i,j) = 0 

For full derivation see 3.30 in Skript
=#


function matrix_multiplication(x::Vector{Float64}, N)
    x_grid = reshape(x, N-2, N-2)' #Row major
    result = zeros(Float64, N-2, N-2)
    for i in 1:(N-2)
        for j in 1:(N-2)
            if i == 1 #First row
                if j == 1 #First column
                    result[i,j] = x_grid[i+1, j] + x_grid[i,j+1]- 4 * x_grid[i, j]
                elseif j == (N-2) #Last column
                    result[i,j] = x_grid[i+1, j] + x_grid[i,j-1]- 4 * x_grid[i, j]
                else #In between
                    result[i,j] = x_grid[i+1, j] + x_grid[i, j+1] + x_grid[i, j-1] - 4 * x_grid[i, j]
                end
            elseif i == (N-2) #Last Row
                if j == 1 #First column
                    result[i,j] = x_grid[i-1, j] + x_grid[i,j+1]- 4 * x_grid[i, j]
                elseif j == (N-2)#Last column
                    result[i,j] = x_grid[i-1, j] + x_grid[i,j-1]- 4 * x_grid[i, j]
                else #In between
                    result[i,j] = x_grid[i-1, j] + x_grid[i, j+1] + x_grid[i, j-1] - 4 * x_grid[i, j]
                end
            elseif j == 1 #First column in between first and last row
                result[i,j] = x_grid[i-1, j] + x_grid[i,j+1] + x_grid[i+1,j]- 4 * x_grid[i, j]
            elseif j == (N-2) #Last column in between first and last row
                result[i,j] = x_grid[i-1, j] + x_grid[i,j-1] + x_grid[i+1,j]- 4 * x_grid[i, j]
            else
                result[i,j] = x_grid[i+1, j] + x_grid[i-1, j] + x_grid[i, j+1] + x_grid[i, j-1] - 4 * x_grid[i, j]
            end
        end
    end
    
    return vec(permutedims(result))
end

function initialize_b(N)
    x = range(-π/2, π/2; length=N)
    y = range(-π/2, π/2; length=N)

    ϕ = zeros(N, N)

    #Initialize "arches on all edges"
    for j in 1:N
        ϕ[1, j] = cos(x[j])
        ϕ[end, j] = cos(x[j])
    end
    for i in 1:N
        ϕ[i, 1] = cos(y[i])
        ϕ[i, end] = cos(y[i])
    end


    b = zeros(Float64,(N - 2)^2)
    b_grid = reshape(b, N-2, N-2)' #Row major
    for i in 1:(N-2)
        for j in 1:(N-2)
            if i == 1 #First row
                if j == 1 #First column
                    b_grid[i,j] = -ϕ[1,2] - ϕ[2,1]
                elseif j == (N-2) #Last column
                    b_grid[i,j] = -ϕ[1,N-1] - ϕ[2,N]
                else #In between
                    b_grid[i,j] = -ϕ[1,j+1]
                end
            elseif i == (N-2) #Last Row
                if j == 1 #First column
                    b_grid[i,j] = -ϕ[N-1,1] - ϕ[N,2]
                elseif j == (N-2)#Last column
                    b_grid[i,j] = -ϕ[N-1, N] - ϕ[N, N-1]
                else #In between
                    b_grid[i,j] = -ϕ[N,j+1]
                end
            elseif j == 1 #First column in between first and last row
                b_grid[i,j] = -ϕ[i+1,1]
            elseif j == (N-2) #Last column in between first and last row
                b_grid[i,j] = -ϕ[i+1,N]
            else
                b_grid[i,j] = 0
            end

        end
    end


    return vec(permutedims(b_grid))
end

function conjugate_gradient(A_mul, b, N; tol=1e-5, maxiter=10000)
    x = zeros(length(b))
    r = b - A_mul(x, N)
    g = copy(r)

    for iter in 1:maxiter
        g_g_k = dot(g,g)
        Ar = A_mul(r, N)
        λ_k =  g_g_k / dot(r, Ar)
        x += λ_k * r
        g -= λ_k * Ar
        g_g_k1 = dot(g,g)
        r = g + (g_g_k1/g_g_k) * r

        if sqrt(g_g_k1) < tol
            println("Converged after $iter iterations with residual $(sqrt(g_g_k1))")
            break
        end
    end

    return x
end

for N in 11:10:51
    b = initialize_b(N)


    x = range(-π/2, π/2; length=N)
    y = range(-π/2, π/2; length=N)

    # Reconstruct the full solution grid including boundary conditions
    ϕ_full = zeros(N, N)

    for j in 1:N
    ϕ_full[1, j] = cos(x[j])
    ϕ_full[end, j] = cos(x[j])
    end
    for i in 1:N
        ϕ_full[i, 1] = cos(y[i])
        ϕ_full[i, end] = cos(y[i])
    end

    solution_vec = conjugate_gradient(matrix_multiplication, b, N)
    solution_grid = reshape(solution_vec, N-2, N-2)

    # Insert the computed solution into the interior
    ϕ_full[2:N-1, 2:N-1] .= solution_grid

    # Now plot the full grid
    fig = Figure(size = (800, 600))
    ax = Axis3(fig[1, 1], title = "Conjugate Gradient", xlabel = "x", ylabel = "y", zlabel = "ϕ")
    surface!(ax, x, y, ϕ_full; colormap = :viridis)
    save("Surface$(N).png", fig)
end




