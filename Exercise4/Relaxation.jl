using GLMakie
using LinearAlgebra

function initialize_phi(Nx, Ny)
    #Divide x and y into Nx (Ny) intervalls
    x = range(-π/2, π/2; length=Nx)
    y = range(-π/2, π/2; length=Ny)
    ϕ = zeros(Nx, Ny)

    #Initialize "arches on all edges"
    for j in 1:Ny
        ϕ[1, j] = cos(y[j])
        ϕ[end, j] = cos(y[j])
    end
    for i in 1:Nx
        ϕ[i, 1] = cos(x[i])
        ϕ[i, end] = cos(x[i])
    end

    return ϕ, x, y
end

function gauss_seidel!(ϕ; tol=1e-5)
    Nx, Ny = size(ϕ)
    iter = 0
    δϕ = Inf
    error_history = Float64[]
    while δϕ > tol
        δϕ = 0.0
        for i in 2:Nx-1, j in 2:Ny-1
            old = ϕ[i, j]
            ϕ[i, j] = 0.25 * (ϕ[i+1, j] + ϕ[i-1, j] + ϕ[i, j+1] + ϕ[i, j-1])
            δϕ = max(δϕ, abs(ϕ[i, j] - old))
        end
        push!(error_history, δϕ)
        iter += 1
    end
    return iter, error_history
end

function sor!(ϕ; tol=1e-5)
    Nx, Ny = size(ϕ)
    ω_opt = 2 / (1 + sqrt(1 - (0.5*(cos(π/Nx)+cos(π/Ny)))^2))
    iter = 0
    δϕ = Inf
    error_history = Float64[]
    while δϕ > tol
        δϕ = 0.0
        for i in 2:Nx-1, j in 2:Ny-1
            old = ϕ[i, j]
            res = 0.25 * (ϕ[i+1, j] + ϕ[i-1, j] + ϕ[i, j+1] + ϕ[i, j-1])
            ϕ[i, j] = (1 - ω_opt) * old + ω_opt * res
            δϕ = max(δϕ, abs(ϕ[i, j] - old))
        end
        push!(error_history, δϕ)
        iter += 1
    end
    return iter, error_history
end

function w_n(n, w_prev, Nx, Ny)
    r = (0.5*(cos(π/Nx)+cos(π/Ny)))
    if isapprox(n, 0.0; atol=1e-5)
        return 1
    elseif isapprox(n, 0.5; atol=1e-5)
        return 1/(1-(0.5*(r^2)))
    else
        return 1/(1-(0.25*(r^2)*w_prev))
    end
end

function sor_chebyshev!(ϕ; tol=1e-5)
    Nx, Ny = size(ϕ)
    iter = 0
    w_prev = 0
    n = 0.0
    start = 0
    δϕ = Inf
    error_history = Float64[]
    while δϕ > tol
        δϕ = 0.0
        start = 2
        w = w_n(n, w_prev,Nx, Ny)
        for i in start:2:Nx-1
            for j in 2:Ny-1
                old = ϕ[i, j]
                res = 0.25 * (ϕ[i+1, j] + ϕ[i-1, j] + ϕ[i, j+1] + ϕ[i, j-1])
                ϕ[i, j] = (1 - w) * old + w * res
                δϕ = max(δϕ, abs(ϕ[i, j] - old))
            end
            if start == 2
                start = 3
            else
                start = 2
            end
        end

        n += 0.5
        w_prev = w
        start = 3
        w = w_n(n, w_prev,Nx, Ny)
        for i in start:2:Nx-1
            for j in 2:Ny-1
                old = ϕ[i, j]
                res = 0.25 * (ϕ[i+1, j] + ϕ[i-1, j] + ϕ[i, j+1] + ϕ[i, j-1])
                ϕ[i, j] = (1 - w) * old + w * res
                δϕ = max(δϕ, abs(ϕ[i, j] - old))
            end
            if start == 2
                start = 3
            else
                start = 2
            end
        end

        n += 0.5
        w_prev = w

        push!(error_history, δϕ)
        iter += 1
    end
    return iter, error_history
end

function solve_plot_surface(N)
    ϕ1, x, y = initialize_phi(N, N)
    ϕ2 = deepcopy(ϕ1)
    ϕ3 = deepcopy(ϕ1)
    t1, error_history1 = gauss_seidel!(ϕ1)

    t2, error_history2 = sor_chebyshev!(ϕ2)

    t3, error_history3 = sor!(ϕ3)

    gauss_surface = Figure(size = (800, 600))
    ax = Axis3(gauss_surface[1, 1], title = "Gauss-Seidel result (N=$N)", xlabel = "x", ylabel = "y", zlabel = "ϕ")
    surface!(ax, x, y, ϕ1; colormap = :viridis)

    sor_chebyshev_surface = Figure(size = (800, 600))
    ax = Axis3(sor_chebyshev_surface[1, 1], title = "SOR-Chebyshev result (N=$N)", xlabel = "x", ylabel = "y", zlabel = "ϕ")
    surface!(ax, x, y, ϕ2; colormap = :viridis)

    sor_surface = Figure(size = (800, 600))
    ax = Axis3(sor_surface[1, 1], title = "SOR result (N=$N)", xlabel = "x", ylabel = "y", zlabel = "ϕ")
    surface!(ax, x, y, ϕ3; colormap = :viridis)
    return gauss_surface, sor_chebyshev_surface, sor_surface
end

function solve_plot_error(N)
    ϕ1, x, y = initialize_phi(N, N)
    ϕ2 = deepcopy(ϕ1)
    ϕ3 = deepcopy(ϕ1)

    its1, error_history1 = gauss_seidel!(ϕ1)

    its2, error_history2 = sor_chebyshev!(ϕ2)

    its3, error_history3 = sor!(ϕ3)

    error_fig = Figure(size = (800, 600))
    ax3 = Axis(error_fig[1, 1], xlabel = "Iteration", ylabel = "Error", title = "Error Convergence (N=$N)")
    lines!(ax3, 1:its1, error_history1, label = "Gauss-Seidel Error")
    lines!(ax3, 1:its2, error_history2, label = "SOR-Chebyshev Error")
    lines!(ax3, 1:its3, error_history3, label = "SOR Error")


    axislegend(ax3)
    return error_fig
end

function solve_plot_iterations(ns)
    N_values = Int32[]
    its_values_gauss = Int32[]
    its_values_sor_chebyshev = Int32[]
    its_values_sor = Int32[]

    for N in ns

        ϕ1, x, y = initialize_phi(N, N)
        ϕ2 = deepcopy(ϕ1)
        ϕ3 = deepcopy(ϕ1)

        its1, error_history1 = gauss_seidel!(ϕ1)
        its2, error_history2 = sor_chebyshev!(ϕ2)
        its3, error_history3 = sor!(ϕ3)
        push!(its_values_gauss, its1)
        push!(its_values_sor_chebyshev, its2)
        push!(its_values_sor, its3)
        push!(N_values, N)
    end

    iteration_fig = Figure(size = (800, 600))
    ax3 = Axis(iteration_fig[1, 1], xlabel = "N", ylabel = "Iterations", title = "Iterations to Converge")
    lines!(ax3, N_values, its_values_gauss, label = "Gauss-Seidel")
    lines!(ax3, N_values, its_values_sor_chebyshev, label = "SOR_Chebyshev")
    lines!(ax3, N_values, its_values_sor, label = "SOR")


    axislegend(ax3)
    save("IterationPlots/IterationPlot.png", iteration_fig)
    return iteration_fig
end

function solve_plots(ns)
    if !isdir("GaussSeidelPlots")
        mkdir("GaussSeidelPlots")
    end
    if !isdir("SORPlots")
        mkdir("SORPlots")
    end
    if !isdir("SOR_Chebyshev_Plots")
        mkdir("SOR_Chebyshev_Plots")
    end
    if !isdir("ErrorPlots")
        mkdir("ErrorPlots")
    end
    if !isdir("IterationPlots")
        mkdir("IterationPlots")
    end
    
    for N in ns
        fig_gs, fig_sor_cheb, fig_sor = solve_plot_surface(N)
        fig_er = solve_plot_error(N)
        
        gs_filename = "GaussSeidelPlots/GaussSeidel_surface_N$(N).png"
        sor_chebyshev_filename = "SOR_Chebyshev_Plots/SOR_Chebyshev_N$(N).png"
        sor_filename = "SORPlots/SOR_surface_N$(N).png"
        er_filename = "ErrorPlots/Error_N$(N).png"
        
        save(gs_filename, fig_gs)
        save(sor_chebyshev_filename, fig_sor_cheb)
        save(sor_filename, fig_sor)
        save(er_filename, fig_er)
    end
    iterations_fig = solve_plot_iterations(ns)
end

solve_plots(11:10:51)
