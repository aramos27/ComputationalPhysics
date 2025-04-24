using LinearAlgebra
using GLMakie
using JLD2
using Colors  
GLMakie.activate!(; float = true, focus_on_show = true)

const kg = 0.2 #Gravitational Interaction
const km = 20.0 #Magnetic Interaction
const γ = 0.5 #Friction
const R = 1.0 #Circle radius
const h = 0.2 #Separation between planes
const τ = 0.001  # time step 
const max_steps = 10^10


# Magnetic pole positions (placed on a circle of radius R)
const magnets = [
    R * [cos(π/2), sin(π/2)],
    R * [cos(7π/6), sin(7π/6)],
    R * [cos(11π/6), sin(11π/6)]
]

#Compute acceleration based on equation (1)
function compute_acceleration(r::Vector{Float64}, v::Vector{Float64})
    a = -kg * r - γ * v
    for r_i in magnets
        d = r - r_i
        a -= km * d / ((dot(d, d) + h^2)^(3/2))
    end
    return a
end

#Compute next step for r and v. 
function rk4_step(r::Vector{Float64}, v::Vector{Float64})
    k1r = τ * v
    k1v = τ * compute_acceleration(r, v)

    k2r = τ * (v + 0.5 * k1v)
    k2v = τ * compute_acceleration(r + 0.5 * k1r, v + 0.5 * k1v)

    k3r = τ * (v + 0.5 * k2v)
    k3v = τ * compute_acceleration(r + 0.5 * k2r, v + 0.5 * k2v)

    k4r = τ * (v + k3v)
    k4v = compute_acceleration(r + k3r, v + k3v)

    r_next = r + (1/6.0)*(k1r + 2*k2r + 2*k3r + k4r)
    v_next = v + (1/6.0)*(k1v + 2*k2v + 2*k3v + k4v)

    return r_next, v_next
end

#Simulates pendulum starting from an initial position assuming initial velocity is always 0
function simulate_pendulum(r0::Vector{Float64})
    r = r0
    v = [0.0, 0.0]

    for step in 1:max_steps
        r_new, v_new = rk4_step(r, v)

        if norm(v_new) < 1e-6 && norm(r_new - r) < 1e-6 #Early termination condition
            return r_new  # final resting position
        end

        r, v = r_new, v_new
    end
    return r  # if it didn't settle, return current position
end

#Finds nearest magnet at endpoint and assings value: {1,2,3} according to which magnet it was.
function classify_endpoint(r::Vector{Float64})
    min_dist = 10
    endpoint_index = 0
    for i in 1:3
        distance = norm(r - magnets[i])
        if distance <= min_dist
            min_dist = distance
            endpoint_index = i 
        end
    end
    if min_dist >0.05
        println("error")
    end
    return endpoint_index
end

#Simulate pendulum with different starting positions on a 2D grid. 
function grid_simulation(gridLength::Int)
    grid = range(-2, 2, gridLength)
    classification = zeros(Int, gridLength,gridLength)
    total = gridLength^2
    count = 0
    for (i, x0) in enumerate(grid)
        for (j, y0) in enumerate(grid)
            r_final = simulate_pendulum([x0,y0])
            classification[i,j] = classify_endpoint(r_final)
            count += 1
            if count % 1000 == 0
                println("Progress: $((count/total)*100) %")
            end
        end
    end
    return classification
end

#Plots results with three different colors. 
function plot_classification_grid(classification::Array{Int,2})
    
    value_to_color = Dict(
    1 => parse(Colorant, "#FF6F00"),  # Vibrant metallic orange (amber/orange anodized aluminum)
    2 => parse(Colorant, "#3B9C9C"),  # Deep metallic teal
    3 => parse(Colorant, "#8A2BE2")   # Electric indigo (bold and complementary)
)

    color_data = [value_to_color[x] for x in classification]
    color_matrix = reshape(color_data, size(classification))

    # Create corresponding x and y axes based on grid range
    gridLength = size(classification, 1)
    x = range(-2, 2, length=gridLength)
    y = range(-2, 2, length=gridLength)

    f = Figure()
    ax = Axis(f[1, 1], title="Classification Grid", xlabel="x", ylabel="y")
    heatmap!(ax, x, y, color_matrix)

    return f
end

#Plots trajectory from a given initial position and initial velocity. 
function simulate_and_plot(initial_r::Vector{Float64}, initial_v::Vector{Float64})
    r = copy(initial_r)
    v = copy(initial_v)
    trajectory = [r]

    for i in 1:max_steps
        r_new, v_new = rk4_step(r, v)
        push!(trajectory, r_new)

        # Termination condition: velocity and position change are both very small
        if norm(v_new) < 1e-6 && norm(r_new - r) < 1e-6
            println("Terminated early at step $i")
            break
        end

        r, v = r_new, v_new
    end

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="x", ylabel="y", title="Magnetic Pendulum Trajectory")
    limits!(ax, -2, 2, -2, 2)
    xs = [p[1] for p in trajectory]
    ys = [p[2] for p in trajectory]
    lines!(ax, xs, ys, linewidth=2)

    # Plot magnets
    for m in magnets
        scatter!(ax, [m[1]], [m[2]], color=:red, marker=:circle, markersize=10)
    end

    return fig
end

#Animates trajectory from a given initial position and inital velocity. 
function animate_trajectory(initial_r::Vector{Float64}, initial_v::Vector{Float64}, filename::String = "pendulum_animation.mp4")
    r = copy(initial_r)
    v = copy(initial_v)
    xs = Float64[]
    ys = Float64[]

    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "Magnetic Pendulum Trajectory")
    limits!(ax, -2, 2, -2, 2)

    lineplot = lines!(ax, xs, ys, linewidth = 2)

    # Plot magnets
    for m in magnets
        scatter!(ax, [m[1]], [m[2]], color = :red, marker = :circle, markersize = 12)
    end

    # Dot for current position
    dot = scatter!(ax, [r[1]], [r[2]], color = :black, markersize = 10)

    # Setup video recorder manually
    Makie.record(fig, filename; framerate = 60) do io
        for frame in 1:max_steps
            r_new, v_new = rk4_step(r, v)
            push!(xs, r_new[1])
            push!(ys, r_new[2])

            lineplot[1] = xs
            lineplot[2] = ys

            dot[1] = [r_new[1]]
            dot[2] = [r_new[2]]

            recordframe!(io)  # Save this frame

            if norm(v_new) < 1e-6 && norm(r_new - r) < 1e-6
                println("Terminated early at frame $frame")
                break  
            end

            r, v = r_new, v_new
        end
    end
end

using Base.Threads
using Base.Threads: Atomic, atomic_add!

function grid_simulation_multithread(gridLength::Int)
    grid = range(-2, 2, gridLength)
    classification = Matrix{Int}(undef, gridLength, gridLength)

    total = gridLength^2
    counter = Atomic{Int}(0)

    @threads for idx in 1:total
        i = fld(idx - 1, gridLength) + 1  # row (y)
        j = mod(idx - 1, gridLength) + 1  # column (x)

        r0 = [grid[j], grid[i]]
        r_final = simulate_pendulum(r0)
        classification[j, i] = classify_endpoint(r_final)

        # Update progress
        current = atomic_add!(counter, 1)
        if current % (total ÷ 100) == 0
            println("Progress: $(round(current / total * 100; digits=1))%")
        end
    end

    return classification
end

function main()
    discretization_size = 10000
    classification = grid_simulation_multithread(discretization_size)
    
    @save "classification_$discretization_size.jld2" classification
    fig = plot_classification_grid(classification)
    save("Fractal_$discretization_size.png", fig)
end

main()
#initial_r = [2.0, 2.0]
#initial_v = [0.0,0.0]
#fig = simulate_and_plot(initial_r, initial_v)
#save("Plots.png", fig)

#animate_trajectory(initial_r, initial_v, "pendulum_animation1.mp4")






