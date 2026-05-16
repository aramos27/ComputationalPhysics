using GLMakie
GLMakie.activate!(; float = true, focus_on_show = true)

const accuracy = 10^(-5)

function chebyshev_accel(half_step, r, ω)
    if half_step == 1
        return 1 / (1 - (1/2) * (r^2))
    end
    return 1 / (1 - (1/4) * (r^2) * ω) 
end

function relaxation(N; sor = false)
    ϕ = zeros((N, N))
    for i in 1:N
        x = cos(pi * (i/N) - pi/2)
        ϕ[1,i] = x
        ϕ[i,1] = x
        ϕ[N,i] = x
        ϕ[i,N] = x
    end

    step = 0
    errors = zeros(0)

    ω = 1
    chebyshev_accel_radius = cos(pi / (N - 1))
    while true
        step += 1
        δ = 0
        for y in 2:N-1
            for x in 2:N-1
                if sor && (x + y) % 2 == step % 2
                    continue
                end
                old = ϕ[x,y]
                ϕ[x,y] = ω * (ϕ[x-1,y] + ϕ[x+1,y] + ϕ[x,y-1] + ϕ[x,y+1])/4 + (1-ω) * old
                δ = max(δ, abs(ϕ[x,y] - old))
            end
        end
        if sor
            ω = chebyshev_accel(step, chebyshev_accel_radius, ω)
        end
        append!(errors, δ)
        if δ < accuracy
            break
        end
    end
    println("Finished $N*$N grid in $step steps, sor: $sor")
    fig = Figure()
    ax = Axis(fig[1, 1])
    heatmap!(ax, ϕ)
    ax = Axis(fig[2, 1], xlabel = "Steps", ylabel = "δΦ")
    lines!(ax, errors)
    display(fig)
    wait(fig.scene)
end

for N in 11:10:51
    relaxation(N)
end

for N in 11:10:51
    relaxation(N, sor = true)
end
