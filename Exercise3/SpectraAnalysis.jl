using FFTW
using GLMakie
GLMakie.activate!(; float = true, focus_on_show = true)

const γ = 0.5
const k = 1
const Ω = 2/3
const T_0 = (2*π)/Ω
const τ = T_0/200

const total_periods = 1000
const ignore_periods = 100
const steps = Int(total_periods * T_0 / τ)

function compute_acceleration(x::Float64, v::Float64, t::Float64, Q::Float64)
    return -k * sin(x) - γ * v + Q * sin(Ω * t)
end

function rk4_step(x::Float64, v::Float64, t::Float64, Q::Float64)
    k1x = τ * v
    k1v = τ * compute_acceleration(x, v, t, Q)

    k2x = τ * (v + 0.5 * k1v)
    k2v = τ * compute_acceleration(x + 0.5 * k1x, v + 0.5 * k1v, t + 0.5 * τ, Q)

    k3x = τ * (v + 0.5 * k2v)
    k3v = τ * compute_acceleration(x + 0.5 * k2x, v + 0.5 * k2v, t + 0.5 * τ, Q)

    k4x = τ * (v + k3v)
    k4v = τ * compute_acceleration(x + k3x, v + k3v, t + τ, Q)

    x_next = x + (1/6.0) * (k1x + 2*k2x + 2*k3x + k4x)
    v_next = v + (1/6.0) * (k1v + 2*k2v + 2*k3v + k4v)
    t_next = t + τ

    x_next = mod(x_next + π, 2π) - π

    return x_next, v_next, t_next
end


function simulate_for_Q(Q::Float64)
    #Initial conditions. 
    x = 1.0
    v = 0.0
    t = 0.0

    x_prev = x

    collected = Float64[]

    for step in 1:steps
        x, v, t = rk4_step(x, v, t, Q)

        # Begin recording at "steady state"
        if t > ignore_periods * T_0
            # Detect crossing from x > 0 to x < 0
            if x_prev > 0 && x < 0
                push!(collected, v)
            end
        end

        x_prev = x
    end

    return collected
end


function stroboscopic_plot(Q::Float64, sampling::Bool)
    x = 1.0
    v = 0.0
    t = 0.0

    x_vals = Float64[]
    v_vals = Float64[]

    for step in 1:steps
        x, v, t = rk4_step(x, v, t, Q)
        
        # Collect only after transient
        if t > ignore_periods * T_0
            # Stroboscopic sampling at multiples of T0
            if sampling
                if isapprox(mod(t, T_0), 0.0, atol=τ/2)
                    push!(x_vals, x)
                    push!(v_vals, v)
                end
            else 
                push!(x_vals, x)
                push!(v_vals, v)
            end
            
        end
    end

    fig = GLMakie.Figure()
    ax = Axis(fig[1, 1]; xlabel="x", ylabel="dx/dt", title="Stroboscopic Plot at Q = $Q")
    #ax = Axis(fig[1, 1];xlabel="x",ylabel="dx/dt",title="Stroboscopic Plot at Q = $Q",limits = ((-4, 4), (-4, 4)))
    scatter!(ax, x_vals, v_vals; markersize=2)
    return fig
    
end


function bifurcation_plot()
    Q_vals = 0.0:0.001:2.0
    Q_plot = Float64[]
    V_plot = Float64[]

    for Q in Q_vals
        velocities = simulate_for_Q(Q)
        for v in velocities
            push!(Q_plot, Q)
            push!(V_plot, v)
        end
    end

    fig = GLMakie.Figure()
    ax = Axis(fig[1, 1]; xlabel="Q", ylabel="v (at x-crossing)", title="Bifurcation Diagram")
    scatter!(ax, Q_plot, V_plot; markersize=1)
    return fig
end

function power_spectrum(Q::Float64)
    steps_b = 2^16

    x = 1.0
    v = 0.0
    t = 0.0

    x_values = Float64[]

    for step in 1:steps_b
        x, v, t = rk4_step(x, v, t, Q)
        push!(x_values, x)
    end

    X = fft(x_values)

    #the mean ”power” of the xj
    P = abs.(X).^2 / steps_b

    f0 = Ω / (2π)
    fs = 1 / τ
    
    freqs = fs * (0:(steps_b ÷ 2)) ./ steps_b

    # Extract positive frequencies
    logP = log.(P[1:(steps_b ÷ 2 + 1)])
    freqs_plot = freqs[1:(steps_b ÷ 2 + 1)]

    fig = Figure(resolution=(800, 500))
    ax = Axis(fig[1, 1], xlabel="Frequency f", ylabel="log P_L(f)", title="Power Spectrum Q = $Q")
    lines!(ax, freqs_plot, logP)
    xlims!(ax, 0, 6f0)
    
    return fig
end

fig1 = power_spectrum(0.5)
fig2 = power_spectrum(1.07)
fig3 = power_spectrum(1.2)
save("PowerSpectrum0.5.png", fig1)
save("PowerSpectrum1.07.png", fig2)
save("PowerSpectrum1.2.png", fig3)

#=
bifurcation = bifurcation_plot()
save("Bifurcation.png", bifurcation)
=#

#=
sampling = false
fig1 = stroboscopic_plot(0.5, sampling)
fig2 = stroboscopic_plot(1.07, sampling)
fig3 = stroboscopic_plot(1.2, sampling)

save("Two_Pi_Omega.png", fig1)
save("Four_Pi_Omega.png", fig2)
save("Aperiodic.png", fig3)
=#


