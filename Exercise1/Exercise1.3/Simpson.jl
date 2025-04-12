using GLMakie

# Define the function to integrate
f(x) = sin(x)

# Simpson's rule implementation
function simpson_rule(f, a, b, N)
    s = 0
    s_even = 0
    h = (b - a) / (N-1) # N-1 is the number of intervals (N is number of points)
    x = range(a, b, length=N)

    #Accoutn for case N is even
    if N % 2 == 0
        s_even += (h/12)*((8*f(x[N-1])) + (5*(f(x[N]))) - f(x[N-2]))
        N -= 1
    end

    #Boundary values
    s += f(x[1]) + f(x[N])

    #Iterate over even x_k
    for i in 2:2:N-1
        s += 4 * f(x[i])
    end

    #Iterate over odd x_k
    for i in 3:2:N-2
        s += 2 * f(x[i])
    end
    return ((h / 3) * s) + s_even
end

# Error analysis setup
a = 0                  # Lower bound of the integral
b = π/2                # Upper bound of the integral
exact = 1.0            # Exact analytical result of ∫₀^{π/2} sin(x) dx = 1
Ns = 3:1:1000          # Range of N values (number of sample points, N-1 intervals)

# Compute absolute errors for each N using Simpson's rule
errors = [abs(simpson_rule(f, a, b, N) - exact) for N in Ns]

# Plotting with GLMakie
fig = Figure(size=(800, 500))  # Create a new figure with size 800x500 pixels

# Create an axis on the figure:
# - Placed at row 1, column 1
# - Label the x-axis and y-axis
# - Set the y-axis to log scale to better visualize small errors
ax = Axis(fig[1, 1],
    xlabel="N (number of intervals)",
    ylabel="Absolute Error",
    yscale=log10
)

# Plot a line of the errors
lines!(ax, Ns, errors, color=:blue, label="Error")

# Overlay scatter points on top of the line for visibility
scatter!(ax, Ns, errors, color=:blue)

# Add a legend to the axis using the label from the lines! plot
axislegend(ax)

# Save the figure to a PNG file
save("error.png", fig)

# Display the figure interactively
display(fig)