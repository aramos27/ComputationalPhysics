using GLMakie
GLMakie.activate!(; float = true, focus_on_show = true)

#=
From the numbers provided in the exercise we can reduce the speed system to cells per second.
The possible speed values reduce to 0,1,2,3,4,5.
=#

# Define the Car struct
mutable struct Car
    position::Int  # Index on the track
    speed::Int  # Cells per time step
    distanceNextCar::Int       # Cells to next car
end

# Function to initialize the track with cars randomly
function init_cells(L::Int, r::Float64)
    cells = Vector{Union{Nothing, Car}}(undef, L)
    fill!(cells, nothing)

    for i in 1:L
        if rand() < r
            speed = 1  
            cells[i] = Car(i, speed, 0)
        end
    end

    compute_gaps(cells)

    return cells
end

#Initialize cells array with all cars at the left hand side (red light)
function init_cells_traffic_light(L::Int, r::Float64)
    cells = Vector{Union{Nothing, Car}}(undef, L)
    fill!(cells, nothing)
    total_cars = trunc(Int, L * r)
    for i in 1:L
        if i <= total_cars
            speed = 1  
            cells[i] = Car(i, speed, 0)
        end
    end

    compute_gaps(cells)

    return cells
end

#Function to compute distanceNextCar for each Car 
function compute_gaps(cells)
    L = length(cells)
    count = 0
    nextCar = nothing
    lastCarIndex = nothing
    for i in L:-1:1
        if (cells[i] !== nothing)
            if lastCarIndex === nothing
                lastCarIndex = i
            end
            if (nextCar !== nothing)
                cells[i].distanceNextCar = nextCar - i         
            end                 
            nextCar = i
        end
    end

    #Account for last car (boundary conditions)
    cells[lastCarIndex].distanceNextCar = (L-lastCarIndex) + nextCar
end

#Moves cars by adjusting their velocity and their position
function update_cells!(cells, p)
    compute_gaps(cells)

    L=length(cells)

    #Compute new speeds
    for i in 1:L
        if(cells[i] !== nothing)

            #Acceleration condition
            if(cells[i].speed < 5)
                cells[i].speed += 1
            end

            #Breaking condition
            if(cells[i].distanceNextCar <= cells[i].speed)
                cells[i].speed = cells[i].distanceNextCar - 1
            end

            #Dawdling condition
            if (rand() < p && cells[i].speed > 0)
                cells[i].speed -= 1
            end
        end
    end


    #Compute new positions under the assumption of no collisions 
    for i in 1:L
        if (cells[i] !== nothing)
            new_position = (cells[i].position + cells[i].speed) % L
            if new_position == 0
                cells[i].position = L
            else
                cells[i].position = new_position
            end
        end
    end

    #Move cars within the array under the assumption of no collisions
    for i in 1:L
        if cells[i] !== nothing
            car = cells[i]
            cells[i] = nothing  # Clear old position
            cells[car.position] = car  # Place in new position
        end
    end
end

#Function to count cars to ensure overwriting does not occurr. 
function count_cars(cells)
    count(car -> car !== nothing, cells)
end


function record_traffic(cells, Nt, p)
    Nc = length(cells)
    data = zeros(Bool, Nc, Nt)

    for t in 1:Nt
        for i in 1:Nc
            data[i, t] = cells[i] !== nothing
        end
        update_cells!(cells, p)
    end
    return data
end

function plot_traffic(data)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Position", ylabel="Time step", title="Traffic Space-Time Diagram")
    heatmap!(ax, data;)
    return fig
end

# Example usage:
Nc = 1000 #Track length
Nt = 1000 #Duration of simulation [s]
ρ = 0.1 #Probability of cell being occupied at initialization
p = 0.15 # Dawdling probability

cells = init_cells(Nc, ρ) #Random positions at beginning
#cells = init_cells_traffic_light(Nc, ρ) #Start at stoplight

filename = "traffic_rho$(ρ)_p$(p).png"

data = record_traffic(cells, Nt, p)
fig = plot_traffic(data)
#save(filename, fig)
display(fig)
wait(fig.scene)