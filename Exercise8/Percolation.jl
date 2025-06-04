using GLMakie
include("Visualize.jl")

function make_labels(L)
    return collect(1:ceil(Int,(((L*L)/2 + 1))))
end


function create_grid(L, p)
    grid = Int.(rand(L, L) .< p)
    properLabels = make_labels(L)

    return grid, properLabels
end

function decideLabel(neighbors, next_label)
    count = 0
    label1 = 0
    label2 = 0

    for k in 1:4
        if !isnothing(neighbors[k])
            if neighbors[k] != 0
                if count < 1
                    label1 = neighbors[k]
                else
                    label2 = neighbors[k]
                end
                count +=1
            end
        end
    end

    if count == 0
        return next_label, false
    elseif count == 1
        return label1, false
    else
        if label1 == label2
            return label1, false
        else 
            return min(label1, label2), true, label1, label2
        end
    end        
end

#Finds "root" label
function find(labels, x)
    while labels[x] != x
        labels[x] = labels[labels[x]]  
        x = labels[x]
    end
    return x
end

#Reassigns root
function union!(labels, x, y)
    root_x = find(labels, x)
    root_y = find(labels, y)
    labels[root_x] = root_y
end

function assign_labels(grid, properLabels)
    side_length = size(grid, 1)
    labelGrid = zeros(Int, side_length, side_length)

    nextLabel = 1

    for i in 1:side_length
        for j in 1:side_length
            if grid[i,j] != 0 #Position is occupied

                up, down, right, left = 0 ,0 ,0 ,0
                up    = i > 1           ? labelGrid[i-1, j] : nothing
                down  = i < side_length ? labelGrid[i+1, j] : nothing
                left  = j > 1           ? labelGrid[i, j-1] : nothing
                right = j < side_length ? labelGrid[i, j+1] : nothing
                
                neighbors = [up,right,down,left]

                output = decideLabel(neighbors, nextLabel)
                label = output[1]
                update_needed = output[2]

                labelGrid[i,j] = label

                if label == nextLabel
                    nextLabel += 1
                end

                if update_needed
                    union!(properLabels,max(output[3], output[4]), min(output[3], output[4]))
                end
                    
            end
        end
    end
    return labelGrid
end

function merge_clusters!(labelGrid, properLabels)
    side_length = size(labelGrid, 1)

    for i in 1:side_length
        for j in 1:side_length
            if labelGrid[i,j] != 0
                labelGrid[i,j] = find(properLabels, labelGrid[i,j])
            end
        end
    end
end

function check_percolating(labelGrid)
    side_length = size(labelGrid, 1)
    for i in 1:side_length
        temp = labelGrid[1,i]
        temp2 = labelGrid[i,1]
        if temp != 0
            for j in 1:side_length
                if temp == labelGrid[side_length, j]
                    return true
                end
            end
        end
        if temp2 != 0
            for j in 1:side_length
                if temp2 == labelGrid[j, side_length]
                    return true
                end
            end
        end
    end
    return false
end

function percolating_probability(L, p, N)
    count = 0
    for i in 1:N
        grid, properLabels = create_grid(10, p)

        labelGrid = assign_labels(grid, properLabels)
        merge_clusters!(labelGrid, properLabels)
        flag = check_percolating(labelGrid)

        if flag
            count += 1
        end
    end
    return (count/N)*100
end

function plot_probability(L_range, N)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], 
              title = "Percolation Probability vs p",
              xlabel = "p",
              ylabel = "Percolation Probability")

    for L in L_range
        p_values = 0.0:0.01:1.0
        probabilities = [percolating_probability(L, p, N) for p in p_values]
        lines!(ax, p_values, probabilities, label = "L = $L")
    end

    axislegend(ax)
    return fig
end

#fig = plot_probability(50:50:400, 100000)
#save("plot.png", fig)


