#Code done by extensive reasoning with ChatGPT

using GLMakie
using ColorTypes
GLMakie.activate!()
include("Visualize.jl")
include("Percolation.jl")

function interactive_labeling_demo(L, p)
    grid, properLabels = create_grid(L, p)
    labelGrid = zeros(Int, L, L)
    properLabels_step = copy(properLabels)
    nextLabel = Observable(1)

    current_i = Observable(1)
    current_j = Observable(1)
    current_phase = Observable(:assign_labels)
    labelGrid_obs = Observable(deepcopy(labelGrid))
    labels_array_obs = Observable(copy(properLabels_step))

    label_color_map = Dict{Int, Colorant}()  # Persistent label-to-color mapping

    fig = Figure(size = (1000, 800))

    ax_occ = Axis(fig[1, 1], title = "Occupation Matrix", aspect = DataAspect())
    ax_lab = Axis(fig[1, 2], title = "Label Matrix", aspect = DataAspect())
    ax_array = Axis(fig[2, 1:2], title = "Label Array", aspect = DataAspect())

    function draw_matrix!(ax, mat; red_highlight = Set(), blue_highlight = Set())
        rows, cols = size(mat)
        highlight_union = union(red_highlight, blue_highlight)

        # First pass: draw all non-highlight cells
        for i in 1:rows, j in 1:cols
            if (i, j) in highlight_union
                continue
            end
            val = mat[i, j]
            rect_x = j - 1
            rect_y = rows - i

            if !haskey(label_color_map, val) && val != 0
                label_color_map[val] = rand(RGBf)
            end

            color = val == 0 ? :gray : label_color_map[val]

            poly!(
                ax,
                [rect_x, rect_x + 1, rect_x + 1, rect_x, rect_x],
                [rect_y, rect_y, rect_y + 1, rect_y + 1, rect_y],
                color = color,
                strokewidth = 2,
                strokecolor = :black
            )

            text!(
                ax,
                string(val),
                position = (rect_x + 0.5, rect_y + 0.5),
                align = (:center, :center),
                color = :black,
                fontsize = 14
            )
        end

        # Second pass: draw blue-highlighted cells
        for (i, j) in blue_highlight
            val = mat[i, j]
            rect_x = j - 1
            rect_y = rows - i

            if !haskey(label_color_map, val) && val != 0
                label_color_map[val] = rand(RGBf)
            end

            color = val == 0 ? :gray : label_color_map[val]

            poly!(
                ax,
                [rect_x, rect_x + 1, rect_x + 1, rect_x, rect_x],
                [rect_y, rect_y, rect_y + 1, rect_y + 1, rect_y],
                color = color,
                strokewidth = 2,
                strokecolor = :blue
            )

            text!(
                ax,
                string(val),
                position = (rect_x + 0.5, rect_y + 0.5),
                align = (:center, :center),
                color = :black,
                fontsize = 14
            )
        end

        # Third pass: draw red-highlighted cells
        for (i, j) in red_highlight
            val = mat[i, j]
            rect_x = j - 1
            rect_y = rows - i

            if !haskey(label_color_map, val) && val != 0
                label_color_map[val] = rand(RGBf)
            end

            color = val == 0 ? :gray : label_color_map[val]

            poly!(
                ax,
                [rect_x, rect_x + 1, rect_x + 1, rect_x, rect_x],
                [rect_y, rect_y, rect_y + 1, rect_y + 1, rect_y],
                color = color,
                strokewidth = 2,
                strokecolor = :red
            )

            text!(
                ax,
                string(val),
                position = (rect_x + 0.5, rect_y + 0.5),
                align = (:center, :center),
                color = :black,
                fontsize = 14
            )
        end

        hidespines!(ax)
        hidedecorations!(ax)
    end

    function draw_array!(ax, arr)
        n = length(arr)
        for j in 1:n
            val = arr[j]
            rect_x = j - 1
            rect_y = 0
            poly!(
                ax,
                [rect_x, rect_x + 1, rect_x + 1, rect_x, rect_x],
                [rect_y, rect_y, rect_y + 1, rect_y + 1, rect_y],
                color = :white,
                strokewidth = 2,
                strokecolor = :black
            )
            text!(
                ax,
                string(val),
                position = (rect_x + 0.5, rect_y + 0.5),
                align = (:center, :center),
                color = :black,
                fontsize = 14
            )
        end
        hidespines!(ax)
        hidedecorations!(ax)
    end

    btn = Button(fig[3, 1:2], label = "Next Step")

    on(btn.clicks) do _
        i, j = current_i[], current_j[]
        neighbors_coords = Set{Tuple{Int,Int}}()

        if current_phase[] == :assign_labels
            if i > L
                current_phase[] = :merge
                current_i[] = 1
                current_j[] = 1
                return
            end

            if grid[i, j] != 0
                up    = i > 1 ? labelGrid[i - 1, j] : nothing
                down  = i < L ? labelGrid[i + 1, j] : nothing
                left  = j > 1 ? labelGrid[i, j - 1] : nothing
                right = j < L ? labelGrid[i, j + 1] : nothing
                neighbors = [up, right, down, left]

                coords = [
                    i > 1 ? (i - 1, j) : nothing,
                    i < L ? (i + 1, j) : nothing,
                    j > 1 ? (i, j - 1) : nothing,
                    j < L ? (i, j + 1) : nothing
                ]
                neighbors_coords = Set(filter(!isnothing, coords))

                output = decideLabel(neighbors, nextLabel[])
                label = output[1]
                update_needed = output[2]

                labelGrid[i, j] = label

                if label == nextLabel[]
                    nextLabel[] += 1
                end

                if update_needed
                    union!(properLabels_step, max(output[3], output[4]), min(output[3], output[4]))
                end

                labelGrid_obs[] = deepcopy(labelGrid)
                labels_array_obs[] = copy(properLabels_step)
            end

            if j == L
                current_i[] += 1
                current_j[] = 1
            else
                current_j[] += 1
            end

        elseif current_phase[] == :merge
            if i > L
                btn.label[] = "Done"
                return
            end

            if labelGrid[i, j] != 0
                labelGrid[i, j] = find(properLabels_step, labelGrid[i, j])
                labelGrid_obs[] = deepcopy(labelGrid)
                labels_array_obs[] = copy(properLabels_step)
            end

            if j == L
                current_i[] += 1
                current_j[] = 1
            else
                current_j[] += 1
            end
        end

        try
            empty!(ax_occ)
            draw_matrix!(ax_occ, grid;
                red_highlight = Set([(i, j)]),
                blue_highlight = neighbors_coords
            )
        catch e
            @error "Error drawing occupation matrix" exception = e
        end

        try
            empty!(ax_lab)
            draw_matrix!(ax_lab, labelGrid_obs[];
                red_highlight = Set([(i, j)]),
                blue_highlight = neighbors_coords
            )
        catch e
            @error "Error drawing label matrix" exception = e
        end

        try
            empty!(ax_array)
            draw_array!(ax_array, labels_array_obs[])
        catch e
            @error "Error drawing label array" exception = e
        end
    end

    # Initial draw
    draw_matrix!(ax_occ, grid)
    draw_matrix!(ax_lab, labelGrid_obs[])
    draw_array!(ax_array, labels_array_obs[])

    fig
end

# Run standalone
if abspath(PROGRAM_FILE) == @__FILE__
    fig = interactive_labeling_demo(10, 0.5)
    display(fig)
    wait(fig.scene)
end