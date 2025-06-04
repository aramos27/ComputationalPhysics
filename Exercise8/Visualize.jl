#Code done with extensive reasoning with ChatGPT
using GLMakie

function visualize_matrix_grid(mat::AbstractMatrix{Int})
    unique_vals = sort(collect(Set(mat)))
    val_to_color = Dict(val => i for (i, val) in enumerate(unique_vals))
    n_colors = length(unique_vals)
    color_palette = cgrad(:Set3_12, n_colors; categorical=true).colors

    rows, cols = size(mat)
    fig = Figure(size = (600, 600))
    ax = Axis(
        fig[1, 1],
        aspect = DataAspect(),
        title = "Matrix Grid",
        xticksvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xlabelvisible = false,
        ylabelvisible = false
    )

    # Draw each cell as a rectangle with matching color
    for i in 1:rows, j in 1:cols
        val = mat[i, j]
        color_idx = val_to_color[val]
        rect_x = j - 1
        rect_y = rows - i  # Flip y to match heatmap style
        poly!(
            ax,
            [rect_x, rect_x + 1, rect_x + 1, rect_x],
            [rect_y, rect_y, rect_y + 1, rect_y + 1],
            color = color_palette[color_idx],
            strokewidth = 1,
            strokecolor = :black
        )
        text!(
            ax,
            string(val),
            position = (rect_x + 0.5, rect_y + 0.5),
            align = (:center, :center),
            color = :black,
            fontsize = 16
        )
    end

    hidespines!(ax)
    return fig
end

function visualize_multiple_matrices(mat1::AbstractMatrix{Int}, mat2::AbstractMatrix{Int}, mat3::AbstractMatrix{Int})
    # Combine all values to build a shared color map
    all_vals = vcat(mat1[:], mat2[:], mat3[:])
    unique_vals = sort(collect(Set(all_vals)))
    val_to_color = Dict(val => i for (i, val) in enumerate(unique_vals))
    n_colors = length(unique_vals)
    color_palette = cgrad(:Set3_12, n_colors; categorical=true).colors

    fig = Figure(size = (900, 300))  # Adjust width for 3 matrices

    mats = [mat1, mat2, mat3]
    titles = ["Occupation Matrix", "Label Matrix", "Merged Labels"]

    for (idx, mat) in enumerate(mats)
        rows, cols = size(mat)
        ax = Axis(
            fig[1, idx],
            aspect = DataAspect(),
            title = titles[idx],
            xticksvisible = false,
            yticksvisible = false,
            xticklabelsvisible = false,
            yticklabelsvisible = false,
            xgridvisible = false,
            ygridvisible = false,
            xlabelvisible = false,
            ylabelvisible = false
        )

        for i in 1:rows, j in 1:cols
            val = mat[i, j]
            color_idx = val_to_color[val]
            rect_x = j - 1
            rect_y = rows - i
            poly!(
                ax,
                [rect_x, rect_x + 1, rect_x + 1, rect_x],
                [rect_y, rect_y, rect_y + 1, rect_y + 1],
                color = color_palette[color_idx],
                strokewidth = 1,
                strokecolor = :black
            )
            text!(
                ax,
                string(val),
                position = (rect_x + 0.5, rect_y + 0.5),
                align = (:center, :center),
                color = :black,
                fontsize = 16
            )
        end

        hidespines!(ax)
    end

    return fig
end
