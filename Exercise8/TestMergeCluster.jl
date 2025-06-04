include("Percolation.jl")

grid, properLabels = create_grid(10, 0.6)

labelGrid = assign_labels(grid, properLabels)
labelGrid_copy = copy(labelGrid)
merge_clusters!(labelGrid, properLabels)

fig = visualize_multiple_matrices(grid, labelGrid_copy, labelGrid)
save("Three_Matrices_Grid.png", fig)

