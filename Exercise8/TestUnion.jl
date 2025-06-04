function make_labels(n)
    return collect(1:n)  # Each label is its own root
end

function find(labels, x)
    while labels[x] != x
        labels[x] = labels[labels[x]]  # Path compression
        x = labels[x]
    end
    return x
end

function union!(labels, x, y)
    root_x = find(labels, x)
    root_y = find(labels, y)
    labels[root_x] = root_y
end

labels = make_labels(10)

union!(labels, 2, 3)   # merge label 2 into label 3
println(labels)
union!(labels, 3, 5)   # merge label 3 into label 5
println(labels)
union!(labels, 4, 5)   # merge label 4 into label 5
println(labels)

println("Label roots:")
for i in 1:6
    println("Label $i → Root ", find(labels, i))
end

println("Internal label array: ", labels)

