using Test

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

@testset "decideLabel tests" begin
    # Format: @test decideLabel(input_array, next_label) == (expected_label, expected_merge)

    @test decideLabel([0,2,0,4], 1) == (2, true, 2, 4)
    @test decideLabel([0,0,1,0], 1) == (1, false)
    @test decideLabel([0,0,0,0], 1) == (1, false)
    @test decideLabel([0,4,0,3], 1) == (3, true, 4, 3)
    @test decideLabel([5,0,0,1], 1) == (1, true, 5, 1)
    @test decideLabel([2,0,0,0], 1) == (2, false)
    @test decideLabel([nothing,nothing,0,0],1) == (1,false)
    @test decideLabel([nothing,0,0,0],1) == (1,false)
    @test decideLabel([3,nothing,0,4],1) == (3,true, 3, 4)

    # Repeat with different next_label to ensure next_label usage correctness
    @test decideLabel([0,2,0,4], 2) == (2, true, 2, 4)
    @test decideLabel([0,0,1,0], 2) == (1, false)
    @test decideLabel([0,0,0,0], 2) == (2, false)
    @test decideLabel([0,4,0,3], 2) == (3, true, 4, 3)
    @test decideLabel([5,0,0,1], 2) == (1, true, 5, 1)
    @test decideLabel([2,0,0,0], 2) == (2, false)
    @test decideLabel([nothing,nothing,0,0],2) == (2,false)
    @test decideLabel([nothing,0,0,0],2) == (2,false)
    @test decideLabel([3,nothing,0,4],2) == (3,true, 3, 4)
end