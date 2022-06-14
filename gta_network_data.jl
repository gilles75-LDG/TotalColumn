using DataFrames
using CSV


function unravel(thing)
    tmp = zeros( Int( sum( [ length(i) for i in thing ] ) ) )
    counter=1
    for i in thing
        for j in i
            tmp[counter] = j
            counter+=1
        end
    end
    thing = tmp
    return thing
end


function unravel_strs(thing)
    # tmp = zeros( Int( sum( [ length(i) for i in thing ] ) ) )
    # tmp = Array{Any,Int( sum( [ length(i) for i in thing ] ) )}
    tmp = []

    # counter=1
    for i in thing
        # println(i)
        for j in i
            println(j)
            # tmp[counter] = String(j)
            push!(tmp,j)
            # counter+=1
        end
    end
    thing = tmp
    return thing
end



#define text strings for things
nc_pre = "/home/lawson/Data/EM27SUN/nc_files/"

# dates = []

variables  = ["xch4", "xco2", "xco", "flag","year","day",
            "hour","xluft","tins","xluft_error","o2_7885_sg",
            "pout","pins"]

instruments=["ta","tb","tc","td"]
