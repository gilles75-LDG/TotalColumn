using DataFrames
using CSV
using Glob

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
nc_pre = "/home/lawson/HDD/Data/em27_ncfiles/"

utsg_files = Glob.glob("t*/UTSG/*.nc",nc_pre)
utsc_files = Glob.glob("t*/UTSC/*.nc",nc_pre)
utm_files = Glob.glob("t*/UTM/*.nc",nc_pre)
dow_files = Glob.glob("t*/ECCC_Downsview/*.nc",nc_pre)



# dates = []

variables  = ["xch4", "xch4_error", "xco2", "xco2_error", "xco", "xco_error",
              "xh2o", "xh2o_error", "flag","year","day",
              "hour","xluft","tins","xluft_error","o2_7885_sg",
              "pout","pins"]

txt_pre = "/home/lawson/HDD/Data/em27_txt_files/"

for file in utsg_files
    cd(txt_pre*"UTSG")
    file_name_head = split(file,"/")[end][1:10]
    println(file_name_head)
    for var in variables
        if !(isfile(file_name_head*"_"*var*".txt"))
            run(pipeline(`ncdump -v $(var) $(file)`, "$(file_name_head)_$(var).txt"))
        end
    end
end

for file in utm_files
    cd(txt_pre*"UTM")
    file_name_head = split(file,"/")[end][1:10]
    println(file_name_head)
    for var in variables
        if !(isfile(file_name_head*"_"*var*".txt"))
            run(pipeline(`ncdump -v $(var) $(file)`, "$(file_name_head)_$(var).txt"))
        end
    end
end

for file in utsc_files
    cd(txt_pre*"UTSC")
    file_name_head = split(file,"/")[end][1:10]
    println(file_name_head)
    for var in variables
        if !(isfile(file_name_head*"_"*var*".txt"))
            run(pipeline(`ncdump -v $(var) $(file)`, "$(file_name_head)_$(var).txt"))
        end
    end
end

for file in dow_files
    cd(txt_pre*"ECCC_Downsview")
    file_name_head = split(file,"/")[end][1:10]
    println(file_name_head)
    for var in variables
        if !(isfile(file_name_head*"_"*var*".txt"))
            run(pipeline(`ncdump -v $(var) $(file)`, "$(file_name_head)_$(var).txt"))
        end
    end
end


# instruments=["ta","tb","tc","td"]
#
# for date in dates
#     for ins in instruments
#         for var in variables
#             cd(nc_pre)
#             nc_path =  ins * date * "_" *  date * ".private.nc"
#             if !(isfile(nc_path))
#                 run(`scp carbon:/export/data/em27/analysis2020/$(ins)/daily/$(date)/*.private.nc ./`)
#             end
#             if !(isfile(nc_pre*ins*date*"_"*var*".txt"))
#                 cd(nc_pre)
#                 run(pipeline(`ncdump -v $(var) $(nc_path)`, "txt/$(ins)$(date)_$(var).txt"))
#             end
#
#         end
#     end
# end
