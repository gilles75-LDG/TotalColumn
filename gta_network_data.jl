using DataFrames
using CSV
using Glob
using Dates
using Plots

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

# for file in utsg_files
#     cd(txt_pre*"UTSG")
#     file_name_head = split(file,"/")[end][1:10]
#     println(file_name_head)
#     for var in variables
#         if !(isfile(file_name_head*"_"*var*".txt"))
#             run(pipeline(`ncdump -v $(var) $(file)`, "$(file_name_head)_$(var).txt"))
#         end
#     end
# end
#
# for file in utm_files
#     cd(txt_pre*"UTM")
#     file_name_head = split(file,"/")[end][1:10]
#     println(file_name_head)
#     for var in variables
#         if !(isfile(file_name_head*"_"*var*".txt"))
#             run(pipeline(`ncdump -v $(var) $(file)`, "$(file_name_head)_$(var).txt"))
#         end
#     end
# end
#
# for file in utsc_files
#     cd(txt_pre*"UTSC")
#     file_name_head = split(file,"/")[end][1:10]
#     println(file_name_head)
#     for var in variables
#         if !(isfile(file_name_head*"_"*var*".txt"))
#             run(pipeline(`ncdump -v $(var) $(file)`, "$(file_name_head)_$(var).txt"))
#         end
#     end
# end
#
# for file in dow_files
#     cd(txt_pre*"ECCC_Downsview")
#     file_name_head = split(file,"/")[end][1:10]
#     println(file_name_head)
#     for var in variables
#         if !(isfile(file_name_head*"_"*var*".txt"))
#             run(pipeline(`ncdump -v $(var) $(file)`, "$(file_name_head)_$(var).txt"))
#         end
#     end
# end



utsg_data = DataFrame(ins = [], month = [])
vars = []
# dates =
#loop to read in ta nc_txt data
for var in variables
    var_data = []
    println(var)
    for file in utsg_files
        # lines = readlines(nc_pre*"txt/ta$(date)_$(var).txt")
        #
        # lines = lines[ta_lin_num:end-1]
        # for line in lines
        #    # println(line)
        #    sstrs = split(line,",")
        #    substrs = Float64[]
        #    for t in sstrs
        #        if contains(t,"=")
        #            push!(substrs, parse(Float64,split(t,"=")[2]))
        #        elseif contains(t,";")
        #            push!(substrs, parse(Float64,split(t,";")[1]))
        #        elseif t == " "
        #        elseif t == ""
        #            push!(substrs,missing)
        #        else
        #            push!(substrs, parse(Float64,t))
        #        end
        #    end
        #    append!(var_data,substrs)
        # end
        file_name_head = split(file,"/")[end][1:10]
        substrs = Float64[]
        nc_txt = read(txt_pre*"UTSG/"*file_name_head*"_$(var).txt",String)
        sp = split(nc_txt,"data:") #creates array of strings with data and extra bits
        things = split(sp[2],"\n")
        inst = split(file,"/")[end][1:2]
        month = file[end-14:end-13]
        tmp_df = DataFrame(ins = [inst], month = [month])
        #
        # println(things)
        # println(file_name_head)
        for j in things
            z = split(j,",")
            for t in z
                # println(t)
               if contains(t,"=")
                   # println(split(split(t,"=")[2],";")[1] )
                   push!(substrs, parse(Float64,split(split(t,"=")[2],";")[1]))
                   append!(utsg_data,tmp_df)
               elseif contains(t,";")
                   push!(substrs, parse(Float64,split(t,";")[1]))
                   append!(utsg_data,tmp_df)
               elseif t == " "
               elseif t == ""
                   # push!(substrs,missing)
               elseif t == "}"
               else
                   try
                       push!(substrs, parse(Float64,t))
                       append!(utsg_data,tmp_df)
                   catch
                       println(t, " ", var)
                       push!(substrs, -99.99)
                       append!(utsg_data,tmp_df)
                       # throw(DomainError(t))
                   end
               end
            end
           end
        push!(var_data,substrs)
        # push!(ta_data,var_data,cols="$var")
        # println(var)
        # println(var_data)
    end
    push!(vars,var_data)
end

#making a dataframe from the data
utsg_data2 = DataFrame(variables[1]=>unravel(vars[1]))
for i in 2:length(variables)
    va_nm = variables[i]
    utsg_data2[!,"$va_nm"] = unravel(vars[i])
end

utsg_data = first(utsg_data,length(eachrow(utsg_data2)))
utsg_data = hcat(utsg_data,utsg_data2)

# utsg_data = utsg_data2
utsg_data.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(utsg_data) ]
# data.sec = [ floor((r.hour - r.min/60)*60)    for r in eachrow(data) ]
utsg_data.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(utsg_data) ]
# utsg_data.month = [f[end-14:end-13] for f in utsg_files]



utsg_data.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(utsg_data)]
dic = Dict([("ta",:red),("tb",:orange),("tc",:blue),("td",:lightblue)])
utsg_data.color = [dic[r] for r in utsg_data.ins]
scatter(utsg_data.dt,utsg_data.xch4,msw=0,color=utsg_data.color)

utsg_data_filt = utsg_data[utsg_data.flag .== 0,:]

utsg_data_filt  = utsg_data_filt[Dates.DateTime(2019) .< utsg_data_filt.dt .< Dates.DateTime(2020),:]

utsg_ta = utsg_data_filt[utsg_data_filt.ins .== "ta", :]
utsg_tb = utsg_data_filt[utsg_data_filt.ins .== "tb", :]
utsg_tc = utsg_data_filt[utsg_data_filt.ins .== "tc", :]
utsg_td = utsg_data[utsg_data.ins .== "td", :]

scatter(utsg_ta.dt, utsg_ta.xch4, color = :red, label = "ta" ,msw = 0,
            markersize = 3)
scatter!(utsg_tb.dt, utsg_tb.xch4, color = :orange, label = "tb", msw = 0, a= .3)
scatter!(utsg_tc.dt, utsg_tc.xch4, color = :green, label = "tc", msw = 0, a= .3)

scatter(utsg_ta.dt, utsg_ta.xluft, color = :red, label = "ta" ,msw = 0,
            markersize = 3)
scatter!(utsg_tb.dt, utsg_tb.xluft, color = :orange, label = "tb", msw = 0, a= .3)
scatter!(utsg_tc.dt, utsg_tc.xluft, color = :green, label = "tc", msw = 0, a= .3)


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
