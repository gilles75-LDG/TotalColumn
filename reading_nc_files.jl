using CSV
using DataFrames
using NetCDF
using Plots
using Glob
using Dates
using Statistics

nc_path = "/home/lawson/Data/EM27SUN/nc_files/ta20220121_20220121.private.nc"
nc_path_pre = "/home/lawson/Research/Data/EM27SUN/nc_files/"

ta_paths = Glob.glob("ta*.nc", nc_path_pre)
tb_paths = Glob.glob("tb*.nc", nc_path_pre)

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

NetCDF.open(nc_path) do nc
    # we now have an NcFile object that can be browsed
    println("Names of the NcFile fields: \n", fieldnames(typeof(nc)))

    for var in nc.vars
        if contains(var.first,"pout")
            println(var.first)
        end
        # println(var)
    end

end

ch4_data = [ncread(path, "xch4") for path in ta_paths
ch4_data = unravel(ch4_data)
time  = [ncread(path, "time") for path in ta_paths]
time = unravel(time)
flags = [ncread(path,"flag") for path in ta_paths]
flags = unravel(flags)
hour_a = [ncread(path,"hour") for path in ta_paths]
hour_a = unravel(hour_a)
year_a = [ncread(path,"year") for path in ta_paths]
year_a = unravel(year_a)
day_a = [ncread(path,"day") for path in ta_paths]
day_a = unravel(day_a)
xluft = [ncread(path,"xluft") for path in ta_paths]
xluft = unravel(xluft)
pout = [ncread(path,"pout") for path in ta_paths]
pout= unravel(pout)


co2_data = [ncread(path, "xco2") for path in ta_paths]
co2_data = unravel(co2_data)
co_data = [ncread(path, "xco") for path in ta_paths]
co_data = unravel(co_data)
o2_cl_tilt = [ncread(path, "o2_7885_ct") for path in ta_paths]
o2_cl_tilt = unravel(o2_cl_tilt)

o2_vsf = [ncread(path,"o2_7885_vsf_o2") for path in ta_paths]
o2_vsf = unravel(o2_vsf)

p_temp = [ncread(path, "prior_temperature") for path in ta_paths]
p_temp = unravel(p_temp)
acs = unravel([ncread(path,"apriori_checksum") for path in ta_paths])
pel = unravel([ncread(path,"prior_effective_latitude") for path in ta_paths])
pmf = unravel([ncread(path, "prior_modfile") for path in ta_paths])
pvf = unravel([ncread(path, "prior_vmrfile") for path in ta_paths])


data = DataFrame(ch4 = ch4_data, time = time, flag = flags, hour = hour_a,
                year = year_a, day = day_a, co2 = co2_data, co_data = co_data,
                xluft = xluft,pres= pout, o2_vsf = o2_vsf, o2_ct = o2_cl_tilt,)
                # acs = acs)

data.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data) ]
# data.sec = [ floor((r.hour - r.min/60)*60)    for r in eachrow(data) ]
data.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data) ]

data.month = [1.0 for i in eachrow(data)]
for (i,val) in enumerate(data.day)
    if val > 31
        data.month[i] = 2.0
        data.day[i] = data.day[i] - 31
    end
end

data.dt = [ Dates.DateTime(r.year,r.month,r.day,floor(r.hour),r.min,r.sec) for r in eachrow(data)]

filtered_data = data[data.flag .== 0, :]

ch4_data_b = [ncread(path, "xch4") for path in tb_paths]
ch4_data_b = unravel(ch4_data_b)
time_b  = [ncread(path, "time") for path in tb_paths]
time_b = unravel(time_b)
flags_b = [ncread(path,"flag") for path in tb_paths]
flags_b = unravel(flags_b)
hour_b = [ncread(path,"hour") for path in tb_paths]
hour_b = unravel(hour_b)
year_b = [ncread(path,"year") for path in tb_paths]
year_b = unravel(year_b)
day_b = [ncread(path,"day") for path in tb_paths]
day_b = unravel(day_b)
co2_data_b = [ncread(path, "xco2") for path in tb_paths]
co2_data_b = unravel(co2_data_b)
co_data_b = [ncread(path, "xco") for path in tb_paths]
co_data_b = unravel(co_data_b)
xluft_b = [ncread(path,"xluft") for path in tb_paths]
xluft_b = unravel(xluft_b)
pout_b = [ncread(path,"pout") for path in tb_paths]
pout_b = unravel(pout_b)
o2_cl_tilt_b = [ncread(path, "o2_7885_ct") for path in tb_paths]
o2_cl_tilt_b = unravel(o2_cl_tilt_b)

o2_vsf_b = [ncread(path,"o2_7885_vsf_o2") for path in tb_paths]
o2_vsf_b = unravel(o2_vsf_b)

p_temp_b = [AbstractChar.(ncread(path,"prior_temperature")) for path in tb_paths]

pel_b = unravel([ncread(path,"prior_effective_latitude") for path in tb_paths])

pmf_b = unravel([ncread(path,"prior_modfile") for path in tb_paths])
pvf_b = unravel([ncread(path, "prior_vmrfile") for path in tb_paths])



data_b = DataFrame(ch4 = ch4_data_b, time = time_b, flag = flags_b,
                    year=year_b,day=day_b,hour=hour_b,
                    co2 = co2_data_b, co_data = co_data_b, xluft = xluft_b,
                     pres= pout_b, o2_vsf = o2_vsf_b, o2_ct = o2_cl_tilt_b)
data_b.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data_b) ]
# data_b.sec = [ floor((r.hour - r.min/60)*60)    for r in eachrow(data_b) ]
data_b.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data_b) ]
data_b.month = [1.0 for i in eachrow(data_b)]
for (i,val) in enumerate(data_b.day)
    if val > 31
        data_b.month[i] = 2.0
        data_b.day[i] = data_b.day[i] - 31
    end
end



data_b.dt = [ Dates.DateTime(r.year,r.month,r.day,floor(r.hour),r.min,r.sec) for r in eachrow(data_b)]


filtered_data_b = data_b[data_b.flag .== 0, :]


d = 14
m = 2
day1=filtered_data[filtered_data.month .== m, :]
day1=day1[day1.day .== d, :]

day1_b = filtered_data_b[filtered_data_b.month .== m, :]
day1_b = day1_b[day1_b.day .== d, :]


date = Date(day1.dt[1])

scatter(day1.dt,day1.ch4,c=:red,msw=0,label="ta CH₄")
scatter!(day1_b.dt,day1_b.ch4,c=:orange,msw=0,xrot=-10,
        label="tb CH₄",ylabel="CH₄ (ppm)",legend=:bottomleft)
title!("$date Coincidental CH₄ Columns")


scatter(day1.dt,day1.o2_vsf,c=:red,msw=0,label="ta O₂ VSF")
scatter!(day1_b.dt,day1_b.o2_vsf,c=:orange,msw=0,xrot=-10,
        label="tb O₂ VSF",ylabel="VMR SF",legend=:bottomleft,ylims = [0.99,1.01])
title!("$date Coincidental O₂ VSF's")


scatter(day1.dt,day1.o2_ct,c=:red,msw=0,label="ta O₂ CT")
scatter!(day1_b.dt,day1_b.o2_ct,c=:orange,msw=0,xrot=-10,
        label="tb O₂ CT",ylabel="CT",legend=:bottomleft,ylims=[-10,-3])
title!("$date Coincidental O₂ Contimum Level Tilt")

#
# scatter(day1.dt,day1.pres,c=:purple,msw=0,label="ta")
# scatter!(day1_b.dt,day1_b.pres,c=:pink,msw=0,label="tb",
#         ylabel="Pout", xrot=-10,legend=:topright,markersize=3)
# title!("2022-01-$d Coindidental Pressure Observations")
# png(pltpath*"pres$d")


#day1 coincidental 2022 days so far
#16,21,28,29, feb 1, 14th

days = [16,21,28,29]
plot1 = plot()
plot2 = plot()
plot3 = plot()

for day in days
    day_of_year = day
    day1=filtered_data[filtered_data.day .== day_of_year, :]
    day1_b = filtered_data_b[filtered_data_b.day .== day_of_year, :]
    plot4=plot()
    scatter(plot4,filtered_data.dt,filtered_data.ch4,c=:red,msw=0)
    scatter!(plot4,filtered_data_b.dt,filtered_data_b.ch4,c=:orange,msw=0)
    pltpath = "/home/lawson/Research/Data/EM27SUN/plots/"
    png(pltpath*"$(day)")

    date = Date(day1.dt[1])

    # scatter(day1.dt,day1.ch4,c=:red,msw=0,label="ta CH₄")
    # scatter!(day1_b.dt,day1_b.ch4,c=:orange,msw=0,xrot=-10,
    #         label="tb CH₄",ylabel="Molecules m⁻²",legend=:bottomleft)
    # title!("$date Coincidental CH₄ Columns")
    #
    # png(pltpath*"$(date)_ch4")

    #bin the data by combining together

    #have some start time
    start_time = day1.dt[1]
    stop_time = day1.dt[end]
    total_min = (stop_time-start_time).value/1000/60
    t_bins = floor(total_min/5)

    #for every n * 5 minutes since then, average the times

    time_bins = DateTime[]
    # counter=1
    ch4_avg = []
    co2_avg = []
    co_avg = []
    ch4_avg_b = []
    co2_avg_b = []
    co_avg_b = []

    for i in 1:t_bins
        tmp = []
        tmp_co2 = []
        tmp_co = []
        for r in eachrow(day1)
            if (start_time + Dates.Minute(5 * i-1) .< r.dt .< start_time +   Dates.Minute(5*i))
                append!(tmp,r.ch4)
                append!(tmp_co2,r.co2)
                append!(tmp_co,r.co_data)
            end
        end

        tmp2 = []
        tmp2_co2 = []
        tmp2_co = []
        for r in eachrow(day1_b)
            if (start_time + Dates.Minute(5 * i-1) .< r.dt .< start_time +   Dates.Minute(5*i))
                append!(tmp2,r.ch4)
                append!(tmp2_co2,r.co2)
                append!(tmp2_co,r.co_data)

            end
        end

        # println(tmp)
        if length(tmp) > 0 && length(tmp2) > 0
            append!(ch4_avg,mean(tmp))
            append!(time_bins,[DateTime(start_time + Dates.Minute(5 * i-1))] )
            append!(ch4_avg_b,mean(tmp2))
            append!(co2_avg,mean(tmp_co2))
            append!(co2_avg_b,mean(tmp2_co2))
            append!(co_avg,mean(tmp_co))
            append!(co_avg_b,mean(tmp2_co))

        end
    end

    # scatter(ch4_avg,ch4_avg_b,xlabel="ta CH₄ molecules m⁻²",
    #     ylabel="tb CH₄ molecules m⁻²",legend=:topleft,label=date)


    scatter!(plot1,ch4_avg,ch4_avg_b,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",legend=:topleft,label=date)
    scatter!(plot2,co2_avg,co2_avg_b,xlabel="ta XCO₂ (ppm)",
        ylabel="tb XCO₂ (ppm)",legend=:topleft,label=date)
    scatter!(plot3,co_avg,co_avg_b,xlabel="ta XCO (ppm)",
        ylabel="tb XCO (ppm)",legend=:topleft,label=date)


end
plot!(plot1,[3.98e19, 4.1e19],[3.98e19, 4.1e19],label="1:1")

title!(plot3,"Coincidental January 2022\nta and tb CO columns, binned")
png(pltpath*"ta_tb_comp_co")
