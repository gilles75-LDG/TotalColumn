using Plots
using CSV
using DataFrames
# using NetCDF
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

#path to netDCF file

#define text strings for things
nc_pre = "/home/lawson/Data/EM27SUN/nc_files/"

# dates = []

variables  = ["xch4", "xco2", "xco", "flag","year","day",
            "hour","xluft","tins","xluft_error","o2_7885_sg",
            "pout","pins"]

instruments=["ta","tb"]

#pick pieces to make sure it works
#
# date = dates[1]
#
# var = variables[1]
#
# ins = instruments[2]

# nc_path =  ins * date * "_" *  date * ".private.nc"

for date in dates
    for ins in instruments
        for var in variables
            cd(nc_pre)
            nc_path =  ins * date * "_" *  date * ".private.nc"
            if !(isfile(nc_path))
                run(`scp carbon:/export/data/em27/analysis2020/$(ins)/daily/$(date)/*.private.nc ./`)
            end
            if !(isfile(nc_pre*ins*date*"_"*var*".txt"))
                cd(nc_pre)
                run(pipeline(`ncdump -v $(var) $(nc_path)`, "txt/$(ins)$(date)_$(var).txt"))
            end

        end
    end
end

# lines = readlines(nc_pre * "txt/$(ins)$(date)_$(var).txt")

#build ta and tb data_frames

# tb_lin_num = 7701
#
# ta_lin_num =  7906

ta_data = DataFrame()
vars = []

#loop to read in ta nc_txt data
for var in variables
    var_data = []
    println(var)
    for date in dates
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

        substrs = Float64[]
        nc_txt = read(nc_pre*"txt/ta$(date)_$(var).txt",String)
        sp = split(nc_txt,"data:") #creates array of strings with data and extra bits
        things = split(sp[2],"\n")
        #
        println(things)
        for j in things
            z = split(j,",")
            for t in z
               if contains(t,"=")
                   push!(substrs, parse(Float64,split(split(t,"=")[2],";")[1]))
               elseif contains(t,";")
                   push!(substrs, parse(Float64,split(t,";")[1]))
               elseif t == " "
               elseif t == ""
                   # push!(substrs,missing)
               elseif t == "}"
               else
                   push!(substrs, parse(Float64,t))
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

vars_b = []

for var in variables
    var_data = []
    println(var)
    for date in dates
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

        substrs = Float64[]
        nc_txt = read(nc_pre*"txt/tb$(date)_$(var).txt",String)
        sp = split(nc_txt,"data:") #creates array of strings with data and extra bits
        things = split(sp[2],"\n")
        #
        println(things)
        for j in things
            z = split(j,",")
            for t in z
               if contains(t,"=")
                   push!(substrs, parse(Float64,split(split(t,"=")[2],";")[1]))
               elseif contains(t,";")
                   push!(substrs, parse(Float64,split(t,";")[1]))
               elseif t == " "
               elseif t == ""
                   # push!(substrs,missing)
               elseif t == "}"
               else
                   push!(substrs, parse(Float64,t))
               end
            end
           end
        push!(var_data,substrs)
        # push!(ta_data,var_data,cols="$var")
        # println(var)
        # println(var_data)
    end
    push!(vars_b,var_data)
end


#making a dataframe from the data
ta_data = DataFrame(variables[1]=>unravel(vars[1]))
for i in 2:length(variables)
    va_nm = variables[i]
    ta_data[!,"$va_nm"] = unravel(vars[i])
end




ta_data.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(ta_data) ]
# data.sec = [ floor((r.hour - r.min/60)*60)    for r in eachrow(data) ]
ta_data.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(ta_data) ]

ta_data.month = [1.0 for i in eachrow(ta_data)]
for (i,val) in enumerate(ta_data.day)
    if 59 > val > 31
        ta_data.month[i] = 2.0
        ta_data.day[i] = ta_data.day[i] - 31
    elseif 60 .< val
        ta_data.month[i] = 3.0
        ta_data.day[i] = ta_data.day[i] - 31 - 28
    end

end


ta_data.dt = [ Dates.DateTime(r.year,r.month,r.day,floor(r.hour),r.min,r.sec) for r in eachrow(ta_data)]




#rinse and repeat for tb
tb_data = DataFrame(variables[1]=>unravel(vars_b[1]))
for i in 2:length(variables)
    va_nm = variables[i]
    tb_data[!,"$va_nm"] = unravel(vars_b[i])
end




tb_data.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(tb_data) ]
# data.sec = [ floor((r.hour - r.min/60)*60)    for r in eachrow(data) ]
tb_data.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(tb_data) ]

tb_data.month = [1.0 for i in eachrow(tb_data)]
for (i,val) in enumerate(tb_data.day)
    if  59 > val > 31
        tb_data.month[i] = 2.0
        tb_data.day[i] = tb_data.day[i] - 31
    elseif val .> 60
        tb_data.month[i] = 3.0
        tb_data.day[i] = tb_data.day[i] - 28  - 31
    end
end


tb_data.dt = [ Dates.DateTime(r.year,r.month,r.day,floor(r.hour),r.min,r.sec) for r in eachrow(tb_data)]

#---
#figure out how to read in spectra from netdcf ...
ta_spec = []

for date in dates
    nc_txt = read(nc_pre*"txt/ta$(date)_spectrum.txt",String)
    sp = split(nc_txt,"data:") #creates array of strings with data and extra bits
    things = split(sp[2],"\n")
    keeps = []
    for t in things
            if contains(t,'\"')
                push!(keeps,split(t,'\"')[2])
            end
    end
    push!(ta_spec,keeps)
end

ta_spec = unravel_strs(ta_spec)

ta_data.spectrum = ta_spec


tb_spec = []

for date in dates
    nc_txt = read(nc_pre*"txt/tb$(date)_spectrum.txt",String)
    sp = split(nc_txt,"data:") #creates array of strings with data and extra bits
    things = split(sp[2],"\n")
    keeps = []
    for t in things
            if contains(t,'\"')
                push!(keeps,split(t,'\"')[2])
            end
    end
    push!(tb_spec,keeps)
end

tb_spec = unravel_strs(tb_spec)

tb_data.spectrum = tb_spec

#-----



filtered_data = ta_data[ta_data.flag .== 0,:]
filtered_data_b = tb_data[tb_data.flag .== 0,:]




d = 17
m = 3

day1=filtered_data[filtered_data.month .== m, :]
day1=day1[day1.day .== d, :]

day1_b = filtered_data_b[filtered_data_b.month .== m, :]
day1_b = day1_b[day1_b.day .== d, :]


date = Date(day1.dt[1])

pltpath = "/home/lawson/Data/EM27SUN/plots/$(date)/"

if ! isdir(pltpath)
    mkdir(pltpath)
end

#------

#Spectra giffing

cd("/home/lawson/Data/EM27SUN/spectra/tb/")


spec1 = day1_b.spectrum[500]

for spec1 in day1_b.spectrum
    if ! (isfile("z"*spec1))
        print(spec1[end-3:end]*",")
    end
end

anim = @animate for spec1 in day1_b.spectrum
    if (isfile("z"*spec1))
        spec_data = CSV.read("z"*spec1,DataFrame;header=3,delim=" ",ignorerepeated=true)
        rename!(spec_data, [strip(n) for n in names(spec_data)])


        l = @layout [a{0.2h}; b]

        plot1 = plot(spec_data.Freq,spec_data.Tc./spec_data.Cont, label = "TC", color =:green)

        plot!(spec_data.Freq, spec_data.Tm./spec_data.Cont, label = "TM", color =:black,
            xlabel = "Wavenumber in cm⁻¹" )

        plot!(plot1, spec_data.Freq, spec_data.o2, label = "O₂", color =:red)
        plot!(plot1, spec_data.Freq, spec_data."0o2", label = "O₂ Cont.", color = :lightgreen)
        plot!(plot1, spec_data.Freq, spec_data.co2, label = "CO₂", color = :yellow)
        plot!(plot1, spec_data.Freq, spec_data.h2o, label = "H₂O", color = :lightblue)
        plot!(plot1, spec_data.Freq, spec_data.hf, label = "HF", color = :grey)
        plot!(plot1, spec_data.Freq, spec_data.other, label = "Other", color = :pink)
        plot!(plot1, spec_data.Freq, spec_data.solar, label = "Solar", color = :orange,
            legend=:bottomright)
        ylabel!(plot1,"Transmittance")


        spec_data.res = (spec_data.Tc - spec_data.Tm)./spec_data.Cont
        #

        plot2 = plot(spec_data.Freq, spec_data.res, ylabel = "Residual", legend = false, linewidth = 0.4,)
                # ylims = (-0.01,0.01) )
                  # ylims =  (-0.08,0.08)) #
        title!(plot2,"$spec1")
        ylims!(plot2, (-0.01,0.01))
        plot3 = plot(plot2, plot1, layout = l, figuresize = (600,400), dpi = 100)
        # file2 = "/home/lawson/Data/EM27/plots/test2.png"
    end
end

gif(anim)
gif(anim,pltpath*"tb_20220128_zo2.gif")


# png(plot3, plotfilename)


#-----


scatter(day1.dt,day1.xch4,c=:red,msw=0,label="ta CH₄")
scatter!(day1_b.dt,day1_b.xch4,c=:orange,msw=0,xrot=-10,
        label="tb CH₄",markersize=3,ylabel="XCH₄ (ppm)",legend=:bottomleft)
plot!(legend=:topleft)
title!("$date Coincidental XCH₄")
png(pltpath*"xch4_$date")


scatter(day1.dt,day1.xco2,c=:red,msw=0,label="ta CO₂")
scatter!(day1_b.dt,day1_b.xco2,c=:orange,msw=0,xrot=-10,
        label="tb CO₂",markersize=3,ylabel="XCO₂ (ppm)",legend=:bottomleft)
plot!(legend=:topleft)
title!("$date Coincidental XCO₂")
png(pltpath*"xco2_$date")

scatter(day1.dt,day1.xco,c=:red,msw=0,label="ta CO")
scatter!(day1_b.dt,day1_b.xco,c=:orange,msw=0,xrot=-10,
        label="tb CO",ylabel="XCO (ppm)",
        markersize=3,legend=:bottomleft)
plot!(legend=:topleft)
title!("$date Coincidental XCO")
png(pltpath*"xco_$date")



scatter(day1.dt,day1.xluft,c=:red,msw=0,label="ta")
scatter!(day1_b.dt,day1_b.xluft,c=:orange,msw=0,xrot=-10,
        label="tb",ylabel="Xluft",legend=:bottomleft,ylims = [0.99,1.01])
title!("$date Coincidental Xluft")
png(pltpath*"xluft_$date")


#
# scatter(day1.dt,day1.o2_ct,c=:red,msw=0,label="ta O₂ CT")
# scatter!(day1_b.dt,day1_b.o2_ct,c=:orange,msw=0,xrot=-10,
#         label="tb O₂ CT",ylabel="CT",legend=:bottomleft,ylims=[-10,-3])
# title!("$date Coincidental O₂ Contimum Level Tilt")


#Create 15 second bins
if day1.dt[end]-day1.dt[1] < day1_b.dt[end]-day1_b.dt[1]
    shorter_day = day1
else
    shorter_day = day1_b
end


t1 = shorter_day.dt[1]
t_f = shorter_day.dt[end]
#counter to combine data later
day1.counter = zeros(length(eachrow(day1)))
day1_b.counter = zeros(length(eachrow(day1_b)))
c=1

while t1 < t_f
    for (i,d) in enumerate(day1.dt)
        if (d > t1) && (d<t1+Dates.Second(20))
            println("Bang ! ",d)
            day1.counter[i] = c
            # day1_b.counter[i]= c
        end

    end
    for (i,d) in enumerate(day1_b.dt)
        if (d > t1) && (d<t1+Dates.Second(20))
            println("Boom ! ",d)
            # day1.counter[i] = c
            day1_b.counter[i]= c
        end

    end
    c=c+1
    t1 = t1 + Dates.Second(20)
    # println(t1)
end

day1 = day1[day1.counter .> 0,:]
day1_b = day1_b[day1_b.counter .> 0 ,:]

combined_data = innerjoin(day1,day1_b; on=:counter,makeunique=true)



scatter(combined_data.dt,combined_data.tins_1,label="tb",ylabel="Inside Temp (°C)",msw=0,xrot=10)
scatter!(combined_data.dt,combined_data.tins,label="ta",msw=0,xrot=10,title="$date Spectrometer Temps")
plot!(legend=:bottom)
png(pltpath*"$(date)_spec_temps")

scatter(combined_data.xch4,combined_data.xch4_1,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,marker_z=combined_data.tins,label=false,
        colorbartitle=" \nta Inside Temp (°C)",rightmargin=5Plots.mm)
title!("$date Coincidental XCH₄")
png(pltpath*"$(date)xch4_ta_tins")

scatter(combined_data.xch4,combined_data.xch4_1,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,marker_z=combined_data.tins_1,label=false,
        colorbartitle=" \ntb Inside Temp (°C)",rightmargin=5Plots.mm)
title!("$date Coincidental XCH₄")
png(pltpath*"$(date)xch4_tb_tins")

scatter(combined_data.xch4,combined_data.xch4_1,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,
        marker_z=combined_data.tins_1.-combined_data.tins,label=false,
        colorbartitle=" \ntb-ta Δ Inside Temp (°C)",rightmargin=5Plots.mm)
title!("$date Coincidental XCH₄")
png(pltpath*"$(date)xch4_tinsdiff")

scatter(combined_data.xch4,combined_data.xch4_1,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,
        marker_z=combined_data.xluft_1.-combined_data.xluft,
        colorbartitle=" \n(tb-ta) ΔXLuft",label=false,rightmargin=5Plots.mm)
title!("$(date) Coincidental XCH₄")
png(pltpath*"$(date)xch4_xluftDelta")


scatter(combined_data.xch4,combined_data.xch4_1,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,
        marker_z=combined_data.xluft_error_1,
        colorbartitle=" \ntb Xluft Error",label=false,rightmargin=5Plots.mm)
title!("$(date) Coincidental XCH₄")
png(pltpath*"$(date)xch4_tb_xluft_err")

hodl~


#-------------------------------------------------------------------------------


filtered_data = ta_data[ta_data.flag .== 0,:]
rename!(filtered_data, [n*"_ta" for n in names(filtered_data)])
filtered_data_b = tb_data[tb_data.flag .== 0,:]
rename!(filtered_data_b, [n*"_tb" for n in names(filtered_data_b)])


#Create 15 second bins
if filtered_data.dt_ta[end]-filtered_data.dt_ta[1] < filtered_data_b.dt_tb[end]-filtered_data_b.dt_tb[1]
    shorter_day = filtered_data
else
    shorter_day = filtered_data_b
end


t1 = shorter_day.dt_tb[1]
t_f = shorter_day.dt_tb[end]
#counter to combine data later
filtered_data.counter = zeros(length(eachrow(filtered_data)))
filtered_data_b.counter = zeros(length(eachrow(filtered_data_b)))
c=1

while t1 < t_f
    for (i,d) in enumerate(filtered_data.dt_ta)
        if (d > t1) && (d<t1+Dates.Second(20))
            println("Bang ! ",d)
            filtered_data.counter[i] = c
            # day1_b.counter[i]= c
        end

    end
    for (i,d) in enumerate(filtered_data_b.dt_tb)
        if (d > t1) && (d<t1+Dates.Second(20))
            println("Boom ! ",d)
            # day1.counter[i] = c
            filtered_data_b.counter[i]= c
        end

    end
    c=c+1
    t1 = t1 + Dates.Second(20)
    # println(t1)
end

filtered_data = filtered_data[filtered_data.counter .> 0,:]
filtered_data_b = filtered_data_b[filtered_data_b.counter .> 0,:]

long_data = innerjoin(filtered_data,filtered_data_b;on=:counter)

long_data.xch4_diff = long_data.xch4_ta - long_data.xch4_tb
long_data.xco2_diff = long_data.xco2_ta - long_data.xco2_tb
long_data.xco_diff = long_data.xco_ta - long_data.xco_tb


scatter(f_long_data.xch4_diff,label="ta-tb 20 second coincidental",ylabel="XCH₄ ta-tb",xlabel="Observation #")
hline!([mean(f_long_data.xch4_diff)],label="Average Offset")
xch4_off = mean(f_long_data.xch4_diff)
title!("XCH₄ ta vs. tb offset\n$xch4_off")
png(pltpath*"xch4_off")


scatter(f_long_data.xco_diff,label="ta-tb 20 second coincidental",ylabel="XCO ta-tb",xlabel="Observation #")
hline!([mean(f_long_data.xco_diff)],label="Average Offset")
xco_off = mean(f_long_data.xco_diff)
title!("XCO ta vs. tb offset\n$xco_off")
png(pltpath*"xco_off")


scatter(f_long_data.xco2_diff,label="ta-tb 20 second coincidental",ylabel="XCO₂ ta-tb",xlabel="Observation #")
hline!([mean(f_long_data.xco2_diff)],label="Average Offset")
xco2_off = mean(f_long_data.xco2_diff)
title!("XCO₂ ta vs. tb offset\n$xco2_off")
png(pltpath*"xco2_off")


short_data = long_data[long_data.day .> 25,:]

pltpath = "/home/lawson/Data/EM27SUN/plots/"
scatter(long_data.dt_tb,long_data.tins_tb,label="tb",ylabel="Inside Temp (°C)",msw=0,xrot=10)
scatter!(long_data.dt_ta,long_data.tins_ta,label="ta",msw=0,xrot=10,title="2022 01/02 Spectrometer Temps")
plot!(legend=:bottom)
png(pltpath*"all_spec_temps")

scatter(long_data.xch4_ta,long_data.xch4_tb,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,marker_z=long_data.tins_ta,label=false,
        colorbartitle=" \nta Inside Temp (°C)",rightmargin=5Plots.mm)
title!("2022 01/02 Coincidental XCH₄")
png(pltpath*"all_xch4_ta_tins")

scatter(long_data.xch4_ta,long_data.xch4_tb,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,marker_z=long_data.tins_tb,label=false,
        colorbartitle=" \ntb Inside Temp (°C)",rightmargin=5Plots.mm)
title!("2022 01/02 Coincidental XCH₄")
png(pltpath*"all_xch4_tb_tins")




scatter(long_data.xch4_ta,long_data.xch4_tb,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,
        marker_z=long_data.tins_tb.-long_data.tins_ta,label=false,
        colorbartitle=" \ntb-ta Δ Inside Temp (°C)",rightmargin=5Plots.mm)
title!("2022 01/02 Coincidental XCH₄")
png(pltpath*"all_xch4_tinsdiff")

scatter(long_data.xch4_ta,long_data.xch4_tb,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,
        marker_z=long_data.xluft_tb.-long_data.xluft_ta,
        colorbartitle=" \n(tb-ta) ΔXLuft",label=false,rightmargin=5Plots.mm)
title!("2022 01/02 Coincidental XCH₄")
png(pltpath*"all_xch4_xluftDelta")


scatter(long_data.xch4_ta,long_data.xch4_tb,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,
        marker_z=long_data.xluft_error_tb,
        colorbartitle=" \n \ntb Xluft Error",label=false,rightmargin=8Plots.mm)
title!("2022 01/02 Coincidental XCH₄")
png(pltpath*"all_xch4_tb_xluft_err")


scatter(long_data.xch4_ta,long_data.xch4_tb,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,
        marker_z=long_data.xluft_error_ta,
        colorbartitle=" \n \nta Xluft Error",label=false,rightmargin=8Plots.mm)
title!("2022 01/02 Coincidental XCH₄")
png(pltpath*"all_xch4_ta_xluft_err")


scatter(long_data.xco2_ta,long_data.xco2_tb,xlabel="ta XCO₂ (ppm)",
        ylabel="tb XCO₂ (ppm)",msw=0,
        marker_z=long_data.xluft_error_tb,
        colorbartitle=" \n \ntb Xluft Error",label=false,rightmargin=8Plots.mm)
title!("2022 01/02 Coincidental XCO₂")
png(pltpath*"all_xco2_tb_xluft_err")


scatter(long_data.xco_ta,long_data.xco_tb,xlabel="ta XCO (ppm)",
        ylabel="tb XCO (ppm)",msw=0,
        marker_z=long_data.xluft_error_tb,
        colorbartitle=" \n \ntb Xluft Error",label=false,rightmargin=8Plots.mm)
title!("2022 01/02 Coincidental XCO")
png(pltpath*"all_xco_tb_xluft_err")

scatter(long_data.xco_ta,long_data.xco_tb,xlabel="ta XCO (ppm)",
        ylabel="tb XCO (ppm)",msw=0,
        marker_z=long_data.xluft_error_ta,
        colorbartitle=" \n \nta Xluft Error",
        label=false,a=.3,rightmargin=8Plots.mm)
title!("2022 01/02 Coincidental XCO")
png(pltpath*"all_xco_ta_xluft_err")

scatter(long_data.xluft_ta,long_data.xluft_tb,xlabel="ta Xluft",
        ylabel="tb Xluft",msw=0,
        marker_z=long_data.xluft_error_tb,
        colorbartitle=" \n \ntb Xluft Error",label=false,rightmargin=8Plots.mm)
title!("2022 01/02 Coincidental Xluft")
png(pltpath*"all_xluft_tbxlufterr")



scatter(long_data.xluft_ta,long_data.xluft_tb,xlabel="ta Xluft",
        ylabel="tb Xluft",msw=0,
        marker_z=long_data.xluft_error_ta,
        colorbartitle=" \n \nta Xluft Error",label=false,rightmargin=8Plots.mm)
title!("2022 01/02 Coincidental Xluft")
png(pltpath*"all_xluft_taxlufterr")

f_long_data = long_data[long_data.xluft_error_ta .< 0.003,:]
f_long_data  = f_long_data[f_long_data.xluft_error_tb .< 0.003, :]


scatter(f_long_data.xch4_ta,f_long_data.xch4_tb,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,
        marker_z=f_long_data.xluft_error_ta,
        colorbartitle=" \n \nta Xluft Error",label=false,rightmargin=8Plots.mm)
title!("2022 01/02 Xluft Error filtered Coincidental XCH₄")
png(pltpath*"all_xch4_ta_xluft_err_filt")


scatter(f_long_data.xco2_ta,f_long_data.xco2_tb,xlabel="ta XCO₂ (ppm)",
        ylabel="tb XCO₂ (ppm)",msw=0,
        marker_z=f_long_data.xluft_error_tb,
        colorbartitle=" \n \ntb Xluft Error",label=false,rightmargin=8Plots.mm)
title!("2022 01/02 Xluft Error Filtered Coincidental XCO₂")
png(pltpath*"all_xco2_tb_xluft_err_filt")


scatter(f_long_data.xco_ta,f_long_data.xco_tb,xlabel="ta XCO (ppm)",
        ylabel="tb XCO (ppm)",msw=0,
        marker_z=f_long_data.xluft_error_tb,
        colorbartitle=" \n \ntb Xluft Error",label=false,rightmargin=8Plots.mm)
title!("2022 01/02 Xluft Error Filtered Coincidental XCO")
png(pltpath*"all_xco_tb_xluft_err_filt")






scatter(long_data.xch4_ta,long_data.xch4_tb,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,
        c=:red,
        # marker_z=long_data.xluft_error_ta,
        # colorbartitle=" \n \nta Xluft Error",
        label="EGI Filtered")

scatter!(f_long_data.xch4_ta,f_long_data.xch4_tb,xlabel="ta XCH₄ (ppm)",
        ylabel="tb XCH₄ (ppm)",msw=0,c=:blue,markersize=3,
        # marker_z=f_long_data.xluft_error_ta,
        # colorbartitle=" \n \nta Xluft Error",
        label="+Xluft filter",legend=:topleft)
png(pltpath*"xch4_filt")

scatter(long_data.xco2_ta,long_data.xco2_tb,xlabel="ta XCO₂ (ppm)",
        ylabel="tb XCO₂ (ppm)",msw=0,
        c=:red,
        # marker_z=long_data.xluft_error_ta,
        # colorbartitle=" \n \nta Xluft Error",
        label="EGI Filtered")

scatter!(f_long_data.xco2_ta,f_long_data.xco2_tb,
        msw=0,c=:blue,markersize=3,
        # marker_z=f_long_data.xluft_error_ta,
        # colorbartitle=" \n \nta Xluft Error",
        label="+Xluft filter",legend=:topleft)
png(pltpath*"xco2_filt")



scatter(long_data.xco_ta,long_data.xco_tb,xlabel="ta XCO (ppm)",
        ylabel="tb XCO (ppm)",msw=0,
        c=:red,
        # marker_z=long_data.xluft_error_ta,
        # colorbartitle=" \n \nta Xluft Error",
        label="EGI Filtered")

scatter!(f_long_data.xco_ta,f_long_data.xco_tb,
        msw=0,c=:blue,markersize=3,
        # marker_z=f_long_data.xluft_error_ta,
        # colorbartitle=" \n \nta Xluft Error",
        label="+Xluft filter",legend=:topleft)
png(pltpath*"xco_filt")


t1 = f_long_data.dt_ta[1]

ta_avgs=[]
tb_avgs=[]

while t1 < f_long_data.dt_ta[end]
    ta_line_avg_xch4 = []
    tb_line_avg_xch4 = []
    ta_line_avg_xco2 = []
    tb_line_avg_xco2 = []
    ta_line_avg_xco = []
    tb_line_avg_xco = []
    ta_line_avg_sg = []
    tb_line_avg_sg = []
    for r in eachrow(f_long_data)
        if t1 < r.dt_ta < t1+Dates.Minute(10)
            push!(ta_line_avg_xch4,r.xch4_ta)
            push!(tb_line_avg_xch4,r.xch4_tb)
            push!(ta_line_avg_xco2,r.xco2_ta)
            push!(tb_line_avg_xco2,r.xco2_tb)
            push!(ta_line_avg_xco,r.xco_ta)
            push!(tb_line_avg_xco,r.xco_tb)
            push!(ta_line_avg_sg,r.o2_7885_sg_ta)
            push!(tb_line_avg_sg,r.o2_7885_sg_tb)
        end
    end
        if length(tb_line_avg_xco) > 0
            println("Boop : ", t1)
            push!(ta_avgs,[mean(ta_line_avg_xch4),mean(ta_line_avg_xco2),mean(ta_line_avg_xco),mean(ta_line_avg_sg)])
            push!(tb_avgs,[mean(tb_line_avg_xch4),mean(tb_line_avg_xco2),mean(tb_line_avg_xco),mean(tb_line_avg_sg)])
        end
    t1 = t1 + Dates.Minute(10)
end


ta_xch4_avgs = [ta_avgs[n][1] for n in 1:length(ta_avgs)]
ta_xco2_avgs = [ta_avgs[n][2] for n in 1:length(ta_avgs)]
ta_xco_avgs = [ta_avgs[n][3] for n in 1:length(ta_avgs)]
ta_sg_avgs = [ta_avgs[n][4] for n in 1:length(ta_avgs)]
tb_xch4_avgs = [tb_avgs[n][1] for n in 1:length(tb_avgs)]
tb_xco2_avgs = [tb_avgs[n][2] for n in 1:length(tb_avgs)]
tb_xco_avgs = [tb_avgs[n][3] for n in 1:length(tb_avgs)]
tb_sg_avgs = [tb_avgs[n][4] for n in 1:length(tb_avgs)]


scatter(ta_xch4_avgs-tb_xch4_avgs,msw=0,marker_z=tb_sg_avgs,
        colorbartitle=" \ntb O₂ SG",rightmargin=5Plots.mm)
hline!([median(ta_xch4_avgs-tb_xch4_avgs)],label="Median")
hline!([mean(ta_xch4_avgs-tb_xch4_avgs)],label="Mean",ylabel="ta - tb XCH₄ (ppm)",
            xlabel="Obs #")
xch4_med = median(ta_xch4_avgs-tb_xch4_avgs)
xch4_mean = mean(ta_xch4_avgs-tb_xch4_avgs)
title!("10 minute averages, ta-tb XCH₄\nMedian=$(xch4_med)\nMean=$(xch4_mean)")
png(pltpath*"10minxch4_tb_sg")


scatter(ta_xco2_avgs-tb_xco2_avgs)
hline!([median(ta_xco2_avgs-tb_xco2_avgs)],label="Median")
hline!([mean(ta_xco2_avgs-tb_xco2_avgs)],label="Mean",
        ylabel="XCO₂ (ppm)",xlabel="Obs #")
xco2_med = median(ta_xco2_avgs-tb_xco2_avgs)
xco2_mean = mean(ta_xco2_avgs-tb_xco2_avgs)
title!("10 minute averages, ta-tb XCO₂\nMedian=$(xco2_med)\nMean=$(xco2_mean)")
png(pltpath*"10minxco2")


scatter(ta_xco_avgs-tb_xco_avgs,msw=0,marker_z=ta_sg_avgs,
        colorbartitle=" \nta O₂ SG",rightmargin=5Plots.mm)
hline!([median(ta_xco_avgs-tb_xco_avgs)],label="Median")
hline!([mean(ta_xco_avgs-tb_xco_avgs)],label="Mean",ylabel="ta-tb XCO (ppm)",
    xlabel="Obs #")
xco_med = median(ta_xco_avgs-tb_xco_avgs)
xco_mean = mean(ta_xco_avgs-tb_xco_avgs)
title!("10 minute averages, ta-tb XCO\nMedian=$(xco_med)\nMean=$(xco_mean)")
png(pltpath*"10minxco_ta_sg")

scatter(ta_xch4_avgs,tb_xch4_avgs,xlabel="ta XCH₄",
            ylabel="tb XCH₄",label="Data",msw=0)
plot!([1.875, 1.92], [1.875,1.92],label="1:1")
scatter!(ta_xch4_avgs,tb_xch4_avgs.+xch4_med,xlabel="ta XCH₄",
            ylabel="tb XCH₄",label="10 Min Median Corrected Data",
            msw=0,legend=:topleft)
title!("2022 Jan-March Coincidental\nXCH₄ 10 Minute Averages")
png(pltpath*"10minCorrectedXCH4")


scatter(ta_xco2_avgs,tb_xco2_avgs,xlabel="ta XCO₂",
            ylabel="tb XCO₂",label="Data",msw=0)
plot!([minimum(ta_xco2_avgs), maximum(ta_xco2_avgs)],
    [minimum(ta_xco2_avgs),maximum(ta_xco2_avgs)],label="1:1")
scatter!(ta_xco2_avgs,tb_xco2_avgs.+xco2_med,xlabel="ta XCO₂",
            ylabel="tb XCO₂",label="10 Min Median Corrected Data",
            msw=0,legend=:topleft)
title!("2022 Jan-March Coincidental\nXCO₂ 10 Minute Averages")
png(pltpath*"10minCorrectedXCO2")



scatter(ta_xco_avgs,tb_xco_avgs,xlabel="ta XCO",
            ylabel="tb XCO",label="Data",msw=0)
plot!([minimum(ta_xco_avgs), maximum(ta_xco_avgs)],
    [minimum(ta_xco_avgs),maximum(ta_xco_avgs)],label="1:1")
scatter!(ta_xco_avgs,tb_xco_avgs.+xco_med,xlabel="ta XCO",
            ylabel="tb XCO",label="10 Min Median Corrected Data",
            msw=0,legend=:topleft)
title!("2022 Jan-March Coincidental\nXCO 10 Minute Averages")
png(pltpath*"10minCorrectedXCO")


# hufl


using GLM


fit = lm(@formula(xch4 ~  xch4_1 ), long_data)
long_data.fit = predict(fit,long_data)

rge = range(minimum(long_data.xch4_1),maximum(long_data.xch4_1);step=0.0001)

dte = predict(fit,DataFrame(xch4_1=collect(rge)))

#
#
# lines = readlines(nc_pre*"txt/tb20220116_xch4.txt")_
# lines = lines[7701:end-1]
#
# substrs = Float64[]
# for line in lines
#    sstrs = split(line,",")
#    for t in sstrs
#        if contains(t,"=")
#            push!(substrs, parse(Float64,split(t,"=")[2]))
#        elseif contains(t,";")
#            push!(substrs, parse(Float64,split(t,";")[1]))
#        elseif t == " "
#        else
#            push!(substrs, parse(Float64,t))
#        end
#    end
# end
