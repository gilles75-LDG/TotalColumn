using NCDatasets
using DataFrames
using CSV
using Dates
using Plots
using Glob
using Statistics
using PyCall
using StatsBase

# #path to netCDF's
nc_path = "/home/lawson/Data/EM27SUN/nc_files/te/"
daterng = 202109
daterng = 2022072

nc_paths = glob("te$(daterng)*.nc",nc_path)
# append!(nc_paths,[i for i in glob("te2023*.nc",nc_path)])
# #make a list of NCDatasets from netCDF's
datasets = []
for path in nc_paths
    ds = Dataset(path)
    push!(datasets,ds)
end

# #list of variables to extract to dataframe
variables  = ["xch4", "xco2", "xco", "flag","year","day",
            "hour","xluft","tins","xluft_error","o2_7885_sg",
            "pout","pins","lat","long","column_luft","xluft","column_h2o","xh2o",
            "xch4_aicf","xch4_aicf","solzen","spectrum"]
kernels = ["ak_xch4","ak_xco2","ak_xco","ak_xh2o","ak_xo2","ak_pressure","ak_altitude"]


# #make a list of dataframes to dump NCDataset data into
dataframes = []
aks = []
spectra = []
for ds in datasets
    tmp_df = DataFrame()
    tmp_df2 = Dict()
    for var in variables
        println(var)
        if var != "spectrum"
            tmp_df[!,"$var"] = ds[var][1:end]
        else
            # strs = []
            # for i in 1:length(ds[var])÷21
            #     str = string()
            #     for j in 1:21
            #         str = str*ds[var][j,i]
            #     end
            #     push!(strs,str)
            # end
            # tmp_df[!,"$var"] = strs
            push!(spectra,ds[var])
        end
        
    end

    for ker in kernels
        println(ds[ker])
        tmp_df2["$ker"] = ds[ker][:,:]
    end
    push!(aks,tmp_df2)
    push!(dataframes,tmp_df)
end

for a in 1:length(dataframes)
    t = spectra[a]
    strs = String[]
    for i in 1:length(t)÷21
        str=string()
        for j in 1:21
            str = str*t[j,i]
        end
        println(str)
        push!(strs,str)
    end
    dataframes[a].spectrum = strs
end


# #make it all into single dataframe
data1_df  = dataframes[1]
for df in dataframes[2:end]
    data1_df = vcat(data1_df,df)
end


data1_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data1_df) ]
data1_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data1_df) ]
data1_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day - 1) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data1_df)]
data1_df.ins = ["te" for i in eachrow(data1_df)]
unique!(data1_df,:dt)

data_te = data1_df[data1_df.flag .== 0,:]

scatter(data_te.dt,data_te.xch4,msw=0.1)

# #now do tf data...

nc_path  = "/home/lawson/Data/EM27SUN/nc_files/tf/"
nc_paths = glob("tf$(daterng)*.nc",nc_path)

# #make a list of NCDatasets from netCDF's
datasets = []
for path in nc_paths
    ds = Dataset(path)
    push!(datasets,ds)
end

# # #list of variables to extract to dataframe
# variables  = ["xch4", "xco2", "xco", "flag","year","day",
#             "hour","xluft","tins","xluft_error","o2_7885_sg",
#             "pout","pins","lat","long","column_luft"]


# #make a list of dataframes to dump NCDataset data into

# #make a list of dataframes to dump NCDataset data into
dataframes = []
aks = []
spectra = []
for ds in datasets
    tmp_df = DataFrame()
    tmp_df2 = Dict()
    for var in variables
        println(var)
        if var != "spectrum"
            tmp_df[!,"$var"] = ds[var][1:end]
        else
            # strs = []
            # for i in 1:length(ds[var])÷21
            #     str = string()
            #     for j in 1:21
            #         str = str*ds[var][j,i]
            #     end
            #     push!(strs,str)
            # end
            # tmp_df[!,"$var"] = strs
            push!(spectra,ds[var])
        end
        
    end

    for ker in kernels
        println(ds[ker])
        tmp_df2["$ker"] = ds[ker][:,:]
    end
    push!(aks,tmp_df2)
    push!(dataframes,tmp_df)
end

for a in 1:length(dataframes)
    t = spectra[a]
    strs = String[]
    for i in 1:length(t)÷21
        str=string()
        for j in 1:21
            str = str*t[j,i]
        end
        println(str)
        push!(strs,str)
    end
    dataframes[a].spectrum = strs
end


# #make it all into single dataframe
data1_df  = dataframes[1]
for df in dataframes[2:end]
    data1_df = vcat(data1_df,df)
end


data1_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data1_df) ]
data1_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data1_df) ]
data1_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day-1) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data1_df)]
data1_df.ins = ["tf" for i in eachrow(data1_df)]
unique!(data1_df,:dt)

data_tf = data1_df[data1_df.flag .== 0,:]

scatter!(data_tf.dt,data_tf.xch4,msw=0.1,α=.3)

#---2022
data_te1 = data_te[data_te.dt .< Date(2022,07,27),:]
data_tf1 = data_tf[data_tf.dt .< Date(2022,07,27),:]

data_te2 = data_te[Date(2022,07,28).> data_te.dt .> Date(2022,07,27),:]
data_tf2 = data_tf[Date(2022,07,28).> data_tf.dt .> Date(2022,07,27),:]

data_te3 = data_te[Date(2022,07,30).> data_te.dt .> Date(2022,07,29),:]
data_tf3 = data_tf[Date(2022,07,30).> data_tf.dt .> Date(2022,07,29),:]

#---2021
# data_te1 = data_te[data_te.dt .< Date(2021,09,16), :]
# data_tf1 = data_tf[Date(2021,09,15) .< data_tf.dt .< Date(2021,09,16), :]

# data_te2 = data_te[Date(2021,09,17) .< data_te.dt .< Date(2021,09,18), :]
# data_tf2 = data_tf[Date(2021,09,17) .< data_tf.dt .< Date(2021,09,18), :]

# data_te3 = data_te[Date(2021,09,18) .< data_te.dt .< Date(2021,09,21), :]
# data_tf3 = data_tf[Date(2021,09,18) .< data_tf.dt .< Date(2021,09,21), :]

# data_te4 = data_te[Date(2021,09,21) .< data_te.dt .< Date(2021,09,25), :]
# data_tf4 = data_tf[Date(2021,09,21) .< data_tf.dt .< Date(2021,09,25), :]

# data_te_11 = data_te1[data_te1.dt .> DateTime(2021,9,15,18),:]
# data_tf_11 = data_tf1[data_tf1.dt .> DateTime(2021,9,15,18),:]

scatter(data_te1.dt,data_te1.xch4,msw=0.1,α=1)
scatter!(data_tf1.dt,data_tf1.xch4,msw=0.1,α=.3)

scatter(data_te2.dt,data_te2.xch4,msw=0.1,α=1,label="te")
scatter!(data_tf2.dt,data_tf2.xch4,msw=0.1,α=.3,label="tf")

scatter(data_te3.dt,data_te3.xch4,msw=0.1,α=1)
scatter!(data_tf3.dt,data_tf3.xch4,msw=0.1,α=.3)

# scatter(data_te4.dt,data_te4.xch4,msw=0.1,α=1)
# scatter!(data_tf4.dt,data_tf4.xch4,msw=0.1,α=.3)

#---
# read in the weather data
# weather_path = "/home/lawson/Data/EM27SUN/weather_data/LandfillCampaigns/harry/"
weather_path = "/home/lawson/Data/EM27SUN/weather_data/old_stuff/"
files = glob("$(daterng)*.txt", weather_path)
# files = glob("202109*vaisala.txt", weather_path)

dataframes = []
for f in files
    push!(dataframes,CSV.read(f,DataFrame))
end

# #make it all into single dataframe
data1_df  = dataframes[1]
for df in dataframes[2:end]
    data1_df = vcat(data1_df,df)
end
data1_df.dt = [DateTime(r.UTCDate*r.UTCTime[1:12], "yyyy/mm/ddHH:MM:SS.sss") for r in eachrow(data1_df)]
weather_data = data1_df[data1_df.dt .> Date(2021,09),:]


rename!(weather_data, [strip(n, ' ') for n in names(weather_data)])
rename!(weather_data, [strip(n, ' ') for n in names(weather_data)])



#2022
weather_te = weather_data[weather_data."ID" .== 5555,:]
weather_tf = weather_data[weather_data."ID" .== 3333,:]
weather_te1 = weather_te[Date(2022,07,26) .< weather_te.dt .< Date(2022,07,27),:]
weather_te2 = weather_te[Date(2022,07,27) .< weather_te.dt .< Date(2022,07,28),:]
weather_te3 = weather_te[Date(2022,07,29) .< weather_te.dt .< Date(2022,07,30),:]

weather_tf1 = weather_tf[Date(2022,07,26) .< weather_tf.dt .< Date(2022,07,27),:]
weather_tf2 = weather_tf[Date(2022,07,27) .< weather_tf.dt .< Date(2022,07,28),:]
weather_tf3 = weather_tf[Date(2022,07,29) .< weather_tf.dt .< Date(2022,07,30),:]
weather_tf22 = weather_tf2[DateTime(2022,07,27,12,10) .< weather_tf2.dt .< DateTime(2022,07,27,12,50),:]
weather_tf23 = weather_tf2[DateTime(2022,07,27,12,50) .< weather_tf2.dt .< DateTime(2022,07,27,13,20),:]

#NW site has bad winds (no tripod)

#2021
# weather_te = weather_data[weather_data."ID" .== 999,:]

# weather_te1 = weather_te[weather_te.dt .< Date(2021,09,16), :]

# weather_te2 = weather_te[Date(2021,09,17) .< weather_te.dt .< Date(2021,09,18), :]

# weather_te3 = weather_te[Date(2021,09,18) .< weather_te.dt .< Date(2021,09,21), :]

# weather_te4 = weather_te[Date(2021,09,21) .< weather_te.dt .< Date(2021,09,25), :]
# weather_te_11 = weather_te1[weather_te1.dt .> DateTime(2021,9,15,18),:]


scatter(weather_te.dt,weather_te.WDIR)
# scatter(weather_te.dt,weather_te.WDIR)



using PyPlot
pygui(true)

windrose = pyimport("windrose")

ax2 = windrose.WindroseAxes.from_ax()
ax2.bar(weather_tf23.WDIR, weather_tf23.WSPD, normed=true, opening = 0.8)
ax2.set_legend()
ax2.set_title("Vaisala Observed Winds\nTwin Creeks 2022-07-27 12:50-13:20")
# show()

gr()

#TODO - calculate 10 minute average concentrations

function x_min_medians(df1, df2, min=10)

    dtmax = Dates.DateTime(0)

    if df1.dt[end] > df2.dt[end]
        dtmax = df1.dt[end]
    else
        dtmax = df2.dt[end]
    end

    dtmin = DateTime(0)

    if df1.dt[1] < df2.dt[1]
        dtmin = df1.dt[1]
    else
        dtmin = df2.dt[1]
    end


    dtt = dtmin
    t1_med_10 = []
    t2_med_10 = []

    while dtt < dtmax
        t1d = df1[ dtt .< df1.dt .< dtt + Dates.Minute(min), :]
        t2d = df2[ dtt .< df2.dt .< dtt + Dates.Minute(min), :]
        if length(eachrow(t1d)) > 5
            push!(t1_med_10,[median(t1d.xch4),median(t1d.xco2),median(t1d.xco),
                    dtt,median(t1d.column_luft),median(t1d.column_h2o),median(t1d.xh2o)])
        end
        if length(eachrow(t2d)) > 5
            push!(t2_med_10,[median(t2d.xch4),median(t2d.xco2),median(t2d.xco),
                    dtt,median(t2d.column_luft),median(t2d.column_h2o),median(t2d.xh2o)])
        end
        println(dtt)
        dtt = dtt + Dates.Minute(1)
    end

    t1_df = DataFrame(xch4_t1 = [r[1] for r in t1_med_10],
                        xco2_t1 = [r[2] for r in t1_med_10],
                        xco_t1 = [r[3] for r in t1_med_10],
                        dt = [r[4] for r in t1_med_10],
                        luft_t1 = [r[5] for r in t1_med_10],
                        col_h2o_t1 = [r[6] for r in t1_med_10],
                        xh2o_t1 = [r[7] for r in t1_med_10])

    t2_df = DataFrame(xch4_t2 = [r[1] for r in t2_med_10],
                        xco2_t2 = [r[2] for r in t2_med_10],
                        xco_t2 = [r[3] for r in t2_med_10],
                        dt = [r[4] for r in t2_med_10],
                        luft_t2 = [r[5] for r in t2_med_10],
                        col_h2o_t2 = [r[6] for r in t2_med_10],
                        xh2o_t2 = [r[7] for r in t2_med_10])

    ta_tb = innerjoin(t1_df,t2_df,on=:dt)
    

    return ta_tb

end


function x_min_means(df1, df2, min=10)

    dtmax = Dates.DateTime(0)

    if df1.dt[end] > df2.dt[end]
        dtmax = df1.dt[end]
    else
        dtmax = df2.dt[end]
    end

    dtmin = DateTime(0)

    if df1.dt[1] < df2.dt[1]
        dtmin = df1.dt[1]
    else
        dtmin = df2.dt[1]
    end


    dtt = dtmin
    t1_med_10 = []
    t2_med_10 = []

    while dtt < dtmax
        t1d = df1[ dtt .< df1.dt .< dtt + Dates.Minute(min), :]
        t2d = df2[ dtt .< df2.dt .< dtt + Dates.Minute(min), :]
        if length(eachrow(t1d)) > 5
            push!(t1_med_10,[mean(t1d.xch4),mean(t1d.xco2),mean(t1d.xco),
                    dtt,mean(t1d.column_luft),mean(t1d.column_h2o),mean(t1d.xh2o)])
        end
        if length(eachrow(t2d)) > 5
            push!(t2_med_10,[mean(t2d.xch4),mean(t2d.xco2),mean(t2d.xco),
                    dtt,mean(t2d.column_luft),mean(t2d.column_h2o),mean(t2d.xh2o)])
        end
        println(dtt)
        dtt = dtt + Dates.Minute(min)
    end

    t1_df = DataFrame(xch4_t1 = [Float64(r[1]) for r in t1_med_10],
                        xco2_t1 = [Float64(r[2]) for r in t1_med_10],
                        xco_t1 = [Float64(r[3]) for r in t1_med_10],
                        dt = [r[4] for r in t1_med_10],
                        luft_t1 = [Float64(r[5]) for r in t1_med_10],
                        col_h2o_t1 = [Float64(r[6]) for r in t1_med_10],
                        xh2o_t1 = [Float64(r[7]) for r in t1_med_10])

    t2_df = DataFrame(xch4_t2 = [Float64(r[1]) for r in t2_med_10],
                        xco2_t2 = [Float64(r[2]) for r in t2_med_10],
                        xco_t2 = [Float64(r[3]) for r in t2_med_10],
                        dt = [r[4] for r in t2_med_10],
                        luft_t2 = [Float64(r[5]) for r in t2_med_10],
                        col_h2o_t2 = [Float64(r[6]) for r in t2_med_10],
                        xh2o_t2 = [Float64(r[7]) for r in t2_med_10])

    ta_tb = innerjoin(t1_df,t2_df,on=:dt)
    ta_tb.Δxch4 = Float64.(ta_tb.xch4_t2 .- ta_tb.xch4_t1)

    return ta_tb

end



data_10min = x_min_means(data_te, data_tf,2)

#TODO - 10 minute average wind speed and directions?

# with rigerous wind uncertainties
function x_min_winds(weather_tf, min=10, time_offset = -10)
    metD=[]
    for dt in data_10min.dt 
        shrt = weather_tf[ dt + Minute(time_offset) .<= weather_tf.dt .< dt + Minute(time_offset) + Dates.Minute(min), :]
        met_data = DataFrame()

        if length(eachrow(shrt)) > 2
            met_data.wx = [i.WSPD .* cosd(270 - i.WDIR) for i in eachrow(shrt)]
            # met_data.wy = [i.ws .* sind(i.wvec) for i in eachrow(met_data)]

            met_data.wy = [i.WSPD .* sind( 270 - i.WDIR) for i in eachrow(shrt)]
            std_x = std(met_data.wx)
            std_y = std(met_data.wy)
            mean_x = mean(met_data.wx)
            mean_y = mean(met_data.wy)

            wind_speed = sqrt(mean_x^2 + mean_y^2)
            wind_dir = ( 270 - atand(mean_y,mean_x)) % 360

            std_wind_speed = wind_speed * sqrt((std_x * mean_x^2)/(mean_x^2+mean_y^2) + (std_y * mean_x^2)/(mean_x^2+mean_y^2) )
            std_wind_dir = wind_dir * sqrt( (mean_x^2 * std_x^2)/(mean_x^2+mean_y^2)^2 + (mean_y^2 * std_y^2)/(mean_x^2+mean_y^2)^2 )
            dtas = [wind_speed,wind_dir,std_wind_speed,std_wind_dir,mean_x,mean_y,std_x,std_y,dt]
            push!(metD,dtas)
        else
            dtas = [0,0,0,0,0,0,0,0,dt]
            push!(metD,dtas)
        end

    end
    return DataFrame(wspd = [r[1] for r in metD],
                     wdir = [r[2] for r in metD],
                     wspd_std = [r[3] for r in metD],
                     wdir_std = [r[4] for r in metD],
                     mean_u = [r[5] for r in metD],
                     mean_v = [r[6] for r in metD],
                     u_std = [r[7] for r in metD],
                     v_std = [r[8] for r in metD],
                     dt = [r[9] for r in metD])
end

wind_10min = x_min_winds(weather_tf,5,-10)

#---
#Split into individual days
#2022
data_10min1 = data_10min[Date(2022,7,27) .> data_10min.dt .> Date(2022,7,26),:]
data_10min2 = data_10min[Date(2022,7,28) .> data_10min.dt .> Date(2022,7,27),:]
data_10min3 = data_10min[Date(2022,7,30) .> data_10min.dt .> Date(2022,7,29),:]

wind_10_min1 = wind_10min[Date(2022,07,27) .> wind_10min.dt .> Date(2022,7,26),:]
wind_10_min2 = wind_10min[Date(2022,07,28) .> wind_10min.dt .> Date(2022,7,27),:]
wind_10_min3 = wind_10min[Date(2022,07,30) .> wind_10min.dt .> Date(2022,7,28),:]


#2021
# data_10min1 = data_10min[data_10min.dt .< Date(2021,9,16),:]
# data_10min2 = data_10min[Date(2021,09,18) .> data_10min.dt .> Date(2021,09,17),:]
# data_10min3 = data_10min[Date(2021,9,21) .> data_10min.dt .> Date(2021,9,18),:]
# data_10min4 = data_10min[Date(2021,9,25) .> data_10min.dt .> Date(2021,9,21),:]

# wind_10_min1 = wind_10min[wind_10min.dt .< Date(2021,9,16),:]
# wind_10_min2 = wind_10min[Date(2021,09,18) .> wind_10min.dt .> Date(2021,09,17),:]
# wind_10_min3 = wind_10min[Date(2021,9,21) .> wind_10min.dt .> Date(2021,9,18),:]
# wind_10_min4 = wind_10min[Date(2021,9,25) .> wind_10min.dt .> Date(2021,9,21),:]



#---
# wm_data1 = wm_data[Date(2022,07,27) .> wm_data.dt .> Date(2022,07,26), : ]

#TODO - split along geometric axis based on coordinates of instruments...

#day 1, 2, & 3 te -  42.977563, -81.878397 
#day 1 tf - 42.975244,-81.864879
#day 2 tf -  42.978983,  -81.866603
#day 3 tf -  42.975441, -81.864877

#active face center-ish 42.972624° -81.873350°
using Geodesy

points = DataFrame(lat_te = [42.977563,42.977563,42.977563],
                   lon_te  = [-81.878397,-81.878397,-81.878397],
                    lon_tf = [-81.864879,-81.866603,-81.864877],
                    lat_tf = [42.975244,42.978983,42.975441],
                    lat_WM = [42.972624,42.972624,42.972624],
                    lon_WM = [-81.873350,-81.873350,-81.873350])

points.lla_te = [LLA(r.lat_te,r.lon_te,0) for r in eachrow(points)]
points.lla_tf = [LLA(r.lat_tf,r.lon_tf,0) for r in eachrow(points)]
points.lla_WM = [LLA(r.lat_WM,r.lon_WM,0) for r in eachrow(points)]

utm_from_lla = UTMfromLLA(17, true, wgs84) #Zone 17 for Toronto?

points.coords_te = map(utm_from_lla, points.lla_te)
points.coords_tf = map(utm_from_lla, points.lla_tf)
points.coords_WM = map(utm_from_lla, points.lla_WM)

points.x_te = [row.coords_te.x - 630000 for row in eachrow(points)]
points.y_te = [row.coords_te.y - 4830000 for row in eachrow(points)]

points.x_tf = [row.coords_tf.x - 630000 for row in eachrow(points)]
points.y_tf = [row.coords_tf.y - 4830000 for row in eachrow(points)]

points.x_WM = [row.coords_WM.x - 630000 for row in eachrow(points)]
points.y_WM = [row.coords_WM.y - 4830000 for row in eachrow(points)]

# points.angle_tf = [90 - atand((row.y_tf - row.y_WM) / (row.x_tf - row.x_WM)) + 180 for row in eachrow(points) ]
points.angle_tf = [270 .- atand((row.y_tf - row.y_WM) / (row.x_tf - row.x_WM)) for row in eachrow(points) ]
# points.angle_te = [180 + atand((row.y_te - row.y_WM) / (row.x_te - row.x_WM)) for row in eachrow(points) ]
points.angle_te = [90 .- atand((row.y_te - row.y_WM) / (row.x_te - row.x_WM)) for row in eachrow(points) ]
# points.angle_tf_w = 270 .- points.angle_tf
# points.angle_te_w = 270 .- points.angle_te

#along axis wind strength...
# wind_10_min1.tf_wpar =
wind_10_min1.tf_wpar = [cosd(points.angle_tf[1] .- i.wdir ) .* i.wspd  for i in eachrow(wind_10_min1)]
wind_10_min1.te_wpar = [cosd(points.angle_te[1] .- i.wdir) .* i.wspd for i in eachrow(wind_10_min1)]

wind_10_min2.tf_wpar = [cosd(points.angle_tf[2] .- i.wdir) .* i.wspd for i in eachrow(wind_10_min2)]
wind_10_min2.te_wpar = [cosd(points.angle_te[2].- i.wdir) .* i.wspd for i in eachrow(wind_10_min2)]

wind_10_min3.tf_wpar = [cosd( points.angle_tf[3] .- i.wdir) .* i.wspd for i in eachrow(wind_10_min3)]
wind_10_min3.te_wpar = [cosd( points.angle_te[3] .- i.wdir) .* i.wspd for i in eachrow(wind_10_min3)]


#---
#---
combined_data1 = innerjoin(data_10min1,wind_10_min1,on=:dt)
combined_data1 = combined_data1[combined_data1.wspd .!= 0, :]
# combined_data1 = combined_data1[(combined_data1.tf_wpar ./ combined_data1.te_wpar) .> -100,:]
combined_data2 = innerjoin(data_10min2,wind_10_min2,on=:dt)
for (i,v) in enumerate(eachrow(combined_data2))
    if v.wspd == 0
        combined_data2.wspd[i] = 0.01
    end
end 
combined_data2 = combined_data2[combined_data2.wspd .!= 0, :]
combined_data3 = innerjoin(data_10min3,wind_10_min3,on=:dt)
combined_data3 = combined_data3[combined_data3.wspd .!= 0, :]

#---

# wind shift?

using LsqFit

x = Float64.(combined_data3.wdir)
y = Float64.(combined_data3.xch4_t1)
x1 = Float64.(combined_data2.wdir)
y1 = Float64.(combined_data2.xch4_t2)

# xy = hcat(x, y)

function oneD_Gaussian(x, p)
    amplitude, μ, σ, offset = p
    g = offset .+ amplitude .* exp.( -1/2 .* ((x .- μ) ./ σ).^2 )
    return g[:]
end


p0 = Float64.([0.08, 250, 30, 1.9])
p1 = Float64.([0.02, 143, 10, 1.9])

# Noisy data
# data_noisy = data + 0.2 * randn(size(data))

fit = LsqFit.curve_fit(oneD_Gaussian, x, y, p0)

newx = 200:340
newx1 = 110:200

#try weighting?
wt = [1 / ((maximum(combined_data3.xch4_t2) - i)^2) for i in combined_data3.xch4_t2]
wt2 = [1 / ((maximum(combined_data3.tf_wpar) - i)^2) for i in combined_data3.tf_wpar]
wt1 = [1 / (.2*(maximum(combined_data2.xch4_t1) - i)^2) for i in combined_data2.xch4_t1]
fit = LsqFit.curve_fit(oneD_Gaussian, x, y, wt, p0)
fit2 = LsqFit.curve_fit(oneD_Gaussian, x, y, wt2, p0)
fit1 = LsqFit.curve_fit(oneD_Gaussian, x1, y1, wt1, p1)
histogram2d(combined_data3.wdir,combined_data3.xch4_t2,show_empty_bins=true,bins=4)
vline!([points.angle_tf[1]])


scatter(combined_data3.wdir,combined_data3.xch4_t2,label="3 minute average tf",marker_z=combined_data3.tf_wpar)
vline!([points.angle_tf[3]])

println(fit1.param)

newy = oneD_Gaussian(newx, fit.param)
oneD_Gaussian(combined_data3.wdir, fit4.param) == oneD_Gaussian(combined_data3.wdir,fit.param)
newyy .- combined_data3.wdir
newy1 = oneD_Gaussian(newx1, fit1.param)

scatter!(combined_data3.wdir,newyy)


scatter!(combined_data2.wdir,combined_data2.xch4_t1)
#---
scatter(combined_data3.tf_wpar .* combined_data3.wspd, combined_data3.xch4_t2 .- combined_data3.xch4_t1)

scatter(wind_10_min1.te_wpar, data_10min1.xch4_t1)

#---

scatter(points.x_te,points.y_te,label="te")
scatter!(points.x_tf,points.y_tf,label="tf")
scatter!(points.x_WM,points.y_WM,label="Active Face")

#---
pltpath = "/home/lawson/Data/EM27SUN/plots/Landfills/"
#2022-07-26/2021-09-15
scatter(data_10min1.dt.+Minute(10),data_10min1.xch4_t1,label="te - downwind\n10min average",
        c=:blue)
scatter!(data_10min1.dt.+Minute(10),data_10min1.xch4_t2,label="tf - upwind\n10min average", 
        # legend=:topleft,
        c=:orange
        )
scatter!(data_te1.dt,data_te1.xch4,msw=0,markersize=3,α=.3,label="te - all data",c=:blue)
scatter!(data_tf1.dt,data_tf1.xch4,msw=0,markersize=3,α=.3,label="tf - all data",c=:orange)
# ch4_diff = round(percentile(data_tf1.xch4,95) - percentile(data_te1.xch4,0.5),digits=3)
ch4_diff = round(1.915-1.908,digits=3)


# plot!(xlims=Dates.value.([DateTime(2022,07,26,13,40),DateTime(2022,07,26,23,45)]))
# plot!(title="Twin Creeks 2022-07-26",ylabel="XCH₄ (ppm)\n10-min median",xrot=-11)
plot!(title="Petrolia 2021-09-15",ylabel="XCH₄ (ppm)",xrot=-11)

hline!([1.908, 1.915],label="ΔXCH₄\n($(ch4_diff)ppm)",c=:black,style=:dot,legend=:bottomleft)
# hline!([percentile(data_tf1.xch4,95) , percentile(data_te1.xch4,0.5)],label="ΔXCH₄\n($(ch4_diff)ppm)",
# c=:black,style=:dot,legend=:right)

png(pltpath*"Petrolia20210915_xch4")
# png(pltpath*"TwinCreeks20220926_xch4")

#2022-07-27 / 2021-09-17
scatter(data_10min2.dt.+Minute(10),data_10min2.xch4_t1,label="te - 10min average",c=:blue)
scatter!(data_10min2.dt.+Minute(10),data_10min2.xch4_t2,label="tf - 10min average",c=:orange 
        # legend=:topleft
        )
scatter(data_te2.dt,data_te2.xch4,msw=0,markersize=3,α=.3,label="te - all data",c=:blue)
scatter!(data_tf2.dt,data_tf2.xch4,msw=0,markersize=3,α=.3,label="tf - all data",c=:orange)
vline!([DateTime(2022,07,27,12,10),DateTime(2022,07,27,12,50)],label="Winds from SE")                
vline!([DateTime(2022,07,27,12,55),DateTime(2022,07,27,13,25)],label="Winds from S")                
# plot!(title="Petrolia 2021-09-17",ylabel="XCH₄ (ppm)",xrot=-11)
plot!(title="Twin Creeks 2022-07-27",ylabel="XCH₄ (ppm)",xrot=-11)
# ch4_diff = round(percentile(data_te2.xch4,99.5) - percentile(data_te2.xch4,.5),digits=3)

data_tf22 = data_tf2[ DateTime(2022,07,27,12,10) .< data_tf2.dt .< DateTime(2022,07,27,12,50), :]
data_te22 = data_te2[ DateTime(2022,07,27,12,10) .< data_te2.dt .< DateTime(2022,07,27,12,50), :]

ch4_diff = abs(round(percentile(data_te22.xch4,99.5) - percentile(data_tf22.xch4,50),digits=3))
# hline!([percentile(data_te2.xch4,99.5),percentile(data_te2.xch4,.5)])
hline!([percentile(data_te22.xch4,99.5),percentile(data_tf22.xch4,50)],
c=:black,linestyle=:dash,label="ΔXCH₄\n($(ch4_diff)ppm)",legend=:topright)
# png(pltpath*"Petrolia20210917")
png(pltpath*"TwinCreeks20220727")

#2022-07-29 / 2021-09-20
scatter(data_10min3.dt,data_10min3.xch4_t1,label="te - 10min average",c=:blue)
scatter!(data_10min3.dt,data_10min3.xch4_t2,label="tf - 10min average", 
        # legend=:topleft
        c=:orange)

scatter!(data_te3.dt,data_te3.xch4,msw=0,markersize=3,α=.3,label="te - all data",c=:blue)
scatter!(data_tf3.dt,data_tf3.xch4,msw=0,markersize=3,α=.3,label="tf - all data",c=:orange)
        
plot!(title="Twin Creeks 2022-07-29",ylabel="XCH₄ (ppm)",xrot=-11)
# plot!(title="Petrolia 2021-09-20",ylabel="XCH₄ (ppm)",xrot=-11)
ch4_diff = round(percentile(data_tf3.xch4,97) - percentile(data_te3.xch4,0.5),digits=3)
hline!([percentile(data_tf3.xch4,97), percentile(data_te3.xch4,0.5)],label="ΔXCH₄\n($(ch4_diff)ppm",
linestyle=:dot,c=:black,xlims=Dates.value.([DateTime(2022,07,29,11,40),DateTime(2022,07,30,1)]),legend=:right)
# png(pltpath*"Petrolia20210920_xch4")
png(pltpath*"TwinCreeks20220729_xch4")

#2021-09-24
scatter(data_10min4.dt.+Minute(10),data_10min4.xch4_t1,label="te - 10min average",c=:blue)
scatter!(data_10min4.dt.+Minute(10),data_10min4.xch4_t2,label="tf - 10min average",c=:orange 
        # legend=:topleft
        )
scatter!(data_te4.dt,data_te4.xch4,msw=0,markersize=3,α=.3,label="te - all data",c=:blue)
scatter!(data_tf4.dt,data_tf4.xch4,msw=0,markersize=3,α=.3,label="tf - all data",c=:orange)
                
plot!(title="Twin Creeks 2021-09-24",ylabel="XCH₄ (ppm)",xrot=-11)
# ch4_diff = round(percentile(data_te2.xch4,99.5) - percentile(data_te2.xch4,.5),digits=3)
ch4_diff = abs(round(percentile(data_tf4.xch4,5) - percentile(data_te4.xch4,95),digits=3))
# hline!([percentile(data_te2.xch4,99.5),percentile(data_te2.xch4,.5)])
hline!([percentile(data_tf4.xch4,5),percentile(data_te4.xch4,95)],
c=:black,style=:dot,label="ΔXCH₄\n($(ch4_diff)ppm)",legend=:topleft)
png(pltpath*"TwinCreeks20210924_xch4")


#---
#XCO₂

#2022-07-29
scatter(data_10min3.dt,data_10min3.xco2_t1,label="te")
scatter!(data_10min3.dt,data_10min3.xco2_t2,label="tf", 
        # legend=:topleft
        )
plot!(title="Twin Creeks 2022-07-29",ylabel="XCO₂ (ppm)",xrot=-11)
png(pltpath*"TwinCreeks20220729_xco2")
#---
# scatter!(data_10min3.dt,data_10min3.xch4_t2 .- data_10min3.xch4_t1,label="tf-te")

#---

scatter(data_10min1.dt,data_10min1.xch4_t1,label="te")
scatter!(data_10min1.dt,data_10min1.xch4_t2,label="tf")

scatter!(twinx(),wind_10_min1.dt,wind_10_min1.tf_wpar,msw=0,marker_size=4,α=0.5,c=:blue)

#---
combined_data1 = innerjoin(data_10min1,wind_10_min1,on=:dt)
combined_data1 = combined_data1[combined_data1.wspd .!= 0, :]
combined_data2 = innerjoin(data_10min2,wind_10_min2,on=:dt)
combined_data2 = combined_data2[combined_data2.wspd .!= 0, :]
combined_data3 = innerjoin(data_10min3,wind_10_min3,on=:dt)
combined_data3 = combined_data3[combined_data3.wspd .!= 0, :]
combined_data4 = innerjoin(data_10min4,wind_10_min4,on=:dt)
combined_data4 = combined_data4[combined_data4.wspd .!= 0, :]

#---

scatter(combined_data1.mean_u,combined_data1.xch4_t2.-combined_data1.xch4_t1)

scatter(combined_data1.dt,combined_data1.xch4_t2.-combined_data1.xch4_t1)
scatter!(combined_data1.dt,abs.(combined_data1.tf_wpar),msw=0,marker_size=4,α=0.5,c=:red)


pltpath = "/home/lawson/Data/EM27SUN/plots/Landfills/wind_shift/"

for i in -60:5:60
    wind = x_min_winds(weather_tf3,15,i)
    wind.tf_wpar = [cosd(j.wdir .- points.angle_tf[1]) * j.wspd for j in eachrow(wind)]
    combined_data_x = innerjoin(data_10min3,wind,on=:dt)
    combined_data_x = combined_data_x[combined_data_x.wspd .> 0,:]
    scatter(combined_data_x.mean_v, combined_data_x.xch4_t2 .- combined_data_x.xch4_t1,
            title = "$(i) minute wind offset")
    png(pltpath*"$(i)_min_wind_offset")
end
    

filt_data = combined_data1[combined_data1.tf_wpar .> 0, :]
using GLM

scatter(filt_data.tf_wpar,filt_data.Δxch4)

lr1 = lm(@formula(Δxch4 ~  tf_wpar),filt_data)
pred = predict(lr1,filt_data)
plot!(filt_data.tf_wpar,pred)
        


#TODO - calculate max/min upwind/downwind concentrations.

#Implement Gaussian plume model?

#TODO - estimate emissions

#calculate distance between active face and tf. 

#Twin Creeks
x_dist = points.x_tf[3] - points.x_WM[3]
y_dist = points.y_tf[3] - points.y_WM[3]
dist = sqrt(x_dist^2 + y_dist^2)
dist = dist


#Petrolia


#------


dist = 720
α1 = 0.22
σ_y = α1 * dist * (1 + 0.0001*dist)^(-.5) #horizontal dispersion
Ueff = 2.6
Ueff = combined_data1.tf_wpar
# Δxch4 = (combined_data3.xch4_t2.-combined_data3.xch4_t1) .* 1e-6 #ppm avg peak of histogram
Δxch4 = (combined_data1.xch4_t2.-combined_data1.xch4_t1) .* 1e-6 #ppm avg peak of histogram
# Δxch4 = (.007 - 0.000) * 1e-6 #in ppm
ch4_ratio = 16.043 / 6.022143e23 #grams per molecule
# column_abundance = maximum(combined_data1.luft_t2.-combined_data1.col_h2o_t2) #molecules per cm^-2, need to subtract water column
column_abundance = mean(combined_data1.luft_t1.-combined_data1.col_h2o_t1) #molecules per cm^-2, need to subtract water column..
column_abundance_molecules_m =column_abundance*100^2 #molecules per cm^-2 
ΔΩ = Δxch4 * ch4_ratio * column_abundance_molecules_m  # we need g/m²
Q₀ = sqrt(2*pi) * σ_y .* Ueff .* ΔΩ #g/s
Q_kgh = Q₀/1000 *3600 # (1/100000)
Q_kgd = Q_kgh*24 /1.2# (1/100000)


histogram(Q_kgd)

~
histogram(combined_data3.xch4_t2.-combined_data3.xch4_t1,bins=25,label="20220729",α=.4)
histogram!(combined_data2.xch4_t1.-combined_data2.xch4_t2,bins=5,label="20220727",α=.3)
histogram!(combined_data1.xch4_t2.-combined_data1.xch4_t1,bins=20,label="20220726",α=.2)


histogram(combined_data3.wspd,bins=25)


data_te3 = data_te3[data_te3.xch4 .> percentile(data_te3.xch4, 0.02),:]

histogram(data_tf3.xch4,α=0.85,label="tf - North",c=:orange)
histogram!(data_te3.xch4,label="te - South",c=:blue,
xlabel="XCH₄",ylabel="# of observations",α=.35)
# plot!(title="20220729 Twin Creeks")
plot!(title="2021-09-20 Petrolia")
# png(pltpath*"20220729_xch4_hist")
png(pltpath*"20210920_xch4_hist")

histogram(data_te2.xch4,label="te - downwind",c=:blue,
xlabel="XCH₄",ylabel="# of observations",bins=20)
histogram!(data_tf2.xch4,α=0.85,label="tf - upwind",c=:orange,)
# plot!(title="20220727 Twin Creeks")
plot!(title="2021-09-17 Petrolia")
# png(pltpath*"20220727_xch4_hist")
png(pltpath*"20210917_xch4_hist")


histogram(data_tf1.xch4,α=1,label="tf - upwind",c=:orange)
histogram!(data_te1.xch4,label="te - downwindish",c=:blue,
xlabel="XCH₄",ylabel="# of observations",α=.45)
plot!(title="2021-09-15 Petrolia")
# png(pltpath*"20220726_xch4_hist")
png(pltpath*"20210915_xch4_hist")

histogram(data_te4.xch4,label="te - downwind",c=:blue,
xlabel="XCH₄",ylabel="# of observations",α=1)
histogram!(data_tf4.xch4,α=.85,label="tf - upwind",c=:orange)

plot!(title="2021-09-24 Twin Creeks")
# png(pltpath*"20220726_xch4_hist")
png(pltpath*"20210924_xch4_hist")

~hodl
# nc_tb_path = "/home/lawson/HDD/Data/em27_ncfiles/tb/cumulative/"
# nc_tb_path2 = "/home/lawson/HDD/Data/em27_ncfiles/tb/UTSG/"
# nc_tb_path3 = "/home/lawson/HDD/Data/em27_ncfiles/tb/UTSC/"
# nc_paths = glob("tb*.nc",nc_tb_path)
# pths2 = glob("tb*.nc",nc_tb_path2)
# append!(nc_paths,append!(pths2,glob("tb*.nc",nc_tb_path3)))
# #make a list of NCDatasets from netCDF's
# datasets = []
# for path in nc_paths
#     ds = Dataset(path)
#     push!(datasets,ds)
# end
# dataframes = []
# for ds in datasets
#     tmp_df = DataFrame()
#     for var in variables
#         tmp_df[!,"$var"] = ds[var][1:end]
#     end
#     push!(dataframes,tmp_df)
# end
# data2_df  = dataframes[1]
# for df in dataframes[2:end]
#     data2_df = vcat(data2_df,df)
# end


# data2_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data2_df) ]
# data2_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data2_df) ]
# data2_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data2_df)]
# data2_df.ins = ["tb" for i in eachrow(data2_df)]
# unique!(data2_df, :dt)

# data1_df = vcat(data1_df,data2_df)

# #now do the same for tc

# nc_tc_path = "/home/lawson/HDD/Data/em27_ncfiles/tc/cumulative/"
# nc_tc_path2 = "/home/lawson/HDD/Data/em27_ncfiles/tc/UTM/"
# nc_tc_path3 = "/home/lawson/HDD/Data/em27_ncfiles/tc/ECCC_Downsview/"
# nc_paths = glob("tc*.nc",nc_tc_path)

# pths2 = glob("tc*.nc",nc_tc_path2)
# append!(pths2,glob("tc*.nc",nc_tc_path3))
# append!(nc_paths,pths2)
# #make a list of NCDatasets from netCDF's
# datasets = []
# for path in nc_paths
#     ds = Dataset(path)
#     push!(datasets,ds)
# end
# dataframes = []
# for ds in datasets
#     tmp_df = DataFrame()
#     for var in variables
#         tmp_df[!,"$var"] = ds[var][1:end]
#     end
#     push!(dataframes,tmp_df)
# end
# data3_df  = dataframes[1]
# for df in dataframes[2:end]
#     data3_df = vcat(data3_df,df)
# end


# data3_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data3_df) ]
# data3_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data3_df) ]
# data3_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data3_df)]
# data3_df.ins = ["tc" for i in eachrow(data3_df)]
# unique!(data3_df, :dt)

# data1_df = vcat(data1_df,data3_df)

# #likewise for td ...

# nc_td_path = "/home/lawson/HDD/Data/em27_ncfiles/td/cumulative/"
# nc_td_path2 = "/home/lawson/HDD/Data/em27_ncfiles/td/ECCC_Downsview/"
# nc_paths = glob("td*.nc",nc_td_path)
# pths2 = glob("td*.nc",nc_td_path2)
# append!(nc_paths,pths2)
# #make a list of NCDatasets from netCDF's
# datasets = []
# for path in nc_paths
#     ds = Dataset(path)
#     push!(datasets,ds)
# end
# dataframes = []
# for ds in datasets
#     tmp_df = DataFrame()
#     for var in variables
#         tmp_df[!,"$var"] = ds[var][1:end]
#     end
#     push!(dataframes,tmp_df)
# end
# data4_df  = dataframes[1]
# for df in dataframes[2:end]
#     data4_df = vcat(data4_df,df)
# end


# data4_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data4_df) ]
# data4_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data4_df) ]
# data4_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data4_df)]
# data4_df.ins = ["td" for i in eachrow(data4_df)]
# unique!(data4_df, :dt)

# data1_df = vcat(data1_df,data4_df)


# utsg_data = data1_df[-79.4 .< data1_df.long .< -79.398, :]
# utsg_data = utsg_data[utsg_data.flag .== 0 , :]

# #split UTSC

# utsc_data = data1_df[(43.78 .< data1_df.lat .< 43.785) .&
#                         (-79.19 .< data1_df.long .< -79.18),  :]
# utsc_data = utsc_data[utsc_data.flag .== 0 , :]


# #split UTM data

# utm_data = data1_df[-79.67 .< data1_df.long .< -79.66, :]
# utm_data = utm_data[utm_data.flag .== 0, : ]

# #split ECCC Downsview data
# dow_data = data1_df[ (-79.47 .< data1_df.long .< -79.465) .& (43.77 .< data1_df.lat .< 43.785) ,:]
# dow_data = dow_data[(dow_data.flag .== 0) .| (dow_data.flag .== 18),:]


# #colour by instrument for plotting ease...
# dic = Dict([("ta",:red),("tb",:orange),("tc",:blue),("td",:lightblue)])
# utsg_data.color = [dic[r] for r in utsg_data.ins]

# ta_data = utsg_data[utsg_data.ins .== "ta", :]
# tb_data = utsg_data[utsg_data.ins .== "tb", :]
# tc_data = utsg_data[utsg_data.ins .== "tc", :]
# td_data = utsg_data[utsg_data.ins .== "td", :]


# #plot xch4
# scatter(ta_data.dt,ta_data.xch4,msw=0,c=ta_data.color,label="ta",
#         xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# scatter!(tb_data.dt,tb_data.xch4,msw=0,c=tb_data.color,label="tb")
# scatter!(tc_data.dt,tc_data.xch4,msw=0,c=tc_data.color,label="tc")
# scatter!(td_data.dt,td_data.xch4,msw=0,c=td_data.color,label="td")
# plot!(xlabel="Date",ylabel="XCH₄ (ppm)",title="UTSG XCH₄ Observations",
#             legend=:bottomright)
# pltpath = "/home/lawson/HDD/Data/"
# png(pltpath*"UTSGXCH4")

# #plotxco2
# scatter(ta_data.dt,ta_data.xco2,msw=0,c=ta_data.color,label="ta",
#             xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# scatter!(tb_data.dt,tb_data.xco2,msw=0,c=tb_data.color,label="tb")
# scatter!(tc_data.dt,tc_data.xco2,msw=0,c=tc_data.color,label="tc")
# scatter!(td_data.dt,td_data.xco2,msw=0,c=td_data.color,label="td")
# plot!(xlabel="Date",ylabel="XCO₂ (ppm)",title="UTSG XCO₂ Observations",
#             legend=:bottomright)
# png(pltpath*"UTSGXCO2")


# #plotxco
# scatter(ta_data.dt,ta_data.xco,msw=0,c=ta_data.color,label="ta",
#         xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# scatter!(tb_data.dt,tb_data.xco,msw=0,c=tb_data.color,label="tb")
# scatter!(tc_data.dt,tc_data.xco,msw=0,c=tc_data.color,label="tc")
# scatter!(td_data.dt,td_data.xco,msw=0,c=td_data.color,label="td")
# plot!(xlabel="Date",ylabel="XCO (ppm)",title="UTSG XCO Observations",
#             legend=:topleft)
# png(pltpath*"UTSGXCO")



# #Plots for UTSC
# #colour by instrument for plotting ease...
# dic = Dict([("ta",:red),("tb",:orange),("tc",:blue),("td",:lightblue)])
# utsc_data.color = [dic[r] for r in utsc_data.ins]

# ta_data = utsc_data[utsc_data.ins .== "ta", :]
# tb_data = utsc_data[utsc_data.ins .== "tb", :]
# tc_data = utsc_data[utsc_data.ins .== "tc", :]
# td_data = utsc_data[utsc_data.ins .== "td", :]


# #plot xch4
# # scatter(ta_data.dt,ta_data.xch4,msw=0,c=ta_data.color,label="ta")
# scatter(tb_data.dt,tb_data.xch4,msw=0,c=tb_data.color,label="tb",
#             xlims = (minimum(data1_df.dt),maximum(data1_df.dt)) )
# # scatter!(tc_data.dt,tc_data.xch4,msw=0,c=tc_data.color,label="tc")
# # scatter!(td_data.dt,td_data.xch4,msw=0,c=td_data.color,label="td")
# plot!(xlabel="Date",ylabel="XCH₄ (ppm)",title="UTSC XCH₄ Observations",
#             legend=:topleft)
# pltpath = "/home/lawson/HDD/Data/"
# png(pltpath*"UTSCXCH4")

# #plotxco2
# # scatter(ta_data.dt,ta_data.xco2,msw=0,c=ta_data.color,label="ta")
# scatter(tb_data.dt,tb_data.xco2,msw=0,c=tb_data.color,label="tb",
#     xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# # scatter!(tc_data.dt,tc_data.xco2,msw=0,c=tc_data.color,label="tc")
# # scatter!(td_data.dt,td_data.xco2,msw=0,c=td_data.color,label="td")
# plot!(xlabel="Date",ylabel="XCO₂ (ppm)",title="UTSC XCO₂ Observations",
#             legend=:left)
# png(pltpath*"UTSCXCO2")


# #plotxco
# # scatter(ta_data.dt,ta_data.xco,msw=0,c=ta_data.color,label="ta")
# scatter(tb_data.dt,tb_data.xco,msw=0,c=tb_data.color,label="tb",
#             xlims = (minimum(data1_df.dt),maximum(data1_df.dt)) )
# # scatter!(tc_data.dt,tc_data.xco,msw=0,c=tc_data.color,label="tc")
# # scatter!(td_data.dt,td_data.xco,msw=0,c=td_data.color,label="td")
# plot!(xlabel="Date",ylabel="XCO (ppm)",title="UTSC XCO Observations",
#             legend=:left)
# png(pltpath*"UTSCXCO")


# #Plots for UTM
# #colour by instrument for plotting ease...
# utm_data.color = [dic[r] for r in utm_data.ins]

# ta_data = utm_data[utm_data.ins .== "ta", :]
# tb_data = utm_data[utm_data.ins .== "tb", :]
# tc_data = utm_data[utm_data.ins .== "tc", :]
# td_data = utm_data[utm_data.ins .== "td", :]


# #plot xch4
# # scatter(ta_data.dt,ta_data.xch4,msw=0,c=ta_data.color,label="ta")
# # scatter(tb_data.dt,tb_data.xch4,msw=0,c=tb_data.color,label="tb")
# scatter(tc_data.dt,tc_data.xch4,msw=0,c=tc_data.color,label="tc",
#         xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# # scatter!(td_data.dt,td_data.xch4,msw=0,c=td_data.color,label="td")
# plot!(xlabel="Date",ylabel="XCH₄ (ppm)",title="UTM XCH₄ Observations",
#             legend=:topleft)
# pltpath = "/home/lawson/HDD/Data/"
# png(pltpath*"UTM_XCH4")

# #plotxco2
# # scatter(ta_data.dt,ta_data.xco2,msw=0,c=ta_data.color,label="ta")
# # scatter(tb_data.dt,tb_data.xco2,msw=0,c=tb_data.color,label="tb")
# scatter(tc_data.dt,tc_data.xco2,msw=0,c=tc_data.color,label="tc",
#             xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# # scatter!(td_data.dt,td_data.xco2,msw=0,c=td_data.color,label="td")
# plot!(xlabel="Date",ylabel="XCO₂ (ppm)",title="UTM XCO₂ Observations",
#             legend=:left)
# png(pltpath*"UTM_XCO2")


# #plotxco
# # scatter(ta_data.dt,ta_data.xco,msw=0,c=ta_data.color,label="ta")
# # scatter(tb_data.dt,tb_data.xco,msw=0,c=tb_data.color,label="tb")
# scatter(tc_data.dt,tc_data.xco,msw=0,c=tc_data.color,label="tc",
#             xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# # scatter!(td_data.dt,td_data.xco,msw=0,c=td_data.color,label="td")
# plot!(xlabel="Date",ylabel="XCO (ppm)",title="UTM XCO Observations",
#             legend=:left)
# png(pltpath*"UTM_XCO")



# #Plots for ECCC Downsview
# #colour by instrument for plotting ease...
# dow_data.color = [dic[r] for r in dow_data.ins]

# ta_data = dow_data[dow_data.ins .== "ta", :]
# tb_data = dow_data[dow_data.ins .== "tb", :]
# tc_data = dow_data[dow_data.ins .== "tc", :]
# td_data = dow_data[dow_data.ins .== "td", :]


# #plot xch4
# # scatter(ta_data.dt,ta_data.xch4,msw=0,c=ta_data.color,label="ta")
# # scatter(tb_data.dt,tb_data.xch4,msw=0,c=tb_data.color,label="tb")
# scatter(td_data.dt,td_data.xch4,msw=0,c=td_data.color,label="td",
#             xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# scatter!(tc_data.dt,tc_data.xch4,msw=0,c=tc_data.color,label="tc")
# plot!(xlabel="Date",ylabel="XCH₄ (ppm)",title="Downsview XCH₄ Observations",
#             legend=:left)
# pltpath = "/home/lawson/HDD/Data/"
# png(pltpath*"DOW_XCH4")

# #plotxco2
# # scatter(ta_data.dt,ta_data.xco2,msw=0,c=ta_data.color,label="ta")
# # scatter(tb_data.dt,tb_data.xco2,msw=0,c=tb_data.color,label="tb")
# scatter(td_data.dt,td_data.xco2,msw=0,c=td_data.color,label="td",
#     xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# scatter!(tc_data.dt,tc_data.xco2,msw=0,c=tc_data.color,label="tc")
# plot!(xlabel="Date",ylabel="XCO₂ (ppm)",title="Downsview XCO₂ Observations",
#             legend=:left)
# png(pltpath*"DOW_XCO2")


# #plotxco
# # scatter(ta_data.dt,ta_data.xco,msw=0,c=ta_data.color,label="ta")
# # scatter(tb_data.dt,tb_data.xco,msw=0,c=tb_data.color,label="tb")
# scatter(td_data.dt,td_data.xco,msw=0,c=td_data.color,label="td",
#         xlims = (minimum(data1_df.dt),maximum(data1_df.dt)) )
# scatter!(tc_data.dt,tc_data.xco,msw=0,c=tc_data.color,label="tc")
# plot!(xlabel="Date",ylabel="XCO (ppm)",title="UTM XCO Observations",
#             legend=:topleft)
# png(pltpath*"DOW_XCO")



# scatter(utsg_data.dt,utsg_data.xch4,msw=0.001,c=utsg_data.color)


# CSV.write("/home/lawson/HDD/Data/EM27/utsg.csv",utsg_data)
# CSV.write("/home/lawson/HDD/Data/EM27/utsc.csv",utsc_data)
# CSV.write("/home/lawson/HDD/Data/EM27/utm_data.csv",utm_data)
# CSV.write("/home/lawson/HDD/Data/EM27/dow_data.csv",dow_data)
