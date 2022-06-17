using NCDatasets
using DataFrames
using CSV
using Dates
using Plots
using Glob

#path to netCDF's
nc_path = "/home/lawson/HDD/Data/em27_ncfiles/ta/cumulative/"
nc_pth2 = "/home/lawson/HDD/Data/em27_ncfiles/ta/UTSG/"
#all .nc files in the directory
nc_paths = glob("ta*.nc",nc_path)
append!(nc_paths,glob("ta*.nc",nc_pth2))
#make a list of NCDatasets from netCDF's
datasets = []
for path in nc_paths
    ds = Dataset(path)
    push!(datasets,ds)
end

#list of variables to extract to dataframe
variables  = ["xch4", "xco2", "xco", "flag","year","day",
            "hour","xluft","tins","xluft_error","o2_7885_sg",
            "pout","pins","lat","long"]

#make a list of dataframes to dump NCDataset data into
dataframes = []
for ds in datasets
    tmp_df = DataFrame()
    for var in variables
        tmp_df[!,"$var"] = ds[var][1:end]
    end
    push!(dataframes,tmp_df)
end

#make it all into single dataframe
data1_df  = dataframes[1]
for df in dataframes[2:end]
    data1_df = vcat(data1_df,df)
end


data1_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data1_df) ]
data1_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data1_df) ]
data1_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data1_df)]
data1_df.ins = ["ta" for i in eachrow(data1_df)]
unique!(data1_df,:dt)



#now do tb data...
nc_tb_path = "/home/lawson/HDD/Data/em27_ncfiles/tb/cumulative/"
nc_tb_path2 = "/home/lawson/HDD/Data/em27_ncfiles/tb/UTSG/"
nc_tb_path3 = "/home/lawson/HDD/Data/em27_ncfiles/tb/UTSC/"
nc_paths = glob("tb*.nc",nc_tb_path)
pths2 = glob("tb*.nc",nc_tb_path2)
append!(nc_paths,append!(pths2,glob("tb*.nc",nc_tb_path3)))
#make a list of NCDatasets from netCDF's
datasets = []
for path in nc_paths
    ds = Dataset(path)
    push!(datasets,ds)
end
dataframes = []
for ds in datasets
    tmp_df = DataFrame()
    for var in variables
        tmp_df[!,"$var"] = ds[var][1:end]
    end
    push!(dataframes,tmp_df)
end
data2_df  = dataframes[1]
for df in dataframes[2:end]
    data2_df = vcat(data2_df,df)
end


data2_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data2_df) ]
data2_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data2_df) ]
data2_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data2_df)]
data2_df.ins = ["tb" for i in eachrow(data2_df)]
unique!(data2_df, :dt)

data1_df = vcat(data1_df,data2_df)

#now do the same for tc

nc_tc_path = "/home/lawson/HDD/Data/em27_ncfiles/tc/cumulative/"
nc_tc_path2 = "/home/lawson/HDD/Data/em27_ncfiles/tc/UTM/"
nc_tc_path3 = "/home/lawson/HDD/Data/em27_ncfiles/tc/ECCC_Downsview/"
nc_paths = glob("tc*.nc",nc_tc_path)

pths2 = glob("tc*.nc",nc_tc_path2)
append!(pths2,glob("tc*.nc",nc_tc_path3))
append!(nc_paths,pths2)
#make a list of NCDatasets from netCDF's
datasets = []
for path in nc_paths
    ds = Dataset(path)
    push!(datasets,ds)
end
dataframes = []
for ds in datasets
    tmp_df = DataFrame()
    for var in variables
        tmp_df[!,"$var"] = ds[var][1:end]
    end
    push!(dataframes,tmp_df)
end
data3_df  = dataframes[1]
for df in dataframes[2:end]
    data3_df = vcat(data3_df,df)
end


data3_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data3_df) ]
data3_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data3_df) ]
data3_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data3_df)]
data3_df.ins = ["tc" for i in eachrow(data3_df)]
unique!(data3_df, :dt)

data1_df = vcat(data1_df,data3_df)

#likewise for td ...

nc_td_path = "/home/lawson/HDD/Data/em27_ncfiles/td/cumulative/"
nc_td_path2 = "/home/lawson/HDD/Data/em27_ncfiles/td/ECCC_Downsview/"
nc_paths = glob("td*.nc",nc_td_path)
pths2 = glob("td*.nc",nc_td_path2)
append!(nc_paths,pths2)
#make a list of NCDatasets from netCDF's
datasets = []
for path in nc_paths
    ds = Dataset(path)
    push!(datasets,ds)
end
dataframes = []
for ds in datasets
    tmp_df = DataFrame()
    for var in variables
        tmp_df[!,"$var"] = ds[var][1:end]
    end
    push!(dataframes,tmp_df)
end
data4_df  = dataframes[1]
for df in dataframes[2:end]
    data4_df = vcat(data4_df,df)
end


data4_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data4_df) ]
data4_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data4_df) ]
data4_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data4_df)]
data4_df.ins = ["td" for i in eachrow(data4_df)]
unique!(data4_df, :dt)

data1_df = vcat(data1_df,data4_df)


utsg_data = data1_df[-79.4 .< data1_df.long .< -79.398, :]
utsg_data = utsg_data[utsg_data.flag .== 0 , :]

#split UTSC

utsc_data = data1_df[(43.78 .< data1_df.lat .< 43.785) .&
                        (-79.19 .< data1_df.long .< -79.18),  :]
utsc_data = utsc_data[utsc_data.flag .== 0 , :]


#split UTM data

utm_data = data1_df[-79.67 .< data1_df.long .< -79.66, :]
utm_data = utm_data[utm_data.flag .== 0, : ]

#split ECCC Downsview data
dow_data = data1_df[ (-79.47 .< data1_df.long .< -79.465) .& (43.77 .< data1_df.lat .< 43.785) ,:]
dow_data = dow_data[(dow_data.flag .== 0) .| (dow_data.flag .== 18),:]


#colour by instrument for plotting ease...
dic = Dict([("ta",:red),("tb",:orange),("tc",:blue),("td",:lightblue)])
utsg_data.color = [dic[r] for r in utsg_data.ins]

ta_data = utsg_data[utsg_data.ins .== "ta", :]
tb_data = utsg_data[utsg_data.ins .== "tb", :]
tc_data = utsg_data[utsg_data.ins .== "tc", :]
td_data = utsg_data[utsg_data.ins .== "td", :]


#plot xch4
scatter(ta_data.dt,ta_data.xch4,msw=0,c=ta_data.color,label="ta",
        xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
scatter!(tb_data.dt,tb_data.xch4,msw=0,c=tb_data.color,label="tb")
scatter!(tc_data.dt,tc_data.xch4,msw=0,c=tc_data.color,label="tc")
scatter!(td_data.dt,td_data.xch4,msw=0,c=td_data.color,label="td")
plot!(xlabel="Date",ylabel="XCH₄ (ppm)",title="UTSG XCH₄ Observations",
            legend=:bottomright)
pltpath = "/home/lawson/HDD/Data/"
png(pltpath*"UTSGXCH4")

#plotxco2
scatter(ta_data.dt,ta_data.xco2,msw=0,c=ta_data.color,label="ta",
            xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
scatter!(tb_data.dt,tb_data.xco2,msw=0,c=tb_data.color,label="tb")
scatter!(tc_data.dt,tc_data.xco2,msw=0,c=tc_data.color,label="tc")
scatter!(td_data.dt,td_data.xco2,msw=0,c=td_data.color,label="td")
plot!(xlabel="Date",ylabel="XCO₂ (ppm)",title="UTSG XCO₂ Observations",
            legend=:bottomright)
png(pltpath*"UTSGXCO2")


#plotxco
scatter(ta_data.dt,ta_data.xco,msw=0,c=ta_data.color,label="ta",
        xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
scatter!(tb_data.dt,tb_data.xco,msw=0,c=tb_data.color,label="tb")
scatter!(tc_data.dt,tc_data.xco,msw=0,c=tc_data.color,label="tc")
scatter!(td_data.dt,td_data.xco,msw=0,c=td_data.color,label="td")
plot!(xlabel="Date",ylabel="XCO (ppm)",title="UTSG XCO Observations",
            legend=:topleft)
png(pltpath*"UTSGXCO")



#Plots for UTSC
#colour by instrument for plotting ease...
dic = Dict([("ta",:red),("tb",:orange),("tc",:blue),("td",:lightblue)])
utsc_data.color = [dic[r] for r in utsc_data.ins]

ta_data = utsc_data[utsc_data.ins .== "ta", :]
tb_data = utsc_data[utsc_data.ins .== "tb", :]
tc_data = utsc_data[utsc_data.ins .== "tc", :]
td_data = utsc_data[utsc_data.ins .== "td", :]


#plot xch4
# scatter(ta_data.dt,ta_data.xch4,msw=0,c=ta_data.color,label="ta")
scatter(tb_data.dt,tb_data.xch4,msw=0,c=tb_data.color,label="tb",
            xlims = (minimum(data1_df.dt),maximum(data1_df.dt)) )
# scatter!(tc_data.dt,tc_data.xch4,msw=0,c=tc_data.color,label="tc")
# scatter!(td_data.dt,td_data.xch4,msw=0,c=td_data.color,label="td")
plot!(xlabel="Date",ylabel="XCH₄ (ppm)",title="UTSC XCH₄ Observations",
            legend=:topleft)
pltpath = "/home/lawson/HDD/Data/"
png(pltpath*"UTSCXCH4")

#plotxco2
# scatter(ta_data.dt,ta_data.xco2,msw=0,c=ta_data.color,label="ta")
scatter(tb_data.dt,tb_data.xco2,msw=0,c=tb_data.color,label="tb",
    xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# scatter!(tc_data.dt,tc_data.xco2,msw=0,c=tc_data.color,label="tc")
# scatter!(td_data.dt,td_data.xco2,msw=0,c=td_data.color,label="td")
plot!(xlabel="Date",ylabel="XCO₂ (ppm)",title="UTSC XCO₂ Observations",
            legend=:left)
png(pltpath*"UTSCXCO2")


#plotxco
# scatter(ta_data.dt,ta_data.xco,msw=0,c=ta_data.color,label="ta")
scatter(tb_data.dt,tb_data.xco,msw=0,c=tb_data.color,label="tb",
            xlims = (minimum(data1_df.dt),maximum(data1_df.dt)) )
# scatter!(tc_data.dt,tc_data.xco,msw=0,c=tc_data.color,label="tc")
# scatter!(td_data.dt,td_data.xco,msw=0,c=td_data.color,label="td")
plot!(xlabel="Date",ylabel="XCO (ppm)",title="UTSC XCO Observations",
            legend=:left)
png(pltpath*"UTSCXCO")


#Plots for UTM
#colour by instrument for plotting ease...
utm_data.color = [dic[r] for r in utm_data.ins]

ta_data = utm_data[utm_data.ins .== "ta", :]
tb_data = utm_data[utm_data.ins .== "tb", :]
tc_data = utm_data[utm_data.ins .== "tc", :]
td_data = utm_data[utm_data.ins .== "td", :]


#plot xch4
# scatter(ta_data.dt,ta_data.xch4,msw=0,c=ta_data.color,label="ta")
# scatter(tb_data.dt,tb_data.xch4,msw=0,c=tb_data.color,label="tb")
scatter(tc_data.dt,tc_data.xch4,msw=0,c=tc_data.color,label="tc",
        xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# scatter!(td_data.dt,td_data.xch4,msw=0,c=td_data.color,label="td")
plot!(xlabel="Date",ylabel="XCH₄ (ppm)",title="UTM XCH₄ Observations",
            legend=:topleft)
pltpath = "/home/lawson/HDD/Data/"
png(pltpath*"UTM_XCH4")

#plotxco2
# scatter(ta_data.dt,ta_data.xco2,msw=0,c=ta_data.color,label="ta")
# scatter(tb_data.dt,tb_data.xco2,msw=0,c=tb_data.color,label="tb")
scatter(tc_data.dt,tc_data.xco2,msw=0,c=tc_data.color,label="tc",
            xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# scatter!(td_data.dt,td_data.xco2,msw=0,c=td_data.color,label="td")
plot!(xlabel="Date",ylabel="XCO₂ (ppm)",title="UTM XCO₂ Observations",
            legend=:left)
png(pltpath*"UTM_XCO2")


#plotxco
# scatter(ta_data.dt,ta_data.xco,msw=0,c=ta_data.color,label="ta")
# scatter(tb_data.dt,tb_data.xco,msw=0,c=tb_data.color,label="tb")
scatter(tc_data.dt,tc_data.xco,msw=0,c=tc_data.color,label="tc",
            xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
# scatter!(td_data.dt,td_data.xco,msw=0,c=td_data.color,label="td")
plot!(xlabel="Date",ylabel="XCO (ppm)",title="UTM XCO Observations",
            legend=:left)
png(pltpath*"UTM_XCO")



#Plots for ECCC Downsview
#colour by instrument for plotting ease...
dow_data.color = [dic[r] for r in dow_data.ins]

ta_data = dow_data[dow_data.ins .== "ta", :]
tb_data = dow_data[dow_data.ins .== "tb", :]
tc_data = dow_data[dow_data.ins .== "tc", :]
td_data = dow_data[dow_data.ins .== "td", :]


#plot xch4
# scatter(ta_data.dt,ta_data.xch4,msw=0,c=ta_data.color,label="ta")
# scatter(tb_data.dt,tb_data.xch4,msw=0,c=tb_data.color,label="tb")
scatter(td_data.dt,td_data.xch4,msw=0,c=td_data.color,label="td",
            xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
scatter!(tc_data.dt,tc_data.xch4,msw=0,c=tc_data.color,label="tc")
plot!(xlabel="Date",ylabel="XCH₄ (ppm)",title="Downsview XCH₄ Observations",
            legend=:left)
pltpath = "/home/lawson/HDD/Data/"
png(pltpath*"DOW_XCH4")

#plotxco2
# scatter(ta_data.dt,ta_data.xco2,msw=0,c=ta_data.color,label="ta")
# scatter(tb_data.dt,tb_data.xco2,msw=0,c=tb_data.color,label="tb")
scatter(td_data.dt,td_data.xco2,msw=0,c=td_data.color,label="td",
    xlims = (minimum(data1_df.dt),maximum(data1_df.dt)))
scatter!(tc_data.dt,tc_data.xco2,msw=0,c=tc_data.color,label="tc")
plot!(xlabel="Date",ylabel="XCO₂ (ppm)",title="Downsview XCO₂ Observations",
            legend=:left)
png(pltpath*"DOW_XCO2")


#plotxco
# scatter(ta_data.dt,ta_data.xco,msw=0,c=ta_data.color,label="ta")
# scatter(tb_data.dt,tb_data.xco,msw=0,c=tb_data.color,label="tb")
scatter(td_data.dt,td_data.xco,msw=0,c=td_data.color,label="td",
        xlims = (minimum(data1_df.dt),maximum(data1_df.dt)) )
scatter!(tc_data.dt,tc_data.xco,msw=0,c=tc_data.color,label="tc")
plot!(xlabel="Date",ylabel="XCO (ppm)",title="UTM XCO Observations",
            legend=:topleft)
png(pltpath*"DOW_XCO")



scatter(utsg_data.dt,utsg_data.xch4,msw=0.001,c=utsg_data.color)


CSV.write("/home/lawson/HDD/Data/EM27/utsg.csv",utsg_data)
CSV.write("/home/lawson/HDD/Data/EM27/utsc.csv",utsc_data)
CSV.write("/home/lawson/HDD/Data/EM27/utm_data.csv",utm_data)
CSV.write("/home/lawson/HDD/Data/EM27/dow_data.csv",dow_data)
