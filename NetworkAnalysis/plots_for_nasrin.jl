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
nc_path = "/home/lawson/Data/EM27SUN/nc_files/cumulative/"
# daterng = 202109
# daterng = 2022072

nc_paths = glob("ta*.nc",nc_path)
# append!(nc_paths,[i for i in glob("te2023*.nc",nc_path)])
# #make a list of NCDatasets from netCDF's
datasets = []
for path in nc_paths
    ds = Dataset(path)
    push!(datasets,ds)
end

# #list of variables to extract to dataframe
variables  = ["xch4", "xco2", "xco", 
                # "flag",
                "year","day",
            "hour","xluft",
            # "tins",
            "xluft_error",
            # "o2_7885_sg",
            "pout",
            # "pins",
            "lat","long",
            # "column_luft",
            "xluft",
            # "column_h2o",
            "xh2o",
            # "xch4_aicf",
            # "xch4_aicf",
            "solzen",
            # "spectrum"
            ]
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
data1_df.ins = ["ta" for i in eachrow(data1_df)]
unique!(data1_df,:dt)

data_ta = data1_df #[data1_df.flag .== 0,:]

# #now do tf data...

# nc_path  = "/home/lawson/Data/EM27SUN/nc_files/cu"
nc_paths = glob("td*.nc",nc_path)

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
data1_df.ins = ["td" for i in eachrow(data1_df)]
unique!(data1_df,:dt)

data_td = data1_df #[data1_df.flag .== 0,:]



scatter(data_ta.dt,data_ta.xch4,msw=0.1,xrot=-11,markersize=2,α=.3,label="ta")
scatter!(data_td.dt,data_td.xch4,msw=0.1,xrot=-11,markersize=2,α=.3,label="td")

pltpath = "/home/lawson/Data/EM27SUN/plots/"
png(pltpath*"ta_td_timeseries_xch4")

scatter(data_ta.dt,data_ta.xco2,msw=0.1,xrot=-11,markersize=2,α=.3)
scatter!(data_td.dt,data_td.xco2,msw=0.1,xrot=-11,markersize=2,α=.3)

pltpath = "/home/lawson/Data/EM27SUN/plots/"
png(pltpath*"ta_td_timeseries_xco2")
