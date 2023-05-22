using NCDatasets
using DataFrames
using CSV
using Dates
using Plots
using Glob
using Statistics
using PyCall

#---useful functions
function data_in(start, finish, ins)
    # #path to netCDF's
    nc_path = "/home/lawson/Data/EM27SUN/nc_files/$(ins)/"
    nc_paths1 = glob("$(ins)*.nc",nc_path)
    nc_paths = []
    println(nc_paths1)
    for pth in nc_paths1
        println(pth)
        d = Date(split(split(pth, "/")[end],"_")[1][3:end],"yyyymmdd")
        if start <= d <= finish
            push!(nc_paths,pth)
        end
    end
    datasets = []
    for path in nc_paths
        ds = Dataset(path)
        push!(datasets,ds)
    end

    # #list of variables to extract to dataframe
    variables  = ["xch4", "xco2", "xco", "flag","year","day",
                "hour","xluft","tins","xluft_error","o2_7885_sg",
                "pout","pins","lat","long","column_luft","xluft","column_h2o","xh2o",
                "xch4_aicf","xch4_adcf"]


    # #make a list of dataframes to dump NCDataset data into
    dataframes = []
    for ds in datasets
        tmp_df = DataFrame()
        for var in variables
            tmp_df[!,"$var"] = ds[var][1:end]
        end
        push!(dataframes,tmp_df)
    end

    # #make it all into single dataframe
    data1_df  = dataframes[1]
    for df in dataframes[2:end]
        data1_df = vcat(data1_df,df)
    end


    data1_df.min = [floor((r.hour - floor(r.hour)) * 60) for r in eachrow(data1_df) ]
    data1_df.sec = [floor((r.hour - floor(r.hour) - r.min/60) * 3600) for r in eachrow(data1_df) ]
    data1_df.dt = [ Dates.DateTime(r.year) + Dates.Day(r.day - 1) + Dates.Hour(floor(r.hour)) + Dates.Minute(r.min) + Dates.Second(r.sec) for r in eachrow(data1_df)]
    data1_df.ins = ["$(ins)" for i in eachrow(data1_df)]
    unique!(data1_df,:dt)
    return data1_df
end

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
        if ! isempty(t1d)
            push!(t1_med_10,[median(t1d.xch4),median(t1d.xco2),median(t1d.xco),
                    dtt,median(t1d.column_luft),median(t1d.column_h2o),median(t1d.xh2o)])
        end
        if ! isempty(t2d)
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
#----

start = Date(2021)
finish = Date(2022)



data_tg = data_in(start,finish,"tg")
# data_tg = data_tg[data_tg.flag .== 0, :]
data_te = data_in(start,finish,"te")
data_te_filt = data_te[data_te.flag .== 0,:]
data_tg_filt = data_tg[data_tg.flag .== 0,:]

data_td = data_in(start,finish,"td")
data_td_filt = data_td[data_td.flag .== 18,:]
xch4_adcf = -0.0045
data_td_filt = data_td_filt[data_td_filt.xch4_adcf .≈ xch4_adcf,:]
data_tf = data_in(start,finish,"tf")
data_tf_filt = data_tf[data_tf.flag .== 0 ,:]


data_te.day

scatter(data_te_filt.dt,data_te_filt.xch4,msw=0.1,label="Egbert")
scatter!(data_tg_filt.dt,data_tg_filt.xch4,msw=0.1,α=.2,label="UTM")
scatter!(data_td_filt.dt,data_td_filt.xch4,msw=0.1,α=.2,label="Downsview") 
scatter!(data_tf_filt.dt,data_tf_filt.xch4,msw=0.1,α=.2,label="UTSC")
title!("GTA GHG-ON $(Date(start)) - $(Date(finish))")
png("GHG_ON_20230101_20230410")

#Egb vs Others
data_10min_te_tg = x_min_medians(data_te_filt, data_tg_filt,10)
data_10min_te_tf = x_min_medians(data_te_filt, data_tf_filt,10)
data_10min_te_td = x_min_medians(data_te_filt, data_td_filt,10)

#Downsview vs. Others
data_10min_td_tf = x_min_medians(data_td_filt, data_tf_filt,10)
data_10min_td_tg = x_min_medians(data_td_filt, data_tg_filt,10)

#UTSC vs. Others
data_10min_tf_tg = x_min_medians(data_tf_filt, data_tg_filt,10)


#Egbert vs rest (-ta)
scatter(data_10min_te_tg.xch4_t1,data_10min_te_tg.xch4_t2,c=:green,msw=0.1,
    xlabel="te - Egbert", ylabel ="Others",label="tg - UTM")

scatter!(data_10min_te_tf.xch4_t1,data_10min_te_tf.xch4_t2,c=:red,
marker=:star,label="tf - UTSC")

scatter!(data_10min_te_td.xch4_t1,data_10min_te_td.xch4_t2,c=:blue,
marker=:ltraingle,label="td - Downsview")
plot!(title="10-minute coincidental median measurements\nEgbert vs. Others")
png("egb_vs_rest_10min")

#Downsview vs. rest (-ta)
scatter(data_10min_te_td.xch4_t2,data_10min_te_td.xch4_t1,c=:green,msw=0.1,
    xlabel="td - Downsview", ylabel ="Others",label="te - Egbert")

scatter!(data_10min_td_tf.xch4_t1,data_10min_td_tf.xch4_t2,c=:red,
marker=:star,label="tf - UTSC")

scatter!(data_10min_td_tg.xch4_t1,data_10min_td_tg.xch4_t2,c=:blue,
marker=:ltraingle,label="tg - UTM")
plot!(title="10-minute coincidental median measurements\nDownsview vs. Others")
png("dow_vs_rest_10min")

#UTSC vs. rest (-ta)
scatter(data_10min_te_tf.xch4_t2,data_10min_te_tf.xch4_t1,c=:green,msw=0.1,
    xlabel="tf - UTSC", ylabel ="Others",label="te - Egbert")

scatter!(data_10min_td_tf.xch4_t2,data_10min_td_tf.xch4_t1,c=:red,
marker=:star,label="td - Downsview")

scatter!(data_10min_tf_tg.xch4_t1,data_10min_tf_tg.xch4_t2,c=:blue,
marker=:ltraingle,label="tg - UTM")
plot!(title="10-minute coincidental median measurements\nUTSC vs. Others")
png("utsc_vs_rest_10min")

#UTM vs. rest (-ta)
scatter(data_10min_te_tg.xch4_t2,data_10min_te_tg.xch4_t1,c=:green,msw=0.1,
    xlabel="tg - UTM", ylabel ="Others",label="te - Egbert")

scatter!(data_10min_td_tg.xch4_t2,data_10min_td_tg.xch4_t1,c=:red,
marker=:star,label="td - Downsview")

scatter!(data_10min_tf_tg.xch4_t2,data_10min_tf_tg.xch4_t1,c=:blue,
marker=:ltraingle,label="tf - UTSC")
plot!(title="10-minute coincidental median measurements\nUTM vs. Others")
png("utm_vs_rest_10min")



scatter(data_10min.xch4_t1,data_10min.xch4_t2)
scatter(data_10min.xco2_t1,data_10min.xco2_t2,marker_z=dayofyear.(data_10min.dt))

#---
scatter(data_10min.dt,
    (data_10min.xch4_t1./data_10min.xch4_t2)./(data_10min.xco2_t1./data_10min.xco2_t2))

 scatter(data_10min.dt,data_10min.xco2_t1./data_10min.xco2_t2)

data_tg.xluft