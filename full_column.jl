using CSV
using Plots
using DataFrames
using Dates
using Statistics

data_path = "/home/lawson/Research/Data/EM27SUN/processed/data/utsg.csv"
data_path_sc = "/home/lawson/Research/Data/EM27SUN/processed/data/utsc.csv"
data_path_m = "/home/lawson/Research/Data/EM27SUN/processed/data/utm_data.csv"
data_path_dow = "/home/lawson/Research/Data/EM27SUN/processed/data/dow_data.csv"

utsg_data = CSV.read(data_path,DataFrame)
utsc_data = CSV.read(data_path_sc, DataFrame)
utm_data = CSV.read(data_path_m, DataFrame)
dow_data = CSV.read(data_path_dow, DataFrame)

utsg_data.dt2 = utsg_data.dt - Dates.Hour(5)
utsc_data.dt2 = utsc_data.dt - Dates.Hour(5)
utm_data.dt2 = utm_data.dt - Dates.Hour(5)
dow_data.dt2 = dow_data.dt - Dates.Hour(5)

utsg_data.date = [Dates.Date(r) for r in utsg_data.dt2]
utsc_data.date = [Dates.Date(r) for r in utsc_data.dt2]
utm_data.date = [Dates.Date(r) for r in utm_data.dt2]
dow_data.date = [Dates.Date(r) for r in dow_data.dt2]

dfs = [utsg_data,utsc_data,utm_data,dow_data]

dmin = minimum(utsg_data.date)
dmax = maximum(utsg_data.date) + Dates.Day(1)

dayy = dmin

utsg_ta = utsg_data[utsg_data.ins .== "ta",:]
utsg_tb = utsg_data[utsg_data.ins .== "tb",:]
utsg_tc = utsg_data[utsg_data.ins .== "tc",:]
utsg_td = utsg_data[utsg_data.ins .== "td",:]

pltpath = "/home/lawson/Research/Data/EM27SUN/"
#plot UTSG
#XCH₄
scatter(utsg_ta.dt,utsg_ta.xch4,c=:red,msw=0,xlims=(dtmin,dtmax),label="ta",
        tickfontsize = 12, legendfontsize=8,xrot=10)
scatter!(utsg_tb.dt,utsg_tb.xch4,c=:orange,msw=0,label="tb")
scatter!(utsg_tc.dt,utsg_tc.xch4,c=:blue,msw=0,label="tc")
scatter!(utsg_td.dt,utsg_td.xch4,c=:lightblue,msw=0,label="td",legendfontsize=12)
plot!(title="UTSG XCH₄ Time Series", legend=:bottomright,ylabel="XCH₄ (ppm)")
png(pltpath*"xch4_utsg")

#XCO₂
scatter(utsg_ta.dt,utsg_ta.xco2,c=:red,msw=0,xlims=(dtmin,dtmax),label="ta",
        tickfontsize = 12, legendfontsize=8,xrot=10)
scatter!(utsg_tb.dt,utsg_tb.xco2,c=:orange,msw=0,label="tb")
scatter!(utsg_tc.dt,utsg_tc.xco2,c=:blue,msw=0,label="tc")
scatter!(utsg_td.dt,utsg_td.xco2,c=:lightblue,msw=0,label="td",legendfontsize=12)
plot!(title="UTSG XCH₄ Time Series", legend=:bottomright,ylabel="XCH₄ (ppm)")
png(pltpath*"xco2_utsg")



scatter(utsg_ta.dt2,utsg_ta.xch4)


dates = []
for d in utsg_ta.date
    if ! (d in dates)
        push!(dates,d)
    end
end

ta_medians = []
for d in dates
    ddata = utsg_ta[utsg_ta.date .== d,:]
    xch4_med = median(ddata.xch4)
    xco2_med = median(ddata.xco2)
    xco_med = median(ddata.xco)
    push!(ta_medians,[xch4_med,xco2_med,xco_med])
end

xch4_ta_med = [x[1] for x in ta_medians]

dtmin = Dates.DateTime(dmin)
dtmax = Dates.DateTime(dmax)

dtt = dtmin
ta_med_10 = []
tb_med_10 = []
tc_med_10 = []
td_med_10 = []

while dtt < dtmax
    tad = utsg_ta[ dtt .< utsg_ta.dt2 .< dtt + Dates.Minute(10), :]
    tbd = utsg_tb[ dtt .< utsg_tb.dt2 .< dtt + Dates.Minute(10), :]
    tcd = utsg_tc[ dtt .< utsg_tc.dt2 .< dtt + Dates.Minute(10), :]
    tdd = utsg_td[ dtt .< utsg_td.dt2 .< dtt + Dates.Minute(10), :]
    if ! isempty(tad)
        push!(ta_med_10,[median(tad.xch4),median(tad.xco2),median(tad.xco),dtt])
    end
    if ! isempty(tbd)
        push!(tb_med_10,[median(tbd.xch4),median(tbd.xco2),median(tbd.xco),dtt])
    end
    if ! isempty(tcd)
        push!(tc_med_10,[median(tcd.xch4),median(tcd.xco2),median(tcd.xco),dtt])
    end
    if ! isempty(tdd)
        push!(td_med_10,[median(tdd.xch4),median(tdd.xco2),median(tdd.xco),dtt])
    end
    println(dtt)
    dtt = dtt + Dates.Minute(10)
end
ta_med_10
ta_df = DataFrame(xch4_ta = [r[1] for r in ta_med_10],
                    xco2_ta = [r[2] for r in ta_med_10],
                    xco_ta = [r[3] for r in ta_med_10],
                    dt = [r[4] for r in ta_med_10])

tb_df = DataFrame(xch4_tb = [r[1] for r in tb_med_10],
                    xco2_tb = [r[2] for r in tb_med_10],
                    xco_tb = [r[3] for r in tb_med_10],
                    dt = [r[4] for r in tb_med_10])

tc_df = DataFrame(xch4_tc = [r[1] for r in tc_med_10],
                    xco2_tc = [r[2] for r in tc_med_10],
                    xco_tc = [r[3] for r in tc_med_10],
                    dt = [r[4] for r in tc_med_10])

td_df = DataFrame(xch4_td = [r[1] for r in td_med_10],
                    xco2_td = [r[2] for r in td_med_10],
                    xco_td = [r[3] for r in td_med_10],
                    dt = [r[4] for r in td_med_10])


ta_tb = innerjoin(ta_df,tb_df,on=:dt)
ta_tc = innerjoin(ta_df,tc_df,on=:dt)
ta_td = innerjoin(ta_df,td_df,on=:dt)

scatter(ta_tb.xch4_ta,ta_tb.xch4_tb,msw=0,color=:red,label="tb")
scatter!(ta_tc.xch4_ta,ta_tc.xch4_tc,msw=0,c=:blue,label="tc")
scatter!(ta_td.xch4_ta,ta_td.xch4_td,msw=0,c=:lightblue,label="td",legendfontsize=12)
plot!([minimum(ta_tb.xch4_ta),maximum(ta_tb.xch4_ta)],
        [minimum(ta_tb.xch4_ta),maximum(ta_tb.xch4_ta)],label="1:1",c=:black)
title!("UTSG Coincidental XCH₄ Observations",tickfontsize=20,titlefontsize=18,
        legend=:bottomright)

tb_mse = mean(sqrt.((ta_tb.xch4_ta .- ta_tb.xch4_tb).^2))*1000
tc_mse = mean(sqrt.((ta_tc.xch4_ta .- ta_tc.xch4_tc).^2))*1000
td_mse = mean(sqrt.((ta_td.xch4_ta .- ta_td.xch4_td).^2))*1000




scatter(ta_tb.xco2_ta,ta_tb.xco2_tb,msw=0,marker_z=ta_tb.dt,label="tb")
scatter!(ta_tc.xco2_ta,ta_tc.xco2_tc,msw=0,c=:blue,label="tc")
scatter!(ta_td.xco2_ta,ta_td.xco2_td,msw=0,c=:lightblue,label="td",legendfontsize=12)
plot!([minimum(ta_tb.xco2_ta),maximum(ta_tb.xco2_ta)],
        [minimum(ta_tb.xco2_ta),maximum(ta_tb.xco2_ta)],label="1:1",c=:black)
title!("UTSG Coincidental XCH₄ Observations",tickfontsize=20,titlefontsize=18,
        legend=:bottomright)

tb_mse = mean(sqrt.((ta_tb.xco2_ta .- ta_tb.xco2_tb).^2))*1000
tc_mse = mean(sqrt.((ta_tc.xco2_ta .- ta_tc.xco2_tc).^2))*1000
td_mse = mean(sqrt.((ta_td.xco2_ta .- ta_td.xco2_td).^2))*1000



scatter(dates,xch4_ta_med)



utsg_data.day_f = floor(ustg_data)


scatter(utsg_data.dt[1:end],utsg_data.xch4[1:end])
