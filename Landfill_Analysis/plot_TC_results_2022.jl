using Plots
using DataFrames
using Statistics
using StatsBase

data_path = "/home/lawson/Data/EM27SUN/Landfills/TwinCreeksCampaigns.csv"
data = CSV.read(data_path,DataFrame)
data.dt = [DateTime(string(d), "yyyymmdd") for d in data.Date]

tc_data = data[data.Site .== "TwinCreeks",:]
tc_data = tc_data[tc_data.dt .> Date(2022),:]
pet_data = data[data.Site .== "Petrolia", :]

pltpath = "/home/lawson/Data/EM27SUN/plots/Landfills/"
scatter(tc_data.dt,tc_data.Q_kgd,yerr=tc_data.Q_kgd_err)

post_data

post_data = CSV.read("/home/lawson/Data/InSituMethane/Landfills/PetroliaTransects.csv",DataFrame)

tran_data_pet = post_data
post_data = CSV.read("/home/lawson/Data/InSituMethane/Landfills/InversionResult/inversion_info_TwinCreeks20220727_C.csv",DataFrame)
tran_data = post_data
tran_data = tran_data[tran_data.source .!= "1L", :]

scatter(tran_data.endtime,tran_data.post_area,yerr=tran_data.post_area_err,label="LGR Mobile In Situ")
scatter!(tc_data.dt.+Hour(12),tc_data.Q_kgd,yerr=tc_data.Q_kgd_err,marker=:star,
legend=:bottomright,label="EM27/Sun Estimate",ylabel="Estimated Emissions CH₄ (kg⋅day⁻¹)",
title="Twin Creeks Landfill\n2022Emissions Estimates")
hline!([4551*1000/365],label="ECCC GHGRP")
png(pltpath*"combined")

scatter(post_data.endtime,post_data.post_area,yerr=post_data.post_area_err)


scatter(tran_data_pet.starttime, tran_data_pet.post_area, 
yerr = tran_data_pet.post_area_err,label="LGR Mobile In Situ")
scatter!(pet_data.dt.+Hour(12),pet_data.Q_kgd,yerr=pet_data.Q_kgd_err,
label="EM27/Sun Estimate",
ylabel="Estimated Emissions CH₄ (kg⋅day⁻¹)",legend=:top,
title="Petrolia Landfill\n2021 Emissions Estimates")

png(pltpath*"combined_pet")


scatter(post_data1.endtime,post_data1.post_area,yerr=post_data1.post_area_err,label="LGR Mobile In Situ",
ylabel="Estimated Emissions CH₄ (kg⋅day⁻¹)",
title="Twin Creeks Landfill\n2022-07-27 Emissions Estimates",xrot=11)
hline!([12500],label="ECCC GHGRP")
png(pltpath*"insitu20220727")