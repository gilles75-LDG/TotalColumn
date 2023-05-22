using Plots
using DataFrames
using Statistics
using StatsBase

data_path = "/home/lawson/Data/EM27SUN/Landfills/TwinCreeksCampaigns.csv"
data = CSV.read(data_path,DataFrame)
data.dt = [DateTime(string(d), "yyyymmdd") for d in data.Date]

tc_data = data[data.Site .== "TwinCreeks",:]
pet_data = data[data.Site .== "Petrolia", :]

scatter(tc_data.dt,tc_data.Q_kgd,yerr=tc_data.Q_kgd_err)