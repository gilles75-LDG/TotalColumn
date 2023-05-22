using CSV
using DataFrames
using Dates
using Statistics
using Glob
using Plots
using GLM

#Read in the KIT Vaisala pressure data....
kit_path = "/home/lawson/Data/PressureCorrections/KIT_travel_Vaisala/"
kit_data_paths = glob("*RawData*",kit_path)
kit_data = CSV.read(kit_data_paths[1],DataFrame)

for path in kit_data_paths[2:end]
    df = CSV.read(path,DataFrame)
    append!(kit_data,df)
end
dfmt1 = "yyyy-mm-ddTHH:MM:SS"
kit_data.dt = [Dates.DateTime(t[1:19],dfmt1) for t in kit_data.SystemUTC]
# kit_data = kit_data[kit_data.dt .< Dates.Date(2022,08,14),:]


scatter(kit_data.dt,kit_data.Pressure,
        markersize=1,msc=0,label="KIT Vaisala",legend=:topleft,
        ylabel="Pressure (hPa)")


#Read in the digi pressure data...

digi_path = "/home/lawson/Data/PressureCorrections/vaisala2_UofT/Digiquartz/20220805_digi.txt"

digi_data = CSV.read(digi_path,DataFrame;header=0)
dfmt2 = "m/d/y HH:MM:SS"

digi_data.dt = [Dates.DateTime(t[1:18],dfmt2)+Dates.Hour(12) for t in digi_data.Column2]

#Read in the Vaisala #2222 data (goes with Kate1)...

vai2_path = "/home/lawson/Data/PressureCorrections/vaisala2_UofT/"
vai2_paths  = glob("*_vaisala.txt", vai2_path)

vai2_data = CSV.read(vai2_paths[1],DataFrame)

for path in vai2_paths[2:end]
    df = CSV.read(path,DataFrame)
    append!(vai2_data,df)
end

dfmt4 = "yyyy/mm/ddHH:MM:SS"
vai2_data.dt = [DateTime(i.UTCDate * i.UTCTime[1:8], dfmt4) for i in eachrow(vai2_data)]

scatter(kit_data.dt,kit_data.Pressure,
        markersize=1,msw=0,label="KIT Vaisala",legend=:topleft,
        ylabel="Pressure (hPa)")

scatter!(vai2_data.dt,vai2_data.Pout,
         markersize=1,msw=0,label="Vai222")

combi_data = innerjoin(kit_data,vai2_data,on=:dt)

scatter(combi_data.dt,combi_data.Pressure,
        markersize=1,msw=0,label="KIT Vaisala",legend=:topright,
        ylabel="Pressure (hPa)")

scatter!(combi_data.dt,combi_data.Pout,
        markersize=1,msw=0,label="Vai222")

scatter(combi_data.Pressure,combi_data.Pout,msw=0,legend=:topleft,marker_z=combi_data.Sensor_Temperature)

fit1 = lm(@formula(Pressure ~  Pout),combi_data)
pred1 = predict(fit1, combi_data)
