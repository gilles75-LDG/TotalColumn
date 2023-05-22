using CSV
using Plots
using DataFrames


data_path = "/home/lawson/Data/EM27SUN/spectra/tf/20220726/spt/ch4_6002/ztf20220726s0e00a.0142"

spec_data = CSV.read(data_path,DataFrame;header=3,delim=" ",ignorerepeated=true)

rename!(spec_data, [strip(n) for n in names(spec_data)])


l = @layout [a{0.2h}; b]

plot1 = plot(spec_data.Freq, spec_data.Tm./spec_data.Cont, label = "Measured Spectrum", color =:black,
xlabel = "Wavenumber in cm⁻¹" )


# plot!(plot1, spec_data.Freq, spec_data.o2, label = "O₂", color =:red)
# plot!(plot1, spec_data.Freq, spec_data."0o2", label = "O₂ Cont.", color = :lightgreen)
plot!(plot1, spec_data.Freq, spec_data.ch4, label = "CH₄", color = :red)

plot!(plot1, spec_data.Freq, spec_data.co2, label = "CO₂", color = :yellow)
plot!(plot1, spec_data.Freq, spec_data.h2o, label = "H₂O", color = :lightblue)
plot!(plot1, spec_data.Freq, spec_data.hdo, label = "HDO", color = :grey)
plot!(plot1, spec_data.Freq, spec_data.other, label = "Other", color = :pink)
plot!(plot1, spec_data.Freq, spec_data.solar, label = "Solar", color = :orange,
    legend=:bottomright)

plot!(plot1,spec_data.Freq,spec_data.Tc./spec_data.Cont, label = "Calculated Spectrum", color =:green)

ylabel!(plot1,"Transmittance",legend=:bottomleft)


spec_data.res = (spec_data.Tc - spec_data.Tm)./spec_data.Cont
#

plot2 = plot(spec_data.Freq, spec_data.res, ylabel = "Residual", legend = false, linewidth = 0.4,)
        # ylims = (-0.01,0.01) )
            # ylims =  (-0.08,0.08)) #
title!(plot2,"tf 20220726 #0142")
# ylims!(plot2, (-0.01,0.01))
plot3 = plot(plot2, plot1, layout = l, figuresize = (600*3,400*3), dpi = 100)
file2 = "/home/lawson/Data/EM27SUN/spectra/test3.png"
png(plot3,file2)
