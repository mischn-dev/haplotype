ENV["R_HOME"] = "C:\\Program Files\\R\\R-3.6.3"
using DataFrames, Statistics, CSV, ProgressMeter, RCall, Crayons, Dates

# calculate the Kcb and Ke values of barley transpiration and evaporation
# ETc = (Kcb + Ke) * ETo
# Kcb = basal crop coeffient
# Ke = evaporation of soil

# klein altendorf
we1 = CSV.read("/media/michael/Uni/WetterDaten/Evaportranspiration/Klein_Altendorf_2011-2019.csv")
we1 = unique(we1)
# remove the missing dates
we1 = we1[.!(ismissing.(we1[!,:Kanal])),:]
# dates
we1[!,:Kanal] .= Date.(we1[!,:Kanal])

# dikopshof
ew = CSV.read("/media/michael/Uni/Wetterdaten Dikopshof/evapotranspiration_calc.csv", missingstring="NA")
ew = unique(ew)
# remove missing dates
ew = ew[.!(ismissing.(ew[!,:Date])),:]

# needed columns are: Date (Kanal), Windspeed, Humidity_min, precipitation, Temp_mean, ETo
ew2 = ew[!,[1,18,17,5,20,6,14,9,43]]
we2 = we1[!,[1,2,3,4,5,7,8,9,30]]
names!(ew2, names(we2))

# humidity minimum for dicokopshof has to be estimated
ew2[!,:Humidity_min] .= round.(ew2[!,:Humidity_min] .- (100 .- ew2[!,:Humidity_min]) , digits=1)
ew2[ew2[!,:Humidity_min] .< 20,:Humidity_min] .= 18.1
we = vcat(ew2,we2)

# first calculate the Temperature sum from sowing to harvest
sh = CSV.read("/media/michael/Uni/WetterDaten/gerste_saat_ernte.csv", delim=";")
df1 = DateFormat("d.m.Y")
for i in names(sh)[2:5]
                        sh[!,i] = Date.(sh[!,i],df1)
                        sh[!,i][5:24] .+= Dates.Year(2000)
                        sh[!,i][1:4] .+= Dates.Year(1900)
end


# subset years by the sowing and harvest dates
kon = Dict()
org = Dict()

for i in sh[!,:Year]
                        hs = sh[sh[!,:Year] .== i,:]
                        kon[i] = we[.&(we[!,:Kanal] .>= hs[!,2], we[!,:Kanal] .< hs[!,3]),:]
                        org[i] = we[.&(we[!,:Kanal] .>= hs[!,4], we[!,:Kanal] .< hs[!,5]),:]
end

# calculate the aggregated temperature sum for each day after sowing
env = Dict()
env[:kon] = kon
env[:org] = org

for set in collect(keys(env))
        for year in collect(keys(env[set]))
                                env[set][year][!,:TSum] .= 0.0
                                env[set][year][!,:SWC] .= 180.0 # Soil water content

        for i in 1:size(env[set][year],1)
                if i == 1
                         env[set][year][i,:TSum] = env[set][year][i,:Temp_Mean]
                         continue
                end
                if ismissing(env[set][year][i,:Temp_Mean]) # if missing, calculate the mean of the month and use it
                                                        val = mean(env[set][year][.&(Month.(env[set][year][!,:Kanal]) .== Month(env[set][year][i,:Kanal]), .!(ismissing.(env[set][year][!,:Temp_Mean]))),:Temp_Mean])
                                                        env[set][year][i,:TSum] = val + env[set][year][i-1,:TSum]
                else
                # when the values is present, use the mean temp of the day and add it to the previous tempsum
                env[set][year][i,:TSum] =  env[set][year][i,:Temp_Mean] + env[set][year][i-1,:TSum]
                end
        end # DOY

        end # of year
end # Environment

#=  basal crop coefficient (Kcb) is defined as the ratio of the crop evapotranspiration over the reference evapotranspiration (ETc/ETo)
 needed are - mean wind speed; minimum humidity; plant height of mid or late season stage (1 meter after fao )
 Kcb ini = 0.15 | Kcb mid = 1.10 | Kcb end = 0.15
the growth period will be split by the temperature sum in 4 categories -
initial - emergenece and leaf development - up to 190 °C
crop development - Tilling and stem elongation - up to 560°C
mid season - Anthesis and Seed fill - up to 1200
late season - Dough and maturity stage - up to 1600 °C
https://store.msuextension.org/publications/AgandNaturalResources/MT200103AG.pdf page 4

- crop development and lat season will be modelled as values in between start - mid - end


Kcb = Kcb(stage) +[0.04(windSpeed -2) - 0.004 * (minHumidity - 45)] * (plantHeight / 3)^0.3
=#

# to access crop development and late season correct, the length of the stages has to be calculated based on the Tmeperature sum and the bounds by the growth stages
initial = 220
cropdev = 650
mids = 1300
endse = 1650

for set in collect(keys(env))
        for year in collect(keys(env[set]))
                                env[set][year][!,:GrowthStage] .= "none"

        for i in 1:size(env[set][year],1)
                if env[set][year][i,:TSum] < initial
                                                                                env[set][year][i,:GrowthStage] = "initial"
                elseif  env[set][year][i,:TSum] > initial && env[set][year][i,:TSum] < cropdev
                                                                                env[set][year][i,:GrowthStage] = "CropDev"
                elseif  env[set][year][i,:TSum] > cropdev && env[set][year][i,:TSum] < mids
                                                                                env[set][year][i,:GrowthStage] = "Anthesis"
                elseif  env[set][year][i,:TSum] >  mids && env[set][year][i,:TSum] < endse
                                                                                env[set][year][i,:GrowthStage] = "Ripening"
                else
                        env[set][year][i,:GrowthStage] = "Ready"
                end

        end # DOY

        end # of year
end # Environment

# calculate the Kci value for each segement of growth period
# Kcb prev + [(day_number - sum(lengthOfPrevStage)) / stageLength] * (Kcnext - Kcprev)

for set in collect(keys(env))
        for year in collect(keys(env[set]))
                                env[set][year][!,:Kcb] .= 0.0

        for i in 1:size(env[set][year],1)
                if env[set][year][i,:GrowthStage] == "initial"
                                                                                env[set][year][i,:Kcb] = 0.15
                elseif  env[set][year][i,:GrowthStage] == "CropDev"
                                                                                env[set][year][i,:Kcb] = 0.15 + ((i - sum(env[set][year][!,:GrowthStage].=="initial")) /  sum(env[set][year][!,:GrowthStage].=="CropDev") ) * (1.10 -0.15)
                elseif  env[set][year][i,:GrowthStage] == "Anthesis"
                                                                                env[set][year][i,:Kcb] = 1.10
                elseif  env[set][year][i,:GrowthStage] == "Ripening"
                                                                                env[set][year][i,:Kcb] = 1.10 + ((i - (sum(env[set][year][!,:GrowthStage].=="initial") + sum(env[set][year][!,:GrowthStage].=="CropDev") + sum(env[set][year][!,:GrowthStage].=="Anthesis"))) /  sum(env[set][year][!,:GrowthStage].=="Ripening") ) * (0.15 - 1.1)
                else
                        env[set][year][i,:Kcb] = 0.15
                end

        end # DOY

        end # of year
end # Environment

# Kcb = Kcb(stage) +[0.04(windSpeed -2) - 0.004 * (minHumidity - 45)] * (plantHeight / 3)^0.3
# calculate the actual transpiration of the barley plant under soil water capacity

for set in collect(keys(env))
        for year in collect(keys(env[set]))
                                env[set][year][!,:kCB] .= 0.0

        for i in 1:size(env[set][year],1)
                                        if ismissing(env[set][year][i,:Temp_Mean]) # if missing, calculate the mean of the month and use it
                                                                                wind = mean(env[set][year][.&(Month.(env[set][year][!,:Kanal]) .== Month(env[set][year][i,:Kanal]), .!(ismissing.(env[set][year][!,:Windspeed]))),:Windspeed])
                                                                                env[set][year][i,:Windspeed] = wind
                                                                                hum = mean(env[set][year][.&(Month.(env[set][year][!,:Kanal]) .== Month(env[set][year][i,:Kanal]), .!(ismissing.(env[set][year][!,:Humidity_min]))),:Humidity_min])
                                                                                env[set][year][i,:Humidity_min] = hum
                                                                                env[set][year][i,:kCB] =  env[set][year][i,:Kcb] + (0.04 * (wind - 2) - 0.004 * (hum - 45 )) * (1 / 3)^0.3

                                        else
                                        # when the values is present, use the mean temp of the day and add it to the previous tempsum
                                        env[set][year][i,:kCB] =  env[set][year][i,:Kcb] + (0.04 * (env[set][year][i,:Windspeed] - 2) - 0.004 * (env[set][year][i,:Humidity_min] - 45 )) * (1 / 3)^0.3
                                        end
        end # DOY
        end # of year
end # Environment


##### estimating the soil component - evapotranspiration Ke
# Ke is the minimum of 1. the exposed and wetted soil * the maximum value of Kc (Kc max) | Kc max - Kcb
# Ke = min (Kr (Kc max - Kcb), few Kc max)

# calculate Kc max
sh[!,:Kcmax] .= 0.0
sh[!,:Kcmin] .= 0.0

for set in collect(keys(env))
        for year in collect(keys(env[set]))
                                        if !(isempty(env[set][year] ))
                                        kcmax = maximum(vcat((1.2 .+ (0.04 .* (env[set][year][!,:Windspeed] .- 2) - 0.004 .* (env[set][year][!,:Humidity_min] .-45) .* (1 / 3)^0.3 )) , env[set][year][!,:kCB] .+ 0.05))
                                        kcmin = minimum(vcat((1.2 .+ (0.04 .* (env[set][year][!,:Windspeed] .- 2) - 0.004 .* (env[set][year][!,:Humidity_min] .-45) .* (1 / 3)^0.3 )) , env[set][year][!,:kCB] .+ 0.05))

                                                        sh[sh[!,:Year] .== year, :Kcmax] .= kcmax
                                                        sh[sh[!,:Year] .== year, :Kcmin] .= kcmin

                                        end
                                end
end

# calculate Kr - soil evaporation reduction coefficient
TEW = 22.0 # total evaporatable water in the 10 cm top soil [mm]
REW = 9.0 # readily evaporatable water [mm] limit for stage one

# decent precipitation level 5 mm rain in one day !
rain = 3.0 # [mm]
# Kr = 1 when De (cumulative depth of evaporation) < REW
# Kr = (TEW -De) / (TEW - REW)
par = Dict("TEW" => TEW, "REW" => REW, "rain" => rain)

for set in collect(keys(env))
        for year in collect(keys(env[set]))
                                env[set][year][!,:Kr] .= 0.0
                                env[set][year][!,:De] .= 0.0
                                env[set][year][!,:Ke] .= 0.0



        for i in 1:size(env[set][year],1)
                                        if ismissing(env[set][year][i,:precipitation]) # if missing, calculate the mean of the month and use it
                                                                                env[set][year][i,:precipitation] = mean(env[set][year][.&(Month.(env[set][year][!,:Kanal]) .== Month(env[set][year][i,:Kanal]), .!(ismissing.(env[set][year][!,:precipitation]))),:precipitation])
                                        end
                                        if ismissing(env[set][year][i,:ETo]) # if missing, calculate the mean of the month and use it
                                                                                env[set][year][i,:ETo] = mean(env[set][year][.&(Month.(env[set][year][!,:Kanal]) .== Month(env[set][year][i,:Kanal]), .!(ismissing.(env[set][year][!,:ETo]))),:ETo])
                                        end

                                        # calculate De - depletion of end of the previous day
                                        if env[set][year][i,:De] < par["REW"] # is there much water in the topsoil?
                                                                                env[set][year][i,:Kr] = 1.0
                                        elseif env[set][year][i,:De] > par["REW"] # is it already depleted?
                                                                                env[set][year][i,:Kr] = (par["TEW"] - env[set][year][i,:De]) / (par["TEW"] - par["REW"]) # might not work for the first row, but should not be effected
                                        end

                                        env[set][year][i,:Ke] = env[set][year][i,:Kr] * (env[set][year][i,:kCB] + sh[sh[!,:Year] .== year, :Kcmax][1])

                                        # next day - if it has rained -> De = 0
                                        if i != size(env[set][year],1)
                                                                                if env[set][year][i,:precipitation] > par["rain"]
                                                                                                                                env[set][year][i+1,:De] = 0.0
                                                                                else     # no sig rain
                                                                                                                                env[set][year][i+1,:De] = env[set][year][i,:Ke] * env[set][year][i,:ETo] + env[set][year][i,:De]
                                                                                end
                                        end
                                        # adjust the REW value - after shallow rain, this can be evaporated still quickly
                                        if env[set][year][i,:precipitation] > par["rain"] && env[set][year][i,:precipitation] < par["REW"]
                                                                                                                                        par["REW"] = round(env[set][year][i,:precipitation])
                                        elseif env[set][year][i,:precipitation] > par["REW"] && env[set][year][i,:precipitation] < 9.0
                                                                                                                                        par["REW"] = round(env[set][year][i,:precipitation])
                                        elseif env[set][year][i,:precipitation] > 9.0
                                                                                                                                        par["REW"] = 9.0
                                        end

        end # DOY
        end # of year
end # Environment


## exposed and wet soil fraction few
# few = min(1 - fc, fw)
# 1- fc => avg exposed soil fraction not covered [0.01 - 1]
# fw => avg fraction of soil wetted by precipitation [0.01 - 1] - does not paly a role in here - omited

# fc = ((Kcb - Kcmin) / (Kcmax - Kcmin))^(1+0.5 * 1)
# few = 1 - fc

# calculate fc value
for set in collect(keys(env))
        for year in collect(keys(env[set]))
                                env[set][year][!,:fe] .= 0.0

        for i in 1:size(env[set][year],1)
                                         if env[set][year][i,:kCB]  > 0.15
                                                                                env[set][year][i,:fe] = 1- ((env[set][year][i,:kCB] - 0.15) / (sh[sh[!,:Year] .== year, :Kcmax][1] - 0.15))^1.5
                                        else
                                                env[set][year][i,:fe] = 1.0
                                        end
        end
    end
end

# get the minimum value if Ke and fe to get the real KE value
for set in collect(keys(env))
        for year in collect(keys(env[set]))
                                env[set][year][!,:KE] .= minimum.(vcat.(env[set][year][!,:Ke], env[set][year][!,:fe]))
        end
end

# calulate ETc
for set in collect(keys(env))
        for year in collect(keys(env[set]))
                                env[set][year][!,:Etc] .= (env[set][year][!,:kCB] .+ env[set][year][!,:KE]) .* env[set][year][!,:ETo]
        end
end


function summ(set)

        println(mean(set))
        println(maximum(set))
        println(minimum(set))
        println(sum(set))
end



## calculate the water balance
# it is differentiated between the actual ETa  and the potential Evapotranspiration ET* - ETA / ET* = 1 if the soil water > 30 % nFK
# if soil water < 30% nFK => linear drop to 0 ET -> 10% available soil water == 0 ET

reg = DataFrame(hcat(11:30,range(0.01,1, step=0.05)))

for set in collect(keys(env))
        for year in collect(keys(env[set]))

                env[set][year][!,:Drought] .= ""
                env[set][year][!,:SWC] .= 180.0 # Soil water content
                env[set][year][!,:ETa] .= 0.0
                                                for i in 1:size(env[set][year],1)
                                                                                 # calculcate the ETA  - ETc equals to ET* potential evapotranspiration
                                                                        if i != size(env[set][year],1)
                                                                                 if env[set][year][i,:SWC] > env[set][year][1,:SWC] * 0.3
                                                                                                                                         env[set][year][i+1,:SWC] =  env[set][year][i,:SWC] -  env[set][year][i,:Etc] +  env[set][year][i,:precipitation]
                                                                                                                                         env[set][year][i,:Drought] = "none"
                                                                                                                                         env[set][year][i,:ETa] = env[set][year][i,:Etc]
                                                                                                                                         if  env[set][year][i,:SWC] > env[set][year][1,:SWC];  env[set][year][i,:SWC] =  env[set][year][1,:SWC]; end # there cannot be stored more water in the soil than the start value
                                                                                elseif env[set][year][i,:SWC] < env[set][year][1,:SWC] * 0.3  && env[set][year][i,:SWC] > env[set][year][1,:SWC] * 0.1 # soil water content is lower 30% nFK
                                                                                                                                         env[set][year][i,:Drought] = "stressed"
                                                                                                                                         # calculate the % soil water content
                                                                                                                                         ta = round(env[set][year][i,:SWC] / env[set][year][1,:SWC] *100 )
                                                                                                                                         tval = reg[reg[!,:x1].== ta,:x2]
                                                                                                                                         env[set][year][i+1,:SWC] =  env[set][year][i,:SWC] -  env[set][year][i,:Etc] * tval[1] +  env[set][year][i,:precipitation]
                                                                                                                                         env[set][year][i,:ETa] = env[set][year][i,:Etc] * tval[1]
                                                                                                                                         if  env[set][year][i,:SWC] > env[set][year][1,:SWC];  env[set][year][i,:SWC] =  env[set][year][1,:SWC]; end
                                                                                elseif env[set][year][i,:SWC] < env[set][year][1,:SWC] * 0.1
                                                                                                                                        env[set][year][i,:Drought] = "severe stress"
                                                                                                                                        env[set][year][i+1,:SWC] =  env[set][year][i,:SWC] -  env[set][year][i,:Etc] * 0.01 +  env[set][year][i,:precipitation]
                                                                                                                                        env[set][year][i,:ETa] = env[set][year][i,:Etc] * 0.01
                                                                                                                                        if  env[set][year][i,:SWC] > env[set][year][1,:SWC];  env[set][year][i,:SWC] =  env[set][year][1,:SWC]; end
                                                                                end
                                                                        end
                                                        end
        end
end

using Gadfly



k96 = stack(kon[:1996], [:ETo,:Etc , :ETa]);
k98 = stack(kon[:1998], [:ETo,:Etc , :ETa]);
k99 = stack(kon[:1999], [:ETo,:Etc , :ETa]);
k00 = stack(kon[:2000], [:ETo,:Etc , :ETa]);
k01 = stack(kon[:2001], [:ETo,:Etc , :ETa]);
k02 = stack(kon[:2002], [:ETo,:Etc , :ETa]);
k03 = stack(kon[:2003], [:ETo,:Etc , :ETa]);
k04 = stack(kon[:2004], [:ETo,:Etc , :ETa]);
k05 = stack(kon[:2005], [:ETo,:Etc , :ETa]);
k06 = stack(kon[:2006], [:ETo,:Etc , :ETa]);
k07 = stack(kon[:2007], [:ETo,:Etc , :ETa]);
k08 = stack(kon[:2008], [:ETo,:Etc , :ETa]);
k09 = stack(kon[:2009], [:ETo,:Etc , :ETa]);

k11 = stack(kon[:2011], [:ETo,:Etc , :ETa]);
k12 = stack(kon[:2012], [:ETo,:Etc , :ETa]);
k13 = stack(kon[:2013], [:ETo,:Etc , :ETa]);
k14 = stack(kon[:2014], [:ETo,:Etc , :ETa]);
k15 = stack(kon[:2015], [:ETo,:Etc , :ETa]);
k16 = stack(kon[:2016], [:ETo,:Etc , :ETa]);
k17 = stack(kon[:2017], [:ETo,:Etc , :ETa]);
k18 = stack(kon[:2018], [:ETo,:Etc , :ETa]);
k19 = stack(kon[:2019], [:ETo,:Etc , :ETa]);

latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=45pt,
                    minor_label_font="CMU Serif", minor_label_font_size=43pt,
                    key_title_font="CMU Serif", key_title_font_size=41pt,
                    key_label_font="CMU Serif", key_label_font_size=49pt)
Gadfly.push_theme(latex_fonts)


p96 = plot(k96, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("1996"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p98 = plot(k98, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("1998"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p99 = plot(k99, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("1999"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p00 = plot(k00, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2000"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p01 = plot(k01, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2001"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p02 = plot(k02, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2002"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p03 = plot(k03, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2003"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p05 = plot(k05, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2005"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p04 = plot(k04, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2004"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p06 = plot(k06, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2006"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p07 = plot(k07, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2007"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p08 = plot(k08, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2008"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p09 = plot(k09, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2009"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));


p11 = plot(k11, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2011"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p12 = plot(k12, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2012"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p13 = plot(k13, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2013"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p14 = plot(k14, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2014"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p15 = plot(k15, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2015"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p16 = plot(k16, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2016"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p18 = plot(k18, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2018"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p17 = plot(k17, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2017"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));
p19 = plot(k19, x=:Kanal, y=:value, color=:variable, Geom.line, Guide.xlabel(""), Guide.ylabel("ET [mm]"), Guide.colorkey(title=""), Guide.title("2019"), Theme(background_color=colorant"white", default_color=colorant"black"), Coord.cartesian(ymin=0, ymax=8));


pp = vstack(p96,p98,p99,p00,p01,p02,p03,p04,p05,p06,p07,p08,p09,p11,p12,p13,p14,p15,p16,p17,p18);
pp2 = hstack(vstack(p96,p98,p99,p00,p01,p02,p03,p04,p05,p06,p07), vstack(p08,p09,p11,p12,p13,p14,p15,p16,p17,p18,p19));
pp3 = hstack(vstack(p96,p98,p99,p00,p01,p02,p03,p04), vstack(p05,p06,p07,p08,p09,p11,p12), vstack(p13,p14,p15,p16,p17,p18,p19))
using Cairo

draw(PNG("/media/michael/Uni/WetterDaten/ETs1996_2019_2col180nFK.png",20inch, 15inch, dpi=600),pp2)
# publikation
draw(PNG("/media/michael/Uni/Allelfreq/Literature_NG_paper/ETs1996_2019_2col180nFK.png",12inch, 18inch, dpi=200),pp2)
# get the information how long the drought happened
for set in collect(keys(env))
        for year in collect(keys(env[set]))

                                            println(hcat(year,size(env[set][year][.&(env[set][year][!,:Drought] .== "stressed", env[set][year][!,:GrowthStage] .== "Anthesis"),:],1)))
                                        end
                                    end


## save the files
env[:kon][:all] = env[:kon][:2011][1:2,:]

for i in collect(keys(env[:kon]))
                                env[:kon][:all] = vcat(env[:kon][:all], env[:kon][i])
                        end
unique!(env[:kon][:all])
sort!(env[:kon][:all], :Kanal)


env[:org][:all] = env[:org][:2011][1:2,:]

for i in collect(keys(env[:org]))
                                env[:org][:all] = vcat(env[:org][:all], env[:org][i])
                        end
unique!(env[:org][:all])
sort!(env[:org][:all], :Kanal)

CSV.write("/media/michael/Uni/WetterDaten/ETa/Kon_Evapotranspiration96-19_180nFk.csv", env[:kon][:all])
CSV.write("/media/michael/Uni/WetterDaten/ETa/Org_Evapotranspiration96-19_180nFk.csv", env[:org][:all])
