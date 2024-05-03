freq = readDict("D:/Allelfreq/ref_snps_dulicates_removed2020/Dict_Marker")
stacker(freq, "D:/Allelfreq/genoytpes_wgs.txt")
genes = intersect(freq[:F3E1_Haplotypes].Gene_ID, freq[:F12E1_Haplotypes].Gene_ID, freq[:F12E2_Haplotypes].Gene_ID, freq[:F16E1_Haplotypes].Gene_ID, freq[:F16E2_Haplotypes].Gene_ID, freq[:F22E1_Haplotypes].Gene_ID, freq[:F22E2_Haplotypes].Gene_ID, freq[:F23E1_Haplotypes].Gene_ID, freq[:F23E2_Haplotypes].Gene_ID)

v1 = vcat(freq[:F3E1_Haplotypes][indexin(genes, freq[:F3E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F12E1_Haplotypes][indexin(genes, freq[:F12E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F12E2_Haplotypes][indexin(genes, freq[:F12E2_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F16E1_Haplotypes][indexin(genes, freq[:F16E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F16E2_Haplotypes][indexin(genes, freq[:F16E2_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F22E1_Haplotypes][indexin(genes, freq[:F22E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F22E2_Haplotypes][indexin(genes, freq[:F22E2_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F23E1_Haplotypes][indexin(genes, freq[:F23E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F23E2_Haplotypes][indexin(genes, freq[:F23E2_Haplotypes].Gene_ID),:Allelfreq_Isr])

v2 = repeat([3,12,12,16,16,22,22,23,23], inner= length(genes))
v5 = repeat(genes, 9)
v4 = repeat(vcat(repeat(["conventional"], length(genes)), repeat(["organic"], length(genes))), 4)
v4 = vcat(repeat(["conventional"], length(genes)), v4)

v3 = DataFrame(hcat(v1,v2,v5,v4))

v3[!,:x1] .= Float64.(v3[!,:x1])
v3[!,:x2] .= Int64.(v3[!,:x2])
v3[!,:x3] .= String.(v3[!,:x3])
v3[!,:x4] .= String.(v3[!,:x4])
rename!(v3, [:wtAF, :Generation, :Gene, :Env])

plot(v3, x=:Generation, y=:wtAF, color=:Gene, Geom.line, Guide.colorkey(pos=[10,10]))

# create a data frame with mean and sd values for each generation and env
qay = DataFrame(Generation = Int64[], Env = String[], meanAF = Float64[], stdAF = Float64[])

for i in unique(v3.Generation), j in ["organic", "conventional"]
                                                                a1 = mean(v3[.&(v3.Generation.== i, v3.Env.== j),1])
                                                                a2 = std(v3[.&(v3.Generation.== i, v3.Env.== j),1])
                                                                push!(qay, [i,j,a1,a2])
end
qay[1,3] = qay[2,3]
qay[1,4] = qay[2,4]


plot(qay, x=:Generation, y=:meanAF, color=:Env, Geom.line, Guide.colorkey(pos=[10,10]))



@rput(v3)

R"""
library(ggplot2)
o3 = v3[v3$Env=="organic",]
o3 = rbind(v3[v3$Generation==3,], o3)
c3 = v3[v3$Env=="conventional",]

p1 = ggplot(c3, aes(x = Generation, y = wtAF, color=Gene)) + geom_line() + labs(x="Generation", y="Wild Type Allele Frequency") + ylim(0,1) + theme(legend.position = "none")
"""
@rget o3 c3
p1 = plot(o3, x=:Generation, y=:wtAF, color=:Gene, Geom.line, Theme(key_position = :none, background_color=colorant"white"), Guide.xlabel(""), Guide.ylabel("wtAF"), Guide.xticks(label=false), Guide.title("Organic"))
p2 = plot(c3, x=:Generation, y=:wtAF, color=:Gene, Geom.line, Theme(key_position = :none, background_color=colorant"white"), Guide.xlabel(""), Guide.ylabel("wtAF"), Guide.xticks(label=false), Guide.title("Conventional"))

# simulated Data

sim = readDict("D:/Allelfreq/Drift_calculation/Popsize5k")
sim2 = (sim[Symbol("Genetic Drift Simulation 1601627968-2")] .+ sim[Symbol("Genetic Drift Simulation 1601627968-3")] .+ sim[Symbol("Genetic_Drift_Simulation_1601627968-1")] .+
sim[Symbol("Genetic Drift Simulation 1601627968-4")] .+ sim[Symbol("Genetic Drift Simulation 1601627968-5")] .+ sim[Symbol("Genetic Drift Simulation 1601627968-6")] .+
sim[Symbol("Genetic Drift Simulation 1601627968-7")] .+ sim[Symbol("Genetic Drift Simulation 1601627968-8")] .+ sim[Symbol("Genetic Drift Simulation 1601627968-9")] .+
sim[Symbol("Genetic Drift Simulation 1601627968-10")]) ./ 10

# sample 5945 genes from 35.000 genes
set = sample(2:35000,5945;replace = false, ordered=true)
set = vcat([1],set)
sim2 = sim2[!,set]

sim3 = stack(sim2, Not(:Line1))
sim3[!,:Line1] .= Int64.(sim3[!,:Line1])

p3 = plot(sim3, x=:Line1, y=:value, color=:variable, Geom.line, Theme(key_position = :none, background_color=colorant"white"), Guide.xlabel("Generation"), Guide.ylabel("wtAF"), Guide.title("Simulated"), Coord.cartesian(ymin=0, ymax=0.8))
p3 = plot(sim3, x=:Line1, y=:value, color=:variable, Geom.line, Theme(key_position = :none, background_color=colorant"white"), Guide.xlabel("Generation"), Guide.ylabel("wtAF"), Guide.title("Simulated"))

pp = vstack(p1,p2,p3)
draw(PNG("/media/michael-uni/Uni/Allelfreq/Drift_calculation/comparison5Ksim_realData.png", 7inch, 7inch, dpi=300), pp)


## test if the simulated data and the real data of organic and conventinal ist equal by a zero inflated Model

## relevant for the Marker level
simreal = DataFrame(Generation = Int64[], Negbin = Float64[], ZeroInf = Float64[], System = String[])
for i in [3,12,16,22,23], j in ["organic", "conventional"]
                        if i == 3 && j == "organic"
                            continue
                        else
                        vx = v3[.&(v3.Generation.==i, v3.Env.==j),:wtAF]
                        vx = Int64.(floor.(vx .* 100))
                        vy = repeat([1], length(vx))
                        vxy = DataFrame(hcat(vx,vy))

                        sx = Int64.(floor.(sim3[sim3.Line1.==i,:value].* 100))
                        sy = repeat([2], length(sx))
                        sxy = DataFrame(hcat(sx,sy))

                        vs = vcat(vxy,sxy)

                        @rput vs
                        R"""
                         library(MASS)
                         library(pscl)
                         vs[,2] = as.factor(vs[,2])
                         m1 =  summary(zeroinfl(x1 ~ x2, vs, dist = "negbin"))
                         negbino = m1$coefficients$count[2,4]
                         zeroinf = m1$coefficients$zero[2,4]
                         """
                         @rget negbino zeroinf
                         push!(simreal, hcat(i,negbino, zeroinf, j))
                     end

end
CSV.write("D:/Allelfreq/Drift_calculation/MarkerToSim.csv", simreal; delim=",")

## gene levels
freq = readDict("D:/Allelfreq/ref_snps_dulicates_removed2020/Dict_HC_Genes")
stacker(freq, "F:/Uni-Arbeit/Allelfreq/genoytpes_wgs.txt")
genes = intersect(freq[:F3E1_Haplotypes].Gene_ID, freq[:F12E1_Haplotypes].Gene_ID, freq[:F12E2_Haplotypes].Gene_ID, freq[:F16E1_Haplotypes].Gene_ID, freq[:F16E2_Haplotypes].Gene_ID, freq[:F22E1_Haplotypes].Gene_ID, freq[:F22E2_Haplotypes].Gene_ID, freq[:F23E1_Haplotypes].Gene_ID, freq[:F23E2_Haplotypes].Gene_ID)

v1 = vcat(freq[:F3E1_Haplotypes][indexin(genes, freq[:F3E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F12E1_Haplotypes][indexin(genes, freq[:F12E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F12E2_Haplotypes][indexin(genes, freq[:F12E2_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F16E1_Haplotypes][indexin(genes, freq[:F16E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F16E2_Haplotypes][indexin(genes, freq[:F16E2_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F22E1_Haplotypes][indexin(genes, freq[:F22E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F22E2_Haplotypes][indexin(genes, freq[:F22E2_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F23E1_Haplotypes][indexin(genes, freq[:F23E1_Haplotypes].Gene_ID),:Allelfreq_Isr],
freq[:F23E2_Haplotypes][indexin(genes, freq[:F23E2_Haplotypes].Gene_ID),:Allelfreq_Isr])

v2 = repeat([3,12,12,16,16,22,22,23,23], inner= length(genes))
v5 = repeat(genes, 9)
v4 = repeat(vcat(repeat(["conventional"], length(genes)), repeat(["organic"], length(genes))), 4)
v4 = vcat(repeat(["conventional"], length(genes)), v4)

v3 = DataFrame(hcat(v1,v2,v5,v4))

v3[!,:x1] .= Float64.(v3[!,:x1])
v3[!,:x2] .= Int64.(v3[!,:x2])
v3[!,:x3] .= String.(v3[!,:x3])
v3[!,:x4] .= String.(v3[!,:x4])
rename!(v3, [:wtAF, :Generation, :Gene, :Env])

# sim data for all genes of F3 starter Haplotypes
f3 =  Float64.(freq[:F3E1_Haplotypes][indexin(genes, freq[:F3E1_Haplotypes].Gene_ID),:Allelfreq_Isr])
f3[1] = 0.1
ngen = 21
n=2
npop = 10000
ntime = 1000
nloci = length(f3)

npop = npop * n
drift = Dict()
#Simulation
Threads.@threads for k in collect(1:ntime)
  X = Matrix{Float64}(undef, ngen, nloci)
  X[1,:] = npop .* f3
    for j in collect(1:nloci), i in collect(2:ngen)
        X[i,j] = rand(Binomial(npop, X[i-1,j]/npop),1)[1]
        #=@rput j i X npop
        R"""
        b = rbinom(1,npop,prob=X[i-1,j]/npop)
        """
        @rget b
        X[i,j] = b =#
    end

    # Change table format and plot it
    X = DataFrame(X./npop)
    rename!(X, Symbol.(string.("Haplotype",1:nloci)))
    drift[k] = X
    println("iteration $k finished")
end


## get one value as 95% confidence interval
# rank the values for each replication , locus and generation wise
drift[:conf] = Matrix{Array{Float64,1}}(undef, ngen, nloci)

@showprogress for j in collect(1:nloci), k in collect(2:ngen)
    a = Float64[]
    for i in collect(1:ntime)
        push!(a, drift[i][k,j])
    end
    sort!(a)
    # 95% conf interval means the 11ths and 189th value
    drift[:conf][k,j] = a[Int64.([ntime*.05,ntime*.5,ntime*.95])]
end

drift[:conf] = DataFrame(drift[:conf][2:21,:])

# the goal is to plot the median value with the confidence interval over the generations
di = hcat(collect(4:23), drift[:conf], makeunique=true)
di = stack(di, Not(:x1))
di[!,:median] .= 0.0
di[!,:low] .= 0.0
di[!,:high] .= 0.0

@showprogress for i in 1:size(di,1)
                    di[i,:low] = di.value[i][1]
                    di[i,:median] = di.value[i][2]
                    di[i,:high] = di.value[i][3]
                end

deletecols!(di, :value)
rename!(di, :x1 => :Generation)

## add the real values for the conv. and organic in the plot - define confidence interval (90%)
real = DataFrame(Generation = Int64[], Environment = String[], Median = Float64[], Low = Float64[], High = Float64[])

for i in [3,12,16,22,23], j in ["conventional", "organic"]
    if i == 3 && j == "organic"
        continue
    end
    as = sort(v3[.&(v3.Env .== j, v3.Generation .== i),:wtAF])
    b = as[Int64.(floor.([length(as)*.05,length(as)*.95]))  ]
    d = mean(as)
    c = vcat(i,j,d,b)
    push!(real, c)
end
push!(real, real[1,:])
real[10,2] = "organic"
# add the dirft sim of 0.1 to the real table to facet
dd = di[1:20,:]
push!(dd, [3,"x1_1", 0.1,0.1,0.1])
rename!(dd, names(real))
dd[!,:Environment] .= "drift simulation"
real = vcat(real, dd)
# set low and high in the F3 in orgnaic and conventional to the median value
real[.&(real.Generation .==3, real.Environment .!= "drift simulation"),:Low] .= real[.&(real.Generation .==3, real.Environment .!= "drift simulation"),:Median]
real[.&(real.Generation .==3, real.Environment .!= "drift simulation"),:High] .= real[.&(real.Generation .==3, real.Environment .!= "drift simulation"),:Median]


@rput real
R"""
library(gridExtra)
library(ggplot2)
p1 = ggplot(real,aes(Generation,Median, color=Environment)) + scale_color_brewer(palette = "Set1") + theme_bw() + geom_line( size = 1) + xlab("") + ylab("wild type Allele frequency") + ggtitle("A") + theme(text = element_text(size=20))
p2 = ggplot(real,aes(Generation,Median, color=Environment, fill=Environment))+ theme_bw() + geom_ribbon(aes(ymin=Low, ymax=High, fill=Environment), alpha=0.3) +
    scale_fill_manual(values = c("organic"="olivedrab1", "conventional"="peru", "Drift Simulation"="gray40")) + ylab("wild type Allele frequency") + ggtitle("B") + theme(text = element_text(size=20))
p3 = grid.arrange(p1,p2)
ggsave("D:/Allelfreq/Drift_calculation/F3asStart_driftCalc/10Startvalue_drift_Comparison_Mean_95ConfInt_2.png",p3, width = 12, height = 12,unit="in", dpi=300)
"""



## calc stats
simreal = DataFrame(Generation = Int64[], Negbin = Float64[], ZeroInf = Float64[], System = String[], CompLevel = Symbol[])
@showprogress for i in [12,16,22,23], j in ["organic", "conventional"], k in [:median, :low, :high]

                        vx = v3[.&(v3.Generation.==i, v3.Env.==j),:wtAF]
                        vx = Int64.(floor.(vx .* 100))
                        vy = repeat([1], length(vx))
                        vxy = DataFrame(hcat(vx,vy))

                        sx = Int64.(floor.(di[di.Generation.==i,k].* 100))
                        sy = repeat([2], length(sx))
                        sxy = DataFrame(hcat(sx,sy))

                        vs = vcat(vxy,sxy)

                        @rput vs
                        R"""
                         library(MASS)
                         library(pscl)
                         vs[,2] = as.factor(vs[,2])
                         m1 =  summary(zeroinfl(x1 ~ x2, vs, dist = "negbin"))
                         negbino = m1$coefficients$count[2,4]
                         zeroinf = m1$coefficients$zero[2,4]
                         """
                         @rget negbino zeroinf
                         push!(simreal, hcat(i,negbino, zeroinf, j, k))
                     end

end

CSV.write("F:/Uni-Arbeit/Allelfreq/Drift_calculation/F3asStart_driftCalc/GeneToSim.csv", simreal; delim=",")
