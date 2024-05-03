genolist = "genoytpes_wgs.txt"
vcffile = "INDELcalling_sambamba.vcf"

ENV["R_HOME"]="C:\\Program Files\\R\\R-3.6.3"
using DelimitedFiles, DataFrames, Statistics, CSV, ProgressMeter, HypothesisTests, GLM, Crayons, Gadfly, Compose, CuArrays, RCall, StatsBase
using Cairo
using DataFrames, Statistics, CSV, ProgressMeter, RCall, Crayons

# a function to write Dictionary entries to disk space
function writeDict(disk, ram, how="text")
                             if how == "text"
                                           for i in keys(ram)
                                                            println(i)
                                                            if typeof(ram[i]) == DataFrame
                                                            CSV.write(string(disk,string(i),".txt"),ram[i];delim = '\t')
                                                            else
                                                              writedlm(string(disk,string(i),".txt"),ram[i], '\t')
                                                            end
                                            end
                             elseif how == "csv"
                                     for i in keys(ram)
                                                      println(i)
                                                      if typeof(ram[i]) == DataFrame
                                                      CSV.write(string(disk,string(i),".csv"),ram[i];delim = ',')
                                                    else
                                                      writedlm(string(disk,string(i),".txt"),ram[i], ',')
                                                     end
                                     end

                             end
 end

function readDict(loc, ok = "txt")
                        p = Dict()
                        if ok != "std"
                                            for i in readdir(loc)[occursin.(ok,readdir(loc))]
                                                                 println(i)
                                                                 r = string(loc,"/", i)
                                                                 r2 = CSV.read(r, DataFrame)

                                                                 p[Symbol(split(i,".")[1])] = r2
                                            end
                                            return p
                        else
                                            for i in readdir(loc)[.!(occursin.(r"f[0-9]", readdir(loc)))]
                                                                 println(i)
                                                                 r = string(loc,"/", i)
                                                                 r2 = CSV.read(r, DataFrame)

                                                                 p[Symbol(split(i,".")[1])] = r2
                                            end
                                            return p

                        end
        end

function Allelfreqcalling(vcffile, genolist, minreaddepth=1, minqual=30, maxqual=1000,posE1=1, posE2=2)

                                                                                                          #=function Rread()
                                                                                                                                                                                                                        @rput vcffile
                                                                                                                                                                                                                        R"""
                                                                                                                                                                                                                        library(data.table)
                                                                                                                                                                                                                        vcf = fread(vcffile)
                                                                                                                                                                                                                        """
                                                                                                                                                                                                                        vcf = rcopy(R"""vcf""")
                                                                                                                                                                                                                        return vcf
                                                                                                                                                                                                         end

                                                                                                          vcf2 = Rread()  =#
                                                                                                          vcf2 = @time CSV.read(vcffile, delim="\t", comment="#", header=false)
                                                                                                          glist = CSV.read(genolist, delim="\t")
                                                                                                          nam = glist[!,1]
                                                                                                          pl = ["Chr", "Pos","id", "Ref", "Alt", "Qual", "Filter", "Info", "Format"]
                                                                                                          # set the name of the snp file
                                                                                                          append!(pl, nam)
                                                                                                          rename!(vcf2,Symbol.(pl))

                                                                                                                                                                                                        #= extract Allelfrequency from vcf files from infinte many GENOTYPES
                                                                                                                                                                                                        vcf = CSV.read(vcffile, delim="\t", comment="#")
                                                                                                                                                                                                        the problem here is that it was not possible to link nam with th first 9 colnames by vcat to produce an Array{Any,1}, instead it turned out to be an Array{Any,2}
                                                                                                                                                                                                        extract the information of the column names from an external file provided
                                                                                                                                                                                                        pl = ["Chr", "Pos","id", "Ref", "Alt", "Qual", "Filter", "Info", "Format"]
                                                                                                                                                                                                        set the name of the snp file
                                                                                                                                                                                                        append!(pl, nam)
                                                                                                                                                                                                        names!(vcf,Symbol.(pl))
                                                                                                                                                                                                        =#

                                                                                                                                                                                                        ### Best is to create Dictionaries for each genotype unique and afterwards extract the information requiered from each of them, e.g. readdepth
                                                                                                                                                                                                        # empty dataframe for each genotype initialized
                                                                                                        println("Data is loaded successfully")
                                                                                                        frequency = Dict()

                                                                                                        ## the file for linking to the individual genotype dicts is processed according to the needs of CuArrays - make all strings become integers
                                                                                                        frequency[:General] = Matrix(vcf2[!,[1,2,4,5,6]])

                                                                                                        # transfer the information to a seperate library for each genotype
                                                                                                        # split the individual colum for each library into 3 columns giving info about the call type, REf readcount and ALT readcount
                                                                                                        Threads.@threads  for j in nam
                                                                                                                     p = vcf2[!,Symbol(j)]
                                                                                                                     k = String[]
                                                                                                                     m = Int64[]
                                                                                                                     n = Int64[]

                                                                                                                    for i in 1:size(p,1)
                                                                                                                                         push!(k, string(split(p[i],":")[1]))
                                                                                                                                         push!(m,parse(Int64,string(split(split(p[i],":")[3],",")[1])))
                                                                                                                                         push!(n,parse(Int64,string(split(split(p[i],":")[3],",")[2])))
                                                                                                                     end

                                                                                                                     frequency[Symbol(j)] = hcat(k,m,n)

                                                                                                        end
                                                                                                         println("seperation in Dictionaries completed")

                                                                                                                        ############## experimental step
                                                                                                                        # this step should introduce an correction of the snp calls generated - it is reported that the reference base has a higher change to be alligend than the alt base has
                                                                                                                        # therefore a correction will be introduced where the over all Loc calculated freq is estimated and the corrected to 0.5
                                                                                                                        println("normalisation of the frequency calls is calculated")
                                                                                                                        re = repeat([0.5],inner=length(nam)-2)
                                                                                                                      Threads.@threads  for i in 3:length(nam)

                                                                                                                                             re[i-2] = sum(frequency[Symbol(nam[i])][:,3])/sum(frequency[Symbol(nam[i])][:,2] .+ frequency[Symbol(nam[i])][:,3])

                                                                                                                                            #=  ACTUAL CALLCULATION IS PEROFORMED IN THE HAPLOTYPE CALLING FOR THOSE WITH ONLY ONE SNP IN A GENE
                                                                                                                                            freq[Symbol(nam[i])][:countref] = floor.(freq[Symbol(nam[i])][:countref] .* (.5/correctref))
                                                                                                                                            freq[Symbol(nam[i])][:countalt] = floor.(freq[Symbol(nam[i])][:countalt] .* (.5/correctalt))
                                                                                                                                            # recalculate the readdepth
                                                                                                                                            freq[Symbol(nam[i])][:readcount] = freq[Symbol(nam[i])][:countalt] .+ freq[Symbol(nam[i])][:countref] =#
                                                                                                                                           end
                                                                                                                                           frequency[:correction_alt] = vcat(fill(0,2),re)


                                                                                                # remove homogen snps in between the parents
                                                                                                                                            println("Process filtering..")
                                                                                                                                            # highlight the snps where the parents have low coverage - convert these to a "./."
                                                                                                                      Threads.@threads     for n in nam[[posE1, posE2]]
                                                                                                                                                                        # calulate the readdepth
                                                                                                                                                                        p = frequency[Symbol(n)][:,2] .+ frequency[Symbol(n)][:,3]
                                                                                                                                                                        frequency[Symbol(n)] = hcat(frequency[Symbol(n)],p)
                                                                                                                                                                        frequency[Symbol(n)] = Array(frequency[Symbol(n)])

                                                                                                                                                                        for i in 1:size(frequency[Symbol(n)],1)
                                                                                                                                                                                                             if frequency[Symbol(n)][i,4] < minreaddepth
                                                                                                                                                                                                             frequency[Symbol(n)][i, 1] = "./."
                                                                                                                                                                                                             end
                                                                                                                                                                         end
                                                                                                                                           end
                                                                                                                                            # check which snps are equal in between the parents and HETEROZYGOT - label them to remove these from ALL Dicts
                                                                                                                                            sel = (.&(frequency[Symbol(nam[posE1])][:,1] .!= "./.", frequency[Symbol(nam[posE1])][:,1] .!= frequency[Symbol(nam[posE2])][:,1], frequency[Symbol(nam[posE2])][:,1] .!= "./.", frequency[Symbol(nam[posE1])][:,1] .!= "0/1", frequency[Symbol(nam[posE2])][:,1] .!= "0/1"))
                                                                                                                     Threads.@threads     for j in nam; frequency[Symbol(j)] = frequency[Symbol(j)][sel,:]; end

                                                                                                                                            frequency[:General] = frequency[:General][sel,:] # the same for the General file
                                                                                                                                            sel=nothing


                                                                                                                                        # filter the Quality of the SNP calls - everything lower qual will be removed
                                                                                                                      Threads.@threads  for i in nam
                                                                                                                                                                    frequency[Symbol(i)] = frequency[Symbol(i)][.&(frequency[:General][:,5] .> minqual, frequency[:General][:,5] .< maxqual),:]
                                                                                                                                                     end

                                                                                                                                            frequency[:General] = frequency[:General][.&(frequency[:General][:,5] .> minqual, frequency[:General][:,5] .< maxqual),:] # the same for the General file


                                                                                                                                          #### get rid of the snps with more than two alleles
                                                                                                                                           # create a new dict to store information from all dict allel calls in one table
                                                                                                                                           frequency[:Alleles] = frequency[:General][:,[1,2]]
                                                                                                                                           for i in 1:size(nam,1); frequency[:Alleles] = hcat(frequency[:Alleles],frequency[Symbol(nam[i])][:,1]); end

                                                                                                                                         frequency[:Alleles] = frequency[:Alleles][:,3:size(frequency[:Alleles],2)]
                                                                                                                                         # save information of the number of different alleles in one value per column
                                                                                                                                         kind = Bool[]
                                                                                                                                         for i in 1:size(frequency[:Alleles],1); append!(kind,all(.&(frequency[:Alleles][i,:] .!= "1/2", frequency[:Alleles][i,:] .!= "0/2", frequency[:Alleles][i,:] .!= "2/2", frequency[:Alleles][i,:] .!= "2/3", frequency[:Alleles][i,:] .!= "1/3", frequency[:Alleles][i,:] .!= "0/3"))
                                                                                                                                                 ) ; end
                                                                                                                                         frequency[:Alleles] = kind

                                                                                                                Threads.@threads        for i in nam
                                                                                                                                                     frequency[Symbol(i)] = frequency[Symbol(i)][frequency[:Alleles],:]
                                                                                                                                           end
                                                                                                                                         # remove the information also for the "General" file
                                                                                                                                         frequency[:General] = frequency[:General][frequency[:Alleles],:]

                                                                                                                                         delete!(frequency, :Alleles)
                                                                                                                                         println("Loci with unnormal Allel pattern removed")


                                                                                                                                        freq = frequency

                                                                                                                                        ### when the coverage is lower 5, convert :Alleles to ./. - calculation of readdepth necessary
                                                                                                                                        # parents should have another consideration that is set as standard, so that only the pools are effected by the change of the minreaddepth setting
                                                                                                                                        # the minimum for parents will be set to 5 statically
                                                                                                                                        nam_pool = nam[.!(.|((nam .==nam[posE1]), (nam .== nam[posE2])))]
                                                                                                                                        nam_parents = nam[.|((nam .==nam[posE1]), (nam .== nam[posE2]))]
                                                                                                                                         # progenies
                                                                                                                                      Threads.@threads for n in nam_pool
                                                                                                                                                                        # calulate the readdepth
                                                                                                                                                                        p = freq[Symbol(n)][:,2].+ freq[Symbol(n)][:,3]
                                                                                                                                                                        freq[Symbol(n)] = hcat(freq[Symbol(n)],p)
                                                                                                                                                                        freq[Symbol(n)] = Array(freq[Symbol(n)])
                                                                                                                                                                         #= keep everything in the library and do not remove loci with low coverage at the parents
                                                                                                                                                                         for i in 1:size(freq[Symbol(n)],1)
                                                                                                                                                                                                          if freq[Symbol(n)][i,4] < minreaddepth
                                                                                                                                                                                                          freq[Symbol(n)][i, 1] = "./."
                                                                                                                                                                                                          end
                                                                                                                                                                        end
                                                                                                                                                                        =#
                                                                                                                                                         end


                                                                                                                                         # calculate the allel freq for each progeny
                                                                                                                                      println("prepare empty frames for Allel freq calculation")

                                                                                                                Threads.@threads      for i in 3:length(nam); freq[Symbol(nam[i])]= hcat(freq[Symbol(nam[i])], Float32.(repeat(2:2, size(freq[Symbol(nam[i])],1))), Float32.(repeat(2:2, size(freq[Symbol(nam[i])],1)))); end
                                                                                                                                                        println("calculate Allel frequency for each progeny")
                                                                                                                                                        Threads.@threads for i in 3:length(nam)

                                                                                                                                                            # gpu based calculation - replaced by threaded code
                                                                                                                                                            #=  println(nam[i])
                                                                                                                                                              freq[Symbol(nam[i])] = hcat(freq[:General][:,1:2], freq[Symbol(nam[i])])
                                                                                                                                                              a1 = freq[Symbol(nam[i])][freq[Symbol(nam[posE1])][:,1].=="0/0",:]
                                                                                                                                                              a1[:,7] .= cu(Int64.(a1[:,4])) ./ cu(Int64.(a1[:,6]))
                                                                                                                                                              a1[:,8] .= cu(Int64.(a1[:,5])) ./ cu(Int64.(a1[:,6]))
                                                                                                                                                              #a1 = hcat(a1,fill("g", 6917))
                                                                                                                                                              a2 = freq[Symbol(nam[i])][freq[Symbol(nam[posE1])][:,1].=="1/1",:]
                                                                                                                                                              a2[:,7] .= cu(Int64.(a2[:,5])) ./ cu(Int64.(a2[:,6]))
                                                                                                                                                              a2[:,8] .= cu(Int64.(a2[:,4])) ./ cu(Int64.(a2[:,6]))
                                                                                                                                                              #a2 = hcat(a2,fill("i", 3422))
                                                                                                                                                              a3 = DataFrame(vcat(a1,a2))
                                                                                                                                                              sort!(a3, (:x1,:x2))
                                                                                                                                                              freq[Symbol(nam[i])] = Matrix(a3[:,3:8]) =#

                                                                                                                                                              println(nam[i])
                                                                                                                                                              freq[Symbol(nam[i])] = hcat(freq[:General][:,1:2], freq[Symbol(nam[i])])
                                                                                                                                                              a1 = freq[Symbol(nam[i])][freq[Symbol(nam[posE1])][:,1].=="0/0",:]
                                                                                                                                                              a1[:,7] .= a1[:,4] ./ a1[:,6]
                                                                                                                                                              a1[:,8] .= a1[:,5] ./ a1[:,6]
                                                                                                                                                              #a1 = hcat(a1,fill("g", 6917))
                                                                                                                                                              a2 = freq[Symbol(nam[i])][freq[Symbol(nam[posE1])][:,1].=="1/1",:]
                                                                                                                                                              a2[:,7] .= a2[:,5] ./ a2[:,6]
                                                                                                                                                              a2[:,8] .= a2[:,4] ./ a2[:,6]
                                                                                                                                                              #a2 = hcat(a2,fill("i", 3422))
                                                                                                                                                              a3 = DataFrame(vcat(a1,a2))
                                                                                                                                                              sort!(a3, (:x1,:x2))
                                                                                                                                                              freq[Symbol(nam[i])] = Matrix(a3[:,3:8])

                                                                                                                                                           end

                                                                                                                                    #=  @time   for i in 3:length(nam)
                                                                                                                                                                println(nam[i])
                                                                                                                                                               for j in 1:size(freq[Symbol(nam[1]),],1)

                                                                                                                                                        if freq[Symbol(nam[i])][:,1][j] != "./."
                                                                                                                                                                                                if freq[Symbol(nam[posE1])][j,1] == "0/0"
                                                                                                                                                                                                                freq[Symbol(nam[i])][j,5] = freq[Symbol(nam[i])][j,2] / freq[Symbol(nam[i])][j,4]
                                                                                                                                                                                                                freq[Symbol(nam[i])][j,6] = freq[Symbol(nam[i])][j,3] / freq[Symbol(nam[i])][j,4]
                                                                                                                                                                                                else
                                                                                                                                                                                                                freq[Symbol(nam[i])][j,5] = freq[Symbol(nam[i])][j,3] / freq[Symbol(nam[i])][j,4]
                                                                                                                                                                                                                freq[Symbol(nam[i])][j,6] = freq[Symbol(nam[i])][j,2] / freq[Symbol(nam[i])][j,4]
                                                                                                                                                                                                end
                                                                                                                                                        else
                                                                                                                                                                                                                freq[Symbol(nam[i])][j,5] = NaN
                                                                                                                                                                                                                freq[Symbol(nam[i])][j,6] = NaN
                                                                                                                                                        end
                                                                                                                                                   end
                                                                                                                                   end
                                                                                                                                    =#

                                                                                                                                         ## calculate the allelfreq for both parents
                                                                                                                                        Threads.@threads for i in 3:length(nam)
                                                                                                                                                                  x = Symbol(nam[i])
                                                                                                                                                                  # subset by > -1 to get rid of the NaN values
                                                                                                                                                                  FREQ = sum((freq[Symbol(nam[i])][freq[Symbol(nam[i])][:,5].>-1,6] .* freq[Symbol(nam[i])][freq[Symbol(nam[i])][:,5].>-1,4]) ./ sum(freq[Symbol(nam[i])][freq[Symbol(nam[i])][:,5].>-1,4]))

                                                                                                                                                                 print(Crayon(foreground = :blue, bold = false),"Allelfreq over all loci of $x: $FREQ\n")

                                                                                                                                                                end

                                                                                                                                                 #freq[Symbol(nam[i])][Symbol("Allelfreq",nam[posE1])][j]
                                                                                                                                                 pl = ["Chr", "Pos", "Ref", "Alt", "Qual", "Alleles", "Refcount", "Altcount", "Readcount"]
                                                                                                                                                 pl2 =  [string("Allelfreq_",nam[posE1]), string("Allelfreq_",nam[posE2])]


                                                                                                                                                 # make ever array a dataframe
                                                                                                                                                 Threads.@threads  for i in nam
                                                                                                                                                                                 freq[Symbol(i)] = hcat(freq[:General],freq[Symbol(i)])
                                                                                                                                                                                 freq[Symbol(i)] = DataFrame(freq[Symbol(i)])
                                                                                                                                                                                 if size(freq[Symbol(i)],2) == 9
                                                                                                                                                                                                                 names!(freq[Symbol(i)], Symbol.(pl))
                                                                                                                                                                                 elseif size(freq[Symbol(i)],2) == 11
                                                                                                                                                                                                                 names!(freq[Symbol(i)] ,Symbol.(vcat(pl,pl2)))

                                                                                                                                                                                end
                                                                                                                                                                 end

                                                                                                                                                        return freq

                                                                                                                                                end


#= relies on the output of function Allelfreqcalling
the following function is designed to extract the information of related genes to the detected SNPS with its Names, functions and Go terms.
The output id written to an other dict entry in the Allelfreqcalling output file, so that all information is stored in one single dictionary =#
function genloc(info, freqfile, genolist, genetic=false,  extention=true)
                                              # gff file might not even be needed
                                              println("read file")
                                              label = CSV.read(genolist, delim="\t")[1,1]
                                              infofile = CSV.read(info,delim="\t",header=1)

                            # two possible ways can be gone - genetic map or classic physical map
                            if genetic == true
                                                print(Crayon(foreground = :yellow, bold = false),"Genetic position is used \n")
                                                ni = names(infofile)
                                                inf = Matrix(infofile)
                                                np = names(freqfile[Symbol(label)])
                                                pp = Matrix(freqfile[Symbol(label)])
                                                qerf = Dict()
                    Threads.@threads             for i in unique(pp[:,1])
                                                                            # subset
                                                                            println(i)
                                                                            ofni = inf[inf[:,2].== i ,:]
                                                                            qerf[Threads.threadid()] = pp[pp[:,1].== i , 1:2]
                                                                            # add the columns from the infotable to the snp table
                                                                            qerf[Threads.threadid()] = hcat(qerf[Threads.threadid()], Array{Union{Nothing, Any}}(nothing, size(qerf[Threads.threadid()],1), size(ni,1)))
                                                                            # create a "trash" file containing some value to set wenn no gene is matching
                                                                            pl = [".",".", 0, 0, 0.0, 0, ".", ".", ".", ".", ".", ".", "."]
                                                                            # write pl to all rows
                                                                            for n in 1:size(inf,2); qerf[Threads.threadid()][:,n+2] .= pl[n]; end
                                                                            # in the first iteration, give the column names to the Dict of Allelfreqcalling
                                                                            if haskey(freqfile, :Location)==false; freqfile[:Location] = Array{Union{Nothing, Any}}(nothing,0, 15); end
                                                                            #if haskey(freqfile, :Location)==false; freqfile[:Location] = DataFrame(describe(qerf)[:eltype], names(qerf), 0); end
                                                                            # do the check and write the information to an new DataFrame with position information of the snp and the gene both

                                                                             for j in 1:size(qerf[Threads.threadid()],1), k in 1:size(ofni,1)

                                                                                                              if qerf[Threads.threadid()][j,2] > ofni[k,3] && qerf[Threads.threadid()][j,2] < ofni[k,4]
                                                                                                                       for u in 1:size(inf,2)
                                                                                                                                       qerf[Threads.threadid()][j,u+2] = ofni[k,u]
                                                                                                                       end
                                                                                                                       continue
                                                                                                              end
                                                                                          end
                                               end
                                               println("done")
                                                             # write the information from the loop to the Dict entry
                                                             for i in keys(qerf); freqfile[:Location] = vcat(freqfile[:Location],qerf[i]); end
                                                             # convert array to a dataframe
                                                             freqfile[:Location] = DataFrame(freqfile[:Location])
                                                             rename!(freqfile[:Location], vcat(np[1:2], ni),makeunique=true)
                                                             # rearrange the file according to the llcation file coming from the physical positions
                                                             freqfile[:Location] = freqfile[:Location][:,[1,2,4,5,6,3,7,8,9,10,11,12,13,14,15]]
                                                             println("Allocation of SNPs to the Genes successfully completed")

                            else # physical map
                                              print(Crayon(foreground = :yellow, bold = false),"Physical position is used ")
                                              # the chromosom, start and end have to be splited up
                                              # two new split functions have been created just for this problem don here
                                              function split3(file,sep1=":",sep2="-")
                                                              es = vcat(split.(file,sep1)...)
                                                              es2 = es[collect(2:2:length(es))]
                                                              es21 = vcat(split.(es2,sep2)...)[collect(1:2:length(es))]
                                                              es22 = vcat(split.(es2,sep2)...)[collect(2:2:length(es))]
                                                              em = DataFrame(Array{Union{Nothing, Any}}(nothing, length(es2), 3))
                                                              em[!,1] = es[collect(1:2:length(es))]
                                                              em[!,2] = es21
                                                              em[!,3] = es22
                                                              return em
                                              end

                                              # split the position column into three columns
                                              println("Split Position Column")
                                              inq = split3(infofile[!,2],":","-")
                                              inq[!,:x2] = parse.(Int,inq[!,:x2]) # convert to an integer from string
                                              inq[!,:x3] = parse.(Int,inq[!,:x3])



                                              #= extend the gene boundaries, due to the fact that with GBS data most of the SNPs will be located outside of genes
                                                 the rule for extention follows the same rule as has been applied to the extention of wheat with different extention sizes
                                              =#
                                              if extention == true
                                                      println("Extention of Gene Bounds is Performed \n")
                                                      inq[!, :start] .= 1
                                                      inq[!, :end] .= 1
                                                      Threads.@threads for i in 1:size(inq,1)
                                                                              if i == 1; inq[!,:start][i] = inq[!,:x2][i]; end
                                                                              if i == size(inq,1)
                                                                                                  inq[!,:end][i] = inq[!,:x3][i]
                                                                              else
                                                                                      w = inq[!,:x2][i+1] - inq[!,:x3][i]
                                                                                      if w > 1000
                                                                                              inq[!,:start][i+1] = floor(inq[!,:x2][i+1] - w * 0.45) # convert an float value to an integer
                                                                                              inq[!,:end][i] = floor(inq[!,:x3][i] + w * 0.45)
                                                                                      else
                                                                                              inq[!,:end][i] = inq[!,:x3][i]
                                                                                              inq[!,:start][i+1] = inq[!,:x2][i+1]
                                                                                      end
                                                                              end
                                                                      end
                                                              else
                                                                      inq[!,:start] .= inq[!,:x2]
                                                                      inq[!,:end] .= inq[!,:x3]
                                                              end
                                              # test if there are overlappings - the overlaps results from the original ref data, extention of overlaps is due to the code above not intended
                                              # overlaps might cause problems when trying to locate the snp correctly to a unique gene
                                       ### OVERLAPS will be kept if the following 3 '#' are existing #######
                                       #         inq[:check] = vcat(false,inq[2:size(inq)[1],:start] .> inq[1:size(inq)[1]-1,:end])
                                              #link the intermed file with the infofile again
                                              infofile = hcat(inq, infofile)
                                              # if its a check == false and a low confidence gene, then remove it from data set
                                       #         infofile = infofile[.&(infofile[:check].!=false, infofile[:confidence_class].!="LC_u"),:]

                                       #         deletecols!(infofile, :check)
                                              select!(infofile, Not(:x2))
                                              select!(infofile, Not(:x3))

                                              # link the SNPs to the gene - it only has to be done once, due to the fact that all libraries of SNPs reffer to the same loci
                                              # saved in an additional dictionary
                                              #= procedure:
                                              subset the chromosoms, so that no mismatch with an other chromosom is possible
                                              =#
                                                ni = names(infofile)
                                                inf = Matrix(infofile)
                                                np = names(freqfile[Symbol(label)])
                                                pp = Matrix(freqfile[Symbol(label)])
                                                qerf = Dict()
                  Threads.@threads              for i in unique(infofile[:,1])
                                                                            # subset
                                                                            println(i)
                                                                            ofni = inf[inf[:,1].== i ,:]
                                                                            qerf[Threads.threadid()] = pp[pp[:,1].== i , 1:2]
                                                                            # add the columns from the infotable to the snp table
                                                                            qerf[Threads.threadid()] = hcat(  qerf[Threads.threadid()], Array{Union{Nothing, Any}}(nothing, size(qerf[Threads.threadid()],1), size(ni,1)))
                                                                            # create a "trash" file containing some value to set wenn no gene is matching
                                                                            pl = [".", 0, 0, ".",".",".",".",".",".","." ]
                                                                            # write pl to all rows
                                                                            for n in 1:size(inf,2); qerf[Threads.threadid()][:,n+2] .= pl[n]; end
                                                                            # in the first iteration, give the column names to the Dict of Allelfreqcalling
                                                                            if haskey(freqfile, :Location)==false; freqfile[:Location] = Array{Union{Nothing, Any}}(nothing,0, 12); end
                                                                            #if haskey(freqfile, :Location)==false; freqfile[:Location] = DataFrame(describe(qerf)[:eltype], names(qerf), 0); end
                                                                            # do the check and write the information to an new DataFrame with position information of the snp and the gene both
                                                                             for j in 1:size(qerf[Threads.threadid()],1), k in 1:size(ofni,1)

                                                                                                              if   qerf[Threads.threadid()][j,2] > ofni[k,2] &&   qerf[Threads.threadid()][j,2] < ofni[k,3]
                                                                                                                       for u in 1:size(inf,2)
                                                                                                                                         qerf[Threads.threadid()][j,u+2] = ofni[k,u]
                                                                                                                       end
                                                                                                                       continue
                                                                                                              end
                                                                            end
                                               end

                                                             # write the information from the loop to the Dict entry
                                                             for i in keys(qerf); freqfile[:Location] = vcat(freqfile[:Location],qerf[i]); end
                                                             # convert array to a dataframe
                                                             freqfile[:Location] = DataFrame(freqfile[:Location])
                                                             rename!(freqfile[:Location], vcat(np[1:2], ni))
                                                             println("Allocation of SNPs to the Genes successfully completed")

                                      end # genetic of physical map

                                      print(Crayon(foreground = :blue, bold = false),"Location file is sorted.. \n")
                                      sort!(freqfile[:Location], (:Chr, :Pos))

                                 end # genloc function

######## HAPLOTYPING ############
# create haplotypes of the snps using the gene bound information
# haplotype have been proven to give a higher consistency in the calling of the correct allelfreq
# for the haplotype calculation 2 ways can be applyed, the first is refering to the readdepth given, so higher rd will contribute more to the Frequency
# second way weights all snps equally, not depending if the are unweighted in readdepth
# this options can be set by using weighted = true for first (default) or false for unweighted
function haplotyping(freqfile,genolist,haplotype="gene", weighted=true)
                                                        nam = CSV.read(genolist, delim='\t')[!,:Name]

                                                        # get a unique list of all GO terms without "none" and "." ; lstrip - remove space from front if exists; vcat(split) - get a vector kind list from the GO terms in the Location subdict
                                                        if haplotype == "gene"
                                                                                hap = unique(freqfile[:Location][!,:Gene_ID])
                                                                                search = Symbol("Gene_ID")
                                                        elseif haplotype == "GO"
                                                                                hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...))[occursin.("GO",unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...)))])
                                                                                search = Symbol("GO-IDs")
                                                        elseif haplotype == "Marker"
                                                                                hap = unique(freqfile[:Location][!,:Marker])
                                                                                search = Symbol("Marker")
                                                        elseif haplotype == "PF"
                                                                                hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...))[occursin.("PF",unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...)))])
                                                                                search = Symbol("PFAM-IDs")
                                                        else
                                                        println("wrong haplotype specified, chose 'gene', 'GO' or 'PF' ")

                                                        end

                                                        # create a new dictionary for each progeny
                                                        for k in 3:length(nam); freqfile[Symbol(string(nam[k], "_Haplotypes"))] = DataFrame(Gene_ID=String[], Allelfreq_Isr=Float64[], Allelfreq_Golf=Float64[], SNPcount=Int64[], Readcount=Int64[], MinFreq=Float64[], MaxFreq=Float64[]); end
                                                                             # the haplotype file information coming with info about the gene, function and FRequency of this particular haplotype
                                                                             println("Make Haplotypes using the Bounds calculated with Genloc..")
                                    ###########################################################


                                    # this part only preformes the haplotyping by gene
                                    if haplotype == "gene" || haplotype == "Marker"
                                       for j in 3:length(nam)
                                                          dat = hcat(freqfile[Symbol(nam[j])],freqfile[:Location][!,Symbol(search)])
                                                                             function callsnps(dat,k)
                                                                                                     m = Int32[]
                                                                                                     for i in 1:nrow(dat)
                                                                                                                         if dat[i,:x1] != k
                                                                                                                         append!(m,i)
                                                                                                                         break
                                                                                                                         end
                                                                                                                     end
                                                                                                   if  !(isempty(m))
                                                                                                      dat2 = dat[1:m[1]-1,:]
                                                                                                      deleterows!(dat, 1:m[1]-1) # remove the gene snps that have already been worked on
                                                                                                    else
                                                                                                      dat2 = dat[1:nrow(dat),:]
                                                                                                    end

                                                                                                 return dat2
                                                                            end

                                                                            dat = dat[dat[!,:x1].!= ".",:] # remove not annotated snps from the list, highlighted by a "." in the x1 gene column
                                                                            sort!(dat, (:Chr, :x1))
                                                                            kap = unique(dat[!,:x1])
                                                                            println(nam[j])
                                                                            @showprogress for k in kap
                                                                                                      qwe = callsnps(dat,k)
                                                                                                      qwe = qwe[.!(isnan.(qwe[:,:11])),:]
                                                                                                      if isempty(qwe) == false
                                                                                                             if size(qwe,1) > 1
                                                                                                                         if weighted == true
                                                                                                                                 p =  sum(qwe[!,:Readcount] .* qwe[!,:Allelfreq_Isr])/sum(qwe[!,:Readcount])
                                                                                                                                 q = 1 - sum(qwe[!,:Readcount] .* qwe[!,:Allelfreq_Isr])/sum(qwe[!,:Readcount]) # Golf freq
                                                                                                                                 wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                 mi = minimum(qwe[!,:Allelfreq_Isr]) # min Allelfreq isr of all markers for isr
                                                                                                                                 mx = maximum(qwe[!,:Allelfreq_Isr]) # max Allelfreq of all markers for isr
                                                                                                                                 push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[k,p,q,size(qwe,1),wi, mi, mx])
                                                                                                                                 #append!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],asd)
                                                                                                                         else # use the none weighted version where each snp contributes equally
                                                                                                                                 p = sum(qwe[!,:Allelfreq_Isr]) ./ size(qwe,1) # Isr freq
                                                                                                                                 q = 1 - sum(qwe[!,:Allelfreq_Isr]) ./ size(qwe,1) # Golf freq
                                                                                                                                 wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                 mi = minimum(qwe[!,:Allelfreq_Isr]) # min Allelfreq isr of all markers for isr
                                                                                                                                 mx = maximum(qwe[!,:Allelfreq_Isr]) # max Allelfreq of all markers for isr
                                                                                                                                 push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[k,p,q,size(qwe,1),wi, mi,mx])
                                                                                                                        end
                                                                                                             else
                                                                                                                 # calculate the correction of the Freq Calling based on the correction factor that has been saved in the freqfile dict "correction"

                                                                                                                 # make sure the Alt Alleles will be corrected - it is not directly constructable if Golf or ISR are the Ref Base call
                                                                                                                 if qwe[:Allelfreq_Isr] == qwe[!,:Altcount] ./ qwe[!,:Readcount]
                                                                                                                                         p = qwe[!,:Allelfreq_Isr]
                                                                                                                                         q = 1 - qwe[!,:Allelfreq_Isr][1]
                                                                                                                                         wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                         mi = 0.0
                                                                                                                                         mx = 0.0
                                                                                                                                         push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[k,p[1],q[1],size(qwe,1),wi, mi, mx])
                                                                                                                 else
                                                                                                                                         p = qwe[!,:Allelfreq_Isr]
                                                                                                                                         q = 1 - qwe[!,:Allelfreq_Isr][1]
                                                                                                                                         wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                         mi = 0.0
                                                                                                                                         mx = 0.0
                                                                                                                                         push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[k,p[1],q[1],size(qwe,1),wi,mi,mx])
                                                                                                                 end
                                                                                                             end
                                                                                                         end

                                                                                                  end # of kap

                                                                                end # of nam

                                                        else
                                                        Threads.@threads       for j in 3:length(nam)
                                                                                      for i in hap
                                                                                                           # the alogrhytm is following the list of genes, so all entries in the list and then alloaction the information to the single pools
                                                                                                           qwe = freqfile[Symbol(nam[j])][freqfile[:Location][search].== i,:] # get all snps in the tested gene
                                                                                                           filter!(row -> !isnan(row[:Allelfreq_Isr]),qwe) # if the Freq contains any NaN, remove these
                                                                                                           # calculate the frequency based on which method has been choosen
                                                                                                                if isempty(qwe) == false
                                                                                                                     if size(qwe,1) > 1
                                                                                                                                 if weighted == true
                                                                                                                                         p =  sum(qwe[!,:Readcount] .* qwe[!,:Allelfreq_Isr])/sum(qwe[!,:Readcount])
                                                                                                                                         q = 1 - sum(qwe[!,:Readcount] .* qwe[!,:Allelfreq_Isr])/sum(qwe[!,:Readcount]) # Golf freq
                                                                                                                                         wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                         mi = minimum(qwe[!,:Allelfreq_Isr]) # min Allelfreq isr of all markers for isr
                                                                                                                                         mx = maximum(qwe[!,:Allelfreq_Isr]) # max Allelfreq of all markers for isr
                                                                                                                                         push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[i,p,q,size(qwe,1),wi, mi, mx])
                                                                                                                                         #append!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],asd)
                                                                                                                                 else # use the none weighted version where each snp contributes equally
                                                                                                                                         p = sum(qwe[!,:Allelfreq_Isr]) ./ size(qwe,1) # Isr freq
                                                                                                                                         q = 1 - sum(qwe[!,:Allelfreq_Isr]) ./ size(qwe,1) # Golf freq
                                                                                                                                         wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                         mi = minimum(qwe[!,:Allelfreq_Isr]) # min Allelfreq isr of all markers for isr
                                                                                                                                         mx = maximum(qwe[!,:Allelfreq_Isr]) # max Allelfreq of all markers for isr
                                                                                                                                         push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[i,p,q,size(qwe,1),wi, mi,mx])
                                                                                                                                end
                                                                                                                     else
                                                                                                                         # calculate the correction of the Freq Calling based on the correction factor that has been saved in the freqfile dict "correction"

                                                                                                                         # make sure the Alt Alleles will be corrected - it is not directly constructable if Golf or ISR are the Ref Base call
                                                                                                                         if qwe[:Allelfreq_Isr] == qwe[!,:Altcount] ./ qwe[!,:Readcount]
                                                                                                                                                 p = qwe[!,:Allelfreq_Isr]
                                                                                                                                                 q = 1 - qwe[!,:Allelfreq_Isr][1]
                                                                                                                                                 wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                                 mi = 0.0
                                                                                                                                                 mx = 0.0
                                                                                                                                                 push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[i,p[1],q[1],size(qwe,1),wi, mi, mx])
                                                                                                                         else
                                                                                                                                                 p = qwe[!,:Allelfreq_Isr]
                                                                                                                                                 q = 1 - qwe[!,:Allelfreq_Isr][1]
                                                                                                                                                 wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                                 mi = 0.0
                                                                                                                                                 mx = 0.0
                                                                                                                                                 push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[i,p[1],q[1],size(qwe,1),wi,mi,mx])
                                                                                                                         end
                                                                                                                     end
                                                                                                                  end
                                                                                                    end # of i loop
                                                                                        end # of j loop

                                                        end # of haplotype selection algorythm

                                                        # give an information how many Genes have 1,2,3,4 or more than 4 snps linked to an haplotype
                                                        println("Calculate overview over all Levels of Haplotypes")
                                                                # build an empty table where the information will be written to
                                                                haplotypecheck = DataFrame(SNPperGene=["all","1","2","3","4",">4"])
                                                      Threads.@threads  for j in 3:length(nam)
                                                                              haplotypecheck[Symbol(nam[j])] = [size(freqfile[Symbol(string(nam[j], "_Haplotypes"))],1),size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .==1,:],1),
                                                                              size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .==2,:],1),
                                                                              size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .==3,:],1),
                                                                              size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .==4,:],1),
                                                                              size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .>4,:],1)]

                                                                      end
                                                                      println("Haplotyping Finished")
                                                                      return haplotypecheck
                                 end # of function haplotyping

######## Contig mean ###########
# create summary stats and mean Frequency values for Contigs of a defined size - homolog to haplotyping, but with fixed windows
function contig(freqfile, genolist, contigsize = 10000000)
                                                          nam = CSV.read(genolist, delim='\t')[!,:Name]
                                                          n = contigsize
                                                          # calculate the maximum count of steps
                                                          steps = Int64(round(maximum(freqfile[:Location][!,:Pos]) / contigsize))

                                                          freqfile[:Location][!,:Contigs] .= "."
                                                          freqfile[:Location][!,:Cstart] .= 0
                                                          freqfile[:Location][!,:Cend] .= 0


                                                          println("prepare the Location file - create the contig bounds")
                                                          Threads.@threads for i in unique(freqfile[:Location][!,:Chr])
                                                                                                  set = freqfile[:Location][freqfile[:Location][!,:Chr] .== i,:]
                                                                                                  set[!,:start] .= 0
                                                                                                  set[!,:end] .= 0

                                                                                                  for j in 1:steps
                                                                                                                 set[.&(set[!,:Pos].> n* (j-1), set[!,:Pos].< n*j),:Contigs] .= "$(i)Contig$j"
                                                                                                                 set[.&(set[!,:Pos].> n* (j-1), set[!,:Pos].< n*j),:start] .= n * (j-1)
                                                                                                                 set[.&(set[!,:Pos].> n* (j-1), set[!,:Pos].< n*j),:end] .= n * j

                                                                                                  end
                                                                                                  freqfile[:Location][freqfile[:Location][!,:Chr] .== i,:Contigs] .= set[!,:Contigs]
                                                                                                  freqfile[:Location][freqfile[:Location][!,:Chr] .== i,:Cstart] .= set[!,:start]
                                                                                                  freqfile[:Location][freqfile[:Location][!,:Chr] .== i,:Cend] .= set[!,:end]

                                                              end

                                                          # merge by Contigs
                                                          for k in 3:length(nam); freqfile[Symbol(string(nam[k], "_Contigs"))] = DataFrame(Gene_ID=String[], Allelfreq_Isr=Float64[], Allelfreq_Golf=Float64[], SNPcount=Int64[], Readcount=Int64[], MinFreq=Float64[], MaxFreq=Float64[]); end

                                                            println("calculate the Contig Haplotypes and save it to a new Dict entry")
                                                            Threads.@threads for j in 3:length(nam)
                                                                                 for k in unique(freqfile[:Location][!,:Contigs])

                                                                                                        qwe = freqfile[Symbol(nam[j])][freqfile[:Location][!,:Contigs].==k,:]
                                                                                                        qwe = qwe[.!(isnan.(qwe[:,:11])),:]
                                                                                                        if isempty(qwe) == false
                                                                                                                                p =  sum(qwe[!,:Readcount] .* qwe[!,:Allelfreq_Isr])/sum(qwe[!,:Readcount])
                                                                                                                                q = 1 - sum(qwe[!,:Readcount] .* qwe[!,:Allelfreq_Isr])/sum(qwe[!,:Readcount]) # Golf freq
                                                                                                                                wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                mi = minimum(qwe[!,:Allelfreq_Isr]) # min Allelfreq isr of all markers for isr
                                                                                                                                mx = maximum(qwe[!,:Allelfreq_Isr]) # max Allelfreq of all markers for isr
                                                                                                                                push!(freqfile[Symbol(string(nam[j], "_Contigs"))],[k,p,q,size(qwe,1),wi, mi, mx])
                                                                                                        end
                                                                                                end
                                                            end
        end # of contig function


########## Test for multi copy variants in the parents ##############
#= The idea behind this is to extract the loci which might be cpoy variants - based on the haplotyp level
The expected coverage per SNP is known, so one can estimate using the av coverage the number of copies of the gene in the genotype
This information has to be linked to the progenies, because the Allelfrequency can be shifted sigificantly from its aspected value if multi copy variants are present
The sets worked on are the first two of the genollist file
covD = coverage of the sample test, can be used with a number or "av" for building an average value for each library
gen is the second generation to print to show in the final table
minSNPcount = the minimum number of accepted SNPs for the Haplotype
=#
function MCVdet(freqfile, genolist, haplotype="gene", covD=10, gen1=3, gen2=23, minSNPcount=3)

                                     ########## first convert the Parents into Haplotype files with information for the readdepth and the snp count
                                                  nam = CSV.read(genolist, delim='\t')[!,:Name][1:2]

                                                  # get a unique list of all GO terms without "none" and "." ; lstrip - remove space from front if exists; vcat(split) - get a vector kind list from the GO terms in the Location subdict
                                                  if haplotype == "gene"
                                                                          hap = unique(freqfile[:Location][!,:Gene_ID])
                                                                          search = Symbol("Gene_ID")
                                                  elseif haplotype == "GO"
                                                                          hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...))[occursin.("GO",unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...)))])
                                                                          search = Symbol("GO-IDs")
                                                  elseif haplotype == "Marker"
                                                                          hap = unique(freqfile[:Location][!,:Marker])
                                                                          search = Symbol("Marker")
                                                  elseif haplotype == "PF"
                                                                          hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...))[occursin.("PF",unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...)))])
                                                                          search = Symbol("PFAM-IDs")
                                                  else
                                                  println("wrong haplotype specified, chose 'gene', 'GO' or 'PF' ")

                                                  end

                                                  # create a new dictionary for each progeny
                                                  for k in 1:length(nam); freqfile[Symbol(string(nam[k], "_Haplotypes"))] = DataFrame(Gene_ID=String[], SNPcount=Int64[], Readcount=Int64[]); end
                                                                       # the haplotype file information coming with info about the gene, function and FRequency of this particular haplotype
                                                                       println("Make Haplotypes using the Bounds calculated with Genloc..")

                                                 if haplotype == "gene" || haplotype == "Marker"
                                                                           for j in 1:length(nam)
                                                                                             dat = hcat(freqfile[Symbol(nam[j])],freqfile[:Location][!,Symbol(search)])
                                                                                                                function callsnps(dat,k)
                                                                                                                                        m = Int32[]
                                                                                                                                        for i in 1:nrow(dat)
                                                                                                                                                            if dat[i,:x1] != k
                                                                                                                                                            append!(m,i)
                                                                                                                                                            break
                                                                                                                                                            end
                                                                                                                                                        end
                                                                                                                                      if  !(isempty(m))
                                                                                                                                         dat2 = dat[1:m[1]-1,:]
                                                                                                                                         deleterows!(dat, 1:m[1]-1) # remove the gene snps that have already been worked on
                                                                                                                                       else
                                                                                                                                         dat2 = dat[1:nrow(dat),:]
                                                                                                                                       end

                                                                                                                                    return dat2
                                                                                                               end

                                                                                                               dat = dat[dat[!,:x1].!= ".",:] # remove not annotated snps from the list, highlighted by a "." in the x1 gene column
                                                                                                               sort!(dat, (:Chr, :x1))
                                                                                                               kap = unique(dat[!,:x1])
                                                                                                               println(nam[j])
                                                                                                               @showprogress for k in kap
                                                                                                                                         qwe = callsnps(dat,k)
                                                                                                                                         if isempty(qwe) == false
                                                                                                                                                                    wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                                                    push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[k,size(qwe,1),wi])
                                                                                                                                            end
                                                                                                             end # of kap
                                                                 end # of nam

                                               else
                                                       Threads.@threads       for j in 1:length(nam)
                                                                                          for i in hap
                                                                                                         # the alogrhytm is following the list of genes, so all entries in the list and then alloaction the information to the single pools
                                                                                                         qwe = freqfile[Symbol(nam[j])][freqfile[:Location][search].== i,:] # get all snps in the tested gene

                                                                                                         # calculate the frequency based on which method has been choosen
                                                                                                         if isempty(qwe) == false

                                                                                                                                  wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                  push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[i,size(qwe,1),wi])
                                                                                                         end
                                                                                            end # of hap
                                                                                end # of nam

                                               end # of haplotype selection algorythm


                   ##### now where the haplotypes are calculated, test if there is a CNV
                  #covD = 10
                  if typeof(covD) != Int64
                                          covP1 = sum(freqfile[Symbol(nam[1])][:,:Readcount]) / size(freqfile[Symbol(nam[1])],1)
                                          covP1 = sum(freqfile[Symbol(nam[2])][:,:Readcount]) / size(freqfile[Symbol(nam[2])],1)
                  else
                                          covP1 = covD
                                          covP2 = covD
                  end

                  # calculate the read count per SNP in the haplotype , real and expected value
                  Threads.@threads for i in 1:size(nam,1)
                            i == 1 ? covd = covP1 : covd = covP2
                            freqfile[Symbol(string(nam[i], "_Haplotypes"))][!,:expRD] = freqfile[Symbol(string(nam[i], "_Haplotypes"))][!,:SNPcount] .* covd
                            freqfile[Symbol(string(nam[i], "_Haplotypes"))][!,:CNVratio] = freqfile[Symbol(string(nam[i], "_Haplotypes"))][!,:Readcount] ./ freqfile[Symbol(string(nam[i], "_Haplotypes"))][!,:SNPcount] .* 0.1

                            freqfile[Symbol(string(nam[i], "_Haplotypes"))][!,:CNVr] = round.(freqfile[Symbol(string(nam[i], "_Haplotypes"))][!,:CNVratio])

                            # get the infromation if here might be a copy number variant present based on the read depth
                            freqfile[Symbol(string(nam[i], "_Haplotypes"))][!,:CNV] .= "n" ::String
                            for o in 1:size(freqfile[Symbol(string(nam[i], "_Haplotypes"))],1)
                                                                                            freqfile[Symbol(string(nam[i], "_Haplotypes"))][o,:CNVratio] > 1.6 ? freqfile[Symbol(string(nam[i], "_Haplotypes"))][o,:CNV] = "yes" : freqfile[Symbol(string(nam[i], "_Haplotypes"))][o,:CNV] = "none"
                            end
                  end

                  # extract the genes with potential CNV for P1 and P2 and merge it into one file
                  qay = join(freqfile[Symbol(string(nam[1], "_Haplotypes"))], freqfile[Symbol(string(nam[2], "_Haplotypes"))], on=:Gene_ID, makeunique=true)
                  # subset the genes with CNV on one side at least
                  qay2 = qay[.|(qay[!,:CNV].=="yes", qay[!,:CNV_1].=="yes"),:]
                  append!(qay2, qay[.&(qay[!,:CNV].=="yes", qay[!,:CNV_1].=="yes"),:])
                  qay=nothing
                  # get only the genes with more than minSNPcount SNPs
                  qay2 = qay2[qay2[!,:SNPcount].>minSNPcount,:]
                  # change the column order and add the genotype name to it
                  qay2 = qay2[!,[1,2,4,3,6,9,12]]
                  rename!(qay2, vcat(names(qay2)[1:3] ,Symbol.(string.(nam[1], "_", names(qay2)[4:5])), Symbol.(string.(nam[2], "_", names(qay2)[4:5]))))

                  # a column for the F3, F23 and the effect and the gene function have to be added as well as the location information
                  # to save time, only the Haplotype data sets will be used for this
                  nam2 = CSV.read(genolist, delim='\t')

                  start = Symbol.(string.(nam2[nam2[!,:Generation].==gen1,:Name], "_Haplotypes"))
                  ek = Symbol.(string.(nam2[.&(nam2[!,:Generation].==gen2, nam2[!,:Env].==1),:Name], "_Haplotypes"))
                  eö = Symbol.(string.(nam2[.&(nam2[!,:Generation].==gen2, nam2[!,:Env].==2),:Name], "_Haplotypes"))

                  # create a column in the qay2 table for each of these datasets
                  qay2[!,Symbol(string("F",gen1))] .= 1.0 ::Float64
                  qay2[!,Symbol(string("F",gen2, "_kon"))] .= 1.0 ::Float64
                  qay2[!,Symbol(string("F",gen2, "_org"))] .= 1.0 ::Float64

                  # calculate the average Allelfreq for ISR in the given generation and env
                  print(Crayon(foreground = :yellow, bold = false),"Get the Allelfreq of ISR for the selected geenrations")
                Threads.@threads for j in qay2[!,:Gene_ID]
                                            ### starter
                                            tablestart = DataFrame(Gene_ID=String[], Allelfreq_Isr=Float64[], Readcount=Int64[])
                                            for i in start; append!(tablestart, freqfile[i][freqfile[i][!,:Gene_ID].==j,[1,2,5]]); end
                                            qay2[qay2[!,:Gene_ID].==j,Symbol(string("F",gen1))] = sum(tablestart[!,:Allelfreq_Isr] .* tablestart[!,:Readcount]) / sum(tablestart[!,:Readcount])

                                            ### ender kon
                                            tableek = DataFrame(Gene_ID=String[], Allelfreq_Isr=Float64[], Readcount=Int64[])
                                            for i in ek; append!(tableek, freqfile[i][freqfile[i][!,:Gene_ID].==j,[1,2,5]]); end
                                            qay2[qay2[!,:Gene_ID].==j, Symbol(string("F",gen2, "_kon"))] = sum(tableek[!,:Allelfreq_Isr] .* tableek[!,:Readcount]) / sum(tableek[!,:Readcount])

                                            #### ender öko
                                            tableeö = DataFrame(Gene_ID=String[], Allelfreq_Isr=Float64[], Readcount=Int64[])
                                            for i in eö; append!(tableeö, freqfile[i][freqfile[i][!,:Gene_ID].==j,[1,2,5]]); end
                                            qay2[qay2[!,:Gene_ID].==j, Symbol(string("F",gen2, "_org"))] = sum(tableeö[!,:Allelfreq_Isr] .* tableeö[!,:Readcount]) / sum(tableeö[!,:Readcount])
                  end

                  qay2[!,:FreqShift_by_CNV_kon] .= qay2[!,Symbol(string("F",gen2, "_kon"))] ./ qay2[!,Symbol(string("F",gen1))]
                  qay2[!,:FreqShift_by_CNV_org] .= qay2[!,Symbol(string("F",gen2, "_org"))] ./ qay2[!,Symbol(string("F",gen1))]

                  # filter out the genes that have been selected by multicpoy advantages -
                  # multi copy ISR > MC Golf + Freqshift > 1 => ISR selection driver
                  # multi copy Golf > MC ISR + Freqshift < 1 => Golf selection driver
                  for i in 1:size(qay2,1)
                                        ### organic system test
                                       if qay2[i,:Golf_CNVr] * 2 .< qay2[i,:Isr_CNVr] && qay2[i,:FreqShift_by_CNV_org] .> 1.5
                                                                                                                              qay2[i,:FreqShift_by_CNV_org] = 2.0
                                      elseif qay2[i,:Golf_CNVr] .> qay2[i,:Isr_CNVr] * 2 && qay2[i,:FreqShift_by_CNV_org] .< 0.7
                                                                                                                              qay2[i,:FreqShift_by_CNV_org] = 1.0
                                      else
                                                                                                                              qay2[i,:FreqShift_by_CNV_org] = 0.0
                                      end

                                      #### konventional system test
                                      if qay2[i,:Golf_CNVr] * 2 .< qay2[i,:Isr_CNVr] && qay2[i,:FreqShift_by_CNV_kon] .> 1.5
                                                                                                                             qay2[i,:FreqShift_by_CNV_kon] = 2.0
                                     elseif qay2[i,:Golf_CNVr] .> qay2[i,:Isr_CNVr] * 2 && qay2[i,:FreqShift_by_CNV_kon] .< 0.7
                                                                                                                             qay2[i,:FreqShift_by_CNV_kon] = 1.0
                                     else
                                                                                                                             qay2[i,:FreqShift_by_CNV_kon] = 0.0
                                     end

                  end
                  # convert the number to a string code
                  qay2[!,:FreqShift_by_CNV_org] = string.(qay2[!,:FreqShift_by_CNV_org])
                  qay2[!,:FreqShift_by_CNV_kon] = string.(qay2[!,:FreqShift_by_CNV_kon])

                  qay2[qay2[!,:FreqShift_by_CNV_org].=="0.0",:FreqShift_by_CNV_org] .= "none"
                  qay2[qay2[!,:FreqShift_by_CNV_kon].=="0.0",:FreqShift_by_CNV_kon] .= "none"

                  qay2[qay2[!,:FreqShift_by_CNV_org].=="1.0",:FreqShift_by_CNV_org] .= nam[1]
                  qay2[qay2[!,:FreqShift_by_CNV_kon].=="1.0",:FreqShift_by_CNV_kon] .= nam[1]

                  qay2[qay2[!,:FreqShift_by_CNV_org].=="2.0",:FreqShift_by_CNV_org] .= nam[2]
                  qay2[qay2[!,:FreqShift_by_CNV_kon].=="2.0",:FreqShift_by_CNV_kon] .= nam[2]


                  # add the position information of the gene  ### use INDEXIN (first in second)
                  qay2[!,:Position] .= freqfile[:Location][indexin(qay2[!,:Gene_ID], freqfile[:Location][!,:Gene_ID]),Symbol("chromosome:start-stop")]
                  # add gene function
                  qay2[!,:Gene_Function] .= freqfile[:Location][indexin(qay2[!,:Gene_ID], freqfile[:Location][!,:Gene_ID]),:description ]

                  return qay2

 end # of function

######## Test for equal variance ###########
#= The datasets should be tested for equal variance using a Negative binomial and Zero Infalted Distribution Model
- because the distribution based on single SNP data can be strongly biased, the test is performed on a level of haplotypes, generated by the function haplotyping
- a dict entry as a table will be created giving the pvalues for the replicates from the neg Bin. Test as well as from a Zero Infl. GLM calculation
- besides the test within the replicates, the generations can be tested against each other as well
- testing on 3 levels:
  # first level = within the relplicates based on a NegBin + ZeroInfl
  # second level = in between the generations
  # thrid level = in between the environmental systems
  second and third level are based on the glm poisson model, because in comparisen to KruskalWallisTest it is less likly to be significant = more strict
  # typer = _Haplotyes or _Contigs
=#
function equaltest(freqfile, genolist, typer="_Haplotypes")
  # for this tool, we need to have the infromation about the years and the evironment conditions of each sample
                          nam = CSV.read(genolist, delim='\t') # Env & Generation = 0 indicate the parents
                          # create a 4th column where env and generation are melted
                          nam[!,:melt] = string.(nam[!,:Env],"_", nam[!,:Generation])
                    ### first level - the the replicates
                          # create an array where to save teh values
                          print(Crayon(foreground = :blue, bold = false),"First level - comparing the replicates of each generation and system for reproducability \n")


                          # the kruskal KruskalWallisTest is much more strict than the glm model based on a poisson distribution
                          first = DataFrame(set = unique(nam[!,:melt]), negbinX2 = 0.0, zeroinfX2 = 0.0, negbinX3 = 0.0, zeroinfX3 = 0.0, negbinX4 = 0.0, zeroinfX4 = 0.0, negbinX5 = 0.0, zeroinfX5 = 0.0)
          for i in unique(nam[!,:melt])
                                                      if i != "0_0" # do not perform this test for the parents
                                                                    m = sum(occursin.(i,nam[!,:melt])) # get the information how many replicates one has
                                                                    nl = nam[nam[!,:melt].==i,:Name] # the names of the libraries to test
                                                                    nl = string.(nl,typer)
                                                                    if m == 1
                                                                                              continue
                                                                    elseif m == 2
                                                                                              t1 = vcat(freqfile[Symbol(nl[1])][!,:Allelfreq_Isr], freqfile[Symbol(nl[2])][!,:Allelfreq_Isr])
                                                                                              t3 = repeat([1], size(freqfile[Symbol(nl[1])],1))
                                                                                              t4 = repeat([2], size(freqfile[Symbol(nl[2])],1))
                                                                                              n3 = DataFrame(hcat(t1,vcat(t3,t4)))
                                                                    elseif m == 3
                                                                                              t1 = vcat(freqfile[Symbol(nl[1])][!,:Allelfreq_Isr], freqfile[Symbol(nl[2])][!,:Allelfreq_Isr], freqfile[Symbol(nl[3])][!,:Allelfreq_Isr])
                                                                                              t2 = repeat([1], size(freqfile[Symbol(nl[1])],1))
                                                                                              t3 = repeat([2], size(freqfile[Symbol(nl[2])],1))
                                                                                              t4 = repeat([3], size(freqfile[Symbol(nl[3])],1))
                                                                                              n3 = DataFrame(hcat(t1,vcat(t2,t3,t4)))
                                                                    elseif m == 4
                                                                                              t1 = vcat(freqfile[Symbol(nl[1])][!,:Allelfreq_Isr], freqfile[Symbol(nl[2])][!,:Allelfreq_Isr],  freqfile[Symbol(nl[3])][!,:Allelfreq_Isr],  freqfile[Symbol(nl[4])][!,:Allelfreq_Isr])
                                                                                              t2 = repeat([1], size(freqfile[Symbol(nl[1])],1))
                                                                                              t3 = repeat([2], size(freqfile[Symbol(nl[2])],1))
                                                                                              t4 = repeat([3], size(freqfile[Symbol(nl[3])],1))
                                                                                              t5 = repeat([4], size(freqfile[Symbol(nl[4])],1))
                                                                                              n3 = DataFrame(hcat(t1,vcat(t2,t3,t4,t4,t5)))
                                                                    elseif m == 5
                                                                                              t1 = vcat(freqfile[Symbol(nl[1])][!,:Allelfreq_Isr], freqfile[Symbol(nl[2])][!,:Allelfreq_Isr],  freqfile[Symbol(nl[3])][!,:Allelfreq_Isr],  freqfile[Symbol(nl[4])][!,:Allelfreq_Isr],  freqfile[Symbol(nl[5])][!,:Allelfreq_Isr])
                                                                                              t2 = repeat([1], size(freqfile[Symbol(nl[1])],1))
                                                                                              t3 = repeat([2], size(freqfile[Symbol(nl[2])],1))
                                                                                              t4 = repeat([3], size(freqfile[Symbol(nl[3])],1))
                                                                                              t5 = repeat([4], size(freqfile[Symbol(nl[4])],1))
                                                                                              t6 = repeat([5], size(freqfile[Symbol(nl[5])],1))
                                                                                              n3 = DataFrame(hcat(t1,vcat(t2,t3,t4,t5,t6)))
                                                                    else println("Name of the first test generation does not work, too many replicates?")
                                                                    end

                                                     n3[!,1] = Int64.(floor.(n3[!,1] .* 100))
                                                     @rput n3
                                                     R"""
                                                     library(MASS)
                                                     library(pscl)
                                                     n3[,2] = as.factor(n3[,2])
                                                     m1 =  summary(zeroinfl(x1~x2, n3, dist = "negbin"))
                                                     negbino = m1$coefficients$count[,4]
                                                     zeroinf = m1$coefficients$zero[,4]
                                                     """
                                                     @rget negbino  zeroinf

                                                     if m == 1
                                                     elseif m == 2
                                                                    first[first[!,:set].==i,:negbinX2] = negbino[2]
                                                                    first[first[!,:set].==i,:zeroinfX2] = zeroinf[2]
                                                     elseif m == 3
                                                                    first[first[!,:set].==i,:negbinX2] = negbino[2]
                                                                    first[first[!,:set].==i,:negbinX3] = negbino[3]
                                                                    first[first[!,:set].==i,:zeroinfX2] = zeroinf[2]
                                                                    first[first[!,:set].==i,:zeroinfX3] = zeroinf[3]
                                                     elseif m == 4
                                                                     first[first[!,:set].==i,:negbinX2] = negbino[2]
                                                                     first[first[!,:set].==i,:negbinX3] = negbino[3]
                                                                     first[first[!,:set].==i,:negbinX4] = negbino[4]
                                                                     first[first[!,:set].==i,:zeroinfX2] = zeroinf[2]
                                                                     first[first[!,:set].==i,:zeroinfX3] = zeroinf[3]
                                                                     first[first[!,:set].==i,:zeroinfX4] = zeroinf[4]
                                                     elseif m == 5
                                                                     first[first[!,:set].==i,:negbinX2] = negbino[2]
                                                                     first[first[!,:set].==i,:negbinX3] = negbino[3]
                                                                     first[first[!,:set].==i,:negbinX4] = negbino[4]
                                                                     first[first[!,:set].==i,:negbinX5] = negbino[5]
                                                                     first[first[!,:set].==i,:zeroinfX2] = zeroinf[2]
                                                                     first[first[!,:set].==i,:zeroinfX3] = zeroinf[3]
                                                                     first[first[!,:set].==i,:zeroinfX4] = zeroinf[4]
                                                                     first[first[!,:set].==i,:zeroinfX5] = zeroinf[5]
                                                     else
                                                         println("Something id wrong in the first step, please have a look in the code")
                                                     end

                                                 end

                              end
            #### second level - between generations
            print(Crayon(foreground = :green, bold = false),"Second level - comparing the generations against each other for each generation against each other within the cultivation system \n")


            eco = nam[nam[!,:Env].==2,:]
            kon = nam[nam[!,:Env].==1,:]
            # create an array where to save teh values - it is a x by x table, writing the values to the table
            es = DataFrame(Array{Float64,2}(undef, size(unique(eco[!,:Generation]),1)*2, size(unique(eco[!,:Generation]),1)))
            rename!(es, Symbol.(string.("F",unique(eco[!,:Generation]))))
            for i in 1:size(es,2); es[!,i] .= 1.0 ; end
            es = hcat(repeat(unique(eco[!,:Generation]),inner = 2),hcat(repeat(["negbin", "zeroinfl"],  length(unique(eco[!,:Generation]))) , es ),makeunique=true)

            ks = DataFrame(Array{Float64,2}(undef, size(unique(kon[!,:Generation]),1)*2, size(unique(kon[!,:Generation]),1)))
            rename!(ks, Symbol.(string.("F",unique(kon[!,:Generation]))))
            for i in 1:size(ks,2); ks[!,i] .= 1.0 ; end
            ks = hcat(repeat(unique(kon[!,:Generation]),inner = 2),hcat(repeat(["negbin", "zeroinfl"],  length(unique(kon[!,:Generation]))) , ks ),makeunique=true)

            # get the names of the libraries of interest
            en = String[]
            eg = Int64[]
            for i in unique(eco[!,:Generation])
                                              append!(en,nam[.&(nam[!,:Env].==2, nam[!,:Generation].== i ),:Name])
                                              append!(eg,nam[.&(nam[!,:Env].==2, nam[!,:Generation].== i ),:Generation])
                                             end
            e = hcat(en,eg)

            kn = String[]
            kg = Int64[]
            for i in unique(kon[!,:Generation])
                                            append!(kn,nam[.&(nam[!,:Env].==1, nam[!,:Generation].== i ),:Name])
                                            append!(kg,nam[.&(nam[!,:Env].==1, nam[!,:Generation].== i ),:Generation])
                                           end
            k = hcat(kn,kg)

            # get all possible combinations for the generations - using R function combn
            @rput e k
            R"""
            if (length(unique(e[,2])) > 1){
                              ekombi = combn(unique(e[,2]),2)
                              }else{ekombi = 0}

            if (length(unique(k[,2])) > 1){
                               kkombi = combn(unique(k[,2]),2)
                               }else{kkombi = 0}
            """
            @rget ekombi kkombi

            combinations = Dict()
            combinations[:ecol] = ekombi
            combinations[:konv] = kkombi

          for j in keys(combinations)
                      # if one is only checking data for one single Environemt, this will lock out the one with less than 2 generations
                      if size(combinations[j]) == ()
                                        print(Crayon(foreground = :red, bold = true),"There is nothing to test for $j \n")
                                        continue
                      end

                     j == :ecol ? p = e : p = k  # check where to select the names from


           for i in 1:size(combinations[j],2)

                      name1 = string.(p[p[:,2].==combinations[j][1,i],1], typer)
                      name2 = string.(p[p[:,2].==combinations[j][2,i],1], typer)
                      # name1
                      if length(name1) == 1
                                                t1 = freqfile[Symbol(name1[1])][!,:Allelfreq_Isr]
                                                t2 = repeat([1], size(freqfile[Symbol(name1[1])],1))
                                                t3 = DataFrame(hcat(t1,t2))
                      elseif length(name1) == 2
                                                t1 = vcat(freqfile[Symbol(name1[1])][!,:Allelfreq_Isr], freqfile[Symbol(name1[2])][!,:Allelfreq_Isr])
                                                t2 = repeat([1], size(freqfile[Symbol(name1[1])],1) + size(freqfile[Symbol(name1[2])],1))
                                                t3 = DataFrame(hcat(t1,t2))
                      elseif length(name1) == 3
                                                t1 = vcat(freqfile[Symbol(name1[1])][!,:Allelfreq_Isr], freqfile[Symbol(name1[2])][!,:Allelfreq_Isr], freqfile[Symbol(name1[3])][!,:Allelfreq_Isr])
                                                t2 = repeat([1], size(freqfile[Symbol(name1[1])],1) + size(freqfile[Symbol(name1[2])],1) + size(freqfile[Symbol(name1[3])],1))
                                                t3 = DataFrame(hcat(t1,t2))
                      elseif length(name1) == 4
                                                t1 = vcat(freqfile[Symbol(name1[1])][!,:Allelfreq_Isr], freqfile[Symbol(name1[2])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[3])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[4])][!,:Allelfreq_Isr])
                                                t2 = repeat([1], size(freqfile[Symbol(name1[1])],1) + size(freqfile[Symbol(name1[2])],1) + size(freqfile[Symbol(name1[3])],1) + size(freqfile[Symbol(name1[4])],1))
                                                t3 = DataFrame(hcat(t1,t2))
                      elseif length(name1) == 5
                                                t1 = vcat(freqfile[Symbol(name1[1])][!,:Allelfreq_Isr], freqfile[Symbol(name1[2])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[3])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[4])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[5])][!,:Allelfreq_Isr])
                                                t2 = repeat([1], size(freqfile[Symbol(name1[1])],1) +  size(freqfile[Symbol(name1[2])],1) + size(freqfile[Symbol(name1[3])],1) + size(freqfile[Symbol(name1[4])],1) + size(freqfile[Symbol(name1[5])],1))
                                                t3 = DataFrame(hcat(t1,t2))
                      else println("Name of the first test generation does not work, too many replicates?")
                      end

                     # name2
                     if length(name2) == 1
                                               m1 = freqfile[Symbol(name2[1])][!,:Allelfreq_Isr]
                                               m2 = repeat([2], size(freqfile[Symbol(name2[1])],1))
                                               m3 = DataFrame(hcat(m1,m2))
                     elseif length(name2) == 2
                                               m1 = vcat(freqfile[Symbol(name2[1])][!,:Allelfreq_Isr], freqfile[Symbol(name2[2])][!,:Allelfreq_Isr])
                                               m2 = repeat([2], size(freqfile[Symbol(name2[1])],1) + size(freqfile[Symbol(name2[2])],1))
                                               m3 = DataFrame(hcat(m1,m2))
                     elseif length(name2) == 3
                                               m1 = vcat(freqfile[Symbol(name2[1])][!,:Allelfreq_Isr], freqfile[Symbol(name2[2])][!,:Allelfreq_Isr], freqfile[Symbol(name2[3])][!,:Allelfreq_Isr])
                                               m2 = repeat([2], size(freqfile[Symbol(name2[1])],1) + size(freqfile[Symbol(name2[2])],1) + size(freqfile[Symbol(name2[3])],1))
                                               m3 = DataFrame(hcat(m1,m2))
                     elseif length(name2) == 4
                                               m1 = vcat(freqfile[Symbol(name2[1])][!,:Allelfreq_Isr], freqfile[Symbol(name2[2])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[3])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[4])][!,:Allelfreq_Isr])
                                               m2 = repeat([2], size(freqfile[Symbol(name2[1])],1) + size(freqfile[Symbol(name2[2])],1) + size(freqfile[Symbol(name2[3])],1) + size(freqfile[Symbol(name2[4])],1))
                                               m3 = DataFrame(hcat(m1,m2))
                     elseif length(name2) == 5
                                               m1 = vcat(freqfile[Symbol(name2[1])][!,:Allelfreq_Isr], freqfile[Symbol(name2[2])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[3])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[4])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[5])][!,:Allelfreq_Isr])
                                               m2 = repeat([2], size(freqfile[Symbol(name2[1])],1) + size(freqfile[Symbol(name2[2])],1) + size(freqfile[Symbol(name2[3])],1) + size(freqfile[Symbol(name2[4])],1) + size(freqfile[Symbol(name2[5])],1))
                                               m3 = DataFrame(hcat(m1,m2))
                     else println("Name of the second test generation does not work, too many replicates?")
                     end

                     # link first abd second generation
                     n3 = vcat(t3,m3)
                     #n3[!,2] = categorical(n3[!,2])
                     n3[!,1] = Int64.(floor.(n3[!,1] .* 100))


                  # choose a model that should be used for calculating the p value
                  print(Crayon(foreground = :yellow, bold = false),"Calculation of the negativ binomial and zero inflated GLM for  $name1 and $name2 \n")

                  @rput n3
                  R"""
                  n3$x2 = as.factor(n3$x2)
                  m1 =  summary(zeroinfl(x1~x2, n3, dist = "negbin"))
                  negbino = m1$coefficients$count[2,4]
                  zeroinf = m1$coefficients$zero[2,4]
                  """
                  @rget negbino  zeroinf
                  # the pvalue has to stored at the correct postion in the cross table
                  #es table is the one to store it in

                  # s1 is the column name, s2 is the generation listed in the first column
                  s1 = Symbol(string("F",combinations[j][:,i][1]))
                  s2 = combinations[j][:,i][2]
                  # set the p value to the correct position in the table
                  if j == :ecol
                       es[.&(es[:,1].==s2, es[:,2].== "negbin"),s1] = negbino
                       es[.&(es[:,1].==s2, es[:,2].== "zeroinfl"),s1] = zeroinf
                   else
                        ks[.&(ks[:,1].==s2, ks[:,2].== "negbin"),s1] = negbino
                         ks[.&(ks[:,1].==s2, ks[:,2].== "zeroinfl"),s1] = zeroinf
                    end

                end # inner loop
      end # outer loop

                # thrid level = in between the environmental systems
                print(Crayon(foreground = :blue, bold = false),"Thrid level - comapring the systems against each other for each generation against each other \n")

                # the table where to store the informations
                bs = DataFrame(Array{Float64,2}(undef, size(unique(kon[!,:Generation]),1)*2, size(unique(eco[!,:Generation]),1)))
                rename!(bs, Symbol.(string.("F",unique(eco[!,:Generation]))))
                for i in 1:size(bs,2); bs[!,i] .= 1.0 ; end
                bs = hcat(string.("F",repeat(unique(kon[!,:Generation]),inner = 2)),hcat(repeat(["negbin", "zeroinfl"],  length(unique(kon[!,:Generation]))) , bs ),makeunique=true)

                # konventional generations are stored in the first column rowwise, while ecological generations are stored in the 2nd to last column

                # get all possible combinations for the generations - using R function combn
                # e and k are selected in the second step already, but have to go to one single array
                for i in unique(k[:,2]), j in unique(e[:,2])

                                name1 = string.(e[e[:,2].== j,1], typer)
                                name2 = string.(k[k[:,2].== i,1], typer)

                                # name1
                                if length(name1) == 1
                                                          t1 = freqfile[Symbol(name1[1])][!,:Allelfreq_Isr]
                                                          t2 = repeat([1], size(freqfile[Symbol(name1[1])],1))
                                                          t3 = DataFrame(hcat(t1,t2))
                                elseif length(name1) == 2
                                                          t1 = vcat(freqfile[Symbol(name1[1])][!,:Allelfreq_Isr], freqfile[Symbol(name1[2])][!,:Allelfreq_Isr])
                                                          t2 = repeat([1], size(freqfile[Symbol(name1[1])],1) + size(freqfile[Symbol(name1[2])],1))
                                                          t3 = DataFrame(hcat(t1,t2))
                                elseif length(name1) == 3
                                                          t1 = vcat(freqfile[Symbol(name1[1])][!,:Allelfreq_Isr], freqfile[Symbol(name1[2])][!,:Allelfreq_Isr], freqfile[Symbol(name1[3])][!,:Allelfreq_Isr])
                                                          t2 = repeat([1], size(freqfile[Symbol(name1[1])],1) + size(freqfile[Symbol(name1[2])],1) + size(freqfile[Symbol(name1[3])],1))
                                                          t3 = DataFrame(hcat(t1,t2))
                                elseif length(name1) == 4
                                                          t1 = vcat(freqfile[Symbol(name1[1])][!,:Allelfreq_Isr], freqfile[Symbol(name1[2])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[3])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[4])][!,:Allelfreq_Isr])
                                                          t2 = repeat([1], size(freqfile[Symbol(name1[1])],1) + size(freqfile[Symbol(name1[2])],1) + size(freqfile[Symbol(name1[3])],1) + size(freqfile[Symbol(name1[4])],1))
                                                          t3 = DataFrame(hcat(t1,t2))
                                elseif length(name1) == 5
                                                          t1 = vcat(freqfile[Symbol(name1[1])][!,:Allelfreq_Isr], freqfile[Symbol(name1[2])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[3])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[4])][!,:Allelfreq_Isr],  freqfile[Symbol(name1[5])][!,:Allelfreq_Isr])
                                                          t2 = repeat([1], size(freqfile[Symbol(name1[1])],1) +  size(freqfile[Symbol(name1[2])],1) + size(freqfile[Symbol(name1[3])],1) + size(freqfile[Symbol(name1[4])],1) + size(freqfile[Symbol(name1[5])],1))
                                                          t3 = DataFrame(hcat(t1,t2))
                                else println("Name of the first test generation does not work, too many replicates?")
                                end

                                # name2
                                if length(name2) == 1
                                                         m1 = freqfile[Symbol(name2[1])][!,:Allelfreq_Isr]
                                                         m2 = repeat([2], size(freqfile[Symbol(name2[1])],1))
                                                         m3 = DataFrame(hcat(m1,m2))
                                elseif length(name2) == 2
                                                         m1 = vcat(freqfile[Symbol(name2[1])][!,:Allelfreq_Isr], freqfile[Symbol(name2[2])][!,:Allelfreq_Isr])
                                                         m2 = repeat([2], size(freqfile[Symbol(name2[1])],1) + size(freqfile[Symbol(name2[2])],1))
                                                         m3 = DataFrame(hcat(m1,m2))
                                elseif length(name2) == 3
                                                         m1 = vcat(freqfile[Symbol(name2[1])][!,:Allelfreq_Isr], freqfile[Symbol(name2[2])][!,:Allelfreq_Isr], freqfile[Symbol(name2[3])][!,:Allelfreq_Isr])
                                                         m2 = repeat([2], size(freqfile[Symbol(name2[1])],1) + size(freqfile[Symbol(name2[2])],1) + size(freqfile[Symbol(name2[3])],1))
                                                         m3 = DataFrame(hcat(m1,m2))
                                elseif length(name2) == 4
                                                         m1 = vcat(freqfile[Symbol(name2[1])][!,:Allelfreq_Isr], freqfile[Symbol(name2[2])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[3])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[4])][!,:Allelfreq_Isr])
                                                         m2 = repeat([2], size(freqfile[Symbol(name2[1])],1) + size(freqfile[Symbol(name2[2])],1) + size(freqfile[Symbol(name2[3])],1) + size(freqfile[Symbol(name2[4])],1))
                                                         m3 = DataFrame(hcat(m1,m2))
                                elseif length(name2) == 5
                                                         m1 = vcat(freqfile[Symbol(name2[1])][!,:Allelfreq_Isr], freqfile[Symbol(name2[2])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[3])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[4])][!,:Allelfreq_Isr],  freqfile[Symbol(name2[5])][!,:Allelfreq_Isr])
                                                         m2 = repeat([2], size(freqfile[Symbol(name2[1])],1) + size(freqfile[Symbol(name2[2])],1) + size(freqfile[Symbol(name2[3])],1) + size(freqfile[Symbol(name2[4])],1) + size(freqfile[Symbol(name2[5])],1))
                                                         m3 = DataFrame(hcat(m1,m2))
                                else println("Name of the second test generation does not work, too many replicates?")
                                end

                                # link first abd second generation
                                n3 = vcat(t3,m3)
                                n3[!,1] = Int64.(floor.(n3[!,1] .* 100))


                                # choose a model that should be used for calculating the p value
                                print(Crayon(foreground = :yellow, bold = false),"Calculation of the generalized linear model based on the negativ binomial and zero inflated distribution for  $name1 and $name2 \n")

                                @rput n3
                                R"""
                                n3$x2 = as.factor(n3$x2)
                                m1 =  summary(zeroinfl(x1~x2, n3, dist = "negbin"))
                                negbino = m1$coefficients$count[2,4]
                                zeroinf = m1$coefficients$zero[2,4]
                                """
                                @rget negbino  zeroinf

                                # the pvalue has to stored at the correct postion in the cross table
                                #es table is the one to store it in
                                # s1 is the column name, s2 is the generation listed in the first column
                                s1 = Symbol(string("F",j))
                                s2 = string("F",i)
                                # set the p value to the correct position in the table
                                bs[.&(bs[:,1].==s2, bs[:,2].== "negbin"),s1] = negbino
                                bs[.&(bs[:,1].==s2, bs[:,2].== "zeroinfl"),s1] = zeroinf

                end

                # create a dict to save all the compariesens of first, second and third level
                res = Dict()
                res[:replicates] = first
                res[:egeneration] = es
                res[:kgeneration] = ks
                res[:system] = bs

                return res
 end # of function

####################################
#= The data points have to be compared for years and the systems
it follows the same pattern as the gen expression analysis - check if a snps allelfreq in a given generation is different to the previous year
and compared to the other cultivating system
 - haplotype - one which level should the compariesen be done - gene, Marker, GO or PF
 - teststat - descide if doing a GLM Poisson model for the generation comparisen or a simple ttest ANOVA
 - aggregate - how should the replicates be treated - get one mean value for each Marker or get one value per replicate added to the list ("mean" / "rep")
=#
function evol(freqfile, genolist, haplotype="gene", transform=false)

                          # A BASE FUNCTION THAT IS NEEDED FOR HAPSTATS AND PVALHAP CALULCATION
                          function callsnps(dat,k)
                                                  m = Int32[]
                                                  for i in 1:nrow(dat)
                                                                      if dat[i,:x1] != k
                                                                      append!(m,i)
                                                                      break
                                                                      end
                                                                  end
                                                if  !(isempty(m))
                                                   dat2 = dat[1:m[1]-1,:]
                                                   deleterows!(dat, 1:m[1]-1) # remove the gene snps that have already been worked on
                                                 else
                                                   dat2 = dat[1:nrow(dat),:]
                                                 end

                                              return dat2
                         end

                        nam = CSV.read(genolist, delim='\t') # Env & Generation = 0 indicate the parents

                        # get a unique list of all GO terms without "none" and "." ; lstrip - remove space from front if exists; vcat(split) - get a vector kind list from the GO terms in the Location subdict
                        if haplotype == "gene"
                                                hap = unique(freqfile[:Location][!,:Gene_ID])
                                                hap = hap[hap.!="."]
                        elseif haplotype == "GO"
                                                hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...))[occursin.("GO",unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...)))])
                                                hap = hap[hap.!="."]
                        elseif haplotype == "Marker"
                                                hap = unique(freqfile[:Location][!,:Marker])
                                                hap = hap[hap.!="."]
                        elseif haplotype == "PF"
                                                hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...))[occursin.("PF",unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...)))])
                                                hap = hap[hap.!="."]

                        else
                        println("wrong haplotype specified, chose 'gene', 'GO' or 'PF' ")

                        end

                        # split the general haplotype description part from the part where the generations are compared to each other
                        #set the column to search in
                        if haplotype == "gene"
                                                search = "Gene_ID"
                        elseif haplotype == "GO"
                                                search = "GO-IDs"
                        elseif haplotype == "Marker"
                                                search = "Marker"
                        elseif haplotype == "PF"
                                                search = "PFAM-IDs"
                        else
                        println("Do you really know what you are doing?!")

                        end

                        # get the info which is the first generation, split by Environment and list all other sample generations
                        pstarte = string.("f" ,nam[.&((nam[!,:Generation] .== minimum(nam[!,:Generation][nam[!,:Generation].>0])), (nam[!,:Env] .== 2)),:Generation],"e2") # get the name of the earliest Population tested for ecological systems
                        pstartk = string.("f" ,nam[.&((nam[!,:Generation] .== minimum(nam[!,:Generation][nam[!,:Generation].>0])), (nam[!,:Env] .== 1)),:Generation],"e1")
                        genE = string.("f" ,unique(nam[.&((nam[!,:Generation].> minimum(nam[!,:Generation][nam[!,:Generation].>0])),(nam[!,:Env] .== 2)),:Generation]),"e2")
                        genK = string.("f" ,unique(nam[.&((nam[!,:Generation].> minimum(nam[!,:Generation][nam[!,:Generation].>0])),(nam[!,:Env] .== 1)),:Generation]),"e1")

                        # creat a dictionary for ecological and conventinal farming system # when there is no starting generation to compare against the dataframes and Dicts have to be arranged differently

                        if isempty(genE) == false && isempty(pstarte) == false
                                        eco = Dict(Symbol("HapStats") =>  DataFrame(Array{Union{Nothing, Any}}(nothing, length(hap),(length(genE)+1)*5)), Symbol.("pvalHap") => DataFrame(Array{Union{Nothing, Any}}(nothing, length(hap),length(genE)*3)))
                                        # the concatination of strings in 2 vectors is not possible with the standard string.() function, therefore a loop with append should help out to get the colnames
                                        colnames = string.(unique(pstarte) ,"_", ["Average", "Var", "Readdepth", "Markercount", "Zerocount"]) # get the starter
                                        for i in 1:length(genE); append!(colnames, string.(genE[i] ,"_", ["Average", "Var", "Readdepth", "Markercount", "Zerocount"])); end
                                        rename!(eco[:HapStats], Symbol.(colnames))
                                        rename!(eco[:pvalHap], Symbol.(string.(repeat(genE,inner=3),repeat(["_negbinManU", "_zeroInfl", "_Stat"],size(genE,1)))))

                        elseif isempty(genE) == false && isempty(pstarte) == true
                                        eco = Dict(Symbol("HapStats") =>  DataFrame(Array{Union{Nothing, Any}}(nothing, length(hap),(length(genE))*5)))
                                        colnames =  Array{String}(undef, 1)
                                        for i in 1:length(genE); append!(colnames, string.(genE[i] ,"_", ["Average", "Var", "Readdepth", "Markercount", "Zerocount"])); end
                                        deleteat!(colnames, 1) # the first entry is #undef und nit usefull
                                        rename!(eco[:HapStats], Symbol.(colnames))
                        else
                                        println("No ecological systems specified")
                        end

                        if isempty(genK) == false && isempty(pstartk) == false
                                            kon = Dict(Symbol("HapStats") =>  DataFrame(Array{Union{Nothing, Any}}(nothing, length(hap),(length(genK)+1)*5)), Symbol.("pvalHap") => DataFrame(Array{Union{Nothing, Any}}(nothing, length(hap),length(genK)*3)))
                                            # the concatination of strings in 2 vectors is not possible with the standard string.() function, therefore a loop with append should help out to get the colnames
                                            colnames = string.(unique(pstartk) ,"_", ["Average", "Var", "Readdepth", "Markercount", "Zerocount"]) # get the starter
                                            for i in 1:length(genK); append!(colnames, string.(genK[i] ,"_", ["Average", "Var", "Readdepth", "Markercount", "Zerocount"])); end
                                            rename!(kon[:HapStats], Symbol.(colnames))
                                            rename!(kon[:pvalHap], Symbol.(string.(repeat(genK,inner=3),repeat(["_negbinManU", "_zeroInfl", "_Stat"],size(genE,1)))))

                        elseif isempty(genK) == false && isempty(pstartk) == true
                                            kon = Dict(Symbol("HapStats") =>  DataFrame(Array{Union{Nothing, Any}}(nothing, length(hap),(length(genK))*5)))
                                            colnames =  Array{String}(undef, 1)
                                            for i in 1:length(genK); append!(colnames, string.(genK[i] ,"_", ["Average", "Var", "Readdepth", "Markercount", "Zerocount"])); end
                                            deleteat!(colnames, 1) # the first entry is #undef und nit usefull

                                            rename!(kon[:HapStats], Symbol.(string.(genK ,"_", ["Average", "Var", "Readdepth", "Markercount", "Zerocount"])))
                        else
                                            println("No konventional Systems specified")
                                            kon = Dict(Symbol("HapStats") =>  DataFrame(Array{Union{Nothing, Any}}(nothing, 0,0)))
                        end
                        # put eco and kon in a Dict to access it more easy in loops
                        comp = Dict()

                        comp[:eco] = eco
                        comp[:kon] = kon
                        # match the haplotypes with the dataframes

                        for i in ["eco","kon"]
                                            comp[Symbol(i)][:HapStats] = hcat(hap,comp[Symbol(i)][:HapStats])
                                            if haskey(comp[Symbol(i)],:pvalHap) == true
                                            comp[Symbol(i)][:pvalHap] = hcat(hap,comp[Symbol(i)][:pvalHap])
                                            end
                                            end

                        #= push the replicates left after the equaltest into one file to calculate the haplotype Statistics - =#

                        # freqfile = freq #  !!!!!!!!!!!!   has to be specified while testing!!!!!!!!!!!!!!!!

        Threads.@threads for i in unique(nam[nam[!,:Generation].>0,:Generation])
                                    for j in 1:2

                                    replications = nam[.&((nam[!,:Env] .== j),(nam[!,:Generation].== i)),:Name]
                                    print(Crayon(foreground = :light_gray, bold = false),"Processing F$i of Environment $j \n")

                                    # create an extra Dict entry for the merged Replicates
                                    try
                                    freqfile[Symbol(string("f",i,"e",j))] = DataFrame(Array{Union{Nothing, Any}}(nothing, 0,ncol(freqfile[Symbol(replications[1])])))
                                    rename!(freqfile[Symbol(string("f",i,"e",j))], names(freqfile[Symbol(replications[1])]))
                                    catch
                                    println("nothing to do for this..")
                                    continue
                                    end

                                    # put all in one file
                                    if size(replications,1) == 1
                                                                      freqfile[Symbol(string("f",i,"e",j))] = freqfile[Symbol(replications[1])]
                                    elseif size(replications,1) == 2
                                                                      freqfile[Symbol(string("f",i,"e",j))] = sort(vcat(freqfile[Symbol(replications[1])], freqfile[Symbol(replications[2])]), (:Chr, :Pos))
                                    elseif size(replications,1) == 3
                                                                      freqfile[Symbol(string("f",i,"e",j))] = sort(vcat(freqfile[Symbol(replications[1])], freqfile[Symbol(replications[2])], freqfile[Symbol(replications[3])]), (:Chr, :Pos))
                                    elseif size(replications,1) == 4
                                                                      freqfile[Symbol(string("f",i,"e",j))] = sort(vcat(freqfile[Symbol(replications[1])], freqfile[Symbol(replications[2])], freqfile[Symbol(replications[3])], freqfile[Symbol(replications[4])]), (:Chr, :Pos))
                                    elseif size(replications,1) == 5
                                                                      freqfile[Symbol(string("f",i,"e",j))] = sort(vcat(freqfile[Symbol(replications[1])], freqfile[Symbol(replications[2])], freqfile[Symbol(replications[3])], freqfile[Symbol(replications[4])], freqfile[Symbol(replications[5])]), (:Chr, :Pos))
                                    else
                                        println("Replicates could not be matched in one file")
                                        break
                                    end

                                    end
                                  end


            Threads.@threads  for j in vcat(pstarte, pstartk, genE, genK)## loop for the haplotypes and the samples

                                                            if occursin("e1",j)
                                                            i = "kon"
                                                            elseif occursin("e2",j)
                                                            i = "eco"
                                                            end

                                                            print(Crayon(foreground = :white, bold = false),"print general stats for $j to table HapStats..\n")
                                                            # create an intermediate file to work with that contains the information of the freq file as well as the corresponding Gene information
                                                            ## the transformation to an array is performed to increase the speed
                                                            # one replicate
                                                            if size(freqfile[Symbol(j)],1) == size(freqfile[:Location],1)
                                                                                dat = hcat(freqfile[Symbol(j)],freqfile[:Location][!,Symbol(search)])
                                                            # two replicates
                                                            elseif size(freqfile[Symbol(j)],1) == size(freqfile[:Location],1) * 2
                                                                                dat = hcat(freqfile[Symbol(j)], sort(vcat(freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                            # three replicates
                                                            elseif size(freqfile[Symbol(j)],1) == size(freqfile[:Location],1) * 3
                                                                                dat = hcat(freqfile[Symbol(j)], sort(vcat(freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                            # four Replicates
                                                            elseif size(freqfile[Symbol(j)],1) == size(freqfile[:Location],1) * 4
                                                                                dat = hcat(freqfile[Symbol(j)], sort(vcat(freqfile[:Location],freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                            else
                                                            println("there is an issue with the count of replicates used in settings aggi=rep, use aggi=mean or start digging")
                                                            #continue
                                                            end

                                                            if any(occursin.(string.(j),string.(names(comp[Symbol(i)][:HapStats])))) == true # checks if there is any "fxY" in the HapStats table - if so the summary entries will be writen in this table

                                                              ####################################################
                                                            dat = dat[dat[!,:x1].!= ".",:] # remove not annotated snps from the list, highlighted by a "." in the x1 gene column
                                                            sort!(dat, (:Chr, :x1))
                                                            kap = unique(dat[!,:x1])

                                                            @showprogress for k in kap

                                                                                               dat2 = callsnps(dat,k)
                                                                                               dat2 = dat2[.!(isnan.(dat2[:,:11])),:] # remove the NaN values from the table
                                                                                                if !(isempty(dat2))
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Average"))] = sum(dat2[:,11] .* dat2[:,9]) / sum(dat2[:,9])
                                                                                                                    # Varianz of the Frequency
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Var"))] = std(dat2[:,11])
                                                                                                                    # Readdepth
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Readdepth"))] = sum(dat2[:,9])
                                                                                                                    # Marker counts
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Markercount"))] = size(dat2,1)
                                                                                                                    # Zero count
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Zerocount"))] = sum(dat2[:,11].== 0)
                                                                                                else
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Average"))] = NaN
                                                                                                                    # Varianz of the Frequency
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Var"))] = NaN
                                                                                                                    # Readdepth
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Readdepth"))] = NaN
                                                                                                                    # Marker counts
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Markercount"))] = NaN
                                                                                                                    # Zero count
                                                                                                                    comp[Symbol(i)][:HapStats][comp[Symbol(i)][:HapStats][!,:x1].== k ,Symbol(string(j,"_Zerocount"))] = NaN
                                                                                                end # if
                                                                                end   # for hap
                                                            end # if any
                                        end # for loop outer



                        ### compare the different generations against each other based on the aggregated and merged replicate sets
                        # using the information defined at the beginning with pstarte, genE, pstartk and genK

                        for i in ["eco", "kon"]    # Generation comparisen based on the presents of the :pvaltab / pstart[e/k]
                                    if haskey(comp[Symbol(i)], :pvalHap)
                                            println("calculate variation in between the generations for $i")
                                            if i == "eco"
                                                        start = pstarte[1]
                                                        test = genE
                                            else
                                                        start = pstartk[1]
                                                        test = genK
                                            end

                                            print(Crayon(foreground = :yellow, bold = false),"Run MannWhitneyUTest calulation in case the SNP sample size per Haplotype is lower 30, otherwise negative binomial + 0 inflated GLM is caluclated, negative Binomial without 0infl or poisson with 0infl \n") # anova model for normal distributed data

                                        for j in test
                                                          println("$j")
                                                          # set the START file against all others are checked against
                                                          if size(freqfile[Symbol(start)],1) == size(freqfile[:Location],1)
                                                                              starter = hcat(freqfile[Symbol(start)],freqfile[:Location][!,Symbol(search)])
                                                          # two replicates
                                                        elseif size(freqfile[Symbol(start)],1) == size(freqfile[:Location],1) * 2
                                                                              starter = hcat(freqfile[Symbol(start)], sort(vcat(freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # three replicates
                                                        elseif size(freqfile[Symbol(start)],1) == size(freqfile[:Location],1) * 3
                                                                              starter = hcat(freqfile[Symbol(start)], sort(vcat(freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # four Replicates
                                                        elseif size(freqfile[Symbol(start)],1) == size(freqfile[:Location],1) * 4
                                                                              starter = hcat(freqfile[Symbol(start)], sort(vcat(freqfile[:Location],freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          else
                                                          println("there is an issue with the count of replicates used in settings aggi=rep, use aggi=mean or start digging")
                                                          #continue
                                                          end

                                                          # set the TEST file
                                                          if size(freqfile[Symbol(j)],1) == size(freqfile[:Location],1)
                                                                              tester = hcat(freqfile[Symbol(j)],freqfile[:Location][!,Symbol(search)])
                                                          # two replicates
                                                          elseif size(freqfile[Symbol(j)],1) == size(freqfile[:Location],1) * 2
                                                                              tester = hcat(freqfile[Symbol(j)], sort(vcat(freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # three replicates
                                                          elseif size(freqfile[Symbol(j)],1) == size(freqfile[:Location],1) * 3
                                                                              tester = hcat(freqfile[Symbol(j)], sort(vcat(freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # four Replicates
                                                          elseif size(freqfile[Symbol(j)],1) == size(freqfile[:Location],1) * 4
                                                                              tester = hcat(freqfile[Symbol(j)], sort(vcat(freqfile[:Location],freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          else
                                                          println("there is an issue with the count of replicates used in settings aggi=rep, use aggi=mean or start digging")
                                                          #continue
                                                          end

                                                          # sort the tester and the started - this is crucial for the alignment with the genes in the hap file
                                                          starter = starter[starter[!,:x1].!= ".",:]
                                                          tester = tester[tester[!,:x1].!= ".",:]
                                                          sort!(starter, (:Chr, :x1))
                                                          sort!(tester, (:Chr, :x1))
                                                          kap = unique(starter[!,:x1])

                @showprogress        for k in kap
                                                           #println(k)
                                                           #= the problem that has to be dealed with is the same as with GLM model, just a bit less complicated
                                                           MannWhitneyUTest allows for unequal size of observations, so one does not have to care about equality
                                                           make sure for each possible combination of replicates for start and test generation the value is calculated correctly=#
                                                           datT = callsnps(tester,k)
                                                           datS = callsnps(starter,k)
                                                           #= only for development
                                                           datT = tester[tester[!,:x1].==k,:]
                                                           datS = starter[starter[!,:x1].==k,:]
                                                           =#

                                                           # if or if not to perform a transformation of the frequency based of the read counts
                                                           if transform == false
                                                                               datT[!,:Allelfreq_Isr] = Float64.(datT[!,:Allelfreq_Isr])  # log.(datT[!,:Readcount]) .* datT[!,:Allelfreq_Isr]
                                                                               datS[!,:Allelfreq_Isr] =  Float64.(datS[!,:Allelfreq_Isr]) # log.(datS[!,:Readcount]) .* datS[!,:Allelfreq_Isr]
                                                          elseif transform == true
                                                                                datT[!,:Allelfreq_Isr] = log.(datT[!,:Readcount]) .* datT[!,:Allelfreq_Isr]
                                                                                datS[!,:Allelfreq_Isr] = log.(datS[!,:Readcount]) .* datS[!,:Allelfreq_Isr]
                                                          end
                                                          # NaN values have to be removed from the list
                                                          datS = datS[.!(isnan.(datS[!,:Allelfreq_Isr])),:]
                                                          datT = datT[.!(isnan.(datT[!,:Allelfreq_Isr])),:]

                                                          # test if one of the sets is empty - if so no next itteration can be started
                                                          if isempty(datS) || isempty(datT)
                                                                                          comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_negbinManU"))] = 5.0
                                                                                          comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_zeroInfl"))] = 5.0
                                                                                          comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_Stat"))] = "MissData"
                                                                                            continue
                                                          end

                                                          # if the sample size is smaller than 30, a nonparmatieric test will be applied - if bigger, a GLM model will be used
                                                          # also , even if the samplesize is bigger 30, it might happen that at one generation, there are only 0 values present - in this case no distribution can be assumed for both sets - this will be only applies if there is one set with completly only 0 values!!
                                                           if size(datS[!,:Allelfreq_Isr],1) < 15 && size(datT[!,:Allelfreq_Isr],1) < 15 || all(datS[!,:Allelfreq_Isr] .== 0) || all(datT[!,:Allelfreq_Isr] .== 0 )
                                                                                     # caluculate the Allelfreq for the whole Haplotype new - use all reads aligned to this region for CALCULATION
                                                                                     # 2.0 indicates a missing value - zero inflated column will be always 2.0 for samples with less than 30 SNPS
                                                                                      #println("MannWhitneyUTest")
                                                                                         comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_negbinManU"))] = pvalue(MannWhitneyUTest(datS[!,:Allelfreq_Isr], datT[!,:Allelfreq_Isr]))
                                                                                         comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_zeroInfl"))] = 2.0
                                                                                         comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_Stat"))] = "MannWhitneyUTest"
                                                            else
                                                                                # create levels where to disect the sets at
                                                                                  datS[!,:set] .= 1
                                                                                  datT[!,:set] .= 2
                                                                                  n3 = vcat(datS, datT)
                                                                                  n3[!,:Allelfreq_Isr] = Int64.(floor.(n3[!,:Allelfreq_Isr] .* 100))
                                                                                  # when there are no 0 values in it, a zero inflated model does not make sense - in this case, a NegativeBinomial Distribution is used only
                                                                                  # if less than 10% of the samples are 0, the distribution will be set to a NegativeBinomial
                                                                                  # with more than 75% zeros in the distribution, negbin zeroinfl will create NaN for the Zero infl model - change to poisson 0infl distribution
                                                                                  if (sum(n3[!,:Allelfreq_Isr] .== 0) / size(n3,1)) < 0.15 || all(datS[!,:Allelfreq_Isr] .!= 0) || all(datT[!,:Allelfreq_Isr] .!= 0 )
                                                                                                                    #println("Julia GLM")
                                                                                                                    categorical!(n3, :set)
                                                                                                                    comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_negbinManU"))] = coeftable(glm(@formula(Allelfreq_Isr ~ set), n3, NegativeBinomial(), LogLink())).cols[4][2]
                                                                                                                    comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_zeroInfl"))] = 3.0
                                                                                                                    comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_Stat"))] = "NegBino"
                                                                                 elseif (sum(n3[!,:Allelfreq_Isr] .== 0) / size(n3,1)) < 0.75
                                                                                                                    #println("ZeroInfl. Negbin")
                                                                                                                    @rput n3
                                                                                                                     # if the data has too many zeros compared to non-0, than a negative binomial model does not work - switch to poisson distribution
                                                                                                                     R"""
                                                                                                                      library(MASS)
                                                                                                                      library(pscl)
                                                                                                                      n3[,12] = as.factor(n3[,12])
                                                                                                                      m1 =  summary(zeroinfl(Allelfreq_Isr ~ set, n3, dist = "negbin"))
                                                                                                                      negbino = m1$coefficients$count[,4]
                                                                                                                      zeroinf = m1$coefficients$zero[,4]
                                                                                                                      """
                                                                                                                      @rget negbino  zeroinf
                                                                                                                      comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_negbinManU"))] = negbino[2]
                                                                                                                      comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_zeroInfl"))] = zeroinf[2]
                                                                                                                      comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_Stat"))] = "NbZi"
                                                                                else
                                                                                                                      #println("ZeroInfl. Poisson")
                                                                                                                try
                                                                                                                  @rput n3
                                                                                                                      R"""
                                                                                                                      library(MASS)
                                                                                                                      library(pscl)
                                                                                                                      n3[,12] = as.factor(n3[,12])
                                                                                                                      m1 =  summary(zeroinfl(Allelfreq_Isr ~ set, n3, dist = "poisson"))
                                                                                                                      negbino = m1$coefficients$count[,4]
                                                                                                                      zeroinf = m1$coefficients$zero[,4]
                                                                                                                      """
                                                                                                                      @rget negbino  zeroinf
                                                                                                                      #print(Crayon(foreground = :red, bold = false),"Poisson Distribution had to be used - too many Zero values! \n")
                                                                                                                      #print(Crayon(foreground = :white, bold = false),"")
                                                                                                                      comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_negbinManU"))] = negbino[2]
                                                                                                                      comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_zeroInfl"))] = zeroinf[2]
                                                                                                                      comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_Stat"))] = "PsZi"

                                                                                                                    catch
                                                                                                                          print(Crayon(foreground = :red, bold = true),"\nHaplotype $k could not be calculated!\n")
                                                                                                                          comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_negbinManU"))] = 5.0
                                                                                                                          comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_zeroInfl"))] = 5.0
                                                                                                                          comp[Symbol(i)][:pvalHap][comp[Symbol(i)][:pvalHap][!,:x1].== k,Symbol(string(j,"_Stat"))] = "none"
                                                                                                                          continue
                                                                                                                    end

                                                                                 end # NegBin or ZeroInfl
                                                         end # of statistical test


                                                      end    # end of for hap
                                                end # end of for test

                                            end # end of haskey

                                        end # end of eco/kon

                                            # select the Haplotypes that have changes significantly over the selected generation - write in a new dictionary entry called results
                                          #=  for i in ["eco", "kon"], j in names(comp[Symbol(i)][:pvalHap][:,2:size(comp[Symbol(i)][:pvalHap],2)])
                                                                # join pvalHap and HapStats for significant markers
                                                                try
                                                                comp[Symbol(i)][Symbol(string("results_", j))] =  join(comp[Symbol(i)][:pvalHap][.&(completecases(comp[Symbol(i)][:pvalHap],Symbol(j)),comp[Symbol(i)][:pvalHap][!,Symbol(j)] .< 0.05),:], comp[Symbol(i)][:HapStats], on= :x1, kind = :inner)
                                                                # add info from the freqfile[:Location] DataFrame
                                                                #comp[Symbol(i)][Symbol(string("results_", j))] = join(comp[Symbol(i)][Symbol(string("results_", j))], freqfile[:Location][[Symbol("Gene_ID"),Symbol("chromosome:start-stop"),Symbol("description")]], on= (:x1, :Gene_ID), kind=:left)
                                                                #rename!(comp[Symbol(i)][Symbol(string("results_", j))], Symbol(j) => :pvalue)
                                                                sort!(comp[Symbol(i)][Symbol(string("results_", j))], Symbol(j))
                                                                catch
                                                                continue
                                                                end
                                                    end =#
                                               for i in ["eco", "kon"]
                                                                      try
                                                                          comp[Symbol(i)][Symbol(string("results"))] = hcat(comp[Symbol(i)][:pvalHap], comp[Symbol(i)][:HapStats], makeunique=true)
                                                                        catch
                                                                          continue
                                                                        end
                                                                      end

                        #eco[:pvalHap][.&(completecases(eco[:pvalHap],:f21e2),eco[:pvalHap][:f21e2] .< 0.1),:]
                        #comp[:eco][:HapStats][comp[:eco][:HapStats][:x1] .== "GO:0006807",:]
                        print(Crayon(foreground = :green),"Haplotype calling has been done and generations have been compared successfully\n")
                        return comp
                    end # of evol function

#=
the following function performes a comaprisen between the two eco systems, konventinel and organic farming
for each system, one generation is entered - the replicates will be merged togehter and the systems are compared on single haplotype level
if kon=NA ,the same generation for organic and konventional will be compared
transformation in regard to the correction of the read depth per snp marker
=#
function scomp(freqfile, genolist, org=23, kon="NA", haplotype="gene", transform=false)

                                    # A BASE FUNCTION THAT IS NEEDED FOR HAPSTATS AND PVALHAP CALULCATION
                                    function callsnps(dat,k)
                                                            m = Int32[]
                                                            for i in 1:nrow(dat)
                                                                                if dat[i,:x1] != k
                                                                                append!(m,i)
                                                                                break
                                                                                end
                                                                            end
                                                          if  !(isempty(m))
                                                             dat2 = dat[1:m[1]-1,:]
                                                             deleterows!(dat, 1:m[1]-1) # remove the gene snps that have already been worked on
                                                           else
                                                             dat2 = dat[1:nrow(dat),:]
                                                           end

                                                        return dat2
                                   end

                                  nam = CSV.read(genolist, delim='\t') # Env & Generation = 0 indicate the parents

                                  # get a unique list of all GO terms without "none" and "." ; lstrip - remove space from front if exists; vcat(split) - get a vector kind list from the GO terms in the Location subdict
                                  if haplotype == "gene"
                                                          hap = unique(freqfile[:Location][!,:Gene_ID])
                                                          hap = hap[hap.!="."]
                                  elseif haplotype == "GO"
                                                          hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...))[occursin.("GO",unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...)))])
                                                          hap = hap[hap.!="."]
                                  elseif haplotype == "Marker"
                                                          hap = unique(freqfile[:Location][!,:Marker])
                                                          hap = hap[hap.!="."]
                                  elseif haplotype == "PF"
                                                          hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...))[occursin.("PF",unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...)))])
                                                          hap = hap[hap.!="."]

                                  else
                                  println("wrong haplotype specified, chose 'gene', 'GO' or 'PF' ")

                                  end

                                  # split the general haplotype description part from the part where the generations are compared to each other
                                  #set the column to search in
                                  if haplotype == "gene"
                                                          search = "Gene_ID"
                                  elseif haplotype == "GO"
                                                          search = "GO-IDs"
                                  elseif haplotype == "Marker"
                                                          search = "Marker"
                                  elseif haplotype == "PF"
                                                          search = "PFAM-IDs"
                                  else
                                  println("Do you really know what you are doing?!")

                                  end

                                  ## get the information which sets to merge
                                  organic = nam[.&(nam[!,:Env].==2, nam[!,:Generation].==org),:Name]
                                  # if there is no 2nd generation specified, the same generation for konventional is choosen
                                  if kon != "NA"
                                  konventinel =  nam[.&(nam[!,:Env].==1, nam[!,:Generation].==kon),:Name]
                                  else
                                  konventinel =  nam[.&(nam[!,:Env].==1, nam[!,:Generation].==org),:Name]
                                  end

                                  # get the replicates merged togehter

                                  # put all in one file
                                  if size(konventinel,1) == 1
                                                                    konv = freqfile[Symbol(konventinel[1])]
                                  elseif size(konventinel,1) == 2
                                                                    konv = sort(vcat(freqfile[Symbol(konventinel[1])], freqfile[Symbol(konventinel[2])]), (:Chr, :Pos))
                                  elseif size(konventinel,1) == 3
                                                                    konv = sort(vcat(freqfile[Symbol(konventinel[1])], freqfile[Symbol(konventinel[2])], freqfile[Symbol(konventinel[3])]), (:Chr, :Pos))
                                  elseif size(konventinel,1) == 4
                                                                    konv = sort(vcat(freqfile[Symbol(konventinel[1])], freqfile[Symbol(konventinel[2])], freqfile[Symbol(konventinel[3])], freqfile[Symbol(konventinel[4])]), (:Chr, :Pos))
                                  elseif size(konventinel,1) == 5
                                                                  konv = sort(vcat(freqfile[Symbol(konventinel[1])], freqfile[Symbol(konventinel[2])], freqfile[Symbol(konventinel[3])], freqfile[Symbol(konventinel[4])], freqfile[Symbol(konventinel[5])]), (:Chr, :Pos))

                                  else
                                      println("Replicates could not be matched in one file for the konventinal system")
                                  end

                                  if size(organic,1) == 1
                                                                  orga = freqfile[Symbol(organic[1])]
                                  elseif size(organic,1) == 2
                                                                  orga = sort(vcat(freqfile[Symbol(organic[1])], freqfile[Symbol(organic[2])]), (:Chr, :Pos))
                                  elseif size(organic,1) == 3
                                                                  orga = sort(vcat(freqfile[Symbol(organic[1])], freqfile[Symbol(organic[2])], freqfile[Symbol(organic[3])]), (:Chr, :Pos))
                                  elseif size(organic,1) == 4
                                                                  orga = sort(vcat(freqfile[Symbol(organic[1])], freqfile[Symbol(organic[2])], freqfile[Symbol(organic[3])], freqfile[Symbol(organic[4])]), (:Chr, :Pos))
                                  elseif size(organic,1) == 5
                                                                  orga = sort(vcat(freqfile[Symbol(organic[1])], freqfile[Symbol(organic[2])], freqfile[Symbol(organic[3])], freqfile[Symbol(organic[4])], freqfile[Symbol(organic[5])]), (:Chr, :Pos))
                                  else
                                    println("Replicates could not be matched in one file for the organic system")
                                  end



                                  # link the location data to the konv and orga files
                                                          ###################### konventional #######################
                                                          if any(occursin.("eLoc",string.(keys(freqfile)))) # when there is already a file with thecorrect size, use this to save the merging time
                                                                              konv = hcat(konv[!,[9,11]],freqfile[:eLoc][!,Symbol(search)])
                                                          elseif size(konv,1) == size(freqfile[:Location],1)
                                                                              konv = hcat(konv[!,[9,11]],freqfile[:Location][!,Symbol(search)])
                                                          # two replicates
                                                          elseif size(konv,1) == size(freqfile[:Location],1) * 2
                                                                              konv = hcat(konv[!,[9,11]], sort(vcat(freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # three replicates
                                                          elseif size(konv,1) == size(freqfile[:Location],1) * 3
                                                                              konv = hcat(konv[!,[9,11]], sort(vcat(freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # four Replicates
                                                          elseif size(konv,1) == size(freqfile[:Location],1) * 4
                                                                              konv = hcat(konv[!,[9,11]], sort(vcat(freqfile[:Location],freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          else
                                                          println("there is an issue with the count of replicates used in settings aggi=rep, use aggi=mean or start digging")
                                                          #continue
                                                          end
                                                          ########################### organic ###################
                                                          if any(occursin.("eLoc",string.(keys(freqfile)))) # when there is already a file with thecorrect size, use this to save the merging time
                                                                              orga = hcat(orga[!,[9,11]],freqfile[:eLoc][!,Symbol(search)])
                                                          elseif size(orga,1) == size(freqfile[:Location],1)
                                                                              orga = hcat(orga[!,[9,11]],freqfile[:Location][!,Symbol(search)])
                                                          # two replicates
                                                          elseif size(orga,1) == size(freqfile[:Location],1) * 2
                                                                              orga = hcat(orga[!,[9,11]], sort(vcat(freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # three replicates
                                                          elseif size(orga,1) == size(freqfile[:Location],1) * 3
                                                                              orga = hcat(orga[!,[9,11]], sort(vcat(freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # four Replicates
                                                          elseif size(orga,1) == size(freqfile[:Location],1) * 4
                                                                              orga = hcat(orga[!,[9,11]], sort(vcat(freqfile[:Location],freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          else
                                                          println("there is an issue with the count of replicates used in settings aggi=rep, use aggi=mean or start digging")
                                                          #continue
                                                          end

                                                          orga = orga[orga[!,:x1].!= ".",:] # remove not annotated snps from the list, highlighted by a "." in the x1 gene column
                                                          konv = konv[konv[!,:x1].!= ".",:]

                                                          kap = unique(konv[!,:x1])

                                                          # sort orga, kap and konv by the haplotyes
                                                          sort!(konv, :x1)
                                                          sort!(orga, :x1)
                                                          sort!(kap)

                                  # set the dataframe where to save the results
                                  result = DataFrame(Array{Union{Nothing, Any}}(nothing, length(kap),14))
                                  # first colum will be kap, 2nd is NegativeBinomial and 3rd is zero inflated - after this the overall stats are following
                                  calnames = [:Haplotype, :NegBin, :Zeroinfl, :StatTest]
                                  rename!(result,vcat(calnames, Symbol.(string.(repeat(["org", "kon"],5),repeat(["Average", "Var", "Readdepth", "Markercount", "Zerocount"],inner=2)))))
                                  # give the columns a type
                                  for i in names(result)[2:14]; result[!,i] .= 0.0 ::Float64 ; end
                                  result[!,:StatTest] = string.(result[!,:StatTest]) # the teststat has to be a string
                                  result[!,:Haplotype] .= kap

                                  @showprogress for k in kap
                                                                                             #println(k)
                                                                                             konv2 = callsnps(konv,k)
                                                                                             #konv2 = konv[konv[!,:x1].==k,:]
                                                                                             konv2 = konv2[.!(isnan.(konv2[:,:2])),:] # remove the NaN values from the table
                                                                                             orga2 = callsnps(orga,k)
                                                                                             #orga2 = orga[orga[!,:x1].==k,:]
                                                                                             orga2 = orga2[.!(isnan.(orga2[:,:2])),:]

                                                                                              # konv
                                                                                              if !(isempty(konv2))
                                                                                                                  result[result[!,:Haplotype].==k,:konAverage] = sum(konv2[:,2] .* konv2[:,1]) / sum(konv2[:,1])
                                                                                                                  # Varianz of the Frequency
                                                                                                                  result[result[!,:Haplotype].==k,:konVar] = std(konv2[:,2])
                                                                                                                  # Readdepth
                                                                                                                  result[result[!,:Haplotype].==k,:konReaddepth] = sum(konv2[:,1])
                                                                                                                  # Marker counts
                                                                                                                  result[result[!,:Haplotype].==k,:konMarkercount] = size(konv2,1)
                                                                                                                  # Zero counts
                                                                                                                  result[result[!,:Haplotype].==k,:konZerocount] = sum(konv2[:,2].== 0)
                                                                                              else
                                                                                                                  result[result[!,:Haplotype].==k,:konAverage] = NaN
                                                                                                                  # Varianz of the Frequency
                                                                                                                  result[result[!,:Haplotype].==k,:konVar] = NaN
                                                                                                                  # Readdepth
                                                                                                                  result[result[!,:Haplotype].==k,:konReaddepth] = NaN
                                                                                                                  # Marker counts
                                                                                                                  result[result[!,:Haplotype].==k,:konMarkercount] = NaN
                                                                                                                  # Zero counts
                                                                                                                  result[result[!,:Haplotype].==k,:konZerocount] = NaN
                                                                                              end # if
                                                                                              # organic
                                                                                              if !(isempty(orga2))
                                                                                                                  result[result[!,:Haplotype].==k,:orgAverage] = sum(orga2[:,2] .* orga2[:,1]) / sum(orga2[:,1])
                                                                                                                  # Varianz of the Frequency
                                                                                                                  result[result[!,:Haplotype].==k,:orgVar] = std(orga2[:,2])
                                                                                                                  # Readdepth
                                                                                                                  result[result[!,:Haplotype].==k,:orgReaddepth] = sum(orga2[:,1])
                                                                                                                  # Marker counts
                                                                                                                  result[result[!,:Haplotype].==k,:orgMarkercount] = size(orga2,1)
                                                                                                                  # Zero counts
                                                                                                                  result[result[!,:Haplotype].==k,:orgZerocount] = sum(orga2[:,2].== 0)
                                                                                              else
                                                                                                                  result[result[!,:Haplotype].==k,:orgAverage] = NaN
                                                                                                                  # Varianz of the Frequency
                                                                                                                  result[result[!,:Haplotype].==k,:orgVar] = NaN
                                                                                                                  # Readdepth
                                                                                                                  result[result[!,:Haplotype].==k,:orgReaddepth] = NaN
                                                                                                                  # Marker counts
                                                                                                                  result[result[!,:Haplotype].==k,:orgMarkercount] = NaN
                                                                                                                  # Zero counts
                                                                                                                  result[result[!,:Haplotype].==k,:orgZerocount] = NaN
                                                                                              end # if


                                                                            # calculate the stats for the pair of data
                                                                            # if or if not to perform a transformation of the frequency based of the read counts
                                                                            if transform == false
                                                                                                jts = Float64.(orga2[!,:Allelfreq_Isr])  # log.(datT[!,:Readcount]) .* datT[!,:Allelfreq_Isr]
                                                                                                sta =  Float64.(konv2[!,:Allelfreq_Isr]) # log.(datS[!,:Readcount]) .* datS[!,:Allelfreq_Isr]
                                                                           elseif transform == true
                                                                                                 jts = log.(orga2[!,:Readcount]) .* orga2[!,:Allelfreq_Isr]
                                                                                                 sta = log.(konv2[!,:Readcount]) .* konv2[!,:Allelfreq_Isr]
                                                                           end
                                                                            # remove the NaN values
                                                                            sta = sta[.!(isnan.(sta))]
                                                                            jts = jts[.!(isnan.(jts))]

                                                                            # test if one of the sets is empty - if so no next itteration can be started
                                                                               if isempty(sta) || isempty(jts)
                                                                                                               result[result[!,:Haplotype].==k,:NegBin] = 5.0
                                                                                                               result[result[!,:Haplotype].==k,:Zeroinfl]   = 5.0
                                                                                                               result[result[!,:Haplotype].==k,:StatTest] = "MissData"
                                                                                                               continue
                                                                               end

                                                                           # if the sample size is smaller than 30, a nonparmatieric test will be applied - if bigger, a GLM model will be used
                                                                                if size(sta,1) < 15 && size(jts,1) < 15 || all(sta .== 0) || all(jts .== 0 )
                                                                                                      # caluculate the Allelfreq for the whole Haplotype new - use all reads aligned to this region for CALCULATION
                                                                                                      # 2.0 indicates a missing value - zero inflated column will be always 2.0 for samples with less than 30 SNPS
                                                                                                       #println("MannWhitneyUTest")
                                                                                                        result[result[!,:Haplotype].==k,:NegBin] = pvalue(MannWhitneyUTest(sta, jts))
                                                                                                        result[result[!,:Haplotype].==k,:Zeroinfl]   = 2.0
                                                                                                        result[result[!,:Haplotype].==k,:StatTest] = "MannWhitneyUTest"
                                                                             else

                                                                                                   t1 = vcat(sta, jts)
                                                                                                   #if isempty(sta) || isempty(jts) ; continue ; end
                                                                                                   t2 = repeat([1], size(sta,1))
                                                                                                   t3 = repeat([2], size(jts,1))
                                                                                                   n3 = DataFrame(hcat(t1, vcat(t2,t3)))
                                                                                                   n3[!,1] = Int64.(floor.(n3[!,1] .* 100))

                                                                                                   # when there are no 0 vaulues in it, a zero inflated model does not make sense - in this case, a NegativeBinomial Distribution is used only
                                                                                                   # if less than 10% of the samples are 0, the distribution will be set to a NegativeBinomial
                                                                                                   # with more than 75% zeros in the distribution, negbin zeroinfl will create NaN for the Zero infl model - change to poisson 0infl distribution
                                                                                                   if (sum(n3[!,:x1] .== 0) / size(n3,1)) < 0.15 || all(sta .!= 0) || all(jts .!= 0 )
                                                                                                                                     #println("Julia GLM")
                                                                                                                                     categorical!(n3, :x2)
                                                                                                                                     result[result[!,:Haplotype].==k,:NegBin] = coeftable(glm(@formula(x1 ~ x2), n3, NegativeBinomial(), LogLink())).cols[4][2]
                                                                                                                                     result[result[!,:Haplotype].==k,:Zeroinfl] = 3.0
                                                                                                                                     result[result[!,:Haplotype].==k,:StatTest] = "NegBino"
                                                                                                  elseif (sum(n3[!,:x1] .== 0) / size(n3,1)) < 0.75
                                                                                                                                     #println("ZeroInfl. Negbin")
                                                                                                                                     @rput n3
                                                                                                                                      # if the data has too many zeros compared to non-0, than a negative binomial model does not work - switch to poisson distribution
                                                                                                                                       R"""
                                                                                                                                       library(MASS)
                                                                                                                                       library(pscl)
                                                                                                                                       n3[,2] = as.factor(n3[,2])
                                                                                                                                       m1 =  summary(zeroinfl(x1~x2, n3, dist = "negbin"))
                                                                                                                                       negbino = m1$coefficients$count[,4]
                                                                                                                                       zeroinf = m1$coefficients$zero[,4]
                                                                                                                                       """
                                                                                                                                       @rget negbino  zeroinf
                                                                                                                                        result[result[!,:Haplotype].==k,:NegBin] = negbino[2]
                                                                                                                                      result[result[!,:Haplotype].==k,:Zeroinfl] = zeroinf[2]
                                                                                                                                      result[result[!,:Haplotype].==k,:StatTest] = "NbZi"
                                                                                                 else
                                                                                                                                       #println("ZeroInfl. Poisson")
                                                                                                                                try
                                                                                                                                       @rput n3
                                                                                                                                       R"""
                                                                                                                                       library(MASS)
                                                                                                                                       library(pscl)
                                                                                                                                       n3[,2] = as.factor(n3[,2])
                                                                                                                                       m1 =  summary(zeroinfl(x1~x2, n3, dist = "poisson"))
                                                                                                                                       negbino = m1$coefficients$count[,4]
                                                                                                                                       zeroinf = m1$coefficients$zero[,4]
                                                                                                                                       """
                                                                                                                                       @rget negbino  zeroinf
                                                                                                                                       #print(Crayon(foreground = :red, bold = false),"Poisson Distribution had to be used - too many Zero values! \n")
                                                                                                                                       #print(Crayon(foreground = :white, bold = false),"")
                                                                                                                                        result[result[!,:Haplotype].==k,:NegBin] = negbino[2]
                                                                                                                                        result[result[!,:Haplotype].==k,:Zeroinfl] = zeroinf[2]
                                                                                                                                        result[result[!,:Haplotype].==k,:StatTest] = "PsZi"

                                                                                                                                    catch
                                                                                                                                          print(Crayon(foreground = :red, bold = true),"\nHaplotype $k could not be calculated!\n")
                                                                                                                                          result[result[!,:Haplotype].==k,:NegBin] = 5.0
                                                                                                                                          result[result[!,:Haplotype].==k,:Zeroinfl] = 5.0
                                                                                                                                          result[result[!,:Haplotype].==k,:StatTest] = "none"
                                                                                                                                          continue
                                                                                                                                    end
                                                                                                end # NegBin or ZeroInfl
                                                                          end # of statistical test
                                                                end   # for hap

                                  return result


 end # of scomp

#=
- haplotypes that are differnt over the the generations but equal between the systems - equivalent to the scomp function above, just within the system - e.g. for environmental factors that effected both systems
- hap. different
=#
function gcomp(freqfile, genolist, env=1, first=22, second=23, haplotype="gene", transform=false)

                                    # A BASE FUNCTION THAT IS NEEDED FOR HAPSTATS AND PVALHAP CALULCATION
                                    function callsnps(dat,k)
                                                            m = Int32[]
                                                            for i in 1:nrow(dat)
                                                                                if dat[i,:x1] != k
                                                                                append!(m,i)
                                                                                break
                                                                                end
                                                                            end
                                                          if  !(isempty(m))
                                                             dat2 = dat[1:m[1]-1,:]
                                                             deleterows!(dat, 1:m[1]-1) # remove the gene snps that have already been worked on
                                                           else
                                                             dat2 = dat[1:nrow(dat),:]
                                                           end

                                                        return dat2
                                   end

                                  nam = CSV.read(genolist, delim='\t') # Env & Generation = 0 indicate the parents

                                  # get a unique list of all GO terms without "none" and "." ; lstrip - remove space from front if exists; vcat(split) - get a vector kind list from the GO terms in the Location subdict
                                  if haplotype == "gene"
                                                          hap = unique(freqfile[:Location][!,:Gene_ID])
                                                          hap = hap[hap.!="."]
                                  elseif haplotype == "GO"
                                                          hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...))[occursin.("GO",unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...)))])
                                                          hap = hap[hap.!="."]
                                  elseif haplotype == "Marker"
                                                          hap = unique(freqfile[:Location][!,:Marker])
                                                          hap = hap[hap.!="."]
                                  elseif haplotype == "PF"
                                                          hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...))[occursin.("PF",unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...)))])
                                                          hap = hap[hap.!="."]

                                  else
                                  println("wrong haplotype specified, chose 'gene', 'GO' or 'PF' ")

                                  end

                                  # split the general haplotype description part from the part where the generations are compared to each other
                                  #set the column to search in
                                  if haplotype == "gene"
                                                          search = "Gene_ID"
                                  elseif haplotype == "GO"
                                                          search = "GO-IDs"
                                  elseif haplotype == "Marker"
                                                          search = "Marker"
                                  elseif haplotype == "PF"
                                                          search = "PFAM-IDs"
                                  else
                                  println("Do you really know what you are doing?!")

                                  end

                                  ## get the information which sets to merge
                                  organic = nam[.&(nam[!,:Env].==env, nam[!,:Generation].==first),:Name]
                                  println("the first test set is $organic")
                                  konventinel =  nam[.&(nam[!,:Env].==env, nam[!,:Generation].==second),:Name]
                                  println("The second test set is $konventinel")

                                  # get the replicates merged togehter

                                  # put all in one file
                                  if size(konventinel,1) == 1
                                                                    konv = freqfile[Symbol(konventinel[1])]
                                  elseif size(konventinel,1) == 2
                                                                    konv = sort(vcat(freqfile[Symbol(konventinel[1])], freqfile[Symbol(konventinel[2])]), (:Chr, :Pos))
                                  elseif size(konventinel,1) == 3
                                                                    konv = sort(vcat(freqfile[Symbol(konventinel[1])], freqfile[Symbol(konventinel[2])], freqfile[Symbol(konventinel[3])]), (:Chr, :Pos))
                                  elseif size(konventinel,1) == 4
                                                                    konv = sort(vcat(freqfile[Symbol(konventinel[1])], freqfile[Symbol(konventinel[2])], freqfile[Symbol(konventinel[3])], freqfile[Symbol(konventinel[4])]), (:Chr, :Pos))
                                  elseif size(konventinel,1) == 5
                                                                  konv = sort(vcat(freqfile[Symbol(konventinel[1])], freqfile[Symbol(konventinel[2])], freqfile[Symbol(konventinel[3])], freqfile[Symbol(konventinel[4])], freqfile[Symbol(konventinel[5])]), (:Chr, :Pos))

                                  else
                                      println("Replicates could not be matched in one file for the konventinal system")
                                  end

                                  if size(organic,1) == 1
                                                                  orga = freqfile[Symbol(organic[1])]
                                  elseif size(organic,1) == 2
                                                                  orga = sort(vcat(freqfile[Symbol(organic[1])], freqfile[Symbol(organic[2])]), (:Chr, :Pos))
                                  elseif size(organic,1) == 3
                                                                  orga = sort(vcat(freqfile[Symbol(organic[1])], freqfile[Symbol(organic[2])], freqfile[Symbol(organic[3])]), (:Chr, :Pos))
                                  elseif size(organic,1) == 4
                                                                  orga = sort(vcat(freqfile[Symbol(organic[1])], freqfile[Symbol(organic[2])], freqfile[Symbol(organic[3])], freqfile[Symbol(organic[4])]), (:Chr, :Pos))
                                  elseif size(organic,1) == 5
                                                                  orga = sort(vcat(freqfile[Symbol(organic[1])], freqfile[Symbol(organic[2])], freqfile[Symbol(organic[3])], freqfile[Symbol(organic[4])], freqfile[Symbol(organic[5])]), (:Chr, :Pos))
                                  else
                                    println("Replicates could not be matched in one file for the organic system")
                                  end



                                  # link the location data to the konv and orga files
                                                          ###################### konventional #######################
                                                          if any(occursin.("eLoc",string.(keys(freqfile)))) # when there is already a file with thecorrect size, use this to save the merging time
                                                                              konv = hcat(konv[!,[9,11]],freqfile[:eLoc][!,Symbol(search)])
                                                          elseif size(konv,1) == size(freqfile[:Location],1)
                                                                              konv = hcat(konv[!,[9,11]],freqfile[:Location][!,Symbol(search)])
                                                          # two replicates
                                                          elseif size(konv,1) == size(freqfile[:Location],1) * 2
                                                                              konv = hcat(konv[!,[9,11]], sort(vcat(freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # three replicates
                                                          elseif size(konv,1) == size(freqfile[:Location],1) * 3
                                                                              konv = hcat(konv[!,[9,11]], sort(vcat(freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # four Replicates
                                                          elseif size(konv,1) == size(freqfile[:Location],1) * 4
                                                                              konv = hcat(konv[!,[9,11]], sort(vcat(freqfile[:Location],freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          else
                                                          println("there is an issue with the count of replicates used in settings aggi=rep, use aggi=mean or start digging")
                                                          #continue
                                                          end
                                                          ########################### organic ###################
                                                          if any(occursin.("eLoc",string.(keys(freqfile)))) # when there is already a file with thecorrect size, use this to save the merging time
                                                                              orga = hcat(orga[!,[9,11]],freqfile[:eLoc][!,Symbol(search)])
                                                          elseif size(orga,1) == size(freqfile[:Location],1)
                                                                              orga = hcat(orga[!,[9,11]],freqfile[:Location][!,Symbol(search)])
                                                          # two replicates
                                                          elseif size(orga,1) == size(freqfile[:Location],1) * 2
                                                                              orga = hcat(orga[!,[9,11]], sort(vcat(freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # three replicates
                                                          elseif size(orga,1) == size(freqfile[:Location],1) * 3
                                                                              orga = hcat(orga[!,[9,11]], sort(vcat(freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          # four Replicates
                                                          elseif size(orga,1) == size(freqfile[:Location],1) * 4
                                                                              orga = hcat(orga[!,[9,11]], sort(vcat(freqfile[:Location],freqfile[:Location], freqfile[:Location], freqfile[:Location]),(:Chr, :Pos))[!,Symbol(search)])

                                                          else
                                                          println("there is an issue with the count of replicates used in settings aggi=rep, use aggi=mean or start digging")
                                                          #continue
                                                          end

                                                          orga = orga[orga[!,:x1].!= ".",:] # remove not annotated snps from the list, highlighted by a "." in the x1 gene column
                                                          konv = konv[konv[!,:x1].!= ".",:]

                                                          kap = unique(konv[!,:x1])

                                                          # sort orga, kap and konv by the haplotyes
                                                          sort!(konv, :x1)
                                                          sort!(orga, :x1)
                                                          sort!(kap)

                                  # set the dataframe where to save the results
                                  result = DataFrame(Array{Union{Nothing, Any}}(nothing, length(kap),14))
                                  # first colum will be kap, 2nd is NegativeBinomial and 3rd is zero inflated - after this the overall stats are following
                                  calnames = [:Haplotype, :NegBin, :Zeroinfl, :StatTest]
                                  rename!(result,vcat(calnames, Symbol.(string.(repeat(["fst", "sec"],5),repeat(["Average", "Var", "Readdepth", "Markercount", "Zerocount"],inner=2)))))
                                  # give the columns a type
                                  for i in names(result)[2:14]; result[!,i] .= 0.0 ::Float64 ; end
                                  result[!,:StatTest] = string.(result[!,:StatTest]) # the teststat has to be a string
                                  result[!,:Haplotype] .= kap

                                  @showprogress for k in kap
                                                                                             #println(k)
                                                                                             konv2 = callsnps(konv,k)
                                                                                             #konv2 = konv[konv[!,:x1].==k,:]
                                                                                             konv2 = konv2[.!(isnan.(konv2[:,:2])),:] # remove the NaN values from the table
                                                                                             orga2 = callsnps(orga,k)
                                                                                             #orga2 = orga[orga[!,:x1].==k,:]
                                                                                             orga2 = orga2[.!(isnan.(orga2[:,:2])),:]

                                                                                              # konv
                                                                                              if !(isempty(konv2))
                                                                                                                  result[result[!,:Haplotype].==k,:secAverage] = sum(konv2[:,2] .* konv2[:,1]) / sum(konv2[:,1])
                                                                                                                  # Varianz of the Frequency
                                                                                                                  result[result[!,:Haplotype].==k,:secVar] = std(konv2[:,2])
                                                                                                                  # Readdepth
                                                                                                                  result[result[!,:Haplotype].==k,:secReaddepth] = sum(konv2[:,1])
                                                                                                                  # Marker counts
                                                                                                                  result[result[!,:Haplotype].==k,:secMarkercount] = size(konv2,1)
                                                                                                                  # Zero counts
                                                                                                                  result[result[!,:Haplotype].==k,:secZerocount] = sum(konv2[:,2].== 0)
                                                                                              else
                                                                                                                  result[result[!,:Haplotype].==k,:secAverage] = NaN
                                                                                                                  # Varianz of the Frequency
                                                                                                                  result[result[!,:Haplotype].==k,:secVar] = NaN
                                                                                                                  # Readdepth
                                                                                                                  result[result[!,:Haplotype].==k,:secReaddepth] = NaN
                                                                                                                  # Marker counts
                                                                                                                  result[result[!,:Haplotype].==k,:secMarkercount] = NaN
                                                                                                                  # Zero counts
                                                                                                                  result[result[!,:Haplotype].==k,:secZerocount] = NaN
                                                                                              end # if
                                                                                              # organic
                                                                                              if !(isempty(orga2))
                                                                                                                  result[result[!,:Haplotype].==k,:fstAverage] = sum(orga2[:,2] .* orga2[:,1]) / sum(orga2[:,1])
                                                                                                                  # Varianz of the Frequency
                                                                                                                  result[result[!,:Haplotype].==k,:fstVar] = std(orga2[:,2])
                                                                                                                  # Readdepth
                                                                                                                  result[result[!,:Haplotype].==k,:fstReaddepth] = sum(orga2[:,1])
                                                                                                                  # Marker counts
                                                                                                                  result[result[!,:Haplotype].==k,:fstMarkercount] = size(orga2,1)
                                                                                                                  # Zero counts
                                                                                                                  result[result[!,:Haplotype].==k,:fstZerocount] = sum(orga2[:,2].== 0)
                                                                                              else
                                                                                                                  result[result[!,:Haplotype].==k,:fstAverage] = NaN
                                                                                                                  # Varianz of the Frequency
                                                                                                                  result[result[!,:Haplotype].==k,:fstVar] = NaN
                                                                                                                  # Readdepth
                                                                                                                  result[result[!,:Haplotype].==k,:fstReaddepth] = NaN
                                                                                                                  # Marker counts
                                                                                                                  result[result[!,:Haplotype].==k,:fstMarkercount] = NaN
                                                                                                                  # Zero counts
                                                                                                                  result[result[!,:Haplotype].==k,:fstZerocount] = NaN
                                                                                              end # if


                                                                            # calculate the stats for the pair of data
                                                                            # if or if not to perform a transformation of the frequency based of the read counts
                                                                            if transform == false
                                                                                                jts = Float64.(orga2[!,:Allelfreq_Isr])  # log.(datT[!,:Readcount]) .* datT[!,:Allelfreq_Isr]
                                                                                                sta =  Float64.(konv2[!,:Allelfreq_Isr]) # log.(datS[!,:Readcount]) .* datS[!,:Allelfreq_Isr]
                                                                           elseif transform == true
                                                                                                 jts = log.(orga2[!,:Readcount]) .* orga2[!,:Allelfreq_Isr]
                                                                                                 sta = log.(konv2[!,:Readcount]) .* konv2[!,:Allelfreq_Isr]
                                                                           end
                                                                            # remove the NaN values
                                                                            sta = sta[.!(isnan.(sta))]
                                                                            jts = jts[.!(isnan.(jts))]

                                                                            # test if one of the sets is empty - if so no next itteration can be started
                                                                               if isempty(sta) || isempty(jts)
                                                                                                               result[result[!,:Haplotype].==k,:NegBin] = 5.0
                                                                                                               result[result[!,:Haplotype].==k,:Zeroinfl]   = 5.0
                                                                                                               result[result[!,:Haplotype].==k,:StatTest] = "MissData"
                                                                                                               continue
                                                                               end

                                                                           # if the sample size is smaller than 30, a nonparmatieric test will be applied - if bigger, a GLM model will be used
                                                                                if size(sta,1) < 15 && size(jts,1) < 15 || all(sta .== 0) || all(jts .== 0 )
                                                                                                      # caluculate the Allelfreq for the whole Haplotype new - use all reads aligned to this region for CALCULATION
                                                                                                      # 2.0 indicates a missing value - zero inflated column will be always 2.0 for samples with less than 30 SNPS
                                                                                                       #println("MannWhitneyUTest")
                                                                                                        result[result[!,:Haplotype].==k,:NegBin] = pvalue(MannWhitneyUTest(sta, jts))
                                                                                                        result[result[!,:Haplotype].==k,:Zeroinfl]   = 2.0
                                                                                                        result[result[!,:Haplotype].==k,:StatTest] = "MannWhitneyUTest"
                                                                             else

                                                                                                   t1 = vcat(sta, jts)
                                                                                                   #if isempty(sta) || isempty(jts) ; continue ; end
                                                                                                   t2 = repeat([1], size(sta,1))
                                                                                                   t3 = repeat([2], size(jts,1))
                                                                                                   n3 = DataFrame(hcat(t1, vcat(t2,t3)))
                                                                                                   n3[!,1] = Int64.(floor.(n3[!,1] .* 100))

                                                                                                   # when there are no 0 vaulues in it, a zero inflated model does not make sense - in this case, a NegativeBinomial Distribution is used only
                                                                                                   # if less than 10% of the samples are 0, the distribution will be set to a NegativeBinomial
                                                                                                   # with more than 75% zeros in the distribution, negbin zeroinfl will create NaN for the Zero infl model - change to poisson 0infl distribution
                                                                                                   if (sum(n3[!,:x1] .== 0) / size(n3,1)) < 0.15 || all(sta .!= 0) || all(jts .!= 0 )
                                                                                                                                     #println("Julia GLM")
                                                                                                                                     categorical!(n3, :x2)
                                                                                                                                     result[result[!,:Haplotype].==k,:NegBin] = coeftable(glm(@formula(x1 ~ x2), n3, NegativeBinomial(), LogLink())).cols[4][2]
                                                                                                                                     result[result[!,:Haplotype].==k,:Zeroinfl] = 3.0
                                                                                                                                     result[result[!,:Haplotype].==k,:StatTest] = "NegBino"
                                                                                                  elseif (sum(n3[!,:x1] .== 0) / size(n3,1)) < 0.75
                                                                                                                                     #println("ZeroInfl. Negbin")
                                                                                                                                     @rput n3
                                                                                                                                      # if the data has too many zeros compared to non-0, than a negative binomial model does not work - switch to poisson distribution
                                                                                                                                       R"""
                                                                                                                                       library(MASS)
                                                                                                                                       library(pscl)
                                                                                                                                       n3[,2] = as.factor(n3[,2])
                                                                                                                                       m1 =  summary(zeroinfl(x1~x2, n3, dist = "negbin"))
                                                                                                                                       negbino = m1$coefficients$count[,4]
                                                                                                                                       zeroinf = m1$coefficients$zero[,4]
                                                                                                                                       """
                                                                                                                                       @rget negbino  zeroinf
                                                                                                                                        result[result[!,:Haplotype].==k,:NegBin] = negbino[2]
                                                                                                                                      result[result[!,:Haplotype].==k,:Zeroinfl] = zeroinf[2]
                                                                                                                                      result[result[!,:Haplotype].==k,:StatTest] = "NbZi"
                                                                                                 else
                                                                                                                                       #println("ZeroInfl. Poisson")
                                                                                                                                try
                                                                                                                                       @rput n3
                                                                                                                                       R"""
                                                                                                                                       library(MASS)
                                                                                                                                       library(pscl)
                                                                                                                                       n3[,2] = as.factor(n3[,2])
                                                                                                                                       m1 =  summary(zeroinfl(x1~x2, n3, dist = "poisson"))
                                                                                                                                       negbino = m1$coefficients$count[,4]
                                                                                                                                       zeroinf = m1$coefficients$zero[,4]
                                                                                                                                       """
                                                                                                                                       @rget negbino  zeroinf
                                                                                                                                       #print(Crayon(foreground = :red, bold = false),"Poisson Distribution had to be used - too many Zero values! \n")
                                                                                                                                       #print(Crayon(foreground = :white, bold = false),"")
                                                                                                                                        result[result[!,:Haplotype].==k,:NegBin] = negbino[2]
                                                                                                                                        result[result[!,:Haplotype].==k,:Zeroinfl] = zeroinf[2]
                                                                                                                                        result[result[!,:Haplotype].==k,:StatTest] = "PsZi"

                                                                                                                                    catch
                                                                                                                                          print(Crayon(foreground = :red, bold = true),"\nHaplotype $k could not be calculated!\n")
                                                                                                                                          result[result[!,:Haplotype].==k,:NegBin] = 5.0
                                                                                                                                          result[result[!,:Haplotype].==k,:Zeroinfl] = 5.0
                                                                                                                                          result[result[!,:Haplotype].==k,:StatTest] = "none"
                                                                                                                                          continue
                                                                                                                                    end
                                                                                                end # NegBin or ZeroInfl
                                                                          end # of statistical test
                                                                end   # for hap

                                  return result


 end # of scomp

# merge the replicates together for the plots
function stacker(freq, genolist, s = "_Haplotypes")
                                                 println("merge the replicates for $s together")
                                                 set = collect(keys(freq))[occursin.(s, String.(collect(keys(freq))))]
                                                 gl = CSV.read(genolist, DataFrame; delim="\t")

                                                 ev = unique(gl[!,:Env])[2:size(unique(gl[!,:Env]),1)]
                                                 gen = unique(gl[!,:Generation])[2:size(unique(gl[!,:Generation]),1)]
                                                 # define the sets that should be matched togehter
                                                 for i in ev , j in gen
                                                                     qay = Symbol.(string.(gl[.&(gl[!,:Env].==i, gl[!,:Generation].==j),:Name], s))
                                                                     println(qay)
                                                                     # calculate the mean allelfreq and the sum of the readcount
                                                                     if size(qay,1) == 2
                                                                         # delete the haplotypes that have no replicate in the other set
                                                                         freq[qay[1]] = freq[qay[1]][ .!(isnothing.(indexin(freq[qay[1]][!,:Gene_ID], freq[qay[2]][!,:Gene_ID]))),:]
                                                                         freq[qay[2]] = freq[qay[2]][ .!(isnothing.(indexin(freq[qay[2]][!,:Gene_ID], freq[qay[1]][!,:Gene_ID]))),:]

                                                                                         a1 = (freq[qay[1]][!,:Allelfreq_Isr] .* freq[qay[1]][!,:Readcount] + freq[qay[2]][!,:Allelfreq_Isr] .* freq[qay[2]][!,:Readcount]) ./ (freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] ) # Allelfreq ISR
                                                                                         a2 = (freq[qay[1]][!,:Allelfreq_Golf] .* freq[qay[1]][!,:Readcount] + freq[qay[2]][!,:Allelfreq_Golf] .* freq[qay[2]][!,:Readcount]) ./ (freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] ) # Allelfreq Golf
                                                                                         a3 = maximum.(vcat.(freq[qay[1]][!,:SNPcount],freq[qay[2]][!,:SNPcount]))
                                                                                         a4 = freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount]
                                                                                         a5 = minimum.(vcat.(freq[qay[1]][!,:MinFreq],freq[qay[2]][!,:MinFreq] ))
                                                                                         a6 = maximum.(vcat.(freq[qay[1]][!,:MaxFreq],freq[qay[2]][!,:MaxFreq] ))
                                                                    elseif size(qay,1) == 3
                                                                                set = vcat(freq[qay[1]][!,:Gene_ID], freq[qay[2]][!,:Gene_ID],freq[qay[3]][!,:Gene_ID])
                                                                                set = countmap(set)
                                                                                # get all haplotypes that are present in all 3
                                                                                df = String[]
                                                                                for i in String.(keys(set))
                                                                                                              if set[i] == 3; push!(df, i); end
                                                                                end
                                                                                freq[qay[1]] = freq[qay[1]][ .!(isnothing.(indexin(freq[qay[1]][!,:Gene_ID], df))),:]
                                                                                freq[qay[2]] = freq[qay[2]][ .!(isnothing.(indexin(freq[qay[2]][!,:Gene_ID], df))),:]
                                                                                freq[qay[3]] = freq[qay[3]][ .!(isnothing.(indexin(freq[qay[3]][!,:Gene_ID], df))),:]



                                                                                      a1 = (freq[qay[1]][!,:Allelfreq_Isr] .* freq[qay[1]][!,:Readcount] + freq[qay[2]][!,:Allelfreq_Isr] .* freq[qay[2]][!,:Readcount] + freq[qay[3]][!,:Allelfreq_Isr] .* freq[qay[3]][!,:Readcount]) ./ (freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] .+ freq[qay[3]][!,:Readcount]) # Allelfreq ISR
                                                                                      a2 = (freq[qay[1]][!,:Allelfreq_Golf] .* freq[qay[1]][!,:Readcount] + freq[qay[2]][!,:Allelfreq_Golf] .* freq[qay[2]][!,:Readcount] + freq[qay[3]][!,:Allelfreq_Golf] .* freq[qay[3]][!,:Readcount]) ./ (freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] .+ freq[qay[3]][!,:Readcount]) # Allelfreq Golf
                                                                                      a3 = maximum.(vcat.(freq[qay[1]][!,:SNPcount],freq[qay[2]][!,:SNPcount], freq[qay[3]][!,:SNPcount]))
                                                                                      a4 = freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] .+ freq[qay[3]][!,:Readcount]
                                                                                      a5 = minimum.(vcat.(freq[qay[1]][!,:MinFreq],freq[qay[2]][!,:MinFreq], freq[qay[3]][!,:MinFreq] ))
                                                                                      a6 = maximum.(vcat.(freq[qay[1]][!,:MaxFreq],freq[qay[2]][!,:MaxFreq], freq[qay[3]][!,:MaxFreq] ))
                                                                     else
                                                                          println("No replicates to merge")
                                                                          continue

                                                                     end

                                                                     a = DataFrame(hcat(freq[qay[1]][!,:Gene_ID], a1,a2,a3,a4,a5,a6), :auto)
                                                                     rename!(a, names(freq[qay[1]]))
                                                                     # push to the Dictionary
                                                                     saver = Symbol(string("F",j,"E",i,s))
                                                                     freq[saver] = a
                                                 end
  end
# stacker(freq, genolist)

# add postion information to the Haplotype ,  col = "Gene_ID" |  col = ["Gene_ID", "Marker"] | col = ["Gene_ID", "GO-IDs"] | col = ["Gene_ID", "PFAM-IDs"] | col = ["Gene_ID", "InterPro-IDs"]
function Position(file, freqfile, col = "Gene_ID")
                                  q3 = unique(freqfile[:Location][!,3:12])
                                  if typeof(col) == String
                                                                        file = innerjoin(file, q3, on=Symbol(col))
                                  elseif typeof(col) == Array{String,1}
                                                                        file = innerjoin(file, q3, on=[Symbol(col[1]) => Symbol(col[2])], makeunique=true)
                                  else
                                  println("""The columns to merge with are not correctly specified, use: \n col = "Gene_ID" |  col = ["Gene_ID", "Marker"] | col = ["Gene_ID", "GO-IDs"] | col = ["Gene_ID", "PFAM-IDs"] | col = ["Gene_ID", "InterPro-IDs"]""")
                                  end

                                  return file

 end
# position( freq[:F23E2_Haplotypes], freq)

#= the following function performes plotting for the Single Genes selected or al list of genes selected =#
function genplot(dict, haplotype,style = "bee", system ="e2", type="Gene_ID", folder = "/home/michael/Schreibtisch/plots/", plotsize=[5inch,5inch], dpi=300)

                                   if typeof(haplotype) == String || length(haplotype) == 1  # for the haplotype, the information is extracted from the freqfile all SNP table to get variation

                                                                         if length(haplotype) == 1; haplotype = haplotype[1];end # this is necessary if a array with only one element is added

                                                                         println( "Print single gene plot")
                                                                         a = string.(collect(keys(dict)))
                                                                         lookinto = Symbol.(a[occursin.(system,a)])
                                                                         # a new dataframe has to be created where all information is stored
                                                                         xplot =   DataFrame(Array{Union{Nothing, Any}}(nothing, 0, size(dict[lookinto[1]],2)+1))
                                                                         rename!(xplot, Symbol.(vcat(string.(names(dict[lookinto[1]])),"gen")))

                                                                         # extract the Marker of this haplotype
                                                                         for i in 1:length(lookinto)

                                                                           if any(occursin.("eLoc", String.(collect(keys(freq)))))
                                                                                              fx = dict[lookinto[i]][dict[:eLoc][!,Symbol(type)].==haplotype,:] # having not to create a double lineup of the Location function makes it much faster
                                                                           elseif size(dict[lookinto[i]],1) == size(dict[:Location],1)
                                                                                              fx = dict[lookinto[i]][dict[:Location][!,Symbol(type)].==haplotype,:]
                                                                           elseif size(dict[lookinto[i]],1) == size(dict[:Location],1) * 2
                                                                                              fx = dict[lookinto[i]][sort(repeat(dict[:Location],2))[!,Symbol(type)].==haplotype,:]
                                                                           elseif size(dict[lookinto[i]],1) == size(dict[:Location],1) * 3
                                                                                              fx = dict[lookinto[i]][sort(repeat(dict[:Location],3))[!,Symbol(type)].==haplotype,:]
                                                                           elseif size(dict[lookinto[i]],1) == size(dict[:Location],1) * 4
                                                                                              fx = dict[lookinto[i]][sort(repeat(dict[:Location],4))[!,Symbol(type)].==haplotype,:]
                                                                           else
                                                                               print("too many replicates, change the script!")
                                                                           end

                                                                                   fx[!,:gen] .= string(lookinto[i])
                                                                                   append!(xplot, fx)
                                                                         end
                                                                         # turn the Lib info into a single value
                                                                         for i in 1:size(xplot, 1); xplot[i,:gen] = SubString(xplot[i,:gen], 2, findfirst.("e", xplot[i,:gen])[1]-1 ) ; end
                                                                         xplot[!,:gen] = parse.(Int64, string.(xplot[!,:gen]))
                                                                         sort!(xplot, :gen)
                                                                         xplot[!,:gen] = string.(xplot[!,:gen]) # this is important for the plot, so that no contimous scale is coming up
                                                                         xplot[!,:Allelfreq_Isr] = xplot[!,:Allelfreq_Isr] .* 100
                                                                         # remove NaN values from Allelfreq_Isr
                                                                         xplot = xplot[.!(isnan.(xplot[!,:Allelfreq_Isr])),:]

                                                                         # the xplot columns that should be plotted has to be float values or Int
                                                                         for j in 1:size(xplot,2)
                                                                                                  # Int64
                                                                                                  if typeof(xplot[1,j]) == Int64
                                                                                                            xplot[!,j] = Int64.(xplot[!,j])
                                                                                                  # Float64
                                                                                                elseif typeof(xplot[1,j]) == Float64
                                                                                                        xplot[!,j] = Float64.(xplot[!,j])
                                                                                                  # String
                                                                                                  else
                                                                                                        xplot[!,j] = String.(xplot[!,j])
                                                                                                  end
                                                                          end

                                                                         if style == "bee"
                                                                         p = plot(xplot, layer(x=:gen, y=:Allelfreq_Isr, color=:Readcount, Geom.beeswarm),
                                                                          layer(x= :gen, y= :Allelfreq_Isr, Geom.boxplot),
                                                                          Guide.xlabel("Generation"), Guide.ylabel("Allelfreq [%]"), Guide.colorkey(title="Read count"), Guide.title(haplotype),
                                                                          Theme(background_color=colorant"white", default_color=colorant"grey"), Coord.cartesian(ymin=0, ymax=100))
                                                                        else
                                                                         p = plot(xplot,x= :gen, y= :Allelfreq_Isr, Geom.boxplot,
                                                                           Guide.xlabel("Generation"), Guide.ylabel("Allelfreq [%]"), Guide.title(haplotype),
                                                                           Theme(background_color=colorant"white", default_color=colorant"grey"), Coord.cartesian(ymin=0, ymax=100),
                                                                           Guide.annotation(compose(context(),  Compose.text(1, 1, "N = $(size(xplot,1))", hleft, vbottom))))
                                                                        end


                                                                          if length(unique(xplot[!,:gen])) == 2
                                                                                                             draw(PNG(string(folder, "/", haplotype,"_", system, ".png"), 5inch, 7.5inch, dpi =  dpi),p)
                                                                          elseif length(unique(xplot[!,:gen])) == 3
                                                                                                             draw(PNG(string(folder, "/", haplotype,"_", system, ".png"), 6.25inch, 7.5inch, dpi =  dpi),p)
                                                                          elseif length(unique(xplot[!,:gen])) == 4
                                                                                                             draw(PNG(string(folder, "/", haplotype,"_", system, ".png"), 7.5inch, 7.5inch, dpi =  dpi),p)
                                                                          elseif length(unique(xplot[!,:gen])) == 5
                                                                                                             draw(PNG(string(folder, "/", haplotype,"_", system, ".png"), 8.75inch, 7.5inch, dpi =  dpi),p)
                                                                          elseif length(unique(xplot[!,:gen])) == 6
                                                                                                             draw(PNG(string(folder, "/", haplotype,"_", system, ".png"), 10inch, 7.5inch, dpi =  dpi),p)
                                                                          elseif length(unique(xplot[!,:gen])) == 7
                                                                                                             draw(PNG(string(folder, "/", haplotype,"_", system, ".png"), 11.25inch, 7.5inch, dpi =  dpi),p)
                                                                          elseif length(unique(xplot[!,:gen])) == 8
                                                                                                             draw(PNG(string(folder, "/", haplotype, "_", system,".png"), 13.5inch, 7.5inch, dpi =  dpi),p)
                                                                          elseif length(unique(xplot[!,:gen])) == 9
                                                                                                             draw(PNG(string(folder, "/", haplotype,"_", system, ".png"), 14.75inch, 7.5inch, dpi =  dpi),p)
                                                                          elseif length(unique(xplot[!,:gen])) == 10
                                                                                                             draw(PNG(string(folder, "/", haplotype,"_", system, ".png"), 16inch, 7.5inch, dpi =  dpi),p)
                                                                          else
                                                                                                             draw(PNG(string(folder, "/", haplotype,"_", system, ".png"), 16inch, 7.5inch, dpi =  dpi),p)
                                                                        end

                                                                        return p

                                             else
                                  # the plots will be stored in an Dict, from which the overall plot will be made and all the single plots
                                                                         println( "Print multi gene plot")
                                                                         plots = Dict()

                                                                         latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=20pt,
                                                                                             minor_label_font="CMU Serif", minor_label_font_size=15pt,
                                                                                             key_title_font="CMU Serif", key_title_font_size=25pt,
                                                                                             key_label_font="CMU Serif", key_label_font_size=20pt)
                                                                         Gadfly.push_theme(latex_fonts)

                                                                         for o in 1:length(haplotype) # create a plot for each gene written in haplotype
                                                                             println(o)
                                                                             a = string.(collect(keys(dict)))
                                                                             lookinto = Symbol.(a[occursin.(system,a)])
                                                                             # a new dataframe has to be created where all information is stored
                                                                             xplot =   DataFrame(Array{Union{Nothing, Any}}(nothing, 0, size(dict[lookinto[1]],2)+1))
                                                                             rename!(xplot, Symbol.(vcat(string.(names(dict[lookinto[1]])),"gen")))

                                                                             # extract the Marker of this haplotype
                                                                             for i in 1:length(lookinto)

                                                                               if any(occursin.("eLoc", String.(collect(keys(freq)))))
                                                                                                  fx = dict[lookinto[i]][dict[:eLoc][!,Symbol(type)].==haplotype[o],:] # having not to create a double lineup of the Location function makes it much faster
                                                                               elseif size(dict[lookinto[i]],1) == size(dict[:Location],1)
                                                                                                  fx = dict[lookinto[i]][dict[:Location][!,Symbol(type)].==haplotype[o],:]
                                                                               elseif size(dict[lookinto[i]],1) == size(dict[:Location],1) * 2
                                                                                                  fx = dict[lookinto[i]][sort(repeat(dict[:Location],2))[!,Symbol(type)].==haplotype[o],:]
                                                                               elseif size(dict[lookinto[i]],1) == size(dict[:Location],1) * 3
                                                                                                  fx = dict[lookinto[i]][sort(repeat(dict[:Location],3))[!,Symbol(type)].==haplotype[o],:]
                                                                               elseif size(dict[lookinto[i]],1) == size(dict[:Location],1) * 4
                                                                                                  fx = dict[lookinto[i]][sort(repeat(dict[:Location],4))[!,Symbol(type)].==haplotype[o],:]
                                                                               else
                                                                                   print("too many replicates, change the script!")
                                                                               end

                                                                                       fx[!,:gen] .= string(lookinto[i])
                                                                                       append!(xplot, fx)
                                                                             end
                                                                             # turn the Lib info into a single value
                                                                             for i in 1:size(xplot, 1); xplot[!,:gen][i] = SubString(xplot[!,:gen][i], 2, findfirst.("e", xplot[!,:gen])[i][1]-1 ) ; end
                                                                             xplot[!,:gen] = parse.(Int64, string.(xplot[!,:gen]))
                                                                             sort!(xplot, :gen)
                                                                             xplot[!,:gen] = string.(xplot[!,:gen])
                                                                             xplot[!,:Allelfreq_Isr] = xplot[!,:Allelfreq_Isr] .* 100
                                                                             # remove NaN values from Allelfreq_Isr
                                                                             xplot = xplot[.!(isnan.(xplot[!,:Allelfreq_Isr])),:]

                                                                             # the xplot columns that should be plotted has to be float values or Int
                                                                             for j in 1:size(xplot,2)
                                                                                                      # Int64
                                                                                                      if typeof(xplot[1,j]) == Int64
                                                                                                                xplot[!,j] = Int64.(xplot[!,j])
                                                                                                      # Float64
                                                                                                    elseif typeof(xplot[1,j]) == Float64
                                                                                                            xplot[!,j] = Float64.(xplot[!,j])
                                                                                                      # String
                                                                                                      else
                                                                                                            xplot[!,j] = String.(xplot[!,j])
                                                                                                      end
                                                                              end

                                                                             if style == "bee"
                                                                             plots[Symbol(string("p",o))] = plot(xplot, layer(x=:gen, y=:Allelfreq_Isr, color=:Readcount, Geom.beeswarm),
                                                                              layer(x= :gen, y= :Allelfreq_Isr, Geom.boxplot()),
                                                                              Guide.xlabel("Generation"), Guide.ylabel("Allelfreq [%]"), Guide.colorkey(title="Read count"), Guide.title(haplotype[o]),
                                                                              Theme(background_color=colorant"white", default_color=colorant"grey"))
                                                                             else
                                                                              plots[Symbol(string("p",o))] = plot(xplot,x= :gen, y= :Allelfreq_Isr, Geom.boxplot,
                                                                                Guide.xlabel("Generation"), Guide.ylabel("Allelfreq [%]"), Guide.title(haplotype[o]),
                                                                                Theme(background_color=colorant"white", default_color=colorant"grey"), Coord.cartesian(ymin=0, ymax=100),
                                                                                Guide.annotation(compose(context(),  Compose.text(1, 1, "N = $(size(xplot,1))", hleft, vbottom))))
                                                                             end

                                                                              #draw(PNG(string(folder, "/", haplotype[o],".png"), plotsize[1], plotsize[2]),plots[Symbol(string("p",o))])

                                                                         end

                                                                         # concat the plots to one plots
                                                                         if length(haplotype) == 2
                                                                                 #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*2, plotsize[2], dpi=dpi),hstack(plots[:p1],plots[:p2]))
                                                                                 p = hstack(plots[:p1],plots[:p2])
                                                                         elseif length(haplotype) == 3
                                                                             #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*3, plotsize[2], dpi=dpi),hstack(plots[:p1],plots[:p2], plots[:p3]))
                                                                                  p = hstack(plots[:p1],plots[:p2], plots[:p3])
                                                                         elseif length(haplotype) == 4
                                                                             #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*2, plotsize[2]*2, dpi=dpi),hstack(vstack(plots[:p1],plots[:p2]), vstack(plots[:p3], plots[:p4])))
                                                                             p =  hstack(vstack(plots[:p1],plots[:p2]), vstack(plots[:p3], plots[:p4]))
                                                                         elseif length(haplotype) == 5
                                                                             #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*3, plotsize[2]*2, dpi=dpi),vstack(hstack(plots[:p1],plots[:p2], plots[:p3]),hstack(plots[:p4], plots[:p5])))
                                                                             p = vstack(hstack(plots[:p1],plots[:p2], plots[:p3]),hstack(plots[:p4], plots[:p5]))
                                                                         elseif length(haplotype) == 6
                                                                             #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*3, plotsize[2]*2, dpi=dpi),vstack(hstack(plots[:p1],plots[:p2], plots[:p3]), hstack(plots[:p4], plots[:p5], plots[:p6])))
                                                                             p = vstack(hstack(plots[:p1],plots[:p2], plots[:p3]), hstack(plots[:p4], plots[:p5], plots[:p6]))
                                                                         elseif length(haplotype) == 7
                                                                             #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*4, plotsize[2]*2, dpi=dpi),vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4]), hstack(plots[:p5], plots[:p6], plots[:p7])))
                                                                             p = vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4]), hstack(plots[:p5], plots[:p6], plots[:p7]))
                                                                         elseif length(haplotype) == 8
                                                                             #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*4, plotsize[2]*2, dpi=dpi),vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4]), hstack(plots[:p5], plots[:p6], plots[:p7], plots[:p8])))
                                                                             p = vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4]), hstack(plots[:p5], plots[:p6], plots[:p7], plots[:p8]))
                                                                         elseif length(haplotype) == 9
                                                                             #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*5, plotsize[2]*2, dpi=dpi),vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4]), hstack(plots[:p5], plots[:p6], plots[:p7], plots[:p8], plots[:p9])))
                                                                             p = vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4]), hstack(plots[:p5], plots[:p6], plots[:p7], plots[:p8], plots[:p9]))
                                                                         elseif length(haplotype) == 10
                                                                             #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*5, plotsize[2]*2, dpi=dpi),vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4], plots[:p5]), hstack(plots[:p6], plots[:p7], plots[:p8], plots[:p9], plots[:p10])))
                                                                             p = vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4], plots[:p5]), hstack(plots[:p6], plots[:p7], plots[:p8], plots[:p9], plots[:p10]))
                                                                         elseif length(haplotype) == 11
                                                                             #draw(PNG(string(folder,"/sumplot.png"), plotsize[1]*6, plotsize[2]*2, dpi=dpi),vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4],plots[:p5], plots[:p6]), hstack(plots[:p7], plots[:p8], plots[:p9], plots[:p10], plots[:p11])))
                                                                             p = vstack(hstack(plots[:p1],plots[:p2], plots[:p3], plots[:p4],plots[:p5], plots[:p6]), hstack(plots[:p7], plots[:p8], plots[:p9], plots[:p10], plots[:p11]))
                                                                         else
                                                                             println("too many Genes for a single plot selected, please choose 7 or less")
                                                                         end
                                                                         return p
                                                             end
                                  end # of function

# use genplot for both sets of environments
function envgenplot(dict, haplotype, style = "bee", type="Gene_ID", folder = "/home/michael/Schreibtisch/plots/", plotsize=[5inch,7.5inch], dpi = 300)
                    println("Plot E1")
                    p1 = genplot(dict, haplotype, style, "e1", type, folder , plotsize, dpi)
                    println("Plot E2")
                    p2 = genplot(dict, haplotype, style, "e2", type, folder , plotsize, dpi)
                    println("Plot both together")
                    p3 = vstack(p1,p2)
                     draw(PNG(string(folder, "/", haplotype,"_both_env", ".png"),plotsize[1]*1.5, plotsize[2]),p3)

  end
# envgenplot(freq, "HORVU1Hr1G009100","box", "Gene_ID", "/media/michael/Arbeit/Allelfrequenz/Gerste/WGS/ref_snps_duplucates_removed/Nutrients/Potassium/", [5inch,8inch])

## a function that is capabile to print the distribution of the snp over the whole chromosom.
# genetic map
function coverplotG(file, set, col=:Readcount,  what=:Allelfreq_Isr ,save=nothing)

                                                                                    if (findfirst("Haplotype", set)===nothing) == false

                                                                                    # check if this has already been done
                                                                                    if !(any(names(file[Symbol(set)]).== :Chr))
                                                                                    file[Symbol(set)] = join(file[Symbol(set)], unique(file[:Location],6)[!,[6,1,7]], on=(:Gene_ID,:Marker))
                                                                                    file[Symbol(set)][!,:Alleles] .= "1/1"
                                                                                    sort!(file[Symbol(set)], (:Chr, :Pos_genetic))
                                                                                    file[Symbol(set)][!,:Pos_genetic] = Float64.(file[Symbol(set)][!,:Pos_genetic])
                                                                                    end
                                                                                    # let all outliers disapear
                                                                                    MEAN = Int64(floor(mean(file[Symbol(set)][!,col])))
                                                                                    STD = Int64(floor(std(file[Symbol(set)][!,col])))
                                                                                    MEDIAN = Int64(floor(median(file[Symbol(set)][!,col])))
                                                                                    AF = round(mean(file[Symbol(set)][!,:Allelfreq_Isr])*100, digits=2)
                                                                                    AF2 = round(sum(file[Symbol(set)][!,:Allelfreq_Isr] .* file[Symbol(set)][!,:Readcount]) / sum(file[Symbol(set)][!,:Readcount])*100,digits=2)

                                                                                    for i in 1:size(file[Symbol(set)],1)
                                                                                                        # if the value is biiger than the std , set it to std + 10
                                                                                                        if file[Symbol(set)][i,col] > STD;  file[Symbol(set)][i,col] = STD+10; end
                                                                                    end
                                                                                    plotx = plot(file[Symbol(set)], x= :Pos_genetic, y=what, color = col, xgroup=:Chr, Geom.subplot_grid(Geom.point), Theme(background_color=colorant"white"), Guide.xlabel("Genetic Position"), Guide.ylabel("Allele Frequency ISR 42-8"),
                                                                                    Guide.annotation(compose(context(), Compose.text(0.1, 0.1, "$(string(col)) - Mean: $MEAN, Median = $MEDIAN, Stdev = $STD \n Allele Frequency - Mean = $AF%, adjusted Mean = $AF2%"))),  Guide.title(set))

                                                                                    if save != nothing
                                                                                    draw(PNG("$save/location$set.png", 10inch, 20inch, dpi=300), plotx)
                                                                                    end

                                                                                    end # haplotyp checker

                                                                                    return plotx

 end

 hap = Dict()
 for i in collect(keys(freq))[occursin.("E", String.(collect(keys(freq))))][[1,3,5,6,7,8,9,10,11,12]]
                                                                           println(i)
                                                                           hap[i] = Position(freq[i], freq)
                                                                         end

pk = vstack(coverplotG(freq, "F3E1_Haplotypes"), coverplotG(freq, "F12E1_Haplotypes"), coverplotG(freq, "F16E1_Haplotypes"), coverplotG(freq, "F22E1_Haplotypes"), coverplotG(freq, "F23E1_Haplotypes"));
po =  vstack(coverplotG(freq, "F3E2_Haplotypes"), coverplotG(freq, "F12E2_Haplotypes"), coverplotG(freq, "F16E2_Haplotypes"), coverplotG(freq, "F22E2_Haplotypes"), coverplotG(freq, "F23E2_Haplotypes"));
pp = hstack(pk,po);
draw(PNG("AllelFreqPlotMarker.png", 32inch, 18inch, dpi=300), pp)



# physical map
function coverplot(file, set, col=:x1,  what=:Allelfreq_Isr ,save=nothing)
                                       # pre step - if one whats to go for plotting the haplotypes, the chr, pos and Alleles column have to be added to the files to plot
                                       if (findfirst("Haplotype", set)===nothing) == false
                                                                                                            latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=22pt,
                                                                                                            minor_label_font="CMU Serif", minor_label_font_size=19pt,
                                                                                                            key_title_font="CMU Serif", key_title_font_size=17pt,
                                                                                                            key_label_font="CMU Serif", key_label_font_size=15pt,
                                                                                                            background_color = "white")
                                                                                                            Gadfly.push_theme(latex_fonts)

                                                                                         if any(names(file[Symbol(set)]).== :x1)
                                                                                         file[Symbol(set)] = join(file[Symbol(set)], unique(file[:Location],6)[[1,2,6]], on=:Gene_ID)
                                                                                         file[Symbol(set)][!,:Alleles] .= "1/1"
                                                                                         sort!(file[Symbol(set)], [:x1, :start])
                                                                                         file[Symbol(set)][!,:Pos] = Float64.(file[Symbol(set)][!,:Pos])
                                                                                         end
                                                                                         # let all outliers disapear
                                                                                         #MEAN = Int64(floor(mean(file[Symbol(set)][!,col])))
                                                                                         #STD = Int64(floor(std(file[Symbol(set)][!,col])))
                                                                                         #MEDIAN = Int64(floor(median(file[Symbol(set)][!,col])))
                                                                                         #AF = round(mean(file[Symbol(set)][!,:Allelfreq_Isr])*100, digits=2)
                                                                                         #AF2 = round(sum(file[Symbol(set)][!,:Allelfreq_Isr] .* file[Symbol(set)][!,:Readcount]) / sum(file[Symbol(set)][!,:Readcount])*100,digits=2)
                                                                                         try
                                                                                           file[Symbol(set)][!,:x1] .= SubString.(file[Symbol(set)][!,:x1],4,5)
                                                                                         catch
                                                                                           nothing
                                                                                         end


                                                                                         # px = plot(file[Symbol(set)], x= :Pos, y=what, color = col, xgroup=:Chr, Geom.subplot_grid(Geom.point), Theme(background_color=colorant"white"), Guide.xlabel("Position"), Guide.ylabel("Allele Frequency ISR 42-8"),
                                                                                         # Guide.annotation(compose(context(), Compose.text(0.1, 0.1, "$(string(col)) - Mean: $MEAN, Median = $MEDIAN, Stdev = $STD \n Allele Frequency - Mean = $AF%, adjusted Mean = $AF2%"))),  Guide.title(set), Guide.xticks(ticks=nothing ))

                                                                                          px = plot(file[Symbol(set)], x= :start, y=what, color = col, ygroup=:x1, yintercept = [0.125], Geom.subplot_grid(Geom.line, Guide.xticks(label=false), Geom.hline(style=:dot)), Theme(background_color=colorant"white"), Guide.xlabel("Physical Position"),
                                                                                          Guide.ylabel("Allele Frequency Parent 2"), Guide.colorkey(pos=[-10,-10]))

                                                                                         if save != nothing
                                                                                         draw(PNG("$save/location$set.png", 10inch, 20inch, dpi=300), px)
                                                                                         end

                                                                                         end # haplotyp checker

                                                                                         return px

      end
#plotx = vstack(coverplotG(freq, "F3E1_Haplotypes"), coverplotG(freq, "F12E1_Haplotypes"), coverplotG(freq, "F12E2_Haplotypes"),coverplotG(freq, "F16E1_Haplotypes"),coverplotG(freq, "F16E2_Haplotypes"),coverplotG(freq, "F22E1_Haplotypes"),coverplotG(freq, "F22E2_Haplotypes"),coverplotG(freq, "F23E1_Haplotypes"),coverplotG(freq, "F23E2_Haplotypes"))
hap = Dict()
for i in collect(keys(freq))[occursin.("E", String.(collect(keys(freq))))][[1,3,5,6,7,8,9,10,11,12]]
                                                                          println(i)
                                                                          hap[i] = Position(freq[i], freq)
                                                                        end
hap[:Location] = freq[:Location]
pk = vstack(coverplot(hap, "F3E1_Haplotypes"), coverplot(hap, "F12E1_Haplotypes"), coverplot(hap, "F16E1_Haplotypes"), coverplot(hap, "F22E1_Haplotypes"), coverplot(hap, "F23E1_Haplotypes"))
po =  vstack(coverplot(hap, "F3E2_Haplotypes"), coverplot(hap, "F12E2_Haplotypes"), coverplot(hap, "F16E2_Haplotypes"), coverplot(hap, "F22E2_Haplotypes"), coverplot(hap, "F23E2_Haplotypes"))
pp = hstack(pk,po);
draw(PNG("F22kon.png", 10inch, 10inch, dpi=300), p2)


function changeplot(set, enter, exit,map="genetic",chrom="all", hair=true)

                                                # only Haplotypes allowed!!
                                                                            map == "genetic" ? pos = :Pos_genetic : pos = :Pos
                                                                            # join the two sets
                                                                            test = join(set[Symbol(enter)], set[Symbol(exit)], on=:Gene_ID, makeunique=true)

                                                                            # calculate a new column with the information on the differnce compared to the first test generation
                                                                            test[!,:diff] = test[!,:Allelfreq_Isr_1] .- test[!,:Allelfreq_Isr]
                                                                            # remove some columns
                                                                            deletecols!(test, [:Allelfreq_Isr, :Allelfreq_Golf, :Allelfreq_Isr_1, :Allelfreq_Golf_1])
                                                                            # join the haplotypes with the location
                                                                            map == "genetic" ?  testp = join(test, unique(set[:Location],6)[:,(1:9)], on=(:Gene_ID,:Marker), makeunique=true) : testp = join(test, unique(set[:Location],6)[:,(1:9)], on=:Gene_ID)
                                                                            # calculate the difference of the Readcount
                                                                            testp[!,:snpDiff] = testp[!,:Readcount_1] ./ testp[!,:Readcount]
                                                                            # create a column with information if its up, down or stable
                                                                            testp[!,:change] .= "equal"
                                                                            testp[testp[!,:diff].< -0.05,:change] .= "decreased"
                                                                            testp[testp[!,:diff].> 0.05 ,:change] .= "increased"

                                                                            #create a dict to loop in  -add color for the snp count for the haplotypes
                                                                            chromosoms = Dict()
                                                                            for i in 1:7
                                                                                chromosoms[Symbol("chr$(i)H")] = testp[testp[!,:Chr].=="chr$(i)H",:]
                                                                                sort!(chromosoms[Symbol("chr$(i)H")], :change)

                                                                                if hair == true
                                                                                    chromosoms[Symbol("PLOTchr$(i)H")] = plot(chromosoms[Symbol("chr$(i)H")], x=pos, y=:diff, color=:change, Geom.hair, Geom.point, Theme(background_color=colorant"white", grid_color=colorant"black"),
                                                                                    Guide.xlabel("Position on $(chromosoms[Symbol("chr$(i)H")][!,:Chr][1])"), Guide.ylabel("ISR freq in Pool"), Scale.color_discrete_manual(colorant"red", colorant"slategray4", colorant"darkgreen"))
                                                                                else
                                                                                    chromosoms[Symbol("PLOTchr$(i)H")] = plot(chromosoms[Symbol("chr$(i)H")], x=pos, y=:diff, color=:change, Geom.point, Theme(background_color=colorant"white", grid_color=colorant"black"),
                                                                                    Guide.xlabel("Position on $(chromosoms[Symbol("chr$(i)H")][!,:Chr][1])"), Guide.ylabel("ISR freq in Pool"), Scale.color_discrete_manual(colorant"red", colorant"slategray4", colorant"darkgreen"))

                                                                                end
                                                                            end
                                                                            title2 = compose(Compose.context(0,0,1w,0.25inch), Compose.text(0.5,2.0, "Variation of $exit to $enter", hcenter, vbottom))
                                                                            plotx = nothing
                                                                            if chrom == "all"
                                                                                                   plotx = vstack(title2, chromosoms[Symbol("PLOTchr1H")],chromosoms[Symbol("PLOTchr2H")],chromosoms[Symbol("PLOTchr3H")],chromosoms[Symbol("PLOTchr4H")],chromosoms[Symbol("PLOTchr5H")],chromosoms[Symbol("PLOTchr6H")],chromosoms[Symbol("PLOTchr7H")])
                                                                            elseif size(chrom,1) == 2
                                                                                                   plotx = vstack(title2, chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")])
                                                                            elseif size(chrom,1) == 1
                                                                                                   plotx = vstack(title2, chromosoms[Symbol("PLOT$(chrom[1])")])
                                                                            elseif size(chrom,1) == 3
                                                                                                   plotx = vstack(title2, chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")], chromosoms[Symbol("PLOT$(chrom[3])")])
                                                                            elseif size(chrom,1) == 4
                                                                                                   plotx = vstack(title2, chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")], chromosoms[Symbol("PLOT$(chrom[3])")], chromosoms[Symbol("PLOT$(chrom[4])")])
                                                                            elseif size(chrom,1) == 5
                                                                                                   plotx = vstack(title2, chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")], chromosoms[Symbol("PLOT$(chrom[3])")], chromosoms[Symbol("PLOT$(chrom[4])")],chromosoms[Symbol("PLOT$(chrom[5])")])
                                                                            elseif size(chrom,1) == 6
                                                                                                   plotx = vstack(title2, chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")], chromosoms[Symbol("PLOT$(chrom[3])")], chromosoms[Symbol("PLOT$(chrom[4])")], chromosoms[Symbol("PLOT$(chrom[5])")], chromosoms[Symbol("PLOT$(chrom[6])")])
                                                                            else
                                                                                                   println("""The "chr" setting was used incorrectly, all Chromosoms are plotted""")
                                                                                                   plotx = hstack(chromosoms[Symbol("PLOTchr1H")],chromosoms[Symbol("PLOTchr2H")],chromosoms[Symbol("PLOTchr3H")],chromosoms[Symbol("PLOTchr4H")],chromosoms[Symbol("PLOTchr5H")],chromosoms[Symbol("PLOTchr6H")],chromosoms[Symbol("PLOTchr7H")])
                                                                            end

                                                                            return plotx
 end

function changeplotT(set, enter, exit,map="genetic",chrom="all", hair=true)

                                                 # only Haplotypes allowed!!
                                                                             map == "genetic" ? pos = :Pos_genetic : pos = :Pos
                                                                             # join the two sets
                                                                             test = join(set[Symbol(enter)], set[Symbol(exit)], on=:Gene_ID, makeunique=true)

                                                                             # calculate a new column with the information on the differnce compared to the first test generation
                                                                             test[!,:diff] = test[!,:Allelfreq_Isr_1] .- test[!,:Allelfreq_Isr]
                                                                             # remove some columns
                                                                             deletecols!(test, [:Allelfreq_Isr, :Allelfreq_Golf, :Allelfreq_Isr_1, :Allelfreq_Golf_1])
                                                                             # join the haplotypes with the location
                                                                             map == "genetic" ?  testp = join(test, unique(set[:Location],6)[:,(1:9)], on=(:Gene_ID,:Marker), makeunique=true) : testp = join(test, unique(set[:Location],6)[:,(1:9)], on=:Gene_ID)
                                                                             # calculate the difference of the Readcount
                                                                             testp[!,:snpDiff] = testp[!,:Readcount_1] ./ testp[!,:Readcount]
                                                                             # create a column with information if its up, down or stable
                                                                             testp[!,:change] .= "equal"
                                                                             testp[testp[!,:diff].< -0.05,:change] .= "decreased"
                                                                             testp[testp[!,:diff].> 0.05 ,:change] .= "increased"

                                                                             #create a dict to loop in  -add color for the snp count for the haplotypes
                                                                             chromosoms = Dict()
                                                                             for i in 1:7
                                                                                 chromosoms[Symbol("chr$(i)H")] = testp[testp[!,:Chr].=="chr$(i)H",:]
                                                                                 sort!(chromosoms[Symbol("chr$(i)H")], :change)

                                                                                 if hair == true
                                                                                     chromosoms[Symbol("PLOTchr$(i)H")] = plot(chromosoms[Symbol("chr$(i)H")], x=pos, y=:diff, color=:change, Geom.hair, Geom.point, Theme(background_color=colorant"white", grid_color=colorant"black"),
                                                                                     Guide.xlabel("Position on $(chromosoms[Symbol("chr$(i)H")][!,:Chr][1])"), Guide.ylabel("ISR freq in Pool"), Scale.color_discrete_manual(colorant"red", colorant"slategray4", colorant"darkgreen"))
                                                                                 else
                                                                                     chromosoms[Symbol("PLOTchr$(i)H")] = plot(chromosoms[Symbol("chr$(i)H")], x=pos, y=:diff, color=:change, Geom.point, Theme(background_color=colorant"white", grid_color=colorant"black"),
                                                                                     Guide.xlabel("Position on $(chromosoms[Symbol("chr$(i)H")][!,:Chr][1])"), Guide.ylabel("ISR freq in Pool"), Scale.color_discrete_manual(colorant"red", colorant"slategray4", colorant"darkgreen"))

                                                                                 end
                                                                             end
                                                                             #title2 = compose(Compose.context(0,0,1w,0.25inch), Compose.text(0.5,2.0, "Variation of $exit to $enter", hcenter, vbottom))
                                                                             plotx = nothing
                                                                           if chrom == "all"
                                                                                                  plotx = hstack(chromosoms[Symbol("PLOTchr1H")],chromosoms[Symbol("PLOTchr2H")],chromosoms[Symbol("PLOTchr3H")],chromosoms[Symbol("PLOTchr4H")],chromosoms[Symbol("PLOTchr5H")],chromosoms[Symbol("PLOTchr6H")],chromosoms[Symbol("PLOTchr7H")])
                                                                           elseif size(chrom,1) == 2
                                                                                                  plotx = hstack(chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")])
                                                                           elseif size(chrom,1) == 1
                                                                                                  plotx = chromosoms[Symbol("PLOT$(chrom[1])")]
                                                                           elseif size(chrom,1) == 3
                                                                                                  plotx = hstack(chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")], chromosoms[Symbol("PLOT$(chrom[3])")])
                                                                           elseif size(chrom,1) == 4
                                                                                                  plotx = hstack(chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")], chromosoms[Symbol("PLOT$(chrom[3])")], chromosoms[Symbol("PLOT$(chrom[4])")])
                                                                           elseif size(chrom,1) == 5
                                                                                                  plotx = hstack(chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")], chromosoms[Symbol("PLOT$(chrom[3])")], chromosoms[Symbol("PLOT$(chrom[4])")],chromosoms[Symbol("PLOT$(chrom[5])")])
                                                                           elseif size(chrom,1) == 6
                                                                                                  plotx = hstack(chromosoms[Symbol("PLOT$(chrom[1])")],chromosoms[Symbol("PLOT$(chrom[2])")], chromosoms[Symbol("PLOT$(chrom[3])")], chromosoms[Symbol("PLOT$(chrom[4])")], chromosoms[Symbol("PLOT$(chrom[5])")], chromosoms[Symbol("PLOT$(chrom[6])")])
                                                                           else
                                                                                                  println("""The "chr" setting was used incorrectly, all Chromosoms are plotted""")
                                                                                                  plotx = hstack(chromosoms[Symbol("PLOTchr1H")],chromosoms[Symbol("PLOTchr2H")],chromosoms[Symbol("PLOTchr3H")],chromosoms[Symbol("PLOTchr4H")],chromosoms[Symbol("PLOTchr5H")],chromosoms[Symbol("PLOTchr6H")],chromosoms[Symbol("PLOTchr7H")])
                                                                           end
                                                                             return plotx
  end

#cp = vstack(changeplot(freq, "F3P1K1_Haplotypes", "F12P1K1_Haplotypes"),  changeplot(freq, "F3P1K1_Haplotypes", "F23P1K1_Haplotypes"),  changeplot(freq, "F3P1K1_Haplotypes", "F12P1Ö1_Haplotypes"),  changeplot(freq, "F3P1K1_Haplotypes", "F23P1Ö1_Haplotypes"))

latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=20pt,
minor_label_font="CMU Serif", minor_label_font_size=15pt,
key_title_font="CMU Serif", key_title_font_size=15pt,
key_label_font="CMU Serif", key_label_font_size=20pt)
Gadfly.push_theme(latex_fonts)
#=
using Cairo
draw(PNG("$save/location_all.png", 30inch, 20inch, dpi = 400),plotxy)
=#

#=
identifing candidate genes that are responsible and not only assosiated with a Freqschift by LD are hard to find and labor intensive to target down
this function should find potential candidate genes by matching different levels of information together and build a set of potential genes of interest
- starting with the information given in the genetic map files, going on with the physical map and genes and extract information of SNPs inside the coding regions as potential driver of variation
- Marker -> Gene Haplotype -> SNP + Sequence
=#
function investRegio(chr, region, target, pval=0.001, equal=true)

                    #chr = "chr4H" # specify the chromosom to look at
                    #region = [0,5] # specifiy the cM region on the chromosom to look at
                    #target = [22,23] # set the generations that should be different from the first generation
                    #pval = 0.001 # the p value as minimum treshold to accept a differnence
                    #equal = true # define if the Pick is present in both environments or just in one of both
                    ### load all data that is needed to do the analysis
                    #println("Load files..")
                    # genetic map
                    gloc = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/Dict_geentic_map_evol/Location.txt", delim="\t")
                    # physical map
                    ploc = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/Dict_physical_map_evol/Location.txt", delim="\t")
                    # results file of Marker and gene sets
                    mrorg = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/evol_gen_map/eco/results.txt", delim="\t")
                    mrkon = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/evol_gen_map/kon/results.txt", delim="\t")
                    # resluts of the physical map
                    prorg = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/evol_phys_map/eco/results.txt", delim="\t")
                    prkon = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/evol_phys_map/kon/results.txt", delim="\t")
                    # multi copy number variation files
                    cnv = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/Copy_number_variation/CNV_F3_F23.txt", delim="\t")
                    # comparisen between the konventinal and organic farming system
                    sys = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/System_comparisen/PhysicalMap/scomp_org_23_kon_23PhysMap.txt", delim ="\t")
                    # the system comparisen for the Markers might be helpful as well

                    ################################
                    # get the Marker names located in the region of interest
                    #println("Extract the Markers in the region of interest..")
                    gloc = unique(gloc[.&(gloc[!,:Chr].==chr, gloc[!,:Pos_genetic].> region[1], gloc[!,:Pos_genetic].< region[2]),3:15])
                    gm = gloc[!,:Marker][gloc[!,:Marker].!= "."] # the markers are extracted from the region - now process with the results files for organic and konventinal farming
                    # check in the results file for significant variation from the first generation for the generations selected
                    #println("Filter Markers according to the effect size by selected p value..")
                    mrorg = mrorg[indexin(gm, mrorg[!,:x1]),:]
                    mrkon = mrkon[indexin(gm, mrkon[!,:x1]),:]
                    # specify the generations to select
                    m = Symbol[]
                    n = Symbol[]
                    for i in 1:size(target,1); append!(m,names(mrorg)[occursin.(string(target[i]), String.(names(mrorg)))][[1,2,4],:]); end
                    for i in 1:size(target,1); append!(n,names(mrkon)[occursin.(string(target[i]), String.(names(mrkon)))][[1,2,4],:]); end
                    # to filter based on the p value difined - one of zeroinfl or NegBino has to be smaller

                    # reduce the list of Markers based on the statistical effect size estimated - the main interest is in an increasing Allele frequency, therefore only the genes with a higher Freq in the later generations will be considered
                    if size(target,1) == 1
                                                mrorg = mrorg[.&(.|(mrorg[!,m[1]].< pval, mrorg[!,m[2]].< pval), mrorg[!,m[3]].> mrorg[!,:f3e2_Average]),:]
                                                mrkon = mrkon[.&(.|(mrkon[n[1]].< pval, mrkon[n[2]].< pval), mrkon[!,n[3]].> mrkon[!,:f3e1_Average]),:]
                    elseif size(target,1) == 2
                                                mrorg = mrorg[.&(.|(mrorg[!,m[1]].< pval, mrorg[!,m[2]].< pval), .|(mrorg[!,m[4]].< pval, mrorg[!,m[5]].< pval), mrorg[!,m[3]].> mrorg[!,:f3e2_Average], mrorg[!,m[6]].> mrorg[!,:f3e2_Average]),:]
                                                mrkon = mrkon[.&(.|(mrkon[!,n[1]].< pval, mrkon[!,n[2]].< pval), .|(mrkon[!,n[4]].< pval, mrkon[!,n[5]].< pval), mrkon[!,n[3]].> mrkon[!,:f3e1_Average], mrkon[!,n[6]].> mrkon[!,:f3e1_Average]),:]
                    elseif size(target,1) == 3
                                                mrorg = mrorg[.&(.|(mrorg[!,m[1]].< pval, mrorg[!,m[2]].< pval), .|(mrorg[!,m[4]].< pval, mrorg[!,m[5]].< pval), .|(mrorg[!,m[7]].< pval, mrorg[!,m[8]].< pval), mrorg[!,m[3]].> mrorg[!,:f3e1_Average], mrorg[!,m[6]].> mrorg[!,:f3e2_Average], mrorg[!,m[9]].> mrorg[!,:f3e2_Average]),:]
                                                mrkon = mrkon[.&(.|(mrkon[!,n[1]].< pval, mrkon[!,n[2]].< pval), .|(mrkon[!,n[4]].< pval, mrkon[!,n[5]].< pval), .|(mrkon[!,n[7]].< pval, mrkon[!,n[8]].< pval), mrkon[!,n[3]].> mrkon[!,:f3e2_Average], mrkon[!,n[6]].> mrkon[!,:f3e1_Average], mrkon[!,n[9]].> mrkon[!,:f3e1_Average]),:]
                    elseif size(target,1) == 4
                                                mrorg = mrorg[.&(.|(mrorg[!,m[1]].< pval, mrorg[!,m[2]].< pval), .|(mrorg[!,m[4]].< pval, mrorg[!,m[5]].< pval), .|(mrorg[!,m[7]].< pval, mrorg[!,m[8]].< pval), .|(mrorg[!,m[10]].< pval, mrorg[!,m[11]].< pval), mrorg[!,m[3]].> mrorg[!,:f3e2_Average], mrorg[!,m[6]].> mrorg[!,:f3e2_Average], mrorg[!,m[9]].> mrorg[!,:f3e2_Average], mrorg[!,m[12]].> mrorg[!,:f3e2_Average]),:]
                                                mrkon = mrkon[.&(.|(mrkon[!,n[1]].< pval, mrkon[!,n[2]].< pval), .|(mrkon[!,n[4]].< pval, mrkon[!,n[5]].< pval), .|(mrkon[!,n[7]].< pval, mrkon[!,n[8]].< pval), .|(mrkon[!,n[10]].< pval, mrkon[!,n[11]].< pval), mrkon[!,n[3]].> mrkon[!,:f3e1_Average], mrkon[!,n[6]].> mrkon[!,:f3e1_Average], mrkon[!,n[9]].> mrkon[!,:f3e1_Average], mrkon[!,n[12]].> mrkon[!,:f3e1_Average]),:]
                    else
                                                println("This does not work, to many target Generations selected")
                    end
                    #########
                    # get the position of these markers to select the corresponding genes
                        # set by equal if to check for the genes that are equally expressed in the system sor the genes that differ
                    if equal == true
                                    seto = mrorg[indexin(mrkon[!,:x1], mrorg[!,:x1])[.!(isnothing.(indexin(mrkon[!,:x1], mrorg[!,:x1])))],:x1]
                                    setk = mrkon[indexin(mrorg[!,:x1], mrkon[!,:x1])[.!(isnothing.(indexin(mrorg[!,:x1], mrkon[!,:x1])))],:x1]
                                    set  = unique(vcat(seto,setk))
                    else
                                    or = indexin(mrkon[!,:x1], mrorg[!,:x1])[.!(isnothing.(indexin(mrkon[!,:x1], mrorg[!,:x1])))]
                                    kr = indexin(mrorg[!,:x1], mrkon[!,:x1])[.!(isnothing.(indexin(mrorg[!,:x1], mrkon[!,:x1])))]
                                    deleterows!(mrkon, kr)
                                    deleterows!(mrorg, or)
                                    set = vcat(mrkon[!,:x1],mrorg[!,:x1])
                    end

                    MR = gloc[indexin(set,gloc[!,:Marker]),1:5]

                    #println("Extract the Genes related to the Markers with diffential signal for the selected generations")
                    # extract the GENES that are located in the Markers identified
                    Ploc = unique(ploc[!,3:12]) # only the gene and not the SNP info is important
                    # create an empty dataframe to append to
                    ep = Ploc[1:2,:]
                    deleterows!(ep,1:2)
                    # write the Genes to these dataframes  # genes that do not fit completly into the marker will be omited without an additional extention of the marker bounds to both sides
                    for i in 1:size(MR,1); append!(ep,unique(Ploc[.&(Ploc[!,:x1].==MR[i,:Chr_1], Ploc[!,:start].> MR[i,:start]-10^4, Ploc[!,:end].< MR[i,:end]+10^4),:]));end
                    # remove "." genes
                    ep = ep[ep[!,:Gene_ID].!= ".",:]
                    # check the expression for these genes on the specified levels in the results file
                    #println("Filter Markers according to the effect size by selected p value..")
                    prorg = prorg[indexin(ep[!,:Gene_ID], prorg[!,:x1]),:]
                    prkon = prkon[indexin(ep[!,:Gene_ID], prkon[!,:x1]),:]

                    # to filter based on the p value difined - one of zeroinfl or NegBino has to be smaller

                    # reduce the list of Markers based on the statistical effect size estimated - the main interest is in an increasing Allele frequency, therefore only the genes with a higher Freq in the later generations will be considered
                    if size(target,1) == 1
                                                prorg = prorg[.&(.|(prorg[!,m[1]].< pval, prorg[!,m[2]].< pval), prorg[!,m[3]].> prorg[!,:f3e2_Average]),:]
                                                prkon = prkon[.&(.|(prkon[n[1]].< pval, prkon[n[2]].< pval), prkon[!,n[3]].> prkon[!,:f3e1_Average]),:]
                    elseif size(target,1) == 2
                                                prorg = prorg[.&(.|(prorg[!,m[1]].< pval, prorg[!,m[2]].< pval), .|(prorg[!,m[4]].< pval, prorg[!,m[5]].< pval), prorg[!,m[3]].> prorg[!,:f3e2_Average], prorg[!,m[6]].> prorg[!,:f3e2_Average]),:]
                                                prkon = prkon[.&(.|(prkon[!,n[1]].< pval, prkon[!,n[2]].< pval), .|(prkon[!,n[4]].< pval, prkon[!,n[5]].< pval), prkon[!,n[3]].> prkon[!,:f3e1_Average], prkon[!,n[6]].> prkon[!,:f3e1_Average]),:]
                    elseif size(target,1) == 3
                                                prorg = prorg[.&(.|(prorg[!,m[1]].< pval, prorg[!,m[2]].< pval), .|(prorg[!,m[4]].< pval, prorg[!,m[5]].< pval), .|(prorg[!,m[7]].< pval, prorg[!,m[8]].< pval), prorg[!,m[3]].> prorg[!,:f3e1_Average], prorg[!,m[6]].> prorg[!,:f3e2_Average], prorg[!,m[9]].> prorg[!,:f3e2_Average]),:]
                                                prkon = prkon[.&(.|(prkon[!,n[1]].< pval, prkon[!,n[2]].< pval), .|(prkon[!,n[4]].< pval, prkon[!,n[5]].< pval), .|(prkon[!,n[7]].< pval, prkon[!,n[8]].< pval), prkon[!,n[3]].> prkon[!,:f3e2_Average], prkon[!,n[6]].> prkon[!,:f3e1_Average], prkon[!,n[9]].> prkon[!,:f3e1_Average]),:]
                    elseif size(target,1) == 4
                                                prorg = prorg[.&(.|(prorg[!,m[1]].< pval, prorg[!,m[2]].< pval), .|(prorg[!,m[4]].< pval, prorg[!,m[5]].< pval), .|(prorg[!,m[7]].< pval, prorg[!,m[8]].< pval), .|(prorg[!,m[10]].< pval, prorg[!,m[11]].< pval), prorg[!,m[3]].> prorg[!,:f3e2_Average], prorg[!,m[6]].> prorg[!,:f3e2_Average], prorg[!,m[9]].> prorg[!,:f3e2_Average], prorg[!,m[12]].> prorg[!,:f3e2_Average]),:]
                                                prkon = prkon[.&(.|(prkon[!,n[1]].< pval, prkon[!,n[2]].< pval), .|(prkon[!,n[4]].< pval, prkon[!,n[5]].< pval), .|(prkon[!,n[7]].< pval, prkon[!,n[8]].< pval), .|(prkon[!,n[10]].< pval, prkon[!,n[11]].< pval), prkon[!,n[3]].> prkon[!,:f3e1_Average], prkon[!,n[6]].> prkon[!,:f3e1_Average], prkon[!,n[9]].> prkon[!,:f3e1_Average], prkon[!,n[12]].> prkon[!,:f3e1_Average]),:]
                    else
                                                println("This does not work, to many target Generations selected")
                    end

                    ### on the genetic map level, only the marker names where considered - on the level of the genes, the comaprisen of the genes will be performed by the system level file
                    # extract the genes of the system file present in prorg and prkon

                    set = unique(vcat(prkon[!,:x1], prorg[!,:x1]))
                    sys = sys[indexin(set, sys[!,:Haplotype]),:]
                    if equal == true  # filter the genes from the list that are not statistically different fro the environment
                                    sys = sys[.&(sys[!,:NegBin].> 0.1, sys[!,:Zeroinfl].> 0.1),: ]
                    else
                                    sys = sys[.|(sys[!,:NegBin].< pval, sys[!,:Zeroinfl].< pval),: ]
                    end

                    system = join(sys,Ploc, on= :Haplotype => :Gene_ID)
                    s2 = stack(system, [:orgAverage,:konAverage])
                    p = plot(s2, layer(x=:start, y=:value, color=:variable, Geom.point), layer(x=:start, y=:value, Geom.smooth),Guide.xlabel("Gene of System DF"), Guide.ylabel("Allelfreq [%]"), Guide.colorkey(title="System"), Guide.title("Allelfrequency of $chr, genetic Position $(region[1]) to $(region[2]) cM "),
                    Theme(background_color=colorant"white", default_color=colorant"grey"))

                    # join the system and the Location information
                    system = join(sys,Ploc, on= :Haplotype => :Gene_ID)
                    return system, p

 end
# qay = investRegio("chr2H", [70,90], [22,23],0.001, true)

function investRegio2(chr, l1, l2, rangek = 1, rangeo = 1, type="gen", pval=0.001)
                # set type of QTl region - genetic or physical map
                #chr,l1,l2 = "chr1H",0,28
                #pval = 0.001
                #type = "phy"
                if type == "phy"; l1, l2 =  l1 * 10^6, l2 * 10^6; end

                ## read data
                # position information
                ploc = CSV.read("/home/michael/Schreibtisch/Allelfreq/ref_snps_dulicates_removed2020/Dict_HC_Genes/Location.txt", delim="\t")
                gloc = CSV.read("/home/michael/Schreibtisch/Allelfreq/ref_snps_dulicates_removed2020/Dict_Marker/Location.txt", delim="\t")
                loc = unique(ploc[!,[6,9,10,11,12,4]])
                locp = unique(ploc[!,[6,1,4]])
                # result infomration
                #prorg = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/evol_phys_map/eco/results.txt", delim="\t")
                #prkon = CSV.read("/media/michael/Volume/WGS_Population_1/SNP_calling/ref_snps_duplucates_removed/evol_phys_map/kon/results.txt", delim="\t")
                sys = CSV.read("/home/michael/Schreibtisch/Allelfreq/ref_snps_dulicates_removed2020/System_comparison/HC_Gene/S23_S23.txt", delim ="\t")
                # protein code information
                prot = CSV.read("/home/michael/Schreibtisch/Allelfreq/ref_snps_dulicates_removed2020/Protein_comparisen_mafft_mega_besthit_SpliceVaraintsAdded.txt", delim="\t")


                # get the genes of the region selected
                if type == "gen" # in a genetic map scenario, translate the genetic to a physical map and then extrat the genes located in the area
                            # extract the marker based on the position
                            gp = unique(gloc[.&(gloc.Pos_genetic .> l1, gloc.Pos_genetic .< l2, gloc.Chr .== chr),[1,6,7,8]])
                            # sort by position
                            sort!(gp,:Pos_phys)
                            # extract the 3rd smallest and biggest position value for the physical map - not taking the biggest or smallest, because miss aligned markers can mess with the data
                            min = gp[3,4]
                            max = gp[size(gp,1)-2,4]


                            # check for genes located in this region and extract the list of genes
                            genes = unique(ploc[.&(ploc.Chr.==chr, ploc.Pos .> min, ploc.Pos .< max),:Gene_ID])
                            # remove "." character from the gene list
                            genes = genes[genes .!= "."]

                            # extract the genes from the list of candidates
                            set = sys[indexin(genes, sys.Haplotype),:]
                            # plot the genes in the area before filtering
                            system = join(set,locp, on= :Haplotype => :Gene_ID)
                            if size(system,1) > 1000
                            system = system[rand(1:size(system,1),1000),:]; end
                            s2 = stack(system, [:orgAverage,:konAverage])
                            p = plot(s2, layer(x=:start, y=:value, color=:variable, Geom.point), layer(x=:start, y=:value, Geom.smooth),Guide.xlabel("Position"), Guide.ylabel("Allelfreq [%]"), Guide.colorkey(title="System"), Guide.title("Allelfrequency of $chr, genetic Position $l1 to $l2 cM "),
                            Theme(background_color=colorant"white", default_color=colorant"grey"))


                            # set a minimum or maximum for the genes to match this allele freq
                            if rangeo > 1 && rangek < 1
                                                          set = set[.&(set.orgAverage .> rangeo-1, set.konAverage .< rangek),:]
                            elseif rangeo < 1 && rangek > 1
                                                          set = set[.&(set.orgAverage .< rangeo, set.konAverage .> rangek-1),:]
                            elseif rangeo < 1 && rangek < 1
                                                          set = set[.&(set.orgAverage .< rangeo, set.konAverage .< rangek),:]
                            elseif rangeo > 1 && rangek > 1
                                                          set = set[.&(set.orgAverage .> rangeo-1, set.konAverage .> rangek-1),:]
                            end


                            # highlight if the two are equal or differently
                            set[!,:Stat] .= "."
                            set[.&(set.NegBin .< pval, .|(set.Zeroinfl .< pval, set.Zeroinfl .== 2, set.Zeroinfl .== 3)),:Stat] .= "diff"
                            set[.&(set.NegBin .> pval*100, set.Zeroinfl .> pval*100),:Stat] .= "equal"
                            # remove unnecessary columns
                            set = set[!,[1,2,3,5,6,9,10,15]]
                            # merge the protein information for all given genes - remain the genes without nnotation in the list
                            set = join(set, prot[!,1:12], on=[:Haplotype => :Gene], kind=:left)

                            # if not equal, the set df can be sorted by the distance of organic to conventinal
                            set[!,:diff] .= abs.(set.orgAverage .- set.konAverage)
                            #set[!,:RD] .= (set.orgReaddepth .+ set.konReaddepth) ./ 2
                            # filter Genes with less than 100 reads coverage
                            #set = set[set.RD.>100,:]
                            #deletecols!(set, :RD)
                    #=        if equal == true
                                              sort!(set,order(:NegBin,rev=true))
                            elseif equal == false
                                                sort!(set,:NegBin)
                            end
                            =#
                            # remove the genes that have a sequenc einfromation and do not show any sight of variation
                            setp = set[.!(ismissing.(set.ISRStart)),:]
                            setp = setp[setp.Mismatches .> 1,:]
                            #setp = setp[.|(.&(setp.DeletionGolf .> 0, setp.DeletionGolf .!= setp.DeletionISR), .&(setp.DeletionISR .> 0, setp.DeletionGolf .!= setp.DeletionISR)),:]
                            # align the function of the genes to this list
                            setp = join(setp, loc, on=[:Haplotype => :Gene_ID])

                            setm =set[ismissing.(set.ISRStart),:]
                            setm = join(setm, loc, on=[:Haplotype => :Gene_ID])

                elseif type == "phy"
                            # extract the marker based on the position
                            gp = unique(ploc[.&(ploc.Pos .> l1, ploc.Pos .< l2, ploc.Chr .== chr),[6,7,9,10,11,12]])
                            # remove "."
                            genes = gp[gp.Gene_ID .!= ".",:Gene_ID]
                            # extract the genes from the list of candidates
                            set = sys[indexin(genes, sys.Haplotype),:]

                            # plot the genes in the area before filtering
                            system = join(set,locp, on= :Haplotype => :Gene_ID)
                            s2 = stack(system, [:orgAverage,:konAverage])
                            if size(system,1) > 1000
                            system = system[rand(1:size(system,1),1000),:]; end
                            p = plot(s2, layer(x=:start, y=:value, color=:variable, Geom.point), layer(x=:start, y=:value, Geom.smooth),Guide.xlabel("Position"), Guide.ylabel("Allelfreq [%]"), Guide.colorkey(title="System"), Guide.title("Allelfrequency of $chr, genetic Position $l1 to $l2 bp "),
                            Theme(background_color=colorant"white", default_color=colorant"grey"))


                            # set a minimum or maximum for the genes to match this allele freq
                            if rangeo > 1 && rangek < 1
                                                          set = set[.&(set.orgAverage .> rangeo-1, set.konAverage .< rangek),:]
                            elseif rangeo < 1 && rangek > 1
                                                          set = set[.&(set.orgAverage .< rangeo, set.konAverage .> rangek-1),:]
                            elseif rangeo < 1 && rangek < 1
                                                          set = set[.&(set.orgAverage .< rangeo, set.konAverage .< rangek),:]
                            elseif rangeo > 1 && rangek > 1
                                                          set = set[.&(set.orgAverage .> rangeo-1, set.konAverage .> rangek-1),:]
                            end

                            # subset the list by a given pval for either sides
                            set[!,:Stat] .= "."
                            set[.&(set.NegBin .< pval, .|(set.Zeroinfl .< pval, set.Zeroinfl .== 2, set.Zeroinfl .== 3)),:Stat] .= "diff"
                            set[.&(set.NegBin .> pval*100, set.Zeroinfl .> pval*100),:Stat] .= "equal"
                            # remove unnecessary columns
                            set = set[!,[1,2,3,5,6,9,10,15]]
                            # merge the protein information for all given genes - remain the genes without nnotation in the list
                            set = join(set, prot[!,1:12], on=[:Haplotype => :Gene], kind=:left)


                            # if not equal, the set df can be sorted by the distance of organic to conventinal
                            set[!,:diff] .= abs.(set.orgAverage .- set.konAverage)
                            #set[!,:RD] .= (set.orgReaddepth .+ set.konReaddepth) ./ 2
                            # filter Genes with less than 100 reads coverage
                            #set = set[set.RD.>100,:]
                            #deletecols!(set, :RD)

                    #=        if equal == true
                                              sort!(set,order(:NegBin,rev=true))
                            elseif equal == false
                                                sort!(set,:NegBin)
                            end
                            =#

                            # remove the genes that have a sequenc einfromation and do not show any sight of variation
                            setp = set[.!(ismissing.(set.ISRStart)),:]
                            setp = setp[setp.Mismatches .> 1,:]
                            #setp = setp[.|(.&(setp.DeletionGolf .> 0, setp.DeletionGolf .!= setp.DeletionISR), .&(setp.DeletionISR .> 0, setp.DeletionGolf .!= setp.DeletionISR)),:]
                            # align the function of the genes to this list
                            setp = join(setp, loc, on=[:Haplotype => :Gene_ID])

                            setm =set[ismissing.(set.ISRStart),:]
                            setm = join(setm, loc, on=[:Haplotype => :Gene_ID])

                end # end of physical genetic choice

                return setp, setm, p
 end
# qw = investRegio2("chr3H", 120,124, "phy", false, 0.001)
