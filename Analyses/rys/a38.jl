using CSV,DataFrames,Distributions,NGRefTools,StatsBase,StatsPlots
gr()
default(legendfont = 10, guidefont = (12,:bold), tickfont = 10, guide = "x")


pth="/bench/PhD/datasets/A38.csv"

a38df=DataFrame(CSV.read(pth))
a38df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a38df)]
a38df.Group[1:10].="sib"
a38df.Group[11:22].="rys"

sib_measure_dict=Dict{String,Vector{Float64}}()
rys_measure_dict=Dict{String,Vector{Float64}}()

for dict in [sib_measure_dict,rys_measure_dict]
    dict["DCMZ"]=Vector{Float64}()
    dict["VCMZ"]=Vector{Float64}()
    dict["TCMZ"]=Vector{Float64}()
    dict["NucPop"]=Vector{Float64}()
    dict["Sphericity"]=Vector{Float64}()
    dict["Volume"]=Vector{Float64}()
end

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for (gdf,dict) in zip(groupby(a38df, "Group"),[sib_measure_dict,rys_measure_dict])
    for row in eachrow(gdf)
        if row."D/V"=="D"
            push!(dict["DCMZ"],row."PCNA +ve")
            push!(dict["TCMZ"],row."Total PCNA")
            push!(dict["NucPop"],row."Total nuclei")
            push!(dict["Sphericity"],row."Sphericity")
            push!(dict["Volume"],row."Volume")
        else
            push!(dict["VCMZ"],row."PCNA +ve")
            dict["Sphericity"][end]+=row."Sphericity";dict["Sphericity"][end]/=2
            dict["Volume"][end]+=row."Volume";dict["Volume"][end]/=2
        end
    end
end
            
#TEST IF DATA IS BETTER MODELLED NORMALLY OR LOGNORMALLY
for ms in ["DCMZ", "VCMZ", "NucPop", "Sphericity", "Volume"]
    n_ms=[fit(Normal,dict[ms]) for dict in [sib_measure_dict, rys_measure_dict]]
    logn_ms=[fit(LogNormal,dict[ms]) for dict in [sib_measure_dict, rys_measure_dict]]

    normal_lh=sum(logpdf(n_ms[1],sib_measure_dict[ms]))+sum(logpdf(n_ms[2],rys_measure_dict[ms]))
    logn_lh=sum(logpdf(logn_ms[1],sib_measure_dict[ms]))+sum(logpdf(logn_ms[2],rys_measure_dict[ms]))

    println("$ms & $(round(normal_lh,digits=3)) & $(round(logn_lh,digits=3)) & $(round(logn_lh-normal_lh,digits=3))\\\\ \\hline")
end


TCMZ_plt=plot_logn_MTDist([sib_measure_dict["TCMZ"],rys_measure_dict["TCMZ"]],[:green,:darkmagenta],[:circle,:dtriangle],["sib","rys"],"Sectional CMZ Population","Posterior mean lh.")
annotate!([(55,3.8,Plots.text("A",18))])

PopRatio_plt=plot_logn_MTDist([sib_measure_dict["NucPop"]./sib_measure_dict["TCMZ"],rys_measure_dict["NucPop"]./rys_measure_dict["TCMZ"]],[:green,:darkmagenta],[:circle,:dtriangle],["sib","rys"],"Total Sectional Population per CMZ Cell","Posterior mean lh.")
annotate!([(7.5,3.,Plots.text("B",18))])

DCMZ_plt=plot_logn_MTDist([sib_measure_dict["DCMZ"],rys_measure_dict["DCMZ"]],[:green,:darkmagenta],[:circle,:dtriangle],["sib","rys"],"Dorsal CMZ Population","Posterior mean lh.")
annotate!([(35,4.5,Plots.text("C",18))])

log_mean_mass_comparator(rys_measure_dict["DCMZ"],sib_measure_dict["DCMZ"],labels=["rys 5dpf sectional DCMZ population","sib 5dpf sectional DCMZ population"])
log_mean_mass_comparator(rys_measure_dict["VCMZ"],sib_measure_dict["VCMZ"],labels=["rys 5dpf sectional VCMZ population","sib 5dpf sectional VCMZ population"])

VCMZ_plt=plot_logn_MTDist([sib_measure_dict["VCMZ"],rys_measure_dict["VCMZ"]],[:green,:darkmagenta],[:circle,:dtriangle],["sib","rys"],"Ventral CMZ Population","Posterior mean lh.")
annotate!([(15,2.25,Plots.text("D",18))])


Vol_plt=plot_n_MTDist([sib_measure_dict["Volume"],rys_measure_dict["Volume"]],[:green,:darkmagenta],[:circle,:dtriangle],["sib","rys"],"Mean nuclear volume","Posterior mean lh.")
annotate!([(.2,105,Plots.text("E",18))])

Sph_plt=plot_n_MTDist([sib_measure_dict["Sphericity"],rys_measure_dict["Sphericity"]],[:green,:darkmagenta],[:circle,:dtriangle],["sib","rys"],"Mean nuclear sphericity","Posterior mean lh.")
annotate!([(.65,325,Plots.text("F",18))])
plot!(legend=:top)

combined=plot(TCMZ_plt,PopRatio_plt,DCMZ_plt,VCMZ_plt,Vol_plt,Sph_plt,layout=grid(3,2), size=(900,900))

savefig(combined,"/bench/PhD/Thesis/images/rys/nuclearstudy.png")