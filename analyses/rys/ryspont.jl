using CSV,DataFrames,Distributions,NGRefTools,StatsBase,StatsPlots
gr()
default(legendfont = (8), guidefont = (12,:bold), tickfont = (8), guide = "x")


a38pth="/bench/PhD/datasets/A38.csv"
a49pth="/bench/PhD/datasets/A49.csv"

a38df=DataFrame(CSV.read(a38pth))
a38df.Group[1:10].="sib"
a38df.Group[11:22].="rys"
a49df=DataFrame(CSV.read(a49pth))
a49df=a49df[2:end,:]

a38df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a38df)]

sib_measure_dict=Dict{String,Vector{Vector{Float64}}}()
rys_measure_dict=Dict{String,Vector{Vector{Float64}}}()

XCMZ=Vector{Float64}()
XEdU=Vector{Float64}()

for dict in [sib_measure_dict,rys_measure_dict]
    dict["DCMZ"]=Vector{Float64}()
    dict["VCMZ"]=Vector{Float64}()
    dict["TCMZ"]=Vector{Float64}()
    dict["EdU"]=Vector{Float64}()
end

for tdf in groupby(a49df,"Age")
    age=unique(tdf."Age")[1]
    push!(XCMZ,age);push!(XEdU,age)
    xidx=length(XCMZ)
    for sdf in groupby(tdf, "Rys/Sib")
        sr=unique(sdf."Rys/Sib")[1]
        for mdf in groupby(sdf,"D/V/P/WE")
            if unique(mdf."D/V/P/WE")[1] == "D"
                if sr=="sib"
                    push!(sib_measure_dict["DCMZ"],mdf."PCNA +ve")
                    push!(sib_measure_dict["TCMZ"],mdf."PCNA +ve")
                    push!(sib_measure_dict["EdU"],mdf."EdU +ve, PCNA +ve")
                else
                    push!(rys_measure_dict["DCMZ"],mdf."PCNA +ve")
                    push!(rys_measure_dict["TCMZ"],mdf."PCNA +ve")
                    push!(rys_measure_dict["EdU"],mdf."EdU +ve, PCNA +ve")
                end
            elseif unique(mdf."D/V/P/WE")[1] == "V"
                if sr=="sib"
                    push!(sib_measure_dict["VCMZ"],mdf."PCNA +ve")
                    sib_measure_dict["TCMZ"][xidx].+=mdf."PCNA +ve"
                    sib_measure_dict["EdU"][xidx].+=mdf."EdU +ve, PCNA +ve"
                else
                    push!(rys_measure_dict["VCMZ"],mdf."PCNA +ve")
                    rys_measure_dict["TCMZ"][xidx].+=mdf."PCNA +ve"
                    rys_measure_dict["EdU"][xidx].+=mdf."EdU +ve, PCNA +ve"
                end
            end
        end
    end
end

push!(XCMZ,5.)
for sdf in groupby(a38df,"Group")
    xidx=length(XCMZ)
    if unique(sdf."Group")[1]=="sib"
        for mdf in groupby(sdf,"D/V")
            if unique(mdf."D/V")[1] == "D"
                push!(sib_measure_dict["DCMZ"],mdf."PCNA +ve")
                push!(sib_measure_dict["TCMZ"],mdf."PCNA +ve")
            else
                push!(sib_measure_dict["VCMZ"],mdf."PCNA +ve")
                sib_measure_dict["TCMZ"][xidx].+=mdf."PCNA +ve"
            end
        end
    else
        for mdf in groupby(sdf,"D/V")
            if unique(mdf."D/V")[1] == "D"
                push!(rys_measure_dict["DCMZ"],mdf."PCNA +ve")
                push!(rys_measure_dict["TCMZ"],mdf."PCNA +ve")
            else
                push!(rys_measure_dict["VCMZ"],mdf."PCNA +ve")
                rys_measure_dict["TCMZ"][xidx].+=mdf."PCNA +ve"
            end
        end
    end
end

perm=sortperm(XCMZ)
for ms in ["DCMZ","VCMZ","TCMZ"]
    sib_measure_dict[ms]=sib_measure_dict[ms][perm]
    rys_measure_dict[ms]=rys_measure_dict[ms][perm]
end
XCMZ=XCMZ[perm]

sp_logn_mts=[fit(MarginalTDist,log.(sib_measure_dict["TCMZ"][n])) for n in 1:length(XCMZ)]
sp_mean=[exp(mean(mt)) for mt in sp_logn_mts]
sp_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in sp_logn_mts]
sp_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in sp_logn_mts]

rp_logn_mts=[fit(MarginalTDist,log.(rys_measure_dict["TCMZ"][n])) for n in 1:length(XCMZ)]
rp_mean=[exp(mean(mt)) for mt in rp_logn_mts]
rp_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in rp_logn_mts]
rp_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in rp_logn_mts]

log_mean_mass_comparator(rys_measure_dict["TCMZ"][1],sib_measure_dict["TCMZ"][1],labels=["rys $(XCMZ[1]) dpf sectional CMZ population","sib $(XCMZ[1]) dpf sectional CMZ population"])
log_mean_mass_comparator(rys_measure_dict["TCMZ"][end],sib_measure_dict["TCMZ"][end],labels=["rys $(XCMZ[end]) dpf sectional CMZ population","sib $(XCMZ[end]) dpf sectional CMZ population"])

println("At 4dpf, $(mean(sib_measure_dict["EdU"][1]./sib_measure_dict["TCMZ"][1])) ± $(std(sib_measure_dict["EdU"][1]./sib_measure_dict["TCMZ"][1])) of sib PCNA+ve are EdU positive")

println("At 4dpf, $(mean(rys_measure_dict["EdU"][1]./rys_measure_dict["TCMZ"][1])) ± $(std(rys_measure_dict["EdU"][1]./rys_measure_dict["TCMZ"][1])) of rys PCNA+ve are EdU positive")

println("At 10dpf, $(mean(sib_measure_dict["EdU"][end]./sib_measure_dict["TCMZ"][end])) ± $(std(sib_measure_dict["EdU"][end]./sib_measure_dict["TCMZ"][end])) of sib PCNA+ve are EdU positive")


println("At 10dpf, $(mean(rys_measure_dict["EdU"][end]./rys_measure_dict["TCMZ"][end])) ± $(std(rys_measure_dict["EdU"][end]./rys_measure_dict["TCMZ"][end])) of rys PCNA+ve are EdU positive")

se_logn_mts=[fit(MarginalTDist,log.(sib_measure_dict["EdU"][n])) for n in 1:length(XEdU)]
se_mean=[exp(mean(mt)) for mt in se_logn_mts]
se_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in se_logn_mts]
se_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in se_logn_mts]

re_logn_mts=[fit(MarginalTDist,log.(filter(val->!iszero(val),rys_measure_dict["EdU"][n]))) for n in 1:length(XEdU)]
re_mean=[exp(mean(mt)) for mt in re_logn_mts]
re_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in re_logn_mts]
re_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in re_logn_mts]


pcna_plt=plot(XCMZ,sp_mean,ribbon=(sp_lower,sp_upper),color=:green, label="sib mean",xlabel="Age (dpf)",ylabel="PCNA+ve CMZ Cells")
plot!(pcna_plt,XCMZ,rp_mean,ribbon=(rp_lower,rp_upper),color=:darkmagenta, label="rys mean",legend=:none)
scatter!(pcna_plt,vcat([[XCMZ[n] for i in 1:length(sib_measure_dict["TCMZ"][n])] for n in 1:length(XCMZ)]...),vcat(sib_measure_dict["TCMZ"]...),marker=:circle,markersize=:3,color=:green,label="sib data")
scatter!(pcna_plt,vcat([[XCMZ[n] for i in 1:length(rys_measure_dict["TCMZ"][n])] for n in 1:length(XCMZ)]...),vcat(rys_measure_dict["TCMZ"]...),marker=:dtriangle,markersize=:3,color=:darkmagenta,label="rys data")
annotate!([(4.5,175,Plots.text("A",18))])

edu_plt=plot(XEdU,se_mean,ribbon=(se_lower,se_upper),color=:green, label="sib mean",xlabel="Age (dpf)",ylabel="EdU+ve CMZ Cells")
plot!(edu_plt,XEdU,re_mean,ribbon=(re_lower,re_upper),color=:darkmagenta, label="rys EdU")
scatter!(edu_plt,vcat([[XEdU[n] for i in 1:length(sib_measure_dict["EdU"][n])] for n in 1:length(XEdU)]...),vcat(sib_measure_dict["EdU"]...),marker=:circle,markersize=:3,color=:green,label="sib data")
scatter!(edu_plt,vcat([[XEdU[n] for i in 1:length(rys_measure_dict["EdU"][n])] for n in 1:length(XEdU)]...),vcat(rys_measure_dict["EdU"]...),marker=:dtriangle,markersize=:3,color=:darkmagenta,label="rys data")
annotate!([(4.5,150,Plots.text("B",18))])


combined_plt=plot(pcna_plt,edu_plt,layout=grid(1,2),size=(900,500))
savefig(combined_plt,"/bench/PhD/Thesis/images/rys/CMZontogeny.png")