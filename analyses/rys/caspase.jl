using CSV, DataFrames, Plots, Distributions, Measures
import Plots:Plot

gr()

caspase_pth="/bench/PhD/Thesis/datasets/rys caspase.csv"
caspasedf=DataFrame(CSV.read(caspase_pth))

sib_ms_dict=Dict{String,Vector{Vector{Float64}}}()
rys_ms_dict=Dict{String,Vector{Vector{Float64}}}()

X=[4.,6.]

for md in [sib_ms_dict,rys_ms_dict]
    md["CMZ"]=[zeros(0) for i in 1:2]
    md["Border"]=[zeros(0) for i in 1:2]
    md["Central"]=[zeros(0) for i in 1:2]
    md["Total"]=[zeros(0) for i in 1:2]
end

for row in eachrow(caspasedf)
    row.dpf==4 ? (xidx=1) : (xidx=2)
    row.genotype=="sib" ? (md=sib_ms_dict) : (md=rys_ms_dict)
    push!(md["CMZ"][xidx],row.CMZ)
    push!(md["Border"][xidx],row."CMZ/central border region")
    push!(md["Central"][xidx],row."central retina")
    push!(md["Total"][xidx],row.CMZ + row."CMZ/central border region" + row."central retina")
end

plot_halves=Vector{Plot}()
for (n, age) in enumerate(X)
    subplots=Vector{Plot}()
    for (m,ms) in enumerate(["CMZ","Border","Central","Total"])
        s_data=sib_ms_dict[ms][n]
        sl=length(s_data)
        r_data=rys_ms_dict[ms][n]
        rl=length(r_data)
        if m == 1
            subplt=bar(1:1:sl,s_data,barwidths=1,color=:green,xticks=(1:1:sl+rl,vcat(["S" for i in 1:sl],["R" for i in 1:rl])),legend=:none, xlabel=ms, margin=0mm,ylabel="Caspase-3+ cells")
        else
            subplt=bar(1:1:sl,s_data,barwidths=1,color=:green,xticks=(1:1:sl+rl,vcat(["S" for i in 1:sl],["R" for i in 1:rl])),legend=:none, xlabel=ms, leftmargin=-7.5mm, yformatter=_->"")
        end 
        bar!(sl+1:1:sl+rl,r_data,color=:darkmagenta, barwidths=1.)
        push!(subplots,subplt)
    end
    title = plot(title = "$(Int64(age)) dpf", grid = false, showaxis = nothing, bottom_margin = -35Plots.px)
    l=@layout [t{.02h};[a{.25w} b{.25w} c{.25w} d{.25w}]]
    half=plot(title, subplots...,layout=l,link=:y)
    push!(plot_halves,half)
end

combined=plot(plot_halves...,layout=grid(1,2), size=(900,450))

savefig(combined,"/bench/PhD/Thesis/images/rys/caspase.png")

println("Rys CMZ caspase 4dpf: $(mean(rys_ms_dict["CMZ"][1])) ± $(std(rys_ms_dict["CMZ"][1]))")
println("Rys CMZ caspase 6dpf: $(mean(rys_ms_dict["CMZ"][2])) ± $(std(rys_ms_dict["CMZ"][2]))")

println("Rys central caspase 4dpf: $(mean(rys_ms_dict["Central"][1])) ± $(std(rys_ms_dict["Central"][1]))")
println("Rys central caspase 6dpf:  $(mean(rys_ms_dict["Central"][2])) ± $(std(rys_ms_dict["Central"][2]))")
println("Sib central caspase 4dpf: $(mean(sib_ms_dict["Central"][1])) ± $(std(sib_ms_dict["Central"][1]))")
println("Sib central caspase 6dpf: $(mean(sib_ms_dict["Central"][2])) ± $(std(sib_ms_dict["Central"][2]))")