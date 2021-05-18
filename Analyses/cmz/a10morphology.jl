using CSV,DataFrames,Distributions,NGRefTools,StatsBase,Plots
gr()
default(legendfont = (10), guidefont = (11), tickfont = (10))

a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["Retina thickness"]=Vector{Vector{Float64}}()
measure_dict["RPE length"]=Vector{Vector{Float64}}()
measure_dict["RPE thickness"]=Vector{Vector{Float64}}()
measure_dict["Lens diameter"]=Vector{Vector{Float64}}()
measure_dict["ON diameter"]=Vector{Vector{Float64}}()
measure_dict["ONL"]=Vector{Vector{Float64}}()
measure_dict["OPL"]=Vector{Vector{Float64}}()
measure_dict["INL"]=Vector{Vector{Float64}}()
measure_dict["IPL"]=Vector{Vector{Float64}}()
measure_dict["GCL"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    push!(X,Float64(unique(t_df."Time point (d)")[1]))
    rtvec, rplvec, rptvec, ldvec, ondvec, onlvec, oplvec, inlvec, iplvec, gclvec = [zeros(0) for i in 1:10]

    for i_df in groupby(t_df,"BSR")
        length([skipmissing(i_df."Lens diameter")...])>0 && push!(ldvec,mean(skipmissing(i_df."Lens diameter")))

        length([skipmissing(i_df."Retina thickness")...])>0 && push!(rtvec,mean(skipmissing(i_df."Retina thickness")))
        length([skipmissing(i_df."RPE length")...])>0 && push!(rplvec,mean(skipmissing(i_df."RPE length")))
        length([skipmissing(i_df."RPE thickness")...])>0 && push!(rptvec,mean(skipmissing(i_df."RPE thickness")))
        length([skipmissing(i_df."ON diameter")...])>0 && push!(ondvec,mean(skipmissing(i_df."ON diameter")))
        length([skipmissing(i_df.ONL)...])>0 && push!(onlvec,mean(skipmissing(i_df.ONL)))   
        length([skipmissing(i_df.OPL)...])>0 && push!(oplvec,mean(skipmissing(i_df.OPL)))  
        length([skipmissing(i_df.INL)...])>0 && push!(inlvec,mean(skipmissing(i_df.INL)))
        length([skipmissing(i_df.IPL)...])>0 && push!(iplvec,mean(skipmissing(i_df.IPL)))  
        length([skipmissing(i_df.GCL)...])>0 && push!(gclvec,mean(skipmissing(i_df.GCL)))
    end

    push!(measure_dict["Retina thickness"],rtvec)
    push!(measure_dict["RPE length"],rplvec)
    push!(measure_dict["RPE thickness"],rptvec)
    push!(measure_dict["Lens diameter"],ldvec)
    push!(measure_dict["ON diameter"],ondvec)
    push!(measure_dict["ONL"],onlvec)
    push!(measure_dict["OPL"],oplvec)
    push!(measure_dict["INL"],inlvec)
    push!(measure_dict["IPL"],iplvec)
    push!(measure_dict["GCL"],gclvec)
end

#TEST IF DATA IS BETTER MODELLED NORMALLY OR LOGNORMALLY
for ms in ["Retina thickness", "RPE length", "RPE thickness", "Lens diameter", "ON diameter", "ONL", "OPL", "INL", "IPL", "GCL"]
    ms=="ON diameter" ? (idxs=length(X)-1) : (idxs=length(X))

    normal_mts=[fit(Normal,measure_dict[ms][n]) for n in 1:idxs]
    logn_mts=[fit(LogNormal,filter!(i->!iszero(i),measure_dict[ms][n])) for n in 1:idxs]

    normal_lh=sum([sum(logpdf(normal_mts[n],measure_dict[ms][n])) for n in 1:idxs])
    logn_lh=sum([sum(logpdf(logn_mts[n],measure_dict[ms][n])) for n in 1:idxs])

    println("$ms & $(round(normal_lh,digits=3)) & $(round(logn_lh,digits=3)) & $(round(logn_lh-normal_lh,digits=3))\\\\ \\hline")
end

plotticks=[30*i for i in 1:12]

rt_n_mts=[fit(MarginalTDist,measure_dict["Retina thickness"][n]) for n in 1:length(X)]
rtn_mean=[mean(mt) for mt in rt_n_mts]
rtn_lower=[mean(mt)-quantile(mt,.025) for mt in rt_n_mts]
rtn_upper=[quantile(mt,.975)-mean(mt) for mt in rt_n_mts]

rtn_chart=scatter(vcat([[X[n] for i in 1:length(measure_dict["Retina thickness"][n])] for n in 1:length(X)]...),vcat(measure_dict["Retina thickness"]...), marker=:cross, color=:black, markersize=3, label="Retinal thic. data", ylabel="Retinal thickness (μm)", showaxis=:y, xticks=plotticks, xformatter=_->"", legend=:bottomright)
plot!(X, rtn_mean, ribbon=(rtn_lower,rtn_upper), color=:darkmagenta, label="Mean thickness")
annotate!([(8,200,Plots.text("A",18))])

rpl_n_mts=[fit(MarginalTDist,measure_dict["RPE length"][n]) for n in 1:length(X)]
rpln_mean=[mean(mt) for mt in rpl_n_mts]
rpln_lower=[mean(mt)-quantile(mt,.025) for mt in rpl_n_mts]
rpln_upper=[quantile(mt,.975)-mean(mt) for mt in rpl_n_mts]

rpln_chart=scatter(vcat([[X[n] for i in 1:length(measure_dict["RPE length"][n])] for n in 1:length(X)]...),vcat(measure_dict["RPE length"]...), marker=:cross, color=:black, markersize=3, label="RPE length data", ylabel="RPE length (μm)", showaxis=:y, xticks=plotticks, xformatter=_->"", legend=:bottomright)
plot!(X, rpln_mean, ribbon=(rpln_lower,rpln_upper), color=:grey, label="Mean thickness")
annotate!([(8,3000,Plots.text("B",18))])

ld_n_mts=[fit(MarginalTDist,measure_dict["Lens diameter"][n]) for n in 1:length(X)]
ldn_mean=[mean(mt) for mt in ld_n_mts]
ldn_mean[end]=ld_n_mts[end].μ
ldn_mean[end-2]=ld_n_mts[end-2].μ
ldn_lower=[mean(mt)-quantile(mt,.025) for mt in ld_n_mts]
ldn_lower[end]=ld_n_mts[end].μ - (ld_n_mts[end].μ-3*ld_n_mts[end].σ)
ldn_lower[end-2]=ld_n_mts[end-2].μ - (ld_n_mts[end-2].μ-3*ld_n_mts[end-2].σ)
ldn_upper=[quantile(mt,.975)-mean(mt) for mt in ld_n_mts]
ldn_upper[end]=(ld_n_mts[end].μ+3*ld_n_mts[end].σ) - ld_n_mts[end].μ
ldn_upper[end-2]=(ld_n_mts[end-2].μ+3*ld_n_mts[end-2].σ) - ld_n_mts[end-2].μ

ldn_chart=scatter(vcat([[X[n] for i in 1:length(measure_dict["Lens diameter"][n])] for n in 1:length(X)]...),vcat(measure_dict["Lens diameter"]...), marker=:cross, color=:black, markersize=3, label="Lens dia. data", ylabel="Lens diameter (μm)", showaxis=:y, xticks=plotticks, xformatter=_->"", legend=:topleft)
plot!(X, ldn_mean, ribbon=(ldn_lower,ldn_upper), color=:darkgreen, label="Diameter mean")
annotate!([(8,600,Plots.text("C",18))])

on_n_mts=[fit(MarginalTDist,measure_dict["ON diameter"][n]) for n in 1:length(X)-1]
onn_mean=[mean(mt) for mt in on_n_mts]
onn_mean[end]=on_n_mts[end].μ
onn_lower=[mean(mt)-quantile(mt,.025) for mt in on_n_mts]
onn_lower[end]=on_n_mts[end].μ - (on_n_mts[end].μ-3*on_n_mts[end].σ)
onn_upper=[quantile(mt,.975)-mean(mt) for mt in on_n_mts]
onn_upper[end]=(on_n_mts[end].μ+3*on_n_mts[end].σ) - on_n_mts[end].μ

onn_chart=scatter(vcat([[X[n] for i in 1:length(measure_dict["ON diameter"][n])] for n in 1:length(X)-1]...),vcat(measure_dict["ON diameter"]...), marker=:cross, color=:black, markersize=3, label="ON dia. data", ylabel="Optic nerve diameter (μm)", showaxis=:y, xticks=plotticks, xformatter=_->"", legend=:topleft, ylims=[0,200],xlims=[0,360])
plot!(X[1:end-1], onn_mean, ribbon=(onn_lower,onn_upper), color=:darkblue, label="Diameter mean")
annotate!([(8,105,Plots.text("D",18))])
lens!([2, 31], [0, 30], inset = (1, bbox(0.55, 0.1, .45, .8)))

stackdata=hcat(mean.(measure_dict["RPE thickness"]),mean.(measure_dict["ONL"]),mean.(measure_dict["OPL"]),mean.(measure_dict["INL"]),mean.(measure_dict["IPL"]),mean.(measure_dict["GCL"]))

layers=areaplot(X, stackdata, xticks=plotticks, legend=:none, xlabel="Age (dpf)", ylabel="Layer thickness (μm)")
plot!([100,100],[220,240], color=:black)
annotate!([(100,250,Plots.text("GCL",12))])
plot!([130,130],[195,240], color=:black)
annotate!([(130,250,Plots.text("IPL",12))])
plot!([160,160],[130,240], color=:black)
annotate!([(160,250,Plots.text("INL",12))])
plot!([190,190],[100,240], color=:black)
annotate!([(190,250,Plots.text("OPL",12))])
plot!([220,220],[75,240], color=:black)
annotate!([(220,250,Plots.text("ONL",12))])
plot!([250,250],[25,240], color=:black)
annotate!([(250,250,Plots.text("RPE",12))])
annotate!([(8.,175,Plots.text("E",18))])

combined=plot(rtn_chart, rpln_chart, ldn_chart, onn_chart, layers, layout=grid(5,1), size=(1000,1200), link=:x)
savefig(combined, "/bench/PhD/Thesis/images/cmz/morphology.png")