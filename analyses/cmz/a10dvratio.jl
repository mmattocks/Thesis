using BayesianLinearRegression,CSV,DataFrames,Distributions,BioSimpleStochastic,GMC_NS, Serialization, Plots, NGRefTools
a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["Dorsal CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Ventral CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Ratio"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    push!(X,Float64(unique(t_df."Time point (d)")[1]))
    rvec,dvec,vvec = [zeros(0) for i in 1:3]

    for i_df in groupby(t_df,"BSR")
        length([skipmissing(i_df."Dorsal CMZ (#)")...])>0 && mean(skipmissing(i_df."Dorsal CMZ (#)")) >0. && push!(dvec,mean(skipmissing(i_df."Dorsal CMZ (#)")))
        length([skipmissing(i_df."Ventral CMZ (#)")...])>0 && mean(skipmissing(i_df."Ventral CMZ (#)")) >0. && push!(vvec,mean(skipmissing(i_df."Ventral CMZ (#)")))
        length([skipmissing(i_df."Dorsal CMZ (#)")...])>0 && mean(skipmissing(i_df."Dorsal CMZ (#)")) >0. && length([skipmissing(i_df."Ventral CMZ (#)")...])>0 && mean(skipmissing(i_df."Ventral CMZ (#)")) >0. && push!(rvec,mean(skipmissing(i_df."Dorsal CMZ (#)"))/mean(skipmissing(i_df."Ventral CMZ (#)")))
    end

    push!(measure_dict["Ratio"],rvec)
    push!(measure_dict["Dorsal CMZ (#)"],dvec)
    push!(measure_dict["Ventral CMZ (#)"],vvec)
end

d_logn_mts=[fit(MarginalTDist,log.(measure_dict["Dorsal CMZ (#)"][n])) for n in 1:length(X)]
d_mean=[exp(mean(mt)) for mt in d_logn_mts]
d_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in d_logn_mts]
d_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in d_logn_mts]

v_logn_mts=[fit(MarginalTDist,log.(measure_dict["Ventral CMZ (#)"][n])) for n in 1:length(X)]
v_mean=[exp(mean(mt)) for mt in v_logn_mts]
v_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in v_logn_mts]
v_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in v_logn_mts]

r_n_mts=[fit(MarginalTDist,measure_dict["Ratio"][n]) for n in 1:length(X)]
r_mean=[mean(mt) for mt in r_n_mts]
r_lower=[mean(mt)-quantile(mt,.025) for mt in r_n_mts]
r_upper=[quantile(mt,.975)-mean(mt) for mt in r_n_mts]

plotticks=[30*i for i in 1:12]

dv_chart=Plots.plot(X, d_mean, ribbon=(d_lower,d_upper), color=:green, label="Dorsal mean", ylabel="CMZ sectional population size", showaxis=:y, xticks=plotticks, xformatter=_->"")
plot!(X, v_mean, ribbon=(v_lower,v_upper), color=:darkmagenta, label="Ventral mean")
scatter!(vcat([[X[n] for i in 1:length(measure_dict["Dorsal CMZ (#)"][n])] for n in 1:length(X)]...),vcat(measure_dict["Dorsal CMZ (#)"]...), marker=:utriangle, color=:green, markersize=3, markerstrokecolor=:green, label="Dorsal data")
scatter!(vcat([[X[n] for i in 1:length(measure_dict["Ventral CMZ (#)"][n])] for n in 1:length(X)]...),vcat(measure_dict["Ventral CMZ (#)"]...), marker=:dtriangle, color=:darkmagenta, markersize=3, markerstrokecolor=:darkmagenta, label="Ventral data")
annotate!([(1.5,125,Plots.text("A",18))])
lens!([2, 31], [10, 150], inset = (1, bbox(0.3, 0.1, 0.45, 0.6)))

ratio_chart=scatter(vcat([[X[n] for i in 1:length(measure_dict["Ratio"][n])] for n in 1:length(X)]...),vcat(measure_dict["Ratio"]...), marker=:cross, color=:black, markersize=3, markerstrokecolor=:black, label="Ratio data", ylabel="Individual asymmetry ratio", xticks=plotticks)
plot!(X, r_mean, ribbon=(r_lower,r_upper), color=:green, label="Ratio mean")
plot!(X, [1. for i in 1:length(X)], style=:dot, color=:black, label="Even ratio")
annotate!([(1.5,.5,Plots.text("B",18))])


l=@layout [a{.6h}
            b]

combined=Plots.plot(dv_chart,ratio_chart,layout=l, link=:x, size=(900,600))

savefig(combined, "/bench/PhD/Thesis/images/cmz/DVontology.png")