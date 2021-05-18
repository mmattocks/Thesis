using BayesianLinearRegression,CSV,DataFrames,Distributions, GMC_NS, Serialization, Plots, NGRefTools
a20pth="/bench/PhD/Thesis/datasets/a20.csv"

a20df=DataFrame(CSV.read(a20pth))
a20df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a20df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["Nasal CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Temporal CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Ratio"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a20df, "Time point (d)")
    push!(X,Float64(unique(t_df."Time point (d)")[1]))
    rvec,nvec,tvec = [zeros(0) for i in 1:3]

    for i_df in groupby(t_df,"BSR")
        length([skipmissing(i_df."Nasal CMZ (#)")...])>0 && mean(skipmissing(i_df."Nasal CMZ (#)")) >0. && push!(nvec,mean(skipmissing(i_df."Nasal CMZ (#)")))
        length([skipmissing(i_df."Temporal CMZ (#)")...])>0 && mean(skipmissing(i_df."Temporal CMZ (#)")) >0. && push!(tvec,mean(skipmissing(i_df."Temporal CMZ (#)")))
        length([skipmissing(i_df."Nasal CMZ (#)")...])>0 && mean(skipmissing(i_df."Nasal CMZ (#)")) >0. && length([skipmissing(i_df."Temporal CMZ (#)")...])>0 && mean(skipmissing(i_df."Temporal CMZ (#)")) >0. && push!(rvec,mean(skipmissing(i_df."Nasal CMZ (#)"))/mean(skipmissing(i_df."Temporal CMZ (#)")))
    end

    push!(measure_dict["Ratio"],rvec)
    push!(measure_dict["Nasal CMZ (#)"],nvec)
    push!(measure_dict["Temporal CMZ (#)"],tvec)
end

n_logn_mts=[fit(MarginalTDist,log.(measure_dict["Nasal CMZ (#)"][n])) for n in 1:length(X)]
n_mean=[exp(mean(mt)) for mt in n_logn_mts]
n_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in n_logn_mts]
n_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in n_logn_mts]

t_logn_mts=[fit(MarginalTDist,log.(measure_dict["Temporal CMZ (#)"][n])) for n in 1:length(X)]
t_mean=[exp(mean(mt)) for mt in t_logn_mts]
t_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in t_logn_mts]
t_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in t_logn_mts]

r_n_mts=[fit(MarginalTDist,measure_dict["Ratio"][n]) for n in 1:length(X)]
r_mean=[mean(mt) for mt in r_n_mts]
r_lower=[mean(mt)-quantile(mt,.025) for mt in r_n_mts]
r_upper=[quantile(mt,.975)-mean(mt) for mt in r_n_mts]

plotticks=[30*i for i in 1:12]

nt_chart=Plots.plot(X, n_mean, ribbon=(n_lower,n_upper), color=:green, label="Nasal mean", ylabel="CMZ sectional population size", showaxis=:y, xticks=plotticks, xformatter=_->"")
plot!(X, t_mean, ribbon=(t_lower,t_upper), color=:darkmagenta, label="Temporal mean")
scatter!(vcat([[X[n] for i in 1:length(measure_dict["Nasal CMZ (#)"][n])] for n in 1:length(X)]...),vcat(measure_dict["Nasal CMZ (#)"]...), marker=:utriangle, color=:green, markersize=3, markerstrokecolor=:green, label="Nasal data")
scatter!(vcat([[X[n] for i in 1:length(measure_dict["Temporal CMZ (#)"][n])] for n in 1:length(X)]...),vcat(measure_dict["Temporal CMZ (#)"]...), marker=:dtriangle, color=:darkmagenta, markersize=3, markerstrokecolor=:darkmagenta, label="Temporal data")
annotate!([(1.5,125,Plots.text("A",18))])

ratio_chart=scatter(vcat([[X[n] for i in 1:length(measure_dict["Ratio"][n])] for n in 1:length(X)]...),vcat(measure_dict["Ratio"]...), marker=:cross, color=:black, markersize=3, markerstrokecolor=:black, label="Ratio data", ylabel="Individual asymmetry ratio", xlabel="Age (dpf)", xticks=plotticks)
plot!(X, r_mean, ribbon=(r_lower,r_upper), color=:green, label="Ratio mean")
plot!(X, [1. for i in 1:length(X)], style=:dot, color=:black, label="Even ratio")
annotate!([(10,.5,Plots.text("B",18))])


l=@layout [a{.6h}
            b]

combined=Plots.plot(nt_chart,ratio_chart,layout=l, link=:x, size=(900,600))

savefig(combined, "/bench/PhD/Thesis/images/cmz/NTontology.png")