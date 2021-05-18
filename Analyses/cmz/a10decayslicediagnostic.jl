using BayesianLinearRegression,CSV,DataFrames,Distributions,CMZNicheSims,GMC_NS, Serialization, Plots

emultipath="/bench/PhD/NGS_binaries/CNS/A10/emulti_decay"

default(legendfont = (10), guidefont = (12), tickfont = (10))

base_scale=[144.,5.]
base_min=[10.,0.]
time_scale=90
time_min=4.

prior=[LogNormal(log(20),log(2)),Uniform(eps(),.03),LogNormal(log(.9),log(1.6))]


em=deserialize(emultipath*"/ens")

models=["132.93","674.18","575.46","745.16"]
MLE=deserialize(em.path*"/132.93").log_Li
X=em.constants[1]
plotticks=[10*i for i in 1:9]

mapplots=Vector{Plots.Plot}()
mapxs=Vector{Float64}()

for model in models
    mapm=deserialize(em.path*"/$model")

    catdobs=vcat([em.obs[1][t] for t in 1:length(X)]...)
    catvobs=vcat([em.obs[2][t] for t in 1:length(X)]...)

    dXs=vcat([[X[n] for i in 1:length(em.obs[1][n])] for n in 1:length(X)]...)
    vXs=vcat([[X[n] for i in 1:length(em.obs[2][n])] for n in 1:length(X)]...)

    md_mean=mapm.slices[1].disp_mat[:,2]
    md_upper=mapm.slices[1].disp_mat[:,3].-mapm.slices[1].disp_mat[:,2]
    md_lower=mapm.slices[1].disp_mat[:,2].-mapm.slices[1].disp_mat[:,1]

    mv_mean=mapm.slices[2].disp_mat[:,2]
    mv_upper=mapm.slices[2].disp_mat[:,3].-mapm.slices[2].disp_mat[:,2]
    mv_lower=mapm.slices[2].disp_mat[:,2].-mapm.slices[2].disp_mat[:,1]

    mapmd_plt=scatter(dXs,catdobs, marker=:cross, color=:black, markersize=3, showaxis=:y, xformatter=_->"", xticks=X, legend=nothing, title=model, ylabel="Population")
    plot!(mapmd_plt, X, md_mean, ribbon=(md_lower,md_upper), color=:green)

    mapmv_plt=scatter(vXs,catvobs, marker=:cross, color=:black, markersize=3, ylabel="Population", xticks=X, xlabel="Age (dpf)", legend=nothing)
    plot!(mapmv_plt, X, mv_mean, ribbon=(mv_lower,mv_upper), color=:darkmagenta)

    combined_map=Plots.plot(mapmd_plt,mapmv_plt,layout=grid(2,1), size=(300,600),  link=:xy)

    push!(mapplots, combined_map)
    push!(mapxs,mapm.θ[1])
end

mkdes=posterior_kde(em, bivar=[1=>2,1=>3,2=>3])

ctkmarg=contour(mkdes[4].x,mkdes[4].y,transpose(mkdes[4].density), c=:viridis, xlabel="iCT (hr)", ylabel="κ", colorbar=:false, xlims=[0,144], ylims=[0,.03], levels=200, xticks=24.:24.:144.)

ctemarg=contour(mkdes[5].x,mkdes[5].y,transpose(mkdes[5].density), c=:viridis, xlabel="iCT (hr)", ylabel="ϵ", colorbar=:false, xlims=[0,144], ylims=[0,5.], levels=200, xticks=24.:24.:144.)

kemarg=contour(mkdes[6].x,mkdes[6].y,transpose(mkdes[6].density), c=:viridis, xlabel="κ", ylabel="ϵ", colorbar=:false, xlims=[0,.03], ylims=[0,5.], levels=200, xticks=24.:24.:144.)

# clim=(0,.5), legend=(.8,.3)

ymarg=map(x->pdf(prior[1][1],x),mkdes[1].x)

ctmarg=plot(mkdes[1].x, ymarg, color=:darkmagenta, fill=true, fillalpha=.5,xlims=[0,144], xlabel="Initial Cycle Time (hr)", ylabel="Probability density", label="Prior")
plot!(mkdes[1].x,mkdes[1].density,color=:green, fill=true, fillalpha=.5, label="Combined  Posterior")

colors=[:darkgreen,:darkred,:darkblue,:darkmagenta]
for (n, model) in enumerate(models)
    plot!([mapxs[n],mapxs[n]],[0,maximum(mkdes[1].density)],color=colors[n], label=model)
end

combined_layout=@layout[ctmarg{.3h}
                        grid(1,4)
                        ]

pmplot=Plots.plot(ctmarg,mapplots...,layout=combined_layout,size=(1200,900))

savefig(pmplot, "/bench/PhD/Defense/images/polymodality.png")