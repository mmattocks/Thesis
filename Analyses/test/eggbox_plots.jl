using GMC_NS, Distributions, Plots

gr()
default(legendfont = (10), guidefont = (12), tickfont = (10))

prior=Uniform(0,1)
priors=[prior,prior]

box=[-1. 1.
     -1. 1.]

gmc=copy(GMC_DEFAULTS)

e=Eggbox_Ensemble("Eggbox_test", 400, priors, box, gmc...)

x1=0.:.01:1.
x2=0.:.01:1.

f(x1, x2) = begin
        (2 + cos(5π * x1) * cos(5π * x2))^5
    end

function harvest_xs(e)
    x1m=[m.pos[1] for m in e.models]
    x2m=[m.pos[2] for m in e.models]
    x1m.+=1.; x2m.+=1.
    x1m./=2.; x2m./=2
    return x1m, x2m
end

plt=contour(x1,x2,f, xlabel="x1", ylabel="x2", colorbar=:none)

init=deepcopy(plt)
plot!(init, title="Initialisation")
scatter!(init, harvest_xs(e), marker=:cross, markercolor=:black, legend=:none)

converge_ensemble!(e,backup=(true,100), max_iterates=1000, converge_factor=1e-6)

thousand=deepcopy(plt)
plot!(thousand, title="1000 iterates")
scatter!(thousand, harvest_xs(e), marker=:cross, markercolor=:black, legend=:none)

converge_ensemble!(e,backup=(true,100), max_iterates=10000, converge_factor=1e-6)

tenthousand=deepcopy(plt)
plot!(tenthousand, title="10000 iterates")
scatter!(tenthousand, harvest_xs(e), marker=:cross, markercolor=:black, legend=:none)

converge_ensemble!(e,backup=(true,100), converge_factor=1e-12)

converged=deepcopy(plt)
plot!(converged, title="Converged ($(length(e.log_Li)) iterates)")
scatter!(converged, harvest_xs(e), marker=:cross, markercolor=:black, legend=:none)

cmb=plot(init, thousand, tenthousand, converged, layout=grid(2,2), size=(900,900))
savefig(cmb, "/bench/PhD/Thesis/images/GMC/eggbox.png")