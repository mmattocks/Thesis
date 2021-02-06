using Distributions, Serialization, CSV, BayesianLinearRegression, Plots, DataFrames, Measurements, GMC_NS
import Measurements:value,uncertainty

a27path="/bench/PhD/datasets/a27 r.csv"
a27df=DataFrame(CSV.File(a27path))

X=unique(a27df["Age"]).*30.
CR=Dict{Float64,Vector{Float64}}()
mo1=Dict{Float64,Vector{Float64}}()

for x in X
    CR[x]=Vector{Float64}()
    x<90. && (mo1[x]=Vector{Float64}())
end

for row in eachrow(a27df)
    x=row."Age"*30.
    !ismissing(row."LR") && push!(CR[x],row."LR")
    !ismissing(row."D/N1") && !ismissing(row."V/T1") && push!(mo1[x],sum([row."D/N1",row."V/T1"]))
end


gmc=GMC_DEFAULTS
gmc[1]=5
gmc[2]=1.e-15

uds=Vector{Vector{Function}}([[convergence_display],[evidence_display],[info_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display]])

n_models=50

evdict=Dict{String,Measurement}()

max_μ=2000.
min_μ=0.
max_λ=500.
min_λ=1e-6

prior=[Uniform(min_μ,max_μ),Uniform(min_λ,max_λ)]
box=[min_μ max_μ;min_λ max_λ]

cobs=vcat([v for v in values(CR)]...)
mobs=vcat([v for v in values(mo1)]...)

for (md, pstring) in zip([CR,mo1],["CR","mo1"])
    for (k,v) in md
        pfx=pstring*string(k)
        enspth="/bench/PhD/NGS_binaries/GMC_NS/A27/"*pfx
        if isfile(enspth*"/ens")
            e=deserialize(enspth*"/ens")
        else
            e=LogNormal_Ensemble(enspth,n_models,v, prior, box, gmc...)
        end

        pstring in keys(evdict) ? (evdict[pstring]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_factor=1e-6)) : (evdict[pstring]=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_factor=1e-6))
    end
end

for (obs,pstring) in zip([cobs,mobs],["cCR","cmo1"])
    enspth="/bench/PhD/NGS_binaries/GMC_NS/A27/"*pstring
    if isfile(enspth*"/ens")
        e=deserialize(enspth*"/ens")
    else
        e=LogNormal_Ensemble(enspth,n_models,obs, prior, box, gmc...)
    end

    evdict[pstring]=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_factor=1e-6)
end

for ms in ["CR","mo1"]
    cev=evdict["c"*ms]
    sev=evdict[ms]
    ratio=cev-sev
    println("$ms & $(round(cev,digits=3)) & $(round(sev,digits=3)) & $(round(ratio,digits=3)) & $(round(ratio.val/ratio.err,digits=1))\\\\ \\hline")
end
