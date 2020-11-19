using CSV,DataFrames,Distributions,NGRefTools,StatsBase,Plots, Measurements,GMC_NS
gr()

a19_1pth="/bench/PhD/datasets/A19GR1.csv"
a19_2pth="/bench/PhD/datasets/A19GR2.csv"
a19_3pth="/bench/PhD/datasets/A19GR3.csv"

gr1df=DataFrame(CSV.read(a19_1pth))
gr2df=DataFrame(CSV.read(a19_2pth))
gr3df=DataFrame(CSV.read(a19_3pth))

X=sort(unique(gr1df."Pulse age (dpf)"))
measure_dict=Dict{String,Vector{Vector{Float64}}}()
measure_dict["ONL"]=[zeros(0) for i in 1:length(X)]
measure_dict["INL"]=[zeros(0) for i in 1:length(X)]
measure_dict["GCL"]=[zeros(0) for i in 1:length(X)]
measure_dict["Isl"]=[zeros(0) for i in 1:length(X)]
measure_dict["Pax6"]=[zeros(0) for i in 1:length(X)]
measure_dict["Isl_Pax6"]=[zeros(0) for i in 1:length(X)]
measure_dict["INL_Pax6"]=[zeros(0) for i in 1:length(X)]
measure_dict["PKCB"]=[zeros(0) for i in 1:length(X)]
measure_dict["GS"]=[zeros(0) for i in 1:length(X)]
measure_dict["HM"]=[zeros(0) for i in 1:length(X)]
measure_dict["Zpr1"]=[zeros(0) for i in 1:length(X)]

function check_layerfrac(row)
    !ismissing(row."#GCL") && !ismissing(row."#INL") && !ismissing(row."#PRL") && !ismissing(row."Total") ? (return true) : (return false)
end

function gr1_dfcheck(row)
    !ismissing(row."GCL - Isl") && !ismissing(row."GCL - Pax6") && !ismissing(row."GCL - Isl/Pax6") && !ismissing(row."INL - Pax6") ? (return true) : (return false)
end

for df in [gr1df,gr2df,gr3df]
    for row in eachrow(df)
        xidx=findfirst(i->i==row."Pulse age (dpf)",X)
        if df === gr1df
            if check_layerfrac(row)
                push!(measure_dict["GCL"][xidx],row."#GCL"/row."Total")
                push!(measure_dict["INL"][xidx],row."#INL"/row."Total")
                push!(measure_dict["ONL"][xidx],row."#PRL"/row."Total")
            end
            if gr1_dfcheck(row)
                push!(measure_dict["Isl"][xidx],row."GCL - Isl"/row."#GCL")
                push!(measure_dict["Pax6"][xidx],row."GCL - Pax6"/row."#GCL")
                push!(measure_dict["Isl_Pax6"][xidx],row."GCL - Isl/Pax6"/row."#GCL")
                push!(measure_dict["INL_Pax6"][xidx],row."INL - Pax6"/row."#INL")
            end
        else
            push!(measure_dict["GCL"][xidx],row."#GCL"/row."Total")
            push!(measure_dict["INL"][xidx],row."#INL"/row."Total")
            push!(measure_dict["ONL"][xidx],row."#PRL"/row."Total")
            if df === gr2df
                push!(measure_dict["PKCB"][xidx],row."INL-PKCB"/row."#INL")
                push!(measure_dict["GS"][xidx],row."INL-GS"/row."#INL")
                push!(measure_dict["HM"][xidx],row."INL-HM"/row."#INL")
            elseif df === gr3df
                push!(measure_dict["Zpr1"][xidx],row."PRL - Zpr1"/row."#PRL")
            end
        end
    end
end

for ms in keys(measure_dict)
    normal_mts=[fit(Normal,filter!(i->!iszero(i)&&!isnan(i),measure_dict[ms][n])) for n in 1:length(X)]
    logn_mts=[fit(LogNormal,filter!(i->!iszero(i)&&!isnan(i),measure_dict[ms][n])) for n in 1:length(X)]

    normal_lh=sum([sum(logpdf(normal_mts[n],measure_dict[ms][n])) for n in 1:length(X)])
    logn_lh=sum([sum(logpdf(logn_mts[n],measure_dict[ms][n])) for n in 1:length(X)])

    println("$ms & $(round(normal_lh,digits=3)) & $(round(logn_lh,digits=3)) & $(round(logn_lh-normal_lh,digits=3))\\\\ \\hline")
end

gmc=GMC_DEFAULTS
gmc[1]=5
gmc[2]=1.e-15

uds=Vector{Vector{Function}}([[convergence_display],[evidence_display],[info_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display]])

n_models=50

evdict=Dict{String,Measurement}()

max_μ=1.
min_μ=0.
max_λ=4000.
min_λ=1.

prior=[Uniform(min_μ,max_μ),Uniform(min_λ,max_λ)]
box=[min_μ max_μ;min_λ max_λ]

for ms in keys(measure_dict)
    for (n,ovec) in enumerate(measure_dict[ms])
        e_basepth="/bench/PhD/NGS_binaries/BSS/A19/n_"*ms
        x=X[n]
        enspth=e_basepth*"/$x"
        if isfile(enspth*"/ens")
            e=deserialize(enspth*"/ens")
        else
            e=Normal_Ensemble(enspth,n_models,filter(i->i>0,ovec), prior, box, gmc...)
        end

        "n_"*ms in keys(evdict) ? (evdict["n_"*ms] += converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)) :  (evdict["n_"*ms] = converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6))
    end

    for (n,ovec) in enumerate(measure_dict[ms])
        e_basepth="/bench/PhD/NGS_binaries/BSS/A19/ln_"*ms
        x=X[n]
        enspth=e_basepth*"/$x"
        if isfile(enspth*"/ens")
            e=deserialize(enspth*"/ens")
        else
            e=LogNormal_Ensemble(enspth,n_models,filter(i->i>0,ovec), prior, box, gmc...)
        end

        "ln_"*ms in keys(evdict) ? (evdict["ln_"*ms] += converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)) :  (evdict["ln_"*ms] = converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6))
    end
end

for ms in keys(measure_dict)
    nev=evdict["n_"*ms]
    lnev=evdict["ln_"*ms]
    ratio=lnev-nev
    println("$ms & $(round(nev,digits=3)) & $(round(lnev,digits=3)) & $(round(ratio,digits=3)) & $(ratio.val/ratio.err)\\\\ \\hline")
end

for ms in keys(measure_dict)
    if ms in ["INL_Pax6","GCL","INL","Isl","Isl_Pax6","ONL","Zpr1"]
        efunc=Normal_Ensemble
        pf="cn_"
    else
        efunc=LogNormal_Ensemble
        pf="cln_"
    end

    ovec=vcat(measure_dict[ms]...)
    e_basepth="/bench/PhD/NGS_binaries/BSS/A19/"*pf*ms

    if isfile(e_basepth*"/ens")
        e=deserialize(e_basepth*"/ens")
    else
        e=efunc(e_basepth,n_models,filter(i->i>0,ovec), prior, box, gmc...)
    end

    evdict[pf*ms] = converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)
end

for ms in keys(measure_dict)
    if ms in ["INL_Pax6","GCL","INL","Isl","Isl_Pax6","ONL","Zpr1"]
        efunc=Normal_Ensemble
        pf="n_"
    else
        efunc=LogNormal_Ensemble
        pf="ln_"
    end

    cev=evdict["c"*pf*ms]
    sev=evdict[pf*ms]
    ratio=cev-sev
    println("$ms & $(round(cev,digits=3)) & $(round(sev,digits=3)) & $(round(ratio,digits=3)) & $(ratio.val/ratio.err)\\\\ \\hline")
end
