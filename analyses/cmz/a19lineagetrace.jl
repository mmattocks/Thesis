using CSV,DataFrames,Distributions,NGRefTools,StatsBase,Plots, Measurements,GMC_NS,Serialization
gr()

a19_1pth="/bench/PhD/datasets/A19GR1.csv"
a19_2pth="/bench/PhD/datasets/A19GR2.csv"
a19_3pth="/bench/PhD/datasets/A19GR3.csv"

gr1df=DataFrame(CSV.read(a19_1pth))
gr2df=DataFrame(CSV.read(a19_2pth))
gr3df=DataFrame(CSV.read(a19_3pth))

X=sort(unique(gr1df."Pulse age (dpf)"))
t_measure_dict=Dict{String,Vector{Vector{Float64}}}()
d_measure_dict=Dict{String,Vector{Vector{Float64}}}()
v_measure_dict=Dict{String,Vector{Vector{Float64}}}()

mds=[t_measure_dict,d_measure_dict,v_measure_dict]

for md in mds
    md["ONL"]=[zeros(0) for i in 1:length(X)]
    md["INL"]=[zeros(0) for i in 1:length(X)]
    md["GCL"]=[zeros(0) for i in 1:length(X)]
    md["Isl"]=[zeros(0) for i in 1:length(X)]
    md["Pax6"]=[zeros(0) for i in 1:length(X)]
    md["Isl_Pax6"]=[zeros(0) for i in 1:length(X)]
    md["INL_Pax6"]=[zeros(0) for i in 1:length(X)]
    md["PKCB"]=[zeros(0) for i in 1:length(X)]
    md["GS"]=[zeros(0) for i in 1:length(X)]
    md["HM"]=[zeros(0) for i in 1:length(X)]
    md["Zpr1"]=[zeros(0) for i in 1:length(X)]
end

function check_layerfrac(row)
    !ismissing(row."#GCL") && !ismissing(row."#INL") && !ismissing(row."#PRL") && !ismissing(row."Total") ? (return true) : (return false)
end

function gr1_dfcheck(subdf)
    pass=falses(size(subdf,1))
    for (n,row) in enumerate(eachrow(subdf))    
        !ismissing(row."GCL - Isl") && !ismissing(row."GCL - Pax6") && !ismissing(row."GCL - Isl/Pax6") && !ismissing(row."INL - Pax6") ? (pass[n]=true) : (pass[n]=false)
    end
    all(pass) ? (return true) : (return false)
end

for df in [gr1df,gr2df,gr3df]
    for agedf in groupby(df, "Pulse age (dpf)")
        xidx=findfirst(i->i==unique(agedf."Pulse age (dpf)")[1],X)
        for inddf in groupby(agedf, "Individual")
            if df === gr1df
                if gr1_dfcheck(inddf)
                    clean=falses(2)
                    for row in eachrow(inddf)
                        if row."D/V" == "D"
                            push!(d_measure_dict["Isl"][xidx],row."GCL - Isl"/row."#GCL")
                            push!(d_measure_dict["Pax6"][xidx],row."GCL - Pax6"/row."#GCL")
                            push!(d_measure_dict["Isl_Pax6"][xidx],row."GCL - Isl/Pax6"/row."#GCL")
                            push!(d_measure_dict["INL_Pax6"][xidx],row."INL - Pax6"/row."#INL")
                            if check_layerfrac(row)
                                push!(d_measure_dict["GCL"][xidx],row."#GCL"/row."Total")
                                push!(d_measure_dict["INL"][xidx],row."#INL"/row."Total")
                                push!(d_measure_dict["ONL"][xidx],row."#PRL"/row."Total")
                                clean[1]=true
                            end
                        else
                            push!(v_measure_dict["Isl"][xidx],row."GCL - Isl"/row."#GCL")
                            push!(v_measure_dict["Pax6"][xidx],row."GCL - Pax6"/row."#GCL")
                            push!(v_measure_dict["Isl_Pax6"][xidx],row."GCL - Isl/Pax6"/row."#GCL")
                            push!(v_measure_dict["INL_Pax6"][xidx],row."INL - Pax6"/row."#INL")
                            if check_layerfrac(row)
                                push!(v_measure_dict["GCL"][xidx],row."#GCL"/row."Total")
                                push!(v_measure_dict["INL"][xidx],row."#INL"/row."Total")
                                push!(v_measure_dict["ONL"][xidx],row."#PRL"/row."Total")
                                clean[2]=true
                            end
                        end
                    end
                    push!(t_measure_dict["Isl"][xidx],sum(inddf."GCL - Isl")/sum(inddf."#GCL"))
                    push!(t_measure_dict["Pax6"][xidx],sum(inddf."GCL - Pax6")/sum(inddf."#GCL"))
                    push!(t_measure_dict["Isl_Pax6"][xidx],sum(inddf."GCL - Isl/Pax6")/sum(inddf."#GCL"))
                    push!(t_measure_dict["INL_Pax6"][xidx],sum(inddf."INL - Pax6")/sum(inddf."#INL"))
                    if all(clean)
                        push!(t_measure_dict["GCL"][xidx],sum(inddf."#GCL")/sum(inddf."Total"))
                        push!(t_measure_dict["INL"][xidx],sum(inddf."#INL")/sum(inddf."Total"))
                        push!(t_measure_dict["ONL"][xidx],sum(inddf."#PRL")/sum(inddf."Total"))
                    end
                end
            elseif df === gr2df
                for row in eachrow(inddf)
                    if row."D/V" == "D"
                        push!(d_measure_dict["PKCB"][xidx],row."INL-PKCB"/row."#INL")
                        push!(d_measure_dict["GS"][xidx],row."INL-GS"/row."#INL")
                        push!(d_measure_dict["HM"][xidx],row."INL-HM"/row."#INL")
                        push!(d_measure_dict["GCL"][xidx],row."#GCL"/row."Total")
                        push!(d_measure_dict["INL"][xidx],row."#INL"/row."Total")
                        push!(d_measure_dict["ONL"][xidx],row."#PRL"/row."Total")
                    else
                        push!(v_measure_dict["PKCB"][xidx],row."INL-PKCB"/row."#INL")
                        push!(v_measure_dict["GS"][xidx],row."INL-GS"/row."#INL")
                        push!(v_measure_dict["HM"][xidx],row."INL-HM"/row."#INL")
                        push!(v_measure_dict["GCL"][xidx],row."#GCL"/row."Total")
                        push!(v_measure_dict["INL"][xidx],row."#INL"/row."Total")
                        push!(v_measure_dict["ONL"][xidx],row."#PRL"/row."Total")
                    end
                end
                push!(t_measure_dict["PKCB"][xidx],sum(inddf."INL-PKCB")/sum(inddf."#INL"))
                push!(t_measure_dict["GS"][xidx],sum(inddf."INL-GS")/sum(inddf."#INL"))
                push!(t_measure_dict["HM"][xidx],sum(inddf."INL-HM")/sum(inddf."#INL"))
                push!(t_measure_dict["GCL"][xidx],sum(inddf."#GCL")/sum(inddf."Total"))
                push!(t_measure_dict["INL"][xidx],sum(inddf."#INL")/sum(inddf."Total"))
                push!(t_measure_dict["ONL"][xidx],sum(inddf."#PRL")/sum(inddf."Total"))
            elseif df === gr3df
                for row in eachrow(inddf)
                    if row."D/V" == "D"
                        push!(d_measure_dict["Zpr1"][xidx],row."PRL - Zpr1"/row."#PRL")
                        push!(d_measure_dict["GCL"][xidx],row."#GCL"/row."Total")
                        push!(d_measure_dict["INL"][xidx],row."#INL"/row."Total")
                        push!(d_measure_dict["ONL"][xidx],row."#PRL"/row."Total")
                    else
                        push!(v_measure_dict["Zpr1"][xidx],row."PRL - Zpr1"/row."#PRL")
                        push!(v_measure_dict["GCL"][xidx],row."#GCL"/row."Total")
                        push!(v_measure_dict["INL"][xidx],row."#INL"/row."Total")
                        push!(v_measure_dict["ONL"][xidx],row."#PRL"/row."Total")
                    end
                end
                push!(t_measure_dict["Zpr1"][xidx],sum(inddf."PRL - Zpr1")/sum(inddf."#PRL"))
                push!(t_measure_dict["GCL"][xidx],sum(inddf."#GCL")/sum(inddf."Total"))
                push!(t_measure_dict["INL"][xidx],sum(inddf."#INL")/sum(inddf."Total"))
                push!(t_measure_dict["ONL"][xidx],sum(inddf."#PRL")/sum(inddf."Total"))
            end
        end
    end
end

for d in mds
    for (k,v) in d
        for ovec in v
            if any(i->i==0.,ovec)
                ovec[findall(i->i==0.,ovec)].=1e-3
            end
        end
        d[k]=v
    end
end

onl_n_mts=[fit(MarginalTDist,t_measure_dict["ONL"][n]) for n in 1:length(X)]
onl_mean=[mean(mt) for mt in onl_n_mts]
onl_lower=[mean(mt)-quantile(mt,.025) for mt in onl_n_mts]
onl_upper=[quantile(mt,.975)-mean(mt) for mt in onl_n_mts]

inl_n_mts=[fit(MarginalTDist,t_measure_dict["INL"][n]) for n in 1:length(X)]
inl_mean=[mean(mt) for mt in inl_n_mts]
inl_lower=[mean(mt)-quantile(mt,.025) for mt in inl_n_mts]
inl_upper=[quantile(mt,.975)-mean(mt) for mt in inl_n_mts]

gcl_n_mts=[fit(MarginalTDist,t_measure_dict["GCL"][n]) for n in 1:length(X)]
gcl_mean=[mean(mt) for mt in gcl_n_mts]
gcl_lower=[mean(mt)-quantile(mt,.025) for mt in gcl_n_mts]
gcl_upper=[quantile(mt,.975)-mean(mt) for mt in gcl_n_mts]

layers_chart=scatter(vcat([[X[n] for i in 1:length(t_measure_dict["GCL"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["GCL"]...), marker=:utriangle, color=:magenta, markersize=3, label="GCL data", ylabel="Fraction of cohort in layer", xlabel="Age (dpf)", xticks=X, foreground_color_legend=nothing, background_color_legend=nothing, legend=:inside)
plot!(X, gcl_mean, ribbon=(gcl_lower,gcl_upper), color=:magenta, label="GCL mean")
scatter!(vcat([[X[n] for i in 1:length(t_measure_dict["INL"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["INL"]...), marker=:cross, color=:blue, markersize=3, label="INL data")
plot!(X, inl_mean,ribbon=(inl_lower,inl_upper), color=:blue, label="INL mean")
scatter!(vcat([[X[n] for i in 1:length(t_measure_dict["ONL"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["ONL"]...), marker=:dtriangle, color=:green, markersize=3, label="ONL data")
plot!(X, onl_mean,ribbon=(onl_lower,onl_upper), color=:green, label="ONL mean")
annotate!([(15,.7,Plots.text("A",18))])

isl_n_mts=[fit(MarginalTDist,t_measure_dict["Isl"][n]) for n in 1:length(X)]
isl_mean=[mean(mt) for mt in isl_n_mts]
isl_lower=[mean(mt)-quantile(mt,.025) for mt in isl_n_mts]
isl_upper=[quantile(mt,.975)-mean(mt) for mt in isl_n_mts]

islchart=scatter(vcat([[X[n] for i in 1:length(t_measure_dict["Isl"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["Isl"]...), marker=:square, color=:magenta, markersize=3, label="Isl+ data", xticks=X, foreground_color_legend=nothing, background_color_legend=nothing, legend=:bottom,  xformatter=_->"", showaxis=:y)
plot!(X, isl_mean,ribbon=(isl_lower,isl_upper), color=:magenta, label="Isl+ mean")
annotate!([(15,.05,Plots.text("B",18))])

pax6_ln_mts=[fit(MarginalTDist,log.(t_measure_dict["Pax6"][n])) for n in 1:length(X)]
pax6_mean=[exp(mean(mt)) for mt in pax6_ln_mts]
pax6_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in pax6_ln_mts]
pax6_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in pax6_ln_mts]

pax6chart=scatter(vcat([[X[n] for i in 1:length(t_measure_dict["Pax6"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["Pax6"]...), marker=:circle, color=:magenta, markersize=3, label="Pax6+ data", xticks=X, foreground_color_legend=nothing, background_color_legend=nothing, legend=:top,  xformatter=_->"", showaxis=:y, ylabel="Fraction of GCL cohort labelled")
plot!(X, pax6_mean,ribbon=(pax6_lower,pax6_upper), color=:magenta, label="Pax6+ mean", ylims=[0,.6])
annotate!([(15,.45,Plots.text("C(λ)",18))])

islp_n_mts=[fit(MarginalTDist,t_measure_dict["Isl_Pax6"][n]) for n in 1:length(X)]
islp_mean=[mean(mt) for mt in islp_n_mts]
islp_lower=[mean(mt)-quantile(mt,.025) for mt in islp_n_mts]
islp_upper=[quantile(mt,.975)-mean(mt) for mt in islp_n_mts]

islpchart=scatter(vcat([[X[n] for i in 1:length(t_measure_dict["Isl_Pax6"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["Isl_Pax6"]...), marker=:triangle, color=:magenta, markersize=3, label="Isl/Pax6+ data", xticks=X, foreground_color_legend=nothing, background_color_legend=nothing, legend=:bottom, xlabel="Age (dpf)")
plot!(X, islp_mean,ribbon=(islp_lower,islp_upper), color=:magenta, label="Isl/Pax6+ mean")
annotate!([(15,.05,Plots.text("D(†)",18))])

gcl_subplot=plot(islchart, pax6chart, islpchart, layout=grid(3,1), link=:x)

pax6i_n_mts=[fit(MarginalTDist,t_measure_dict["INL_Pax6"][n]) for n in 1:length(X)]
pax6i_mean=[mean(mt) for mt in pax6i_ln_mts]
pax6i_lower=[mean(mt)-quantile(mt,.025) for mt in pax6i_n_mts]
pax6i_upper=[quantile(mt,.975)-mean(mt) for mt in pax6i_n_mts]

inlpax6chart=scatter(vcat([[X[n] for i in 1:length(t_measure_dict["INL_Pax6"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["INL_Pax6"]...), marker=:square, color=:blue, markersize=3, label="Pax6+ data", xticks=X, foreground_color_legend=nothing, background_color_legend=nothing, legend=:bottom, xformatter=_->"", showaxis=:y)
plot!(X, pax6i_mean,ribbon=(pax6i_lower,pax6i_upper), color=:blue, label="Pax6+ mean")
annotate!([(15,.075,Plots.text("E",18))])

pkc_ln_mts=[fit(MarginalTDist,log.(t_measure_dict["PKCB"][n])) for n in 1:length(X)]
pkc_mean=[exp(mean(mt)) for mt in pkc_ln_mts]
pkc_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in pkc_ln_mts]
pkc_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in pkc_ln_mts]

pkcchart=scatter(vcat([[X[n] for i in 1:length(t_measure_dict["PKCB"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["PKCB"]...), marker=:circle, color=:blue, markersize=3, label="PKCβ+ data", xticks=X, foreground_color_legend=nothing, background_color_legend=nothing, legend=:top, xformatter=_->"", showaxis=:y, ylabel="Fraction of INL cohort labelled")
plot!(X, pkc_mean,ribbon=(pkc_lower,pkc_upper), color=:blue, label="PKCβ+ mean", ylims=[0,.2])
annotate!([(15,.18,Plots.text("F(λ†)",18))])

gs_ln_mts=[fit(MarginalTDist,log.(t_measure_dict["GS"][n])) for n in 1:length(X)]
gs_mean=[exp(mean(mt)) for mt in gs_ln_mts]
gs_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in gs_ln_mts]
gs_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in gs_ln_mts]

gschart=scatter(vcat([[X[n] for i in 1:length(t_measure_dict["GS"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["GS"]...), marker=:diamond, color=:blue, markersize=3, label="GS+ data", xticks=X, foreground_color_legend=nothing, background_color_legend=nothing, legend=:topright, xformatter=_->"", showaxis=:y)
plot!(X, gs_mean,ribbon=(gs_lower,gs_upper), color=:blue, label="GS+ mean")
annotate!([(15,.125,Plots.text("G(λ)",18))])

hm_ln_mts=[fit(MarginalTDist,log.(t_measure_dict["HM"][n])) for n in 1:length(X)]
hm_mean=[exp(mean(mt)) for mt in hm_ln_mts]
hm_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in hm_ln_mts]
hm_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in hm_ln_mts]

hmchart=scatter(vcat([[X[n] for i in 1:length(t_measure_dict["HM"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["HM"]...), marker=:cross, color=:blue, markersize=3, label="HM+ data", xticks=X, foreground_color_legend=nothing, background_color_legend=nothing, legend=:topright, xlabel="Age (dpf)")
plot!(X, hm_mean,ribbon=(hm_lower,hm_upper), color=:blue, label="HM+ mean")
annotate!([(15,.13,Plots.text("H(λ)",18))])

inl_subplot=plot(inlpax6chart, pkcchart, gschart, hmchart, layout=grid(4,1), link=:x)

z_ln_mts=[fit(MarginalTDist,log.(t_measure_dict["Zpr1"][n])) for n in 1:length(X)]
z_mean=[exp(mean(mt)) for mt in z_ln_mts]
z_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in z_ln_mts]
z_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in z_ln_mts]

zchart=scatter(vcat([[X[n] for i in 1:length(t_measure_dict["Zpr1"][n])] for n in 1:length(X)]...),vcat(t_measure_dict["Zpr1"]...), marker=:square, color=:blue, markersize=3, label="Zpr1+ data", xticks=X, foreground_color_legend=nothing, background_color_legend=nothing, legend=:bottom, xlabel="Age (dpf)", ylabel="Frac. OPL labelled")
plot!(X, z_mean,ribbon=(z_lower,z_upper), color=:blue, label="HM+ mean", ylims=[0,.55])
annotate!([(15,.15,Plots.text("I(λ†)",18))])

l=@layout [[layers{0.5h};gcl] [inl{0.8h};opl]]
combined=plot(layers_chart,gcl_subplot, inl_subplot,zchart, size=(1200,1200),layout=l)
savefig(combined, "/bench/PhD/Thesis/images/cmz/layercontributions.png")

for ms in keys(t_measure_dict)
    normal_mts=[fit(Normal,filter!(i->!isnan(i),t_measure_dict[ms][n])) for n in 1:length(X)]
    logn_mts=[fit(LogNormal,filter!(i->!isnan(i),t_measure_dict[ms][n])) for n in 1:length(X)]

    normal_lh=sum([sum(logpdf(normal_mts[n],t_measure_dict[ms][n])) for n in 1:length(X)])
    logn_lh=sum([sum(logpdf(logn_mts[n],t_measure_dict[ms][n])) for n in 1:length(X)])

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

for ms in keys(t_measure_dict)
    for (n,ovec) in enumerate(t_measure_dict[ms])
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

    for (n,ovec) in enumerate(t_measure_dict[ms])
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

for ms in keys(t_measure_dict)
    nev=evdict["n_"*ms]
    lnev=evdict["ln_"*ms]
    ratio=lnev-nev
    println("$ms & $(round(nev,digits=3)) & $(round(lnev,digits=3)) & $(round(ratio,digits=3)) & $(round(ratio.val/ratio.err, digits=1))\\\\ \\hline")
end

for ms in keys(t_measure_dict)
    if ms in ["INL_Pax6","INL","Isl","Isl_Pax6","ONL","Zpr1"]
        efunc=Normal_Ensemble
        pf="cn_"
    else
        efunc=LogNormal_Ensemble
        pf="cln_"
    end

    ovec=vcat(t_measure_dict[ms]...)
    e_basepth="/bench/PhD/NGS_binaries/BSS/A19/"*pf*ms

    if isfile(e_basepth*"/ens")
        e=deserialize(e_basepth*"/ens")
    else
        e=efunc(e_basepth,n_models,filter(i->i>0,ovec), prior, box, gmc...)
    end

    evdict[pf*ms] = converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)
end

for ms in keys(t_measure_dict)
    if ms in ["INL_Pax6","INL","Isl","Isl_Pax6","ONL","Zpr1"]
        pf="n_"
    else
        pf="ln_"
    end

    cev=evdict["c"*pf*ms]
    sev=evdict[pf*ms]
    ratio=cev-sev
    println("$ms & $(round(cev,digits=3)) & $(round(sev,digits=3)) & $(round(ratio,digits=3)) & $(round(ratio.val/ratio.err,digits=1))\\\\ \\hline")
end

dv_evdict=Dict{String,Measurement}()

for ms in keys(d_measure_dict)
    if ms in ["INL_Pax6","INL","Isl","Isl_Pax6","ONL","Zpr1"]
        efunc=Normal_Ensemble
        pf="d_cn_"
    else
        efunc=LogNormal_Ensemble
        pf="d_cln_"
    end

    ovec=vcat(d_measure_dict[ms]...)
    e_basepth="/bench/PhD/NGS_binaries/BSS/A19/"*pf*ms

    if isfile(e_basepth*"/ens")
        e=deserialize(e_basepth*"/ens")
    else
        e=efunc(e_basepth,n_models,filter(i->!isnan(i),ovec), prior, box, gmc...)
    end

    dv_evdict[pf*ms] = converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)
end

for ms in keys(v_measure_dict)
    if ms in ["INL_Pax6","INL","Isl","Isl_Pax6","ONL","Zpr1"]
        efunc=Normal_Ensemble
        pf="v_cn_"
    else
        efunc=LogNormal_Ensemble
        pf="v_cln_"
    end

    ovec=vcat(d_measure_dict[ms]...)
    e_basepth="/bench/PhD/NGS_binaries/BSS/A19/"*pf*ms

    if isfile(e_basepth*"/ens")
        e=deserialize(e_basepth*"/ens")
    else
        e=efunc(e_basepth,n_models,filter(i->!isnan(i),ovec), prior, box, gmc...)
    end

    dv_evdict[pf*ms] = converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)
end


for ms in keys(t_measure_dict)
    if ms in ["INL_Pax6","INL","Isl","Isl_Pax6","ONL","Zpr1"]
        pf="n_"
    else
        pf="ln_"
    end

    cev=evdict["c"*pf*ms]
    dev=dv_evdict["d_c"*pf*ms]
    vev=dv_evdict["v_c"*pf*ms]
    sev=dev+vev
    ratio=cev-sev
    println("$ms & $(round(cev,digits=3)) & $(round(sev,digits=3)) & $(round(ratio,digits=3)) & $(round(ratio.val/ratio.err,digits=1))\\\\ \\hline")
end