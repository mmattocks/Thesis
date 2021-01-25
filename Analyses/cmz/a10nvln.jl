using CSV,DataFrames,Distributions,NGRefTools,StatsBase,Plots,GMC_NS,Serialization, Measurements

a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"
pne_pth="/bench/PhD/NGS_binaries/GMC_NS/A10/pop_norm"
plne_pth="/bench/PhD/NGS_binaries/GMC_NS/A10/pop_lognorm"
vne_pth="/bench/PhD/NGS_binaries/GMC_NS/A10/vol_norm"
vlne_pth="/bench/PhD/NGS_binaries/GMC_NS/A10/vol_lognorm"

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["PopEst"]=Vector{Vector{Float64}}()
measure_dict["VolEst"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    push!(X,Float64(unique(t_df."Time point (d)")[1]))
    pevec,ldvec,spvec = [zeros(0) for i in 1:3]

    for i_df in groupby(t_df,"BSR")
        length([skipmissing(i_df."Lens diameter")...])>0 && push!(ldvec,mean(skipmissing(i_df."Lens diameter")))
    end

    for i_df in groupby(t_df,"BSR")
        if length([skipmissing(i_df."CMZ Sum")...])>0 && mean(skipmissing(i_df."CMZ Sum")) > 0
            if length([skipmissing(i_df."Lens diameter")...])>0
                push!(pevec,mean(skipmissing(i_df."CMZ Sum"))*mean(skipmissing(i_df."Lens diameter"))/14.)
            else
                push!(pevec,mean(skipmissing(i_df."CMZ Sum"))*mean(ldvec)/14.)
            end
        end
        if length([skipmissing(i_df."Retina thickness")...])>0 &&
         length([skipmissing(i_df."IPL")...])>0 &&
         length([skipmissing(i_df."OPL")...])>0 &&
         length([skipmissing(i_df."RPE length")...])>0

            rthi=mean(skipmissing(i_df."Retina thickness"))-mean(skipmissing(i_df."OPL"))-mean(skipmissing(i_df."IPL"))
            rpel=mean(skipmissing(i_df."RPE length"))
            if length([skipmissing(i_df."Lens diameter")...])>0
                rcirc=rpel+mean(skipmissing(i_df."Lens diameter"))
            else
                rcirc=rpel+mean(ldvec)
            end

            or=rcirc/2π
            ir=or-.5*rthi

            ov=(4/3)*π*(or^3)
            iv=(4/3)*π*(ir^3)

            push!(spvec,(ov-iv)*(4/5))
        end
    end

    push!(measure_dict["PopEst"],pevec)
    push!(measure_dict["VolEst"],spvec)
end

gmc=GMC_DEFAULTS
gmc[1]=5
gmc[2]=1.e-15

uds=Vector{Vector{Function}}([[convergence_display],[evidence_display],[info_display]])

lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display]])

evdict=Dict{String,Measurement}()

pmax_μ=2e5
pmin_μ=1
pmax_λ=1e-3
pmin_λ=1e-10
vmax_μ=1e10
vmin_μ=100
vmax_λ=1e-10
vmin_λ=1e-20

pn_priors=[Uniform(pmin_μ,pmax_μ),Uniform(pmin_λ,pmax_λ)]
pn_box=[pmin_μ pmax_μ;pmin_λ pmax_λ]
vn_priors=[Uniform(vmin_μ,vmax_μ),Uniform(vmin_λ,vmax_λ)]
vn_box=[vmin_μ vmax_μ;vmin_λ vmax_λ]

n_models=50

for (pth,prior,box,eststring) in zip([pne_pth, vne_pth],[pn_priors, vn_priors],[pn_box,vn_box],["PopEst","VolEst"])
    for (nx,x) in enumerate(X)
        enspth=pth*"/$x"
        if isfile(enspth*"/ens")
            e=deserialize(enspth*"/ens")
        else
            e=Normal_Ensemble(enspth,n_models,filter(i->!iszero(i), measure_dict[eststring][nx]), prior, box, gmc...)
        end

        pth in keys(evdict) ? (evdict[pth] += converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)) :  (evdict[pth] = converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6))
    end
end

ln_pmax_μ=log(pmax_μ^2/sqrt(pmax_μ^2 + inv(pmax_λ)))
ln_pmin_μ=log(pmin_μ^2/sqrt(pmin_μ^2 + inv(pmin_λ)))
ln_pmax_λ=1/log(1+(inv(pmax_λ)/pmax_μ^2))
ln_pmin_λ=1/log(1+(inv(pmin_λ)/pmin_μ^2))

ln_vmax_μ=log(vmax_μ^2/sqrt(vmax_μ^2 + inv(vmax_λ)))
ln_vmin_μ=log(vmin_μ^2/sqrt(vmin_μ^2 + inv(vmin_λ)))
ln_vmax_λ=1/log(1+(inv(vmax_λ)/vmax_μ^2))
ln_vmin_λ=1/log(1+(inv(vmin_λ)/vmin_μ^2))

pln_priors=[Uniform(ln_pmin_μ,ln_pmax_μ),Uniform(ln_pmin_λ,ln_pmax_λ)]
pln_box=[ln_pmin_μ ln_pmax_μ;ln_pmin_λ ln_pmax_λ]
vln_priors=[Uniform(ln_vmin_μ,ln_vmax_μ),Uniform(ln_vmin_λ,ln_vmax_λ)]
vln_box=[ln_vmin_μ ln_vmax_μ;ln_vmin_λ ln_vmax_λ]

for (pth,prior,box,eststring) in zip([plne_pth, vlne_pth],[pn_priors, vn_priors],[pn_box,vn_box],["PopEst","VolEst"])
    for (nx,x) in enumerate(X)
        enspth=pth*"/$x"
        if isfile(enspth*"/ens")
            e=deserialize(enspth*"/ens")
        else
            e=LogNormal_Ensemble(enspth,n_models,filter(i->i>0,measure_dict[eststring][nx]), prior, box, gmc...)
        end

        pth in keys(evdict) ? (evdict[pth] += converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)) :  (evdict[pth] = converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6))
    end
end

pop_evidence_ratio=evdict[plne_pth]-evdict[pne_pth]
vol_evidence_ratio=evdict[vlne_pth]-evdict[vne_pth]

println("CMZ Population & $(round(evdict[pne_pth],digits=3)) & $(round(evdict[plne_pth],digits=3)) & $(round(pop_evidence_ratio,digits=3)) & $(pop_evidence_ratio.val/pop_evidence_ratio.err)\\\\ \\hline")
println("Estimated Retinal Volume & $(round(evdict[vne_pth],digits=3)) & $(round(evdict[vlne_pth],digits=3)) & $(round(vol_evidence_ratio,digits=3)) & $(vol_evidence_ratio.val/vol_evidence_ratio.err)\\\\ \\hline")