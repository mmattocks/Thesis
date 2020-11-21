using CSV,DataFrames,Distributions,NGRefTools,StatsBase,Plots
gr()
default(legendfont = (8), guidefont = (10), tickfont = (8))


a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["PopEst"]=Vector{Vector{Float64}}()
measure_dict["CMZ Sum"]=Vector{Vector{Float64}}()
measure_dict["Dorsal CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Ventral CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["VolEst"]=Vector{Vector{Float64}}()
measure_dict["Retina thickness"]=Vector{Vector{Float64}}()
measure_dict["RPE length"]=Vector{Vector{Float64}}()
measure_dict["Lens diameter"]=Vector{Vector{Float64}}()
measure_dict["ON diameter"]=Vector{Vector{Float64}}()
measure_dict["ONL"]=Vector{Vector{Float64}}()
measure_dict["OPL"]=Vector{Vector{Float64}}()
measure_dict["INL"]=Vector{Vector{Float64}}()
measure_dict["IPL"]=Vector{Vector{Float64}}()
measure_dict["GCL"]=Vector{Vector{Float64}}()
measure_dict["Pop/VolEst"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    push!(X,Float64(unique(t_df."Time point (d)")[1]))
    pevec,tvec,dvec,vvec,vevec, rtvec, rplvec, rptvec, ldvec, ondvec, onlvec, oplvec, inlvec, iplvec, gclvec, pvest = [zeros(0) for i in 1:16]

    for i_df in groupby(t_df,"BSR")
        length([skipmissing(i_df."Lens diameter")...])>0 && push!(ldvec,mean(skipmissing(i_df."Lens diameter")))
    end

    for i_df in groupby(t_df,"BSR")
        ipop=false
        ivol=false
        if length([skipmissing(i_df."CMZ Sum")...])>0
            push!(tvec,mean(skipmissing(i_df."CMZ Sum")))
            if length([skipmissing(i_df."Lens diameter")...])>0
                push!(pevec,mean(skipmissing(i_df."CMZ Sum"))*mean(skipmissing(i_df."Lens diameter"))/14.)
                ipop=true
            else
                push!(pevec,mean(skipmissing(i_df."CMZ Sum"))*mean(ldvec)/14.)
                ipop=true
            end
        end

        length([skipmissing(i_df."Dorsal CMZ (#)")...])>0 && push!(dvec,mean(skipmissing(i_df."Dorsal CMZ (#)")))
        length([skipmissing(i_df."Ventral CMZ (#)")...])>0 && push!(vvec,mean(skipmissing(i_df."Ventral CMZ (#)")))
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

            push!(vevec,(ov-iv)*(4/5))
        end
        length([skipmissing(i_df."Retina thickness")...])>0 && push!(rtvec,mean(skipmissing(i_df."Retina thickness")))
        length([skipmissing(i_df."RPE length")...])>0 && push!(rptvec,mean(skipmissing(i_df."RPE length")))
        length([skipmissing(i_df."ON diameter")...])>0 && push!(ondvec,mean(skipmissing(i_df."ON diameter")))
        length([skipmissing(i_df.ONL)...])>0 && push!(onlvec,mean(skipmissing(i_df.ONL)))   
        length([skipmissing(i_df.OPL)...])>0 && push!(oplvec,mean(skipmissing(i_df.OPL)))  
        length([skipmissing(i_df.INL)...])>0 && push!(inlvec,mean(skipmissing(i_df.INL)))
        length([skipmissing(i_df.IPL)...])>0 && push!(iplvec,mean(skipmissing(i_df.IPL)))  
        length([skipmissing(i_df.GCL)...])>0 && push!(gclvec,mean(skipmissing(i_df.GCL)))
        ivol && ipop && 
        push!(pvest,
        mean(skipmissing(i_df."CMZ Sum"))/
            ((mean(skipmissing(i_df."Retina thickness"))-mean(skipmissing(i_df."OPL"))-mean(skipmissing(i_df."IPL")))
                *mean(skipmissing(i_df."RPE length"))*(π/4)))
    end

    push!(measure_dict["PopEst"],pevec)
    push!(measure_dict["CMZ Sum"],tvec)
    push!(measure_dict["Dorsal CMZ (#)"],dvec)
    push!(measure_dict["Ventral CMZ (#)"],vvec)
    push!(measure_dict["VolEst"],vevec)
    push!(measure_dict["Retina thickness"],rtvec)
    push!(measure_dict["RPE length"],rptvec)
    push!(measure_dict["Lens diameter"],ldvec)
    push!(measure_dict["ON diameter"],ondvec)
    push!(measure_dict["ONL"],onlvec)
    push!(measure_dict["OPL"],oplvec)
    push!(measure_dict["INL"],inlvec)
    push!(measure_dict["IPL"],iplvec)
    push!(measure_dict["GCL"],gclvec)
    push!(measure_dict["Pop/VolEst"],pvest)
end

#TEST IF DATA IS BETTER MODELLED NORMALLY OR LOGNORMALLY
for ms in ["CMZ Sum", "Lens diameter", "PopEst", "RPE length", "Retina thickness", "VolEst"]
    normal_mts=[fit(Normal,measure_dict[ms][n]) for n in 1:length(X)]
    logn_mts=[fit(LogNormal,filter!(i->!iszero(i),measure_dict[ms][n])) for n in 1:length(X)]

    normal_lh=sum([sum(logpdf(normal_mts[n],measure_dict[ms][n])) for n in 1:length(X)])
    logn_lh=sum([sum(logpdf(logn_mts[n],measure_dict[ms][n])) for n in 1:length(X)])

    println("$ms & $(round(normal_lh,digits=3)) & $(round(logn_lh,digits=3)) & $(round(logn_lh-normal_lh,digits=3))\\\\ \\hline")

end

pe_n_mts=[fit(MarginalTDist,measure_dict["PopEst"][n]) for n in 1:length(X)]
pen_lower=[mean(mt)-quantile(mt,.025) for mt in pe_n_mts]
pen_upper=[quantile(mt,.975)-mean(mt) for mt in pe_n_mts]

pe_logn_mts=[fit(MarginalTDist,log.(filter!(i->!iszero(i),measure_dict["PopEst"][n]))) for n in 1:length(X)]
pe_mean=[exp(mean(mt)) for mt in pe_logn_mts]
pe_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in pe_logn_mts]
pe_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in pe_logn_mts]

plotticks=[30*i for i in 1:12]

ann_chart=scatter(vcat([[X[n] for i in 1:length(measure_dict["PopEst"][n])] for n in 1:length(X)]...),vcat(measure_dict["PopEst"]...), marker=:cross, color=:black, markersize=3, label="Pop. data", ylabel="Est. CMZ annulus population", showaxis=:y, xticks=plotticks, xformatter=_->"")
plot!(ann_chart, X, pe_mean, ribbon=(pe_lower,pe_upper), color=:magenta, label="Pop. mean")
lens!([2, 20], [500, 2500], inset = (1, bbox(0.3, 0.1, 0.45, 0.45)))
annotate!([(8,4500,Plots.text("A",18))])

poprate_mean,poprate_lower,poprate_upper=[Vector{Float64}() for i in 1:3]
for (nx, x) in enumerate(X)
    if nx > 1
        d=X[nx]-X[nx-1]
        l,m,u=MTDist_MC_func((a,b)->((exp(a)-exp(b))/d), [log.(filter!(i->!iszero(i),measure_dict["PopEst"][nx])), log.(filter!(i->!iszero(i),measure_dict["PopEst"][nx-1]))], summary=true)
        push!(poprate_mean,m);push!(poprate_lower,m-l);push!(poprate_upper,u-m)
    end
end

poprate_chart=plot(X[2:end],poprate_mean,ribbon=(poprate_lower, poprate_upper), color=:grey, ylabel="Est. CMZ pop. rate of change (cells/d)", label="Simulated mean rate", legend=:topright, showaxis=:y, xticks=plotticks, xformatter=_->"")
annotate!([(8,400,Plots.text("B",18))])


l3,m3,u3=MTDist_MC_func((a,b)->(exp(a)-exp(b))/5., [filter!(!isinf,log.(measure_dict["PopEst"][2])), filter!(!isinf,log.(measure_dict["PopEst"][1]))], summary=true)

ratedist23=MTDist_MC_func((a,b)->(exp(a)-exp(b))/6., [filter!(!isinf,log.(measure_dict["PopEst"][6])), filter!(!isinf,log.(measure_dict["PopEst"][5]))], summary=false)

println("3dpf poprate mean: $(m3) 23dpf poprate mean $(mean(ratedist23))poprate 23dpf quantile : $(quantile(ratedist23,.025))")

ve_logn_mts=[fit(MarginalTDist,log.(filter!(i->!iszero(i),measure_dict["VolEst"][n]))) for n in 1:length(X)]
ve_mean=[exp(mean(mt)) for mt in ve_logn_mts]
ve_lower=[exp(mean(mt))-exp(quantile(mt,.025)) for mt in ve_logn_mts]
ve_upper=[exp(quantile(mt,.975))-exp(mean(mt)) for mt in ve_logn_mts]

vol_chart=scatter(vcat([[X[n] for i in 1:length(measure_dict["VolEst"][n])] for n in 1:length(X)]...),vcat(measure_dict["VolEst"]...), marker=:cross, color=:black, markersize=3, label="Vol. data", ylabel="Est. retinal volume (μm³)",  showaxis=:y, xticks=plotticks, xformatter=_->"", legend=:bottomright)
plot!(vol_chart, X, ve_mean, ribbon=(ve_lower,ve_upper), color=:blue, label="Vol. mean")
lens!([2, 20], [2.5e6, 1e7], inset = (1, bbox(0.12, 0., 0.39, 0.40)))
annotate!([(8,4.5e8,Plots.text("C",18))])

println("$(ccdf(ve_logn_mts[2],log(ve_mean[1]))) % of marginal posterior mass of 5dpf volume estimate is greater than the 3dpf estimated mean")

volrate_mean,volrate_lower,volrate_upper=[Vector{Float64}() for i in 1:3]
for (nx, x) in enumerate(X)
    if nx > 1
        d=X[nx]-X[nx-1]
        l,m,u=MTDist_MC_func((a,b)->(exp(a)-exp(b))/d, [filter!(!isinf,log.(measure_dict["VolEst"][nx])), filter!(!isinf,log.(measure_dict["VolEst"][nx-1]))], summary=true)
        push!(volrate_mean,m);push!(volrate_lower,l);push!(volrate_upper,u)
    end
end
volrate_lower=volrate_mean.-volrate_lower
volrate_upper=volrate_upper.-volrate_mean

volrate_chart=plot(X[2:end],volrate_mean,ribbon=(volrate_lower, volrate_upper), color=:grey, label="Simulated mean rate", legend=:topright, xlabel="Age(dpf)", ylabel="Est. volume rate (μm³d⁻¹)", xticks=plotticks)
annotate!([(8,4.5e6,Plots.text("D",18))])

l=@layout [a{0.3h}
            b{0.2h}
            c{0.3h}
            d]
combined=plot(ann_chart,poprate_chart,vol_chart,volrate_chart,layout=(4,1),size=(1200,1200),link=:x)
savefig(combined, "/bench/PhD/Thesis/images/cmz/CMZoverall.png")

l30,m30,u30=MTDist_MC_func((a,b)->(exp(a)-exp(b))/7., [filter!(!isinf,log.(measure_dict["VolEst"][7])), filter!(!isinf,log.(measure_dict["VolEst"][6]))], summary=true)

ratedist60=MTDist_MC_func((a,b)->(exp(a)-exp(b))/30., [filter!(!isinf,log.(measure_dict["VolEst"][8])), filter!(!isinf,log.(measure_dict["VolEst"][7]))], summary=false)

println("6*30dpf volrate mean: $(6*m30) volrate 60dpf quantile : $(quantile(ratedist60,.01))")