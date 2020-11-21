using BayesianLinearRegression,CSV,DataFrames,Distributions,NGRefTools,StatsBase,Plots,GMC_NS,Serialization, Measurements
import Measurements:value,uncertainty

default(legendfont = (8), guidefont = (10,), tickfont = (8), guide = "x")

a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"
pne_pth="/bench/PhD/NGS_binaries/BSS/A10/pop_norm"
plne_pth="/bench/PhD/NGS_binaries/BSS/A10/pop_lognorm"
vne_pth="/bench/PhD/NGS_binaries/BSS/A10/vol_norm"
vlne_pth="/bench/PhD/NGS_binaries/BSS/A10/vol_lognorm"

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["PopEst"]=Vector{Vector{Float64}}()
measure_dict["VolEst"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    push!(X,Float64(unique(t_df."Time point (d)")[1]))
    pevec,tvec,dvec,vvec,vevec, rtvec, rplvec, rptvec, ldvec, ondvec, onlvec, oplvec, inlvec, iplvec, gclvec = [zeros(0) for i in 1:15]

    for i_df in groupby(t_df,"BSR")
        length([skipmissing(i_df."Lens diameter")...])>0 && push!(ldvec,mean(skipmissing(i_df."Lens diameter")))
    end

    for i_df in groupby(t_df,"BSR")
        if length([skipmissing(i_df."CMZ Sum")...])>0
            push!(tvec,mean(skipmissing(i_df."CMZ Sum")))
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
   
               push!(vevec,(ov-iv)*(4/5))
           end
    end

    push!(measure_dict["PopEst"],filter(p->p>0,pevec))
    push!(measure_dict["VolEst"],vevec)
end

x=measure_dict["PopEst"][1]
uncorrX3d=ones(length(x),1)
corrX3d=hcat(ones(length(x)),x)
y=measure_dict["VolEst"][1]
uncorrm=BayesianLinearRegression.fit!(BayesianLinReg(uncorrX3d,y))
corrm=BayesianLinearRegression.fit!(BayesianLinReg(corrX3d,y))

function blr_plt(x,y,Xs,blrs,colors,labels,q=.975)
    plt=scatter(x,y,marker=:diamond,markersize=8,markerstrokestyle=:bold,markercolor=:black,label="Data",ylabel="Estimated retinal volume",xlabel="Estimated CMZ Annulus Population",legend=:bottomleft)
    for (X,blr,color,label) in zip(Xs, blrs,colors,labels)
        line=zeros(size(X,1))
        ribbon=zeros(size(X,1))
        predictions=BayesianLinearRegression.predict(blr,X)
        for (n,p) in enumerate(predictions)
            normal=Normal(value(p),uncertainty(p))
            ribbon[n]=quantile(normal,q)-mean(normal)
            line[n]=mean((value(p)))
        end
        plot!(plt,x,line,ribbon=ribbon,color=color,label=label)
    end
    return plt
end

p=blr_plt(x,y,[corrX3d,uncorrX3d],[corrm,uncorrm],[:magenta,:green],["Correlated Model", "Uncorrelated Model"])

savefig(p, "/bench/PhD/Thesis/images/cmz/3dcorr.png")
@info "Saved fig."

println("\\begin{tabular}{|l|l|l|l|}")
println("\\hline")
println("{\\bf Age (dpf)} & {\\bf Uncorrelated logZ} & {\\bf Correlated logZ} & {\\bf logZR}\\\\ \\hline")

for (x,xs,ys) in zip(X,measure_dict["PopEst"],measure_dict["VolEst"])
    uncorrX=ones(length(xs),1)
    corrX=hcat(ones(length(xs)),xs)
    uncorr=BayesianLinearRegression.fit!(BayesianLinReg(uncorrX,ys))
    corr=BayesianLinearRegression.fit!(BayesianLinReg(corrX,ys))
    println("$x & {\\bf $(round(logEvidence(uncorr),digits=3))} & $(round(logEvidence(corr),digits=3)) & $(round(logEvidence(uncorr)-logEvidence(corr),digits=3))\\\\ \\hline")
end

println("\\end{tabular}")


        