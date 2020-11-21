using Distributions, Serialization, CSV, BayesianLinearRegression, Plots, DataFrames, Measurements
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

xCR=zeros(0); xmo1=zeros(0)
yCR=zeros(0); ymo1=zeros(0)

for (age,vec) in CR
    println([age for i in 1:length(vec)])
    global xCR=vcat(xCR,[age for i in 1:length(vec)])
    global yCR=vcat(yCR,vec)
end

for (time,vec) in mo1
    global xmo1=vcat(xmo1,[time for i in 1:length(vec)])
    global ymo1=vcat(ymo1,vec)
end

xCR=hcat(ones(length(xCR)),xCR)
xmo1=hcat(ones(length(xmo1)),xmo1)

uncorrCR2=BayesianLinearRegression.fit!(BayesianLinReg(ones(length(yCR),1),yCR))
corrCR2=BayesianLinearRegression.fit!(BayesianLinReg(xCR,yCR))
uncorrCRX=ones(length(X),1)
corrCRX=hcat(ones(length(X)),X)

uncorrmo1=BayesianLinearRegression.fit!(BayesianLinReg(ones(length(ymo1),1),ymo1))
corrmo1=BayesianLinearRegression.fit!(BayesianLinReg(xmo1,ymo1))
uncorrmo1X=ones(length(X[1:2]),1)
corrmo1X=hcat(ones(length(X[1:2])),X[1:2])

function blr_plt(i,x,y,Xs,blrs,colors,labels,q=.975)
    plt=scatter(x,y,marker=:diamond,markersize=8,markerstrokestyle=:bold,markercolor=:black,label="Data",ylabel="Cohort size",xlabel="Age (dpf)",legend=:topright,ylims=[0,maximum(y)+25])
    for (X,blr,color,label) in zip(Xs, blrs,colors,labels)
        line=zeros(size(X,1))
        ribbon=zeros(size(X,1))
        predictions=BayesianLinearRegression.predict(blr,X)
        for (n,p) in enumerate(predictions)
            normal=Normal(value(p),uncertainty(p))
            ribbon[n]=quantile(normal,q)-mean(normal)
            line[n]=mean((value(p)))
        end
        plot!(plt,i,line,ribbon=ribbon,color=color,label=label)
    end
    return plt
end

pCR=blr_plt(X,xCR[:,2],yCR,[corrCRX,uncorrCRX],[corrCR2,uncorrCR2],[:magenta,:green],["Correlated Model - CR", "Uncorrelated Model - CR"])
annotate!([(35,250,Plots.text("A",18))])

pmo1=blr_plt(X[1:2],xmo1[:,2],ymo1,[corrmo1X,uncorrmo1X],[corrmo1,uncorrmo1],[:magenta,:green],["Correlated Model - 1mo", "Uncorrelated Model - 1mo CMZ"])
annotate!([(32,350,Plots.text("B",18))])

combined=plot(pCR,pmo1,layout=grid(2,1),size=(600,800))
savefig(combined,"/bench/PhD/Thesis/images/cmz/a27linreg.png")

println("\\begin{tabular}{|l|l|l|l|}")
println("\\hline")
println("{\\bf Measurement} & {\\bf Uncorrelated logZ} & {\\bf Correlated logZ} & {\\bf logZR}\\\\ \\hline")
println("1dpf Central Remnant & {\\bf $(round(logEvidence(uncorrCR2),digits=3))} & $(round(logEvidence(corrCR2),digits=3)) & $(round(logEvidence(uncorrCR2)-logEvidence(corrCR2),digits=3))\\\\ \\hline")
println("30dpf Cohort & {\\bf $(round(logEvidence(uncorrmo1),digits=3))} & $(round(logEvidence(corrmo1),digits=3)) & $(round(logEvidence(uncorrmo1)-logEvidence(corrmo1),digits=3))\\\\ \\hline")
println("\\end{tabular}")