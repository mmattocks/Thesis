using Distributions, Serialization, CSV, BayesianLinearRegression, Plots, DataFrames, Measurements
import Measurements:value,uncertainty


a25path="/bench/PhD/datasets/a25.csv"
a25df=DataFrame(CSV.File(a25path))

X=Vector{Float64}()
yD=Vector{Vector{Float64}}()
yV=Vector{Vector{Float64}}()

for timedf in groupby(a25df,"Time")
    push!(X,timedf[1,"Time"])
    dvec=Vector{Float64}()
    vvec=Vector{Float64}()
    for row in eachrow(timedf)
        dv=row."D/V"
        if dv == "D"
            push!(dvec, row.EdU/row.PCNA)
        else
            push!(vvec,row.EdU/row.PCNA)
        end
    end

    push!(yD,dvec)
    push!(yV,vvec)
end

sp=sortperm(X)
yD=yD[sp]; yV=yV[sp]
sort!(X)

xYD=zeros(0); xYV=zeros(0)
YD=zeros(0); YV=zeros(0)

for (n,vec) in enumerate(yD)
    global xYD=vcat(xYD,[X[n] for i in 1:length(vec)])
    global YD=vcat(YD,vec)
end

for (n,vec) in enumerate(yV)
    global xYV=vcat(xYV,[X[n] for i in 1:length(vec)])
    global YV=vcat(YV,vec)
end

xYD=hcat(ones(length(xYD)),xYD)
xYV=hcat(ones(length(xYV)),xYV)

xY=vcat(xYD,xYV);Y=vcat(YD,YV)

dm=BayesianLinearRegression.fit!(BayesianLinReg(xYD,YD))
vm=BayesianLinearRegression.fit!(BayesianLinReg(xYV,YV))
tm=BayesianLinearRegression.fit!(BayesianLinReg(xY,Y))

println("Dorsal model logZ: $(logEvidence(dm))")
println("Ventral model logZ: $(logEvidence(vm))")
println("Joint D/V logZ: $(logEvidence(dm) + logEvidence(vm)) ")
println("Combined model logZ: $(logEvidence(tm))")
println("Combined/Joint ratio : $(logEvidence(tm)-(logEvidence(dm) + logEvidence(vm)))")

function blr_plt(x,y,Xs,blrs,colors,labels,q=.975)
    plt=scatter(x,y,marker=:cross,markersize=3,markercolor=:black,label="Data",ylabel="PCNA+ve labelled fraction",xlabel="Pulse time (hr)",legend=:bottomright, ylims=(-0.1,1.))
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

dplt=blr_plt(xYD[:,2],YD,[xYD],[dm],[:green],["Dorsal model"])
annotate!([(1,.9,Plots.text("A",12))])
vplt=blr_plt(xYV[:,2],YV,[xYV],[vm],[:magenta],["Ventral model"])
annotate!([(1,.9,Plots.text("B",12))])
tplt=blr_plt(xY[:,2],Y,[xY],[tm],[:black],["Combined model"])
annotate!([(1,.9,Plots.text("C",12))])

combined_plt=plot(dplt,vplt,tplt,layout=(3,1), size=(600,900))

savefig(combined_plt, "/bench/PhD/Thesis/images/cmz/3ddvlinreg.png")
@info "Saved fig."


println("\\begin{tabular}{|l|l|l|l|}")
println("\\hline")
println("{\\bf Model} & {\\bf Implied \$T_c\$} & {\\bf Implied \$T_s\$} & {\\bf logZ}\\\\ \\hline")

for (name,model) in zip(["Dorsal","Ventral","Combined"],[dm,vm,tm])
    local intercept,slope=posteriorWeights(model)
    println("$name & $(1/slope) & $(intercept*1/slope) & $(round(logEvidence(model),digits=3))\\\\ \\hline")
end

println("\\end{tabular}")
