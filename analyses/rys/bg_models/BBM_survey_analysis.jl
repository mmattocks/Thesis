using BioBackgroundModels, Serialization, Plots, Measures
import Measures:mm
import Plots:Plot

hmm_output = "/bench/PhD/NGS_binaries/BBM/survey_chains"
sample_output = "/bench/PhD/NGS_binaries/BBM/survey_samples"

survey_folders = "/bench/PhD/NGS_binaries/BBM/survey_folders"

chains=deserialize(hmm_output)
sample_dfs = deserialize(sample_output)
training_sets, test_sets = split_obs_sets(sample_dfs)

report_folders=generate_reports(chains, test_sets)
serialize(survey_folders, report_folders)

partplts=Vector{Plot}()

for o in [0,1,2]
    for p in ["exon","intergenic","periexonic"]
        rf=report_folders[p]
        pr=rf.partition_report
        od=pr.orddict
        or=od[o]

        p=="exon" ? (ylbl="Order $o lh"; yfmt=y->y; lm=0mm) : (ylbl=""; yfmt=_->""; lm=-7mm)
        o==2 ? (xlbl=p; xfmt=x->Int64(x);bm=0mm) : (xlbl=""; xfmt=_->""; bm=-7mm)

        plt=scatter(or.converged_K,or.converged_lh, marker=:cross, markersize=15, markercolor=:black, label="Ord. $o model lhs", xlabel=xlbl, ylabel=ylbl, legend=nothing, xformatter=xfmt, yformatter=yfmt, xticks=[1,2,4,6], leftmargin=lm,bottommargin=bm)
        plot!(unique(or.converged_K),[pr.naive_lh for l in 1:length(unique(or.converged_K))], label="Naive lh", color=:darkmagenta,linewidth=2.5)
        push!(partplts,plt)
    end
end
combined=plot(partplts...,layout=grid(3,3),link=:xy, size=(600,600))
savefig(combined,"/bench/PhD/Thesis/images/rys/bghmm.png")