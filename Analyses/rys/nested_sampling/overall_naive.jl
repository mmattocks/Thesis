using BioMotifInference, BioSequences, DataFrames, BioBackgroundModels, CSV, Serialization
import BioMotifInference: IPM_likelihood

include("/bench/PhD/Thesis/Analyses/rys/position_analysis/position_overlap.jl")

sib_pos = "/bench/PhD/danpos_results/pooled/sib.Fnor.smooth.positions.xls"
rys_pos = "/bench/PhD/danpos_results/pooled/rys.Fnor.smooth.positions.xls"

danio_genome_path = "/bench/PhD/seq/GRCz11/GCA_000002035.4_GRCz11_genomic.fna"
danio_gen_index_path = "/bench/PhD/seq/GRCz11/GCA_000002035.4_GRCz11_genomic.fna.fai"
danio_gff_path = "/bench/PhD/seq/GRCz11/Danio_rerio.GRCz11.94.gff3"

refined_folders_path = "/bench/PhD/NGS_binaries/BBM/refined_folders"
selected_hmms = "/bench/PhD/NGS_binaries/BBM/selected_hmms"

sib_df_binary = "/bench/PhD/NGS_binaries/BMI/sib_diff_positions"
rys_df_binary = "/bench/PhD/NGS_binaries/BMI/rys_diff_positions"
combined_df_binary = "/bench/PhD/NGS_binaries/BMI/combined_diff_positions"

@info "Reading from position.xls..."
sib_df=CSV.read(sib_pos)
rys_df=CSV.read(rys_pos)

@info "Mapping..."
map_positions!(sib_df, rys_df)
map_positions!(rys_df, sib_df)

@info "Making differential position dataframes..."
sib_diff_df = deepcopy(sib_df[findall(iszero,sib_df.mapped_pos),:])
sib_diff_df = sib_diff_df[findall(!isequal("MT"), sib_diff_df.chr),:]
rys_diff_df = deepcopy(rys_df[findall(iszero,rys_df.mapped_pos),:])
rys_diff_df = rys_diff_df[findall(!isequal("MT"), rys_diff_df.chr),:]
sib_mapped_df = deepcopy(sib_df[findall(!iszero,sib_df.mapped_pos),:])
sib_mapped_df = sib_mapped_df[findall(x->x>0, sib_mapped_df.start),:]
sib_mapped_df = sib_mapped_df[findall(!isequal("MT"), sib_mapped_df.chr),:]
rys_mapped_df = deepcopy(rys_df[findall(!iszero,rys_df.mapped_pos),:])
rys_mapped_df = rys_mapped_df[findall(x->x>0, rys_mapped_df.start),:]
rys_mapped_df = rys_mapped_df[findall(!isequal("MT"), rys_mapped_df.chr),:]

add_position_sequences!(sib_diff_df, danio_genome_path, danio_gen_index_path)
add_position_sequences!(rys_diff_df, danio_genome_path, danio_gen_index_path)
add_position_sequences!(sib_mapped_df, danio_genome_path, danio_gen_index_path)
add_position_sequences!(rys_mapped_df, danio_genome_path, danio_gen_index_path)

@info "Filtering ambiguous sequences..."
deleterows!(sib_diff_df, [hasambiguity(seq) for seq in sib_diff_df.seq])
deleterows!(rys_diff_df, [hasambiguity(seq) for seq in rys_diff_df.seq])
deleterows!(sib_mapped_df, [hasambiguity(seq) for seq in sib_mapped_df.seq])
deleterows!(rys_mapped_df, [hasambiguity(seq) for seq in rys_mapped_df.seq])

@info "Masking positions by genome partition and strand..."
BioBackgroundModels.add_partition_masks!(sib_diff_df, danio_gff_path, 500, (:chr,:seq,:start))
BioBackgroundModels.add_partition_masks!(rys_diff_df, danio_gff_path, 500, (:chr,:seq,:start))
BioBackgroundModels.add_partition_masks!(sib_mapped_df, danio_gff_path, 500, (:chr,:seq,:start))
BioBackgroundModels.add_partition_masks!(rys_mapped_df, danio_gff_path, 500, (:chr,:seq,:start))

@info "Creating coded observation sets..."
sib_codes = Matrix(transpose(observation_setup(sib_diff_df, order=0, symbol=:seq)))
rys_codes = Matrix(transpose(observation_setup(rys_diff_df, order=0, symbol=:seq)))
sib_mapped_codes = Matrix(transpose(observation_setup(sib_mapped_df, order=0, symbol=:seq)))
rys_mapped_codes = Matrix(transpose(observation_setup(rys_mapped_df, order=0, symbol=:seq)))

@info "Setting up for BBM likelihood calculations..."
BHMM_dict = Dict{String,BHMM}()
refined_folders=deserialize(refined_folders_path)
for (part, folder) in refined_folders
    BHMM_dict[part]=folder.partition_report.best_model[2]
end

sib_lh_matrix = BGHMM_likelihood_calc(sib_diff_df, BHMM_dict, symbol=:seq)
rys_lh_matrix = BGHMM_likelihood_calc(rys_diff_df, BHMM_dict, symbol=:seq)
sibm_lh_matrix= BGHMM_likelihood_calc(sib_mapped_df, BHMM_dict, symbol=:seq)
rysm_lh_matrix = BGHMM_likelihood_calc(rys_mapped_df, BHMM_dict, symbol=:seq)

sib_e_pth = "/bench/PhD/NGS_binaries/BMI/sib_e"
rys_e_pth = "/bench/PhD/NGS_binaries/BMI/rys_e"

sibe=deserialize(sib_e_pth*"/ens")
ryse=deserialize(rys_e_pth*"/ens")

sib_MAP=deserialize(sibe.models[findmax([rec.log_Li for rec in sibe.models])[2]].path)
rys_MAP=deserialize(ryse.models[findmax([rec.log_Li for rec in ryse.models])[2]].path)

sib_bg_lh = IPM_likelihood(sib_MAP.sources, sib_codes, [findfirst(iszero,sib_codes[:,o])-1 for o in 1:size(sib_codes,2)], sib_lh_matrix, falses(size(sib_codes,2),length(sib_MAP.sources)))

sibm_bg_lh = IPM_likelihood(sib_MAP.sources, sib_mapped_codes, [findfirst(iszero,sib_mapped_codes[:,o])-1 for o in 1:size(sib_mapped_codes,2)], sibm_lh_matrix, falses(size(sib_mapped_codes,2),length(sib_MAP.sources)))

rys_bg_lh = IPM_likelihood(rys_MAP.sources, rys_codes, [findfirst(iszero,rys_codes[:,o])-1 for o in 1:size(rys_codes,2)], rys_lh_matrix, falses(size(rys_codes,2),length(rys_MAP.sources)))

rysm_bg_lh = IPM_likelihood(rys_MAP.sources, rys_mapped_codes, [findfirst(iszero,rys_mapped_codes[:,o])-1 for o in 1:size(rys_mapped_codes,2)], rysm_lh_matrix, falses(size(rys_mapped_codes,2),length(rys_MAP.sources)))

println("\\begin{tabular}{|l|l|l|l|}")
println("\\hline")
println("{\\bf Sequence set} & {\\bf \\textit{rys} naive LH} & {\\bf sib naive LH} & {\\bf LH ratio}\\\\ \\hline")
println("Differential & $(round(rys_bg_lh, sigdigits=5)) & $(round(sib_bg_lh, sigdigits=5)) & $(round(rys_bg_lh-sib_bg_lh,sigdigits=5))\\\\ \\hline")
println("Mapped & $(round(rysm_bg_lh, sigdigits=5)) & $(round(sibm_bg_lh, sigdigits=5)) & $(round(rysm_bg_lh-sibm_bg_lh,sigdigits=5))\\\\ \\hline")
println("\\end{tabular}")