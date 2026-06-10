module CampaignLib

using BoundaryIntegral
using LinearAlgebra
using Printf
using Serialization
using TOML

export Campaign, load_campaign, batch_path, v_path, manifest_path, centers_path,
       targets_path, rho_store_path, logs_dir,
       CenterInfo, BatchSpec, enumerate_centers, enumerate_pairs, build_batches,
       write_manifest, read_manifest, write_centers, read_centers,
       pending_batches, load_templates!,
       solve_batch, consolidate, eval_batch, assemble_v, prepare

include("config.jl")
include("manifest.jl")
include("tasks.jl")

end
