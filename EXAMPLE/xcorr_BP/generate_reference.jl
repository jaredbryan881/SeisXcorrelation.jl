using SeisIO, SeisNoise, JLD2

include("../../src/reference.jl")

corrname  = "./dataset/BPnetwork_Jan03_neq_xcorrs.jld2"
refname   = "refXCorrs/BPnet_neq_Jan03_refxcorr.jld2"
convname  = "convData/BPnet_Jan03_rms_neq.jld2"
xcorrSize = 4001

#compute_reference_xcorr(corrname, refname, xcorrSize)
convergence(corrname, refname, convname)
