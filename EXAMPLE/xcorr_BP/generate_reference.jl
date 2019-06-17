using SeisIO, Noise, JLD2

include("../../src/reference.jl")

corrname  = "outputData/BPnetworkxcorr_neq.jld2"
refname   = "refXCorrs/reference_xcorr_neq.jld2"
convname  = "convData/rms_neq_wol.jld2"
xcorrSize = 4001
numcorr   = 46

compute_reference_xcorr(corrname, refname, xcorrSize)
convergence(corrname, refname, convname, numcorr)
