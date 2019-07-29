using SeisIO, SeisNoise, JLD2

include("../../src/reference.jl")

dir = "BPnetwork_2003"
# input file
corrname  = "/Users/jared/SCECintern2019/data/xcorrs/$dir/BPnetwork_2003_xcorrs"
# reference output file
refname   = "/Users/jared/SCECintern2019/data/reference/$dir/BPnetwork_2003_ref.jld2"
# convergence output files
convname_snr  = "/Users/jared/SCECintern2019/data/convergence/$dir/BPnet_2003_snr.jld2"
convname_cc = "/Users/jared/SCECintern2019/data/convergence/$dir/BPnet_2003_cc.jld2"

# generate reference
compute_reference_xcorr(corrname, refname)
# compute convergence of successively long stacks to the reference xcorr
convergence(corrname, refname, convname_cc, metric="cc")
convergence(corrname, refname, convname_snr, metric="snr")
