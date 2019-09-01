__precompile__()
module SeisXcorrelation

using SeisIO, SeisNoise, Dates, FFTW, JLD2, Sockets, ORCA, PlotlyJS, Plots, Distributed

# compute cross-correlation
include("seisxcorr.jl")
# compute reference cross-correlation function
include("reference.jl")
# compute stacking (linear stack, selective stack, ...)
include("pstack.jl")

end # module
