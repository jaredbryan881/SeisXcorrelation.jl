using Printf

export remove_nanandzerocol, append_wtcorr!, append_wtcorr, copy_without_wtcorr


"""
	remove_nancol(A::AbstractArray)

Remove column (i.e. trace) which has NaN.
"""
function remove_nanandzerocol(A::AbstractArray)

	N = size(A, 2)
	nancol = ones(Int64, N)
	for i = 1:N
		if any(isnan.(A[:, i])) || all(iszero, A[:,i])
			# this has NaN in its column
			nancol[i] = 0
		end
	end

	nancol=convert.(Bool, nancol)
	return A[:, nancol], nancol

end

"""
	append_wtcorr!(C::CorrData, Nfreqband::Int)

Append wavelet transformed cc functions into CorrData.

# Arguments
- `C::CorrData`    : CorrData contains raw cc functions
- `Nfreqband::Int` : Number of frequency band
- `dj::Float64`    : Discretization span of frequency. See doc of SeisNoise.cwt
- `figdir::String` : if !isempty(figdir), output debugging figures to figdir.

"""
function append_wtcorr!(C::CorrData, Nfreqband::Int; dj = 1/12, figdir::String = "")
	# process flow:
	# 1. compute frequency band
	# 2. apply wavelet transform
	# 3. reconstruct waveform with each frequency band
	# 4. store to C.misc["wtcorr"]

	T, N = size(C.corr)

	# append wtcorr to store all cc functions after cwt
	C.misc["wtcorr"] = Array{Float32, 3}(undef, T, N, Nfreqband)
	C.misc["freqband"] = Array{Float32, 3}(undef, T, N, Nfreqband)

	dt = 1/C.fs

	tvec = collect(-C.maxlag:1/C.fs:C.maxlag)

	if !isempty(figdir) && !ispath(figdir)
		mkdir(figdir)
	end

	for traceid = 1:N

		tr = C.corr[:,traceid]

		W, sj, freqs, coi = SeisNoise.cwt(tr, dt, C.freqmin, C.freqmax, dj=dj)

		if !isempty(figdir)
			freqs[freqs .== 0.0] .= 1e-6
			Plots.heatmap(tvec, freqs, abs.(W)'[end:-1:1, :],  fill=:viridis,
						xlabel = "Time(s)", ylabel ="Frequency(Hz)") #, yaxis=:log
			cone = 1.0 ./ coi #[Hz]
			Plots.plot!(tvec, cone)
			ylims!((C.freqmin, C.freqmax))
			Plots.savefig(figdir*"/wavelet_scalogram_$(C.name)_$(C.id)_$(traceid).png")
		end

		# reconstruct signals
		# 1. Full reconstruction

		tr_freq_reconstructed = zeros(Float32, T, Nfreqband)
		#===freqband_linspace===#
		freqband = range(minimum(freqs), maximum(freqs), length = Nfreqband+1)
		#===freqband_logspace===#
		#freqband = logfreqband(minimum(freqs), maximum(freqs), Nfreqband)

		for i = 1:Nfreqband
		   freqbandmin = freqband[i]
		   freqbandmax = freqband[i+1]
		   inds = findall(f -> freqbandmin<= f <= freqbandmax, freqs)
		   for ind in inds
		      inv_W = SeisNoise.icwt(W[:,ind],sj[ind],dt, dj=dj)
		      if any(isnan.(inv_W))
				 println("nan found in icwt.")
		         println(inv_W)
		      end
		      tr_freq_reconstructed[:, i] .+= inv_W
		   end
		end

		if !isempty(figdir)

			Plots.plot(bg=:white)
			for i = 1:Nfreqband
			   freqbandmin = freqband[i]
			   freqbandmax = freqband[i+1]

			   # normalize
			   norm_amp = maximum(abs.(tr_freq_reconstructed[:, i]))
			   Plots.plot!(tvec, 2*(i-1) .+ tr_freq_reconstructed[:, i] ./ norm_amp, label=@sprintf("%4.2f-%4.2f", freqbandmin, freqbandmax))
			end
			#plot full reconstruct
			norm_amp_tr = maximum(abs.(tr))

			Plots.plot!(tvec, 2*(Nfreqband+1) .+ tr ./ norm_amp_tr, label="Raw trace")

			pt = Plots.plot!(size=(800, 600),
			   xtickfontsize=12,
			   ytickfontsize=12,
			   legendfontsize=12,
			   legend=:outertopright)

			Plots.savefig(figdir*"/wavelet_cc_$(C.name)_$(C.id)_$(traceid).png")
		end

		# store tr_freq_reconstructed into C.misc["wtcorr"]
		C.misc["wtcorr"][:, traceid, :] = tr_freq_reconstructed

		if traceid == 1
			# store freq band under assumption that all cc functions have same
			# freqband during cwt.
			C.misc["freqband"] = freqband
		end
	end

end
append_wtcorr(C::CorrData, Nfreqband::Int; dj=1/12, figdir::String = "") = (U =
 deepcopy(C);append_wtcorr!(U, Nfreqband, dj=dj, figdir=figdir);return U)


"""
	copy_without_wtcorr(C::CorrData)

copy corrdata without wtcorr.
"""
copy_without_wtcorr(C::CorrData) = (@views U = deepcopy(C);
   delete!(U.misc, "wtcorr"); return U)

"""
	logfreqband(fmin::AbstractFloat, fmax::AbstractFloat, N::Int)

return logalithmically spaced frequency bands used for banding of wavelet transform.
"""
function logfreqband(fmin::AbstractFloat, fmax::AbstractFloat, N::Int)

	if fmax <= fmin
		error("fmax should be larger than fmin.")
	elseif fmin < 0.0
		error("fmin should be larger than zero.")
	elseif fmin == 0.0
		warning("fmin == 0.0. replaced to 1e-6")
		fmin = 1e-6
	end

	x = log10.(range(fmin, stop=fmax, length=N+1))
	h1 = fmax - fmin
	h2 = log10(fmax/fmin)
	h3 = log10(fmax) .- x
	y = fmax .- (h1 .* h3 ./ h2)
	# avoid numerical error
	y[1] = fmin
	y[end] = fmax
	return y
end
