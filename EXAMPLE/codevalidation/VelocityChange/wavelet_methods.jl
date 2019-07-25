"""
Allen and Smith autoregressive lag-1 autocorrelation coefficient.
In an AR(1) model
    x(t) - <x> = γ(x(t-1) - <x>) + α z(t)
where <x> is the process mean, γ and α are process parameters and z(t) is a
Gaussian unit-variance white noise

References
    ----------
    [1] Allen, M. R. and Smith, L. A. Monte Carlo SSA: detecting
        irregular oscillations in the presence of colored noise.
        *Journal of Climate*, **1996**, 9(12), 3373-3404.
        <http://dx.doi.org/10.1175/1520-0442(1996)009<3373:MCSDIO>2.0.CO;2>
    [2] http://www.madsci.org/posts/archives/may97/864012045.Eg.r.html

Original function translated from https://github.com/regeirk/pycwt/blob/master/pycwt/helpers.py
"""
function ar1(x::AbstractArray)
    N = size(x)
    xm = mean(x)
    x = x .- xm

    # estimate the lag zero and one covariance
    c0 = (x' .* x) / N
    c1 = (x[0:N-1]' .* x[1:N]) / (N-1)

    # according to A. Grinsteds' substitutions
    A = c0 * N^2
    B = -c1 .* N .- A - 2*c0 .+ 2*c1 .- c1 * N^2 .+ c0*N
    C = N * (c0 .+ c1*N .- c1)
    D = B^2 - 4*A .* C

    if D > 0
        g = (-B - sqrt(D) / (2*A))
    else
        @warn "Cannot place an upperbound on the unbiased AR(1). Series is
               too short or trend too large."
    end

    # according to Allen & Smith (1996), footnote 4
    μ2 = -1 / N + (2 / N^2) * ((N - g^N) / (1-g) - g * (1 - g^(N-1)) / (1-g)^2)
    c0t = c0 / (1-μ2)
    a = sqrt((1-g^2) * c0t)

    return g, a, μ2
end

"""
Wavelet coherence transform
Modified from original by Congcong Yuan
"""
function wct_modified(y1::Union{Array{Float64,1}, Array{Float32,1}}, y2::Union{Array{Float64,1}, Array{Float32,1}}, dt::Float64, dj::Float64, s0::Int64; J::Int64=109, wavelet::String="morlet", norm::Bool=true, sig::Bool=false, siglvl::Float64=0.95)
    # define wavelet
    if wavelet=="morlet"
        wav = WT.Morlet()
    end

    # calculate std dev of both input signals (std(y1), std(y2))
    std1 = std(y1)
    std2 = std(y2)
    # normalize
    if norm
        y1 = (y1 .- mean(y1)) ./ std1
        y2 = (y2 .- mean(y2)) ./ std2
    end

    # calculate the CWT of the time series, using identical parameters for both calculations
    W1 = cwt(y1, wav, J1=J)
    W2 = cwt(y2, wav, J1=J)

    per, sj, coi = Transforms.caveats(y1, wav)
    freq = 1 ./ per

    scales1 = ones(length(y1)) .* sj
    scales2 = ones(length(y2)) .* sj

    # smooth wavelet spectra before truncating
    S1 = smooth(abs.(W1).^2 ./ scales1, dt, dj, sj) ##########
    S2 = smooth(abs.(W2).^2 ./ scales1, dt, dj, sj)

    # compute wavelet transform coherence
    W12 = W1 .* conj(W2)
    scales = ones(length(y1)) .* sj1

    S12 = wavelet.smooth(W12 ./ scales, dt, dj, sj) ############
    WCT = abs.(S12).^2 ./ (S1 .* S2)
    aWCT = angle.(W12)

    # calculate cross spectrum and its amplitude
    WXS = W12
    WXA = abs.(S12)

    # calculate the significance using a Monte Carlo simulation with 95%
    # confidence as a function of scale
    if sig
        a1, b1, c1 = pycwt.ar1(y1)
        a2, b2, c2 = pycwt.ar1(y2)
        sig = pycwt.wct_significance(a1, a2, dt=dt, dj=dj, s0=s0, J=J, significance_level=siglvl, wavelet=wavelet)
    else
        sig = [0]
    end

    return WXS, WXA, WCT, aWCT, coi, freq, sig
end

"""
Compute wavelet cross spectrum given an array of frequency band minima and maxima.
"""
function wxs_freqbands(cur, ref, t, twindow, freqbands, dj, s0; J::S=NaN, wvn::String="morlet", unwrap::Bool=false, sig::Bool=false, siglvl::Float64=0.95) where S<:Real
    # define sample frequency
    dt = t[2] - t[1]
    fs = 1/dt
    (tmin, tmax) = twindow

    # trim time vector
    t_ind = findall(x->(x≤tmax && x≥tmin), t)

    # trim vectors to correspond to time
    t = t[t_ind]
    cur = cur[t_ind]
    ref = ref[t_ind]

    WXS, WXA, WCT, aWCT, coi, freq, sig = wct_modified(cur, ref, dt, dj, s0, J=J, wavelet=wvn, sig=sig, siglvl=siglvl)
end
