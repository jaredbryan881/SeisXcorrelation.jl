using Plots
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

function smooth(W::AbstractArray, dt::Float64, dj::Float64, scales::Array{Float64,1})
    # The smoothing is performed by using a filter given by the absolute
    # value of the wavelet function at each scale, normalized to have a
    # total weight of unity, according to suggestions by Torrence &
    # Webster (1999) and by Grinsted et al. (2004).
    m, n = size(W)

    # Filter in time
    k = 2π*fftfreq(length(W[1,:]))
    k2 = k.^2
    snorm = scales ./ dt

    # Smoothing by Gaussian window (absolute value of wavelet function)
    # using the convolution theorem: multiplication by Gaussian curve in
    # Fourier domain for each scale, outer product of scale and frequency
    F = exp.(-0.5 .* (snorm.^2) .* k2') # outer product
    smooth = (ifft(F .* fft(W,2), 2))

    T = smooth[:, 1:n] # Remove possible padded region due to FFTW

    # Filter in scale. For the Morlet wavelet, this is simply a boxcar with 0.6 width
    # construct boxcar
    wsize = convert(Int64, round(0.60 / dj * 2))
    if wsize % 2 == 0
        wsize+=1
    end
    halfWin = div(wsize,2)

    # iterate over 'horizontal' and smooth in the 'vertical'
    # this could also be done by adding an axis to the transpose and performing a 2d convolution
    for i=1:size(T,2)
        # pad signal for 'same' padding
        paddedT = vcat(T[1, i]*ones(halfWin), T[:, i], T[end, i]*ones(halfWin))
        # smooth signal
        win = Array{eltype(paddedT), 1}(undef, wsize)
        win.=ones(wsize)/wsize
        smoothT = conv(paddedT, win)
        # trim signal
        T[:, i] = smoothT[2*halfWin+1:end-2*halfWin]
    end

    return T
end

"""
Modified wavelet coherence transform
Modified from original by Congcong Yuan
"""
function wct_modified(y1, y2, dt, dj, s0, J; wavelet::String="morlet", norm::Bool=true, sig::Bool=false, siglvl::Float64=0.95)
    # define wavelet
    if wavelet=="morlet"
        wav = WT.Morlet(6)
    end

    fλ = (4*π) / (wav.σ + sqrt(2 + wav.σ^2))

    # calculate std dev of both input signals (std(y1), std(y2))
    std1 = std(y1)
    std2 = std(y2)
    # normalize
    if norm
        y1 = (y1 .- mean(y1)) ./ std1
        y2 = (y2 .- mean(y2)) ./ std2
    end

    # calculate the CWT of the time series, using identical parameters for both calculations
    W1, sj, freqs, coi = cwt(y1, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)
    W2, sj, freqs, coi = cwt(y2, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)

    # fourier frequencies corresponding to wavelet scales
    freqs = 1 ./ (fλ .* sj)

    scales1 = Array{Float64,2}(undef, size(W1))
    scales2 = Array{Float64,2}(undef, size(W1))
    for i=1:size(scales1, 2)
        scales1[:,i] .= sj
        scales2[:,i] .= sj
    end

    # smooth wavelet spectra before truncating
    S1 = smooth(abs.(W1).^2 ./ scales1, dt, dj, sj)
    S2 = smooth(abs.(W2).^2 ./ scales2, dt, dj, sj)

    # compute wavelet transform coherence
    W12 = W1 .* conj(W2)
    scales = scales1
    S12 = smooth(W12 ./ scales, dt, dj, sj)
    WCT = abs.(S12).^2 ./ (S1 .* S2)
    aWCT = angle.(W12)

    # calculate cross spectrum and its amplitude
    WXS = W12
    WXA = abs.(S12)

    # calculate the significance using a Monte Carlo simulation with 95%
    # confidence as a function of scale
    if sig
        a1, b1, c1 = ar1(y1)
        a2, b2, c2 = ar1(y2)
        sig = wct_significance(a1, a2, dt=dt, dj=dj, s0=s0, J=J, significance_level=siglvl, wavelet=wavelet)
    else
        sig = [0]
    end

    return WXS, WXA, WCT, aWCT, coi, freq, sig
end

"""
Compute wavelet cross spectrum given an array of frequency band minima and maxima.
"""
function wxs_freqbands(cur, ref, t, twindow, freqbands, dj, s0, J; wvn::String="morlet", unwrap::Bool=false, sig::Bool=false, siglvl::Float64=0.95)
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

    WXS, WXA, WCT, aWCT, coi, freq, sig = wct_modified(cur, ref, dt, dj, s0, J, wavelet=wvn, sig=sig, siglvl=siglvl)
end
