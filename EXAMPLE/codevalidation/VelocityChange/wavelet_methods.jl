using Plots, SeisNoise, DTWDT, PlotlyJS
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
        # TODO: implement and test the wct_significance and ar1 functions
        a1, b1, c1 = ar1(y1)
        a2, b2, c2 = ar1(y2)
        sig = wct_significance(a1, a2, dt=dt, dj=dj, s0=s0, J=J, significance_level=siglvl, wavelet=wavelet)
    else
        sig = [0]
    end

    return WXS, WXA, real.(WCT), aWCT, coi, freqs, sig
end

"""
Compute wavelet cross spectrum given an array of frequency band minima and maxima.
"""
function wxs_freqbands(cur, ref, t, twindow, freqbands, dj, s0, J; wvn::String="morlet", unwrapflag::Bool=false, sig::Bool=false, siglvl::Float64=0.95, windowflag::Bool=true, nwindow::Float64=1.5)
    # define sample frequency
    dt = t[2] - t[1]
    fs = 1/dt

    # perform wavelet coherence transform
    WXS, WXA, WCT, aWCT, coi, freqs, sig = wct_modified(cur, ref, dt, dj, s0, J, wavelet=wvn, sig=sig, siglvl=siglvl)

    # do inverse cwt for different frequency bands
    if unwrapflag==true
        phase = unwrap(aWCT, dims=ndims(aWCT))
    else
        phase = aWCT
    end

    (nband, nfreq) = size(freqbands)
    if nband==0 || nfreq!=2
        println("Error: please check if inputs in freqbands is right!")
    else
        dtt = zeros(nband)
        err = zeros(nband)

        # iterate over frequency bands
        for iband=1:nband
            fmin = freqbands[iband, 1]
            fmax = freqbands[iband, 2]
            if (fmax > maximum(freqs)) || (fmax <= fmin)
                println("Error: please ensure columns 1 and 2 are right frequency limits in freqbands!")
            else
                freq_ind = findall(f->(f>=fmin && f<=fmax), freqs)

                iphase = []
                iphase = phase[freq_ind, :]

                if windowflag
                    tmin = twindow[1]
                    tmax = twindow[2]

                    # time checks
                    if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                        println("Error: please input right time limits in the time window!")
                    else
                        # truncate data with the time window
                        t_ind = findall(x->(x≤tmax && x≥tmin), t)
                        wt = t[t_ind]
                    end
                else # using dynamic time window
                    tmin = twindow[1]
                    tmax = twindow[1] + nwindow*(1/fmin + 1/fmax)/2.

                    if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                        println("Error: please input right time limits in the time window!")
                    else
                        # truncate data with the time window
                        t_ind = findall(x->(x≤tmax && x≥tmin), t)
                        wt = t[t_ind]
                    end
                end

                delta_t = zeros(nband, length(t))
                delta_t_err = zeros(nband, length(t))
                # do regression to get slope dt_m between phase delays and frequency band
                for itime=1:length(wt)
                    # get weights
                    # w = 1 ./ WXA[freq_ind,itime]
                    w = 1 ./ WCT[freq_ind, itime]
                    infNaN = findall(x->(isnan.(x) || isinf.(x)), w)
                    if length(infNaN)!=0
                        w[infNaN] .= 1.0
                    end

                    # perform a WLS inversion
                    # This does NOT force the best fit line through the origin
                    model = glm(@formula(Y ~ X),DataFrame(X=freqs[freq_ind]*2π,Y=iphase[:,itime]),Normal(),
                    IdentityLink(),wts=w)

                    delta_t[iband, itime] = coef(model)[2]
                    delta_t_err[iband, itime] = stderror(model)[2]
                end

                # do regression in time to get time shift
                # w2 = 1/np.mean(WXA[freq_ind,:],axis=0)
                w2 = 1 ./ mean(WCT[freq_ind,:], dims=1)
                infNaN =findall(x->(isnan.(x) || isinf.(x)), w2[:])
                if length(infNaN)!=0
                    w2[:, infNaN] .= 1.0
                end
            end

            if length(wt)>2
                if delta_t == nothing
                    continue
                end

                # perform a WLS inversion
                # This DOES force the best fit line through the origin
                model = glm(@formula(Y ~0 + X),DataFrame(X=wt,Y=delta_t[iband, :]),Normal(),
                IdentityLink(),wts=w2[:])
                dtt[iband] = coef(model)[1]
                err[iband] = stderror(model)[1]
            else
                print("not enough points to estimate dv/v")
                dtt[iband]=NaN
                err[iband]=NaN
            end
        end
    end

    return freqbands, dtt.*100, err.*100
end

"""
Compute wavelet cross spectrum given an array of frequencies
"""
function wxs_allfreq(cur, ref, t, twindow, fwindow, dj, s0, J; wvn::String="morlet", unwrapflag::Bool=false, sig::Bool=false, siglvl::Float64=0.95, windowflag::Bool=true, nwindow::Float64=1.5)
    # Part 1 ----- Cross spectrum analysis
    # define sample frequency
    dt = t[2] - t[1]
    fs = 1/dt
    (fmin, fmax) = fwindow

    # perform wavelet coherence transform
    WXS, WXA, WCT, aWCT, coi, freqs, sig = wct_modified(cur, ref, dt, dj, s0, J, wavelet=wvn, sig=sig, siglvl=siglvl)

    # do inverse cwt for different frequency bands
    if unwrapflag==true
        phase = unwrap(aWCT, dims=2)
    else
        phase = aWCT
    end

    delta_t = Array{Float64,2}(undef, size(phase))
    for i=1:size(phase,2)
        delta_t[:, i] = phase[:, i] ./ (2π.*freqs)
    end

    # zero out cone of influence and data outside frequency band
    if (fmax > maximum(freqs) || fmax<=fmin)
        println("Error: please input correct frequency limits in the frequency window!")
    else
        freq_indout = findall(f->(f < fmin || f > fmax), freqs)
        delta_t[freq_indout, :] .= 0.
        freq_indin = findall(f->(f >= fmin && f <= fmax), freqs)
    end

    # Part 2 ------ dv/v
    dvv = zeros(size(freqs))
    err = zeros(size(freqs))

    if windowflag
        (tmin, tmax) = twindow

        # time checks
        if tmin < minimum(twindow) || tmax > maximum(twindow) || tmax <= tmin
            println("Error: please input right time limits in the time window!")
        else
            # trim time vector
            t_ind = findall(x->(x≤tmax && x≥tmin), t)
            # trim vectors to correspond to time
            wt = t[t_ind]
        end

        # perform linear regression
        for ifreq in freq_indin
            if length(wt) > 2
                if delta_t[ifreq, :] == nothing
                    continue
                end

                # use WXA as weight for regression
                w = 1 ./ WCT[ifreq, :]
                infNaN = findall(x->(isnan.(x) || isinf.(x)), w)
                if length(infNaN)!=0
                    w[infNaN] .= 1.0
                end
                # perform linear regression
                model = glm(@formula(Y ~0 + X),DataFrame(X=wt,Y=delta_t[ifreq,:]),Normal(),
                            IdentityLink(),wts=w)

                dvv[ifreq] = coef(model)[1]
                err[ifreq] = stderror(model)[1]
            else
                println("Not enough points to estimate dv/v")
                dvv[ifreq] = NaN
                err[ifreq] = NaN
            end
        end
    else # using dynamic time window
        for ifreq in freq_indin
            tmin = twindow[1]
            tmax = twindow[1] + nwindow*(1.0 / freqs[ifreq])

            # time checks
            if tmin < minimum(twindow) || tmax > maximum(twindow) || tmax <= tmin
                println("Error: please input right time limits in the time window!")
            else
                # trim time vector
                t_ind = findall(x->(x≤tmax && x≥tmin), t)
                # trim vectors to correspond to time
                wt = t[t_ind]
            end

            if length(wt) > 2
                if delta_t[ifreq, :] == nothing
                    continue
                end

                # use WXA as weight for regression
                w = 1 ./ WCT[ifreq, :]
                infNaN = findall(x->(isnan.(x) || isinf.(x)), w)
                if length(infNaN)!=0
                    w[infNaN] .= 1.0
                end
                # perform linear regression
                model = glm(@formula(Y ~0 + X),DataFrame(X=wt,Y=delta_t[ifreq,:]),Normal(),
                            IdentityLink(),wts=w)

                dvv[ifreq] = coef(model)[1]
                err[ifreq] = stderror(model)[1]
            else
                println("Not enough points to estimate dv/v")
                dvv[ifreq] = NaN
                err[ifreq] = NaN
            end
        end
    end

    return freqs, dvv.*100, err.*100
end

function wts_freqbands(cur, ref, t, twindow, freqbands, dj, s0, J; dvmax::Float64=0.1, ndv::Int64=100, wvn::String="morlet", normalize::Bool=true, windowflag::Bool=true, nwindow::Float64=1.5)
    # define sample frequency
    dt = t[2] - t[1]
    fs = 1/dt

    # apply cwt to two traces
    # define wavelet
    if wvn=="morlet"
        wav = WT.Morlet(6)
    end
    # calculate the CWT of the time series, using identical parameters for both calculations
    cwt1, sj, freqs, coi = cwt(cur, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)
    cwt2, sj, freqs, coi = cwt(ref, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)

    (nband, nfreq) = size(freqbands)
    if nband==0 || nfreq != 2
        println("Error: please check if inputs in freqbands is right!")
    else
        dvv = zeros(Float32, nband)
        cc = zeros(Float32, nband)
        cdp = zeros(Float32, nband)
        err = zeros(Float32, nband)

        # iterate over frequency bands
        for iband=1:nband
            (fmin, fmax) = freqbands[iband, :]

            if (fmax > maximum(freqs) || fmax <= fmin)
                println("Error: please ensure columns 1 and 2 are right frequency limits in freqbands!")
            else
                freq_ind = findall(f->(f>=fmin && f<=fmax), freqs)

                # perform icwt
                icwt1 = icwt(cwt1[freq_ind, :], CFW(wav,1/dj), sj[freq_ind], dt=dt, dj=dj)
                icwt2 = icwt(cwt2[freq_ind, :], CFW(wav,1/dj), sj[freq_ind], dt=dt, dj=dj)

                if windowflag
                    tmin = twindow[1]
                    tmax = twindow[2]

                    # time checks
                    if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                        println("Error: please input right time limits in the time window!")
                    else
                        # trim time vector
                        t_ind = findall(x->(x≤tmax && x≥tmin), t)
                    end
                else
                    tmin = twindow[1]
                    tmax = twindow[1] + nwindow*(1/fmin + 1/fmax)/2.

                    # time checks
                    if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                        println("Error: please input right time limits in the time window!")
                    else
                        # trim time vector
                        t_ind = findall(x->(x≤tmax && x≥tmin), t)
                    end
                end
                # trim vectors to correspond to time
                wt = t[t_ind]
                window = collect(1:length(wt))
                rcwt1 = real.(icwt1)
                rcwt2 = real.(icwt2)
                wcwt1 = rcwt1[t_ind]
                wcwt2 = rcwt2[t_ind]

                # normalize both signals, if appropriate
                if normalize
                    ncwt1 = ((wcwt1 .- mean(wcwt1)) ./ std(wcwt1))[:]
                    ncwt2 = ((wcwt2 .- mean(wcwt2)) ./ std(wcwt2))[:]
                else
                    ncwt1 = wcwt1[:]
                    ncwt2 = wcwt2[:]
                end

                # use stretching to extract dv/v
                (dvv[iband], cc[iband], cdp[iband], eps, err[iband], allC) = stretching(ncwt1, ncwt2, wt, window, fmin, fmax, dvmin=-dvmax, dvmax=dvmax)
            end
        end
    end

    return freqbands, dvv, cc, cdp, err
end

function wts_allfreqs(cur, ref, t, twindow, fwindow, dj, s0, J; dvmax::Float64=0.1, ndv::Int64=100, wvn::String="morlet", normalize::Bool=true, windowflag::Bool=true, nwindow::Float64=1.5)
    # define sample frequency
    dt = t[2] - t[1]
    fs = 1/dt
    (fmin, fmax) = fwindow

    # apply cwt to two traces
    # define wavelet
    if wvn=="morlet"
        wav = WT.Morlet(6)
    end
    # calculate the CWT of the time series, using identical parameters for both calculations
    cwt1, sj, freqs, coi = cwt(cur, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)
    cwt2, sj, freqs, coi = cwt(ref, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)

    # extract real part of cwt
    rcwt1 = real.(cwt1)
    rcwt2 = real.(cwt2)

    # zero out cone of influence and data outside frequency band
    if (fmax > maximum(freqs) || fmax<=fmin)
        println("Error: please input correct frequency limits in the frequency window!")
    else
        freq_indout = findall(f->(f < fmin || f > fmax), freqs)
        rcwt1[freq_indout, :] .= 0.
        rcwt2[freq_indout, :] .= 0.
        freq_indin = findall(f->(f >= fmin && f <= fmax), freqs)

        # Use stretching method to extract dvv
        nfreq=length(freqs)

        dvv = zeros(Float32, nfreq)
        cc = zeros(Float32, nfreq)
        cdp = zeros(Float32, nfreq)
        err = zeros(Float32, nfreq)

        for ifreq in freq_indin
            if windowflag
                tmin = twindow[1]
                tmax = twindow[2]

                # time checks
                if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                    println("Error: please input right time limits in the time window!")
                else
                    # trim time vector
                    t_ind = findall(x->(x≤tmax && x≥tmin), t)
                end
            else
                tmin = twindow[1]
                tmax = twindow[1] + nwindow*(1/fmin + 1/fmax)/2.

                # time checks
                if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                    println("Error: please input right time limits in the time window!")
                else
                    # trim time vector
                    t_ind = findall(x->(x≤tmax && x≥tmin), t)
                end
            end

            # prepare time axis and its indices
            wt = t[t_ind]
            it = collect(1:length(wt))
            # prepare windowed data
            wcwt1 = rcwt1[ifreq, t_ind]
            wcwt2 = rcwt2[ifreq, t_ind]

            # normalize both signals if appropriate
            # normalize both signals, if appropriate
            if normalize
                ncwt1 = ((wcwt1 .- mean(wcwt1)) ./ std(wcwt1))[:]
                ncwt2 = ((wcwt2 .- mean(wcwt2)) ./ std(wcwt2))[:]
            else
                ncwt1 = wcwt1[:]
                ncwt2 = wcwt2[:]
            end

            # compute the stretching
            window = collect(1:length(wt))
            # use stretching to extract dv/v
            (dvv[ifreq], cc[ifreq], cdp[ifreq], eps, err[ifreq], allC) = stretching(ncwt1, ncwt2, wt, window, fmin, fmax, dvmin=-dvmax, dvmax=dvmax)
        end
    end

    return freqs, dvv, cc, cdp, err
end

function wtdtw_freqbands(cur, ref, t, twindow, freqbands, maxLag, b, direction, dj, s0, J; wvn::String="morlet", normalize::Bool=true, windowflag::Bool=true, nwindow::Float64=1.5)
    # define sample frequency
    dt = t[2] - t[1]
    fs = 1/dt

    # apply cwt to two traces
    # define wavelet
    if wvn=="morlet"
        wav = WT.Morlet(6)
    end
    # calculate the CWT of the time series, using identical parameters for both calculations
    cwt1, sj, freqs, coi = cwt(cur, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)
    cwt2, sj, freqs, coi = cwt(ref, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)

    (nband, nfreq) = size(freqbands)
    if nband==0 || nfreq != 2
        println("Error: please check if inputs in freqbands is right!")
    else
        dvv = zeros(Float32, nband)
        err = zeros(Float32, nband)

        # iterate over frequency bands
        for iband=1:nband
            (fmin, fmax) = freqbands[iband, :]

            if (fmax > maximum(freqs) || fmax <= fmin)
                println("Error: please ensure columns 1 and 2 are right frequency limits in freqbands!")
            else
                freq_ind = findall(f->(f>=fmin && f<=fmax), freqs)

                # perform icwt
                icwt1 = icwt(cwt1[freq_ind, :], CFW(wav,1/dj), sj[freq_ind], dt=dt, dj=dj)
                icwt2 = icwt(cwt2[freq_ind, :], CFW(wav,1/dj), sj[freq_ind], dt=dt, dj=dj)

                if windowflag
                    tmin = twindow[1]
                    tmax = twindow[2]

                    # time checks
                    if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                        println("Error: please input right time limits in the time window!")
                    else
                        # trim time vector
                        t_ind = findall(x->(x≤tmax && x≥tmin), t)
                    end
                else
                    tmin = twindow[1]
                    tmax = twindow[1] + nwindow*(1/fmin + 1/fmax)/2.

                    # time checks
                    if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                        println("Error: please input right time limits in the time window!")
                    else
                        # trim time vector
                        t_ind = findall(x->(x≤tmax && x≥tmin), t)
                    end
                end
                # trim vectors to correspond to time
                wt = t[t_ind]
                window = collect(1:length(wt))
                rcwt1 = real.(icwt1)
                rcwt2 = real.(icwt2)
                wcwt1 = rcwt1[t_ind]
                wcwt2 = rcwt2[t_ind]

                # normalize both signals, if appropriate
                if normalize
                    ncwt1 = ((wcwt1 .- mean(wcwt1)) ./ std(wcwt1))[:]
                    ncwt2 = ((wcwt2 .- mean(wcwt2)) ./ std(wcwt2))[:]
                else
                    ncwt1 = wcwt1[:]
                    ncwt2 = wcwt2[:]
                end

                (stbarTime, stbar, dist, error) = dtwdt(ncwt1, ncwt2, dt, maxLag=maxLag, b=b, direction=direction)

                # perform linear regression
                model = glm(@formula(Y ~0 + X),DataFrame(X=wt,Y=stbarTime),Normal(),
                            IdentityLink(),wts=ones(length(wt)))

                dvv[iband] = coef(model)[1]
                err[iband] = stderror(model)[1]
            end
        end
    end
    return freqs, dvv, err
end

function wtdtw_allfreqs(cur, ref, t, twindow, fwindow, maxLag, b, direction, dj, s0, J; wvn::String="morlet", normalize::Bool=true, windowflag::Bool=true, nwindow::Float64=1.5)
    # define sample frequency
    dt = t[2] - t[1]
    fs = 1/dt
    (fmin, fmax) = fwindow

    # apply cwt to two traces
    # define wavelet
    if wvn=="morlet"
        wav = WT.Morlet(6)
    end
    # calculate the CWT of the time series, using identical parameters for both calculations
    cwt1, sj, freqs, coi = cwt(cur, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)
    cwt2, sj, freqs, coi = cwt(ref, CFW(wav,1/dj), J1=J, dt=dt, s0=s0)

    # extract real part of cwt
    rcwt1 = real.(cwt1)
    rcwt2 = real.(cwt2)

    # zero out cone of influence and data outside frequency band
    if (fmax > maximum(freqs) || fmax<=fmin)
        println("Error: please input correct frequency limits in the frequency window!")
    else
        freq_indout = findall(f->(f < fmin || f > fmax), freqs)
        rcwt1[freq_indout, :] .= 0.
        rcwt2[freq_indout, :] .= 0.
        freq_indin = findall(f->(f >= fmin && f <= fmax), freqs)

        # Use stretching method to extract dvv
        nfreq=length(freqs)

        dvv = zeros(Float32, nfreq)
        err = zeros(Float32, nfreq)

        for ifreq in freq_indin
            if windowflag
                tmin = twindow[1]
                tmax = twindow[2]

                # time checks
                if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                    println("Error: please input right time limits in the time window!")
                else
                    # trim time vector
                    t_ind = findall(x->(x≤tmax && x≥tmin), t)
                end
            else
                tmin = twindow[1]
                tmax = twindow[1] + nwindow*(1/fmin + 1/fmax)/2.

                # time checks
                if tmin < minimum(t) || tmax > maximum(t) || tmax <= tmin
                    println("Error: please input right time limits in the time window!")
                else
                    # trim time vector
                    t_ind = findall(x->(x≤tmax && x≥tmin), t)
                end
            end

            # prepare time axis and its indices
            wt = t[t_ind]
            it = collect(1:length(wt))
            # prepare windowed data
            wcwt1 = rcwt1[ifreq, t_ind]
            wcwt2 = rcwt2[ifreq, t_ind]

            # normalize both signals if appropriate
            # normalize both signals, if appropriate
            if normalize
                ncwt1 = ((wcwt1 .- mean(wcwt1)) ./ std(wcwt1))[:]
                ncwt2 = ((wcwt2 .- mean(wcwt2)) ./ std(wcwt2))[:]
            else
                ncwt1 = wcwt1[:]
                ncwt2 = wcwt2[:]
            end

            (stbarTime, stbar, dist, error) = dtwdt(ncwt1, ncwt2, dt, maxLag=maxLag, b=b, direction=direction)

            # perform linear regression
            model = glm(@formula(Y ~0 + X),DataFrame(X=wt,Y=stbarTime),Normal(),
                        IdentityLink(),wts=ones(length(wt)))

            dvv[ifreq] = coef(model)[1]
            err[ifreq] = stderror(model)[1]
        end
    end

    return freqs, dvv, err
end
