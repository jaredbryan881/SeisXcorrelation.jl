export compute_fft_c3
"""
    compute_fft(S::CorrData,freqmin::Float64,freqmax::Float64, fs::Float64
                         cc_step::Int, cc_len::Int;
                         time_norm::Union{Bool,String}=false,
                         to_whiten::Bool=false)
Computes windowed fft of cross-correlated ambient noise data.
Cross-correlates data using either cross-correlation, deconvolution or
cross-coherence.
# Arguments
- `C::CorrData`: CorrData structure.
- `freqmin::Float64`: minimun frequency for whitening.
- `freqmax::Float64`: maximum frequency for whitening.
- `fs::Float64`: Sampling rate to downsample `S`.
- `cc_step::Int`: time, in seconds, between successive cross-correlation windows.
- `cc_len::Int`: length of noise data window, in seconds, to cross-correlate.
- `time_norm::Union{Bool,String}`: time domain normalization to perform.
- `to_whiten::Bool`: Apply whitening in frequency domain.
"""
function compute_fft_c3(C::CorrData,freqmin::Float64,freqmax::Float64,fs::Float64,
                       cc_step::Int, cc_len::Int;
                       time_norm::Union{Bool,String}=false,
                       to_whiten::Bool=false)
    # window data
    corr_neg, corr_pos = slide_c3(C, cc_len, cc_step)

    # compute FFT on windowed data for positive and negative windows
    FFT_neg = process_fft(corr_neg, freqmin, freqmax, fs, time_norm=time_norm,
                          to_whiten=to_whiten)
    FFT_pos = process_fft(corr_pos, freqmin, freqmax, fs, time_norm=time_norm,
                          to_whiten=to_whiten)
    FFT = hcat(FFT_neg, FFT_pos)

    return F = FFTData(C.name, C.id,
                       C.loc, C.fs, C.gain, freqmin, freqmax,
                       cc_len, cc_step, to_whiten, time_norm, C.resp,
                       C.misc, C.notes, C.t, FFT)
end
