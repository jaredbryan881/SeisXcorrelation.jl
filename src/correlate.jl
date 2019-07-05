export compute_cc_c3

"""
    compute_cc(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;
               smoothing_half_win::Int=20,
               corr_type::String="cross-correlation" )

"""
function compute_cc_c3(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;
                    smoothing_half_win::Int=20,
                    corr_type::String="cross-correlation")

    N = convert(Int,round(FFT1.cc_len * FFT1.fs)) # number of data points
    comp = FFT1.name[end] * FFT2.name[end]

    # assume perfect overlap of dates for this cross-correlation
    corr = correlate(FFT1.fft, FFT2.fft, N,
                     convert(Int,round(maxlag * FFT1.fs)),
                     corr_type=corr_type)
    rotated = false

    return CorrData(FFT1, FFT2, comp, rotated, corr_type,
                    maxlag, FFT1.t, corr)

end
