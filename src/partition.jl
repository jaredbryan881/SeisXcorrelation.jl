using SeisIO, SeisNoise
export partition

"""
    partition(A::Array{Float64,2}, offset::Int64, size::Int64)

Partition an array into two sub-arrays, centered on the midpoint.

# Arguments
- `A::Array{Float32,2},`    : Input array to be partitioned
- `offset::Int64`    : Number of samples away from the midpoint of the array to start the window.
- `size::Int64`    : Length of the window in number of samples.

# Output
- `A`::Array{Float32,2}    : Partitioned array.
"""
function partition(A::Array{Float32,2}, offset::Int64, size::Int64)
    # midpoint of the array
    zerolag = convert(Int64, (length(A[:,1])+1)/2)
    # define windows
    pos_win = zerolag+offset:zerolag+offset+size
    neg_win = zerolag-offset-size:zerolag-offset
    # stack subarrays
    A = hcat(A[neg_win, :], A[pos_win, :])

    return A
end

"""
    slide(C::CorrData, cc_len::Int64, cc_step::Int64)

Generate equal length sliding windows into a CorrData object for both positive and negative codas.

# Arguments
- `C::CorrData,`    : Input array to be partitioned
- `cc_len::Int64`    : Number of samples away from the midpoint of the array to start the window
- `cc_step::Int64`    : Length of the window in number of samples

# Output
- `corr_neg::Array{Float32,2}`    : Partitioned array of negative-lag data
- `corr_pos::Array{Float32,2}`    : Partitioned array of positive-lag data
"""
function slide_c3(C::CorrData, cc_len::Int64, cc_step::Int64)
    window_samples = Int(cc_len * C.fs)
    t=[1 0; length(C.corr[:, 1]) 0]
    su,eu = SeisIO.t_win(t, C.fs) * 1e-6
    starts = Array(range(su,stop=eu,step=cc_step))
    ends = starts .+ cc_len .- 1. / C.fs
    ind = findlast(x -> x <= eu,ends)
    starts = starts[1:ind]
    ends = ends[1:ind]

    # fill arrays with overlapping windows
    corr_neg = Array{Float64,2}(undef, window_samples,length(starts))
    corr_pos = Array{Float64,2}(undef, window_samples,length(starts))
    s = convert.(Int,round.((hcat(starts,ends) .- su) .* C.fs .+ 1.))
    for ii = 1:length(starts)
        corr_neg[:,ii] = C.corr[s[ii,1]:s[ii,2], 1]
        corr_pos[:,ii] = C.corr[s[ii,1]:s[ii,2], 2]
    end

    return corr_neg, corr_pos
end

"""
    window(A::Array{R,1}, win_len::Int64, win_step::Int64) where R<:Real

Generate equal length sliding windows into an array.

# Arguments
- `A::Array{R,1} where R<:Real,`    : Input array to be windowed
- `win_len::Int64`    : Number of samples in each window
- `win_step::Int64`    : Step size for successive windows

# Output
- `Awin::Array{R,2} where R<:Real`    : Windowed array
"""
function window(A::Array{R,1}, win_len::Int64, win_step::Int64) where R<:Real
    # left-most indices of each window
    minind = 1:win_len:length(A) - win_len

    # number of windows
    N = length(minind)

    # initialize array of size (win_len, nWindows)
    Awin = zeros(win_len, N)
    # fill matrix
    for i=1:N
        Awin[:,i] = A[minind[i]:minind[i]+win_len-1]
    end

    return Awin
end
