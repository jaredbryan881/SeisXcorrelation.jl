export remove_nanandzerocol


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