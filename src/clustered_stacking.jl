using Dates, FileIO, SeisIO, SeisNoise, PlotlyJS, Statistics, Printf
using IterTools
using PyCall
using PyPlot

@pyimport scipy.cluster.hierarchy as hi
@pyimport scipy.spatial.distance as di
@pyimport scipy.signal as sig

include("stacking.jl")

# set up plot
# layout = Layout(width=1200, height=800,
#                 xaxis_title="Lags [s]",
#                 yaxis_title="Days",
#                 font=attr(size=14),
#                 showlegend=false)
#
# p = PlotlyJS.plot([NaN], layout)
#
# distance_threshold = 0.8

"""
    clustered_selective_stacking(data::CorrData, reference::CorrData; threshold::Float64=0.0, slice::Union{Bool, Float64, Array{Float64,1}}=false, metric::String="cc", win_len::Float64=10.0, win_step::Float64=5.0)

Stack the windows in a CorrData object that exceed a correlation-coefficient threshold with respect to a reference.

# Arguments
- `data::CorrData,`    : Input CorrData matrix used to define the reference and each window
- `reference::CorrData`    : Reference cross-correlation function to use in computing the correlation coefficient
- `threshold::Float64,`    : Value of correlation-coefficient below which we zero windows
- `slice::Union{Bool, Float64, Array{Float64,1}}`    : whether to slice the cross-correlations before computing convergence. If so, a single Float64 specifies the
                                                       lag at which to slice to keep just the ballistic wave and an array of Float64 specifies the start and end lag
                                                       to keep only the coda.

# Output
- `stackedData::CorrData,`    : Slectively stacked data
- `cList::Array{Float64,1}`    : Correlation coefficient of each window with respect to the reference
"""
function clustered_selective_stacking(data::CorrData, reference::CorrData; threshold::Float64=0.0, slice::Union{Bool, Float64, Array{Float64,1}}=false, metric::String="cc", coh_win_len::Float64=10.0, coh_win_step::Float64=5.0, filter::Union{Bool, Array{Float64,1}}=false)
    # slice data if given a time window
    # TODO: find a better way to pass arguments to this function that only apply to coh, or only to cc. e.g., win_len has no meaning if metric=="cc"
    #if typeof(slice) != Bool
    ref = copy(reference)
    d_origin = copy(data)

    # time vector
    tvec = -reference.maxlag:1/reference.fs:reference.maxlag
    # find all times within [-slice,slice] time
    if slice ==false
        slice = 0.0
    end

    if typeof(slice) != Bool
        if typeof(slice)==Float64
            # find all times within [-slice,slice] time
            slice_inds = findall(x->(x<=-slice || slice <= x), tvec)
        elseif typeof(slice)==Array{Float64,1}
            # convert startlag/endlag[s] to startlag/windowlength[samples]
            slice_inds = findall(x->(-slice[2]<=x || x<=-slice[1] || slice[1]<=x || x<=slice[2]), tvec)
        else
            println("Please choose an allowable slicing operation. Exiting.")
            exit()
        end
    else
        # default to full reference cross-correlation
        tvec = -reference.maxlag:1/reference.fs:reference.maxlag
        slice_inds = findall(x->(x<=slice ||slice<= x), tvec)
    end

    #-------------------------------------------#
    # Work flow
    # 1. clustered by N groups ( N = 2-15 according to computation power)
    # 2. stacking them within the cluster
    # 3. grid search for best combination of stacking by taking correlation between reference
    # 4. select a set of cc function
    #-------------------------------------------#
    plot=false
    wiener = false
    # 1. clustering cc functions
    # distance is computed by pdist correlation coefficient

    if metric == "cc"

        # before doing clustering, remove incoherent shorttime window ccs by correlation coefficient like canonical selective stacking
        cList = get_cc(d_origin, ref)
        good_fit = findall(x->(x>=threshold), cList)
        d = copy(d_origin)
        #update d::CorrData only with high correlation coefficient
        d.corr = d_origin.corr[:, good_fit]

        pdis = di.pdist(transpose(d.corr[slice_inds,:]), metric="correlation") #'euclidean' 'seuclidean' 'matching' correlation
        #pdis = []
        #numofcc = size(d.corr)[2]
        #println(numofcc)
        # compute pdist array [d(2,1), d(3,1), ..., d(n, 1), d(3,2), ..., d(n, n-1)]
        # @simd for j = 1:numofcc-1
        #     print(j)
        #     @simd for i = j+1:numofcc
        #         # u1 = d.corr[:,i]
        #         # u2 = d.corr[:,j]
        #         # dt = 1/C.fs
        #         # b = 10
        #         # direction = 1
        #         # (_, _, _, error) = dtwdt(u1, u2, dt, maxLag=round(Int64, 0.05*C.maxlag*C.fs), b=b, direction=direction)
        #         dtemp = 1 - cor(d.corr[:,i], d.corr[:,j]) #https://www.researchgate.net/publication/265604603_Minimum_Pearson_Distance_Detection_for_Multilevel_Channels_With_Gain_andor_Offset_Mismatch
        #         push!(pdis, dtemp)
        #     end
        # end

    else
        error("Please specify metric to calculate distance array.")
    end

    Z = hi.linkage(pdis, method="ward") #method="single"

    #normalize distance
    Z[:,3] = Z[:,3]./maximum(Z[:,3])
    #println(Z)
    # this is used for dendrogram plotting
    maxclustnum = 12 # this is equal to the num of grid searching windows

    if plot
        PyPlot.figure(figsize=(12,8))
        fig = hi.dendrogram(Z, p=maxclustnum, truncate_mode="lastp")
        PyPlot.show()
        str = d.id
        PyPlot.savefig("./fig/dendrogram_cc_$str.png", dpi=300)
        PyPlot.close()
    end
    #clustering
    #group = hi.fcluster(Z, distance_threshold, criterion="distance")
    group = hi.fcluster(Z, maxclustnum, criterion="maxclust")

    #println(group)

    #stack each clusters

    #---plotting parameters ---#

    plot_xlag = length(d.corr[:,1]) # xlag of cc function
    plot_maxgroupnum = maxclustnum # max num of plotting group
    plot_maxynum = 15 # max num of traces in y axis
    xlabel = "clusters"
    ylabel = "clustered cc functions"

    if plot;
        layout = Layout(width=1400, height=800,
                xaxis_title=xlabel,
                yaxis_title=ylabel,
                font =attr(size=18),
                showlegend=false,
                title="Clustered cross-correlation function",
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)")
        p = PlotlyJS.plot([NaN], layout);
    end

    Stackedcc = copy(d)
    Stackedcc.corr = zeros(length(d.corr[:,1]), maximum(group))

    numtrace = 0

    #sort by correlation coefficient
    corrcoef = []
    for i = 1:maximum(group)
        idx = findall(x -> (x==i), group)
        tempData = copy(d)
        tempData.corr = tempData.corr[:, idx]
        stackedtemp = stack(tempData, allstack=true)
        Stackedcc.corr[:,i] = stackedtemp.corr[:,1]
        push!(corrcoef, cor(stackedtemp.corr[slice_inds,1], ref.corr[slice_inds,1]))
        #plot cc fnctions
    end

    sortidx = sortperm(corrcoef, rev=true)
    #print("corrcoef is: ")
    #println(corrcoef)

    #plot with sorted corr
    if plot
        traceall = GenericTrace{Dict{Symbol,Any}}[]
        for (i, id) = enumerate(sortidx)
            if id <= plot_maxgroupnum
                #idx = findall(x -> (x==id), group)
                idx = findall(x -> (x==i), group)
                tempData = copy(d)
                tempData.corr = tempData.corr[:, idx]
                stackedtemp = stack(tempData, allstack=true)

                # plot each cc function
                lags = (collect(1:plot_xlag) .+ (i-1)*(1.5 * plot_xlag)) ./ tempData.fs
                for j = 1:minimum([length(idx), plot_maxynum])
                    plot_ccfunc = 0.5 * tempData.corr[:, j] ./ maximum(tempData.corr[:, j]) .+ float(plot_maxynum - j)
                    trace1 = PlotlyJS.scatter(;x=lags, y=plot_ccfunc, mode="lines", line_color="rgb(0,0,0)")
                    #println(trace1)
                    #println(trace1.fields)
                    #PlotlyJS.addtraces!(p, trace1) # this may cause bug
                    push!(traceall, trace1)
                end
                # plot stacked trace
                plot_stackedccfunc = 0.5 * stackedtemp.corr[:, 1] ./ maximum(stackedtemp.corr[:, 1]) .+ float(plot_maxynum)
                #trace2 = PlotlyJS.scatter(;x=lags, y=plot_stackedccfunc, mode="lines+text",
                # line_color="rgb(255,0,0)", name="stacked trace", textposition="top center", text=[@sprintf("%4.3f",corrcoef[id])])
                trace2 = PlotlyJS.scatter(;x=lags, y=plot_stackedccfunc, mode="lines+text",
                 line_color="rgb(255,0,0)", name="stacked trace")
                push!(traceall, trace2)
                #PlotlyJS.addtraces!(p, trace2) # this may cause bug
            end
        end
        p = PlotlyJS.plot(traceall, layout)
        display(p)
        str = d.id
        PlotlyJS.savefig(p, "./fig/clustered_ccfunctions_$str.png")

    end

    if plot
        #plot reference
        lags = collect(1:plot_xlag)./ ref.fs
        plot_ref = 0.5 * ref.corr[:, 1] ./ maximum(ref.corr[:, 1])
        trace1 = PlotlyJS.scatter(;x=lags, y=plot_ref, mode="lines", line_color="rgb(0,0,255)", name="reference")
        p = PlotlyJS.plot(trace1, layout)
        display(p)
        str = d.id
        PlotlyJS.savefig(p, "./fig/refference_$str.png")
    end

    # performing grid search to find best combination

    # init coeff is of all stack with correlation coefficient > threshold
    stackedalltemp = stack(d, allstack=true)

    init_coef = cor(stackedalltemp.corr[slice_inds,1], ref.corr[slice_inds,1])
    maxcrosscoef = init_coef
    bestgroups = collect(1:maximum(group))
    #println(corrcoef)
    totalcalcnum = length(IterTools.subsets(1:maximum(group)))

    for (i, subs) = enumerate(IterTools.subsets(1:maximum(group)))
        if isempty(subs); continue; end

        if mod((totalcalcnum - i), 1000) == 0; println(totalcalcnum - i); end

        tempData = copy(Stackedcc)
        tempData.corr = tempData.corr[:, subs]
        stackedtemp = stack(tempData, allstack=true)
        # evaluate
        trial_crosscoef = cor(stackedtemp.corr[slice_inds,1], ref.corr[slice_inds,1])
        if trial_crosscoef > maxcrosscoef
            # crosscoef is increased, so update info
            #print("coef update: ")
            #println([maxcrosscoef, trial_crosscoef])
            bestgroups = subs
            maxcrosscoef = trial_crosscoef
        end
    end

    println([init_coef, maxcrosscoef])

    # prepare stacked data
    #println(bestgroups)
    stackedData = copy(d)
    stackedData.corr = zeros(size(d.corr)[1],1)
    for mergegroupid = bestgroups
        #println(mergegroupid)
        idx = findall(x -> (x==mergegroupid), group)
        tempData = copy(d)
        tempData.corr = tempData.corr[:, idx]
        stackedtemp = stack(tempData, allstack=true)
        stackedData.corr += stackedtemp.corr
    end

    if wiener
        # apply wiener filter before clustering
        wienerwindowsize = [3,3]
        d.corr = sig.wiener(stackedData.corr, wienerwindowsize)
    end

    return stackedData, init_coef, maxcrosscoef

end
