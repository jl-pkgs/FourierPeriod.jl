using Ipaper: @par


"""
# return
`:date, :freq, :explain, :accu_explain, :sd, :A`
"""
function moving_Fourier(y::AbstractVector, t::AbstractVector;
  Δt=1 / 24, win_half::Int=3, m::Int=3, threshold=0.7)

  dates = Date.(t) |> unique_sort
  n = length(dates)

  res = Vector{Any}(undef, n)
  # @showprogress 
  @inbounds @par for i in 1:n
    ibeg = max(i - win_half, 1)
    iend = min(i + win_half, n)

    time_beg = DateTime(dates[ibeg])
    time_end = DateTime(dates[iend]) + Hour(23)
    inds_win = (time_beg .<= t .<= time_end) |> findall

    _y = @view y[inds_win]
    _t = @view t[inds_win]

    inds_1d = DateTime(dates[i]) .<= _t .<= DateTime(dates[i]) + Hour(23)
    _ypred, info = Fourier(_y, Δt; threshold, m)

    info.date .= dates[i]
    info.sd .= nanstd(_y)
    info.A .= nanrange(@view(_y[inds_1d]))

    res[i] = info[:, Cols(:date, :freq, :explain, :accu_explain, :sd, :A)]
    # push!(res, )
  end
  vcat(res...)
end
