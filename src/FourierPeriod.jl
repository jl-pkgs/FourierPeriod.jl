module FourierPeriod


using FFTW, DataFrames
using Ipaper, Dates
using NaNStatistics: nansum, nancumsum, nanrange

include("Fourier.jl")
include("moving_Fourier.jl")

export moving_Fourier, Fourier


end # module FourierPeriod
