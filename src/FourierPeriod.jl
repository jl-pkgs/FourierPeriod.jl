module FourierPeriod


using FFTW, DataFrames
using Ipaper, Dates
using NaNStatistics: nansum, nancumsum, nanrange, nanstd

include("Fourier.jl")
include("moving_Fourier.jl")
include("Rs_toa.jl")


export moving_Fourier, Fourier

end # module FourierPeriod
