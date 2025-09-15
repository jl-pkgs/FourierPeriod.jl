using FourierPeriod, Test, Dates

include("test-Fourier.jl")

@testset "moving_Fourier" begin
  t = DateTime(2010, 8, 1):Hour(1):DateTime(2010, 8, 10)
  t = t[1:100]
  y = rand(100)
  @test_nowarn moving_Fourier(y, t)
end
