using FourierPeriod, Test, Dates

include("test-Fourier.jl")

@testset "moving_Fourier" begin
  t = DateTime(2010, 8, 1):Hour(1):DateTime(2010, 8, 10)
  t = t[1:100]
  y = rand(100)
  @test_nowarn moving_Fourier(y, t)
end


@testset "Rs_toa" begin
  date = DateTime(2010, 8, 1, 15)
  lat = 30.0

  Rs_toa_mean(date, lat; l=0, r=1) < Rs_toa_inst(date, lat)
  @test Rs_toa_inst(date, lat) ≈ 977.4907861007081
  @test Rs_toa_inst(date, lat) > Rs_toa_mean(date, lat)
  @test Rs_toa_inst(date, lat) - Rs_toa_mean(date, lat) < 2.3

  # Vector
  dates = DateTime(2010, 8, 1):Hour(1):DateTime(2010, 8, 3)
  @test_nowarn Rs_toa_inst(dates, lat)
  @test_nowarn Rs_toa_mean(dates, lat)
end
