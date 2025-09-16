using Dates


# 基本函数：地日距离、赤纬、日落时角
deg2rad(x) = x * π / 180
dr(doy::Int) = 1.0 + 0.033 * cos(2π * doy / 365)
# declination(doy::Int) = deg2rad(23.45) * sin(2π * (284 + doy) / 365)
declination(J::Int) = 0.409 * sin(2pi / 365 * J .- 1.39) # Allen, Eq. 24


function omega_sunset(phi::Float64, delta::Float64)
  x = -tan(phi) * tan(delta)
  if x ≤ -1      # 极昼
    return π
  elseif x ≥ 1   # 极夜
    return 0.0
  else
    return acos(x)
  end
end

# 由 DateTime 得到时角 ω (弧度)
function hour_angle(dt::DateTime)
  tsol = Dates.hour(dt) + Dates.minute(dt) / 60 + Dates.second(dt) / 3600
  return (tsol - 12) * (2π / 24)
end

# 平均 TOA Rs (W/m^2)，闭式解
function Rs_toa_mean(doy::Int, lat::Float64, w1::Float64, w2::Float64; Rg=1367.0)
  δ = declination(doy)
  φ = lat
  dr_ = dr(doy)

  ωs = omega_sunset(φ, δ)
  a = max(w1, -ωs)
  b = min(w2, ωs)
  a >= b && return 0.0

  term = cos(φ) * cos(δ) * (sin(b) - sin(a)) + sin(φ) * sin(δ) * (b - a)
  return Rg * dr_ * term / (w2 - w1)
end


function Rs_toa_mean(dt::DateTime, lat_deg::Float64; Rg=1367.0, l = -0.5, r=0.5)
  lat = deg2rad(lat_deg)
  doy = dayofyear(dt)
  ω = hour_angle(dt)
  Δl = l * 2π / 24   # 半区间弧度宽度
  Δr = r * 2π / 24   # 半区间弧度宽度
  Rs_toa_mean(doy, lat, ω + Δl, ω + Δr; Rg)
end


function Rs_toa_inst(dt::DateTime, lat_deg::Float64; Rg=1367.0)
  doy = dayofyear(dt)
  lat = deg2rad(lat_deg)
  δ = declination(doy)
  ω = hour_angle(dt)

  cosz = cos(lat) * cos(δ) * cos(ω) + sin(lat) * sin(δ)
  max(0, Rg * dr(doy) * cosz)
end


function Rs_toa_mean(dates::AbstractVector{DateTime}, lat_deg::Float64; Rg=1367.0, l=-0.5, r=0.5)
  map(date -> Rs_toa_mean(date, lat_deg; Rg, l, r), dates)
end

function Rs_toa_inst(dates::AbstractVector{DateTime}, lat_deg::Float64; Rg=1367.0)
  map(date -> Rs_toa_inst(date, lat_deg; Rg), dates)
end

export Rs_toa_inst, Rs_toa_mean
