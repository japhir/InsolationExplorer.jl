"""
    insolation(eccentricity, obliquity, lpx;
               longitude = pi/2, latitude = deg2rad(65),
               S0 = 1360.7, H = nothing)

Calculate insolation at the top of the atmosphere.

# Arguments
- `eccentricity::Real` the Earth's eccentricity.
- `obliquity::Real` the Earth's obliquity or tilt in radians.
- `lpx::Real` the longitude of perihelion with respect
   to the moving equinox in radians.

# Keyword Arguments
- `longitude::Real` is the true solar longitude in radians.
- `latitude::Real` is Earth's latitude.
- `S0::Real` the total solar irradiance.
- `H::Real` the sun hour angle in radians.

True solar longitude can be specified as:
- `pi/2`   for the Summer Solstice
- `pi`     for the Autumn Equinox
- `3pi/2`  for the Winter Solstice
- `0`      for the Spring Equinox

The `S0` value defualts to 1360.7 Wm¯² which is the same as the new PMIP4
transient simulations (Otto-Bliesner, 2017, Menviel, 2019). A value of 1368
Wm¯² was used in the La04 webtool, while 1361 Wm¯² was used in DeepMIP (Lunt et
al., 2017, Matthes et al., 2016).

# References

Crucifix, M. (2023). palinsol : A R package to compute Incoming Solar Radiation
(insolation) for palaeoclimate studies. Zenodo. [doi:
10.5281/zenodo.14893715](https://doi.org/10.5281/zenodo.14893715)

# Examples
Calculate 65°N summer insolation at t₀
```jldoctest
julia> insolation(0.016705, 0.4090928042223287, 1.7962486166737615;
                  longitude = pi/2, latitude = deg2rad(65),
                  S0 = 1360.7, H = nothing)
509.9881957773474
```

Apply the function to vectors of input arguments with `insolation.()`
```jldoctest
julia> insolation.([0.016705, 0.016846], [0.40909, 0.40995], [1.796, 1.684])
2-element Vector{Float64}:
 509.9863194415539
 511.3106728099949
```
"""
function insolation(eccentricity, obliquity, lpx;
                    longitude = pi/2,
                    latitude = deg2rad(65), # 65*pi/180,
                    S0 = 1360.7,
                    H = nothing)
    nu = longitude - lpx
    rho = (1 - eccentricity^2) / (1 + eccentricity * cos(nu))

    sindelta = sin(obliquity) * sin(longitude)
    cosdelta = sqrt(1 - sindelta^2)
    sinlatsindelta = sin(latitude) * sindelta
    coslatcosdelta = cos(latitude) * cosdelta

    if isnothing(H)
      cosH0 = min(max(-1, -sinlatsindelta / coslatcosdelta), 1)
      sinH0 = sqrt(1 - cosH0^2)
      H0 = acos(cosH0)
      insol = S0 / (pi * rho^2) * (H0 * sinlatsindelta + coslatcosdelta * sinH0)
    else
      insol = max(0, S0 / (rho^2) * (sinlatsindelta + coslatcosdelta * cos(H)))
    end

    return(insol)
end
