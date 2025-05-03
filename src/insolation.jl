function insolation(eccentricity, obliquity, lpx;
                    longitude = pi/2,
                    latitude = 65 * pi / 180,
                    S0 = 1360.7, # PMIP4 transient Otto-Bliesner 2017, Menviel 2019
                    # 1368, # La04 webtool
                    # 1361, # based on Lunt et al., 2017 DeepMIP => Matthes et al., 2016,
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
