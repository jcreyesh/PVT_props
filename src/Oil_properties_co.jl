# ---- Compresibilidad del aceite - Co ----
# Correlación de de Vásquez, M.E. y Beggs, H.D.
""" Correlación de Vásquez, M.E. y Beggs, H.D. - 1980 \n
co\\_Vasquez\\_Beggs\\_sp(Py, Ty, Psp, Tsp, Pb, API, Ge_gas) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Psp -> Presión en el separador - (psia) \n
Tsp -> Temperatura en el separador - (°F) \n
Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Vasquez_Beggs_sp(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64

    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end

    Gec_gas = Ge_gas * (1 + 5.912e-5*API*Tsp*log10(Psp/114.7))
    rs = c1*Gec_gas*Pb^c2*exp(c3*API/(Ty + 460))
    co = (-1433 + 5*rs + 17.2*Ty - 1180*Gec_gas + 12.61*API)/(Py*10^5)

    return round(co, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de de Vásquez, M.E. y Beggs, H.D.
""" Correlación de Vásquez, M.E. y Beggs, H.D. - 1980 \n
co\\_Vasquez\\_Beggs(Py, Ty, Pb, API, Ge_gas) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Vasquez_Beggs(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64

    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else 
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end

    rs = c1*Ge_gas*Pb^c2*exp(c3*API/(Ty + 460))
    co = (-1433 + 5*rs + 17.2*Ty - 1180*Ge_gas + 12.61*API)/(Py*10^5)

    return round(co, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de de Petrosky, G.E. y Farshad, F.F.
""" Correlación de Petrosky, G.E. y Farshad, F.F. - 1993 \n
co\\_Petrosky\\_Farshad(Py, Ty, Pb, API, Ge_gas) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Petrosky_Farshad(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64

    rs = ((Ge_gas^0.8439) * (Pb/112.727 + 12.34) * 10^(7.916e-4*API^1.541 - 4.561e-5*Ty^1.3911))^1.73184
    co = (1.705e-7*rs^0.69357)*(Ge_gas^0.1885)*(API^0.3272)*(Ty^0.6729)*(Py^-0.5906)

    return round(co, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T. y Schmidt, Z.
""" Correlación de Kartoatmodjo, T. y Schmidt, Z. - 1994 \n
co\\_Kartoatmodjo\\_Schmidt\\_sp(Py, Ty, Psp, Tsp, Pb, API, Ge_gas) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Psp -> Presión en el separador - (psia) \n
Tsp -> Temperatura en el separador - (°F) \n
Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Kartoatmodjo_Schmidt_sp(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64

    Gec_gas = Ge_gas * (1 + 0.1595*API^0.4078*Tsp^-0.2466*log10(Psp/114.7))

    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end

    rs = c1 * (Gec_gas^c2) * Pb^(1/c4) *10^(c3*API/(Ty + 460))

    co = (6.8257*rs^0.5002)*(API^0.3613)*(Ty^0.76606)*(Gec_gas^0.35505)/(Py*10^6)

    return round(co, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T. y Schmidt, Z.
""" Correlación de Kartoatmodjo, T. y Schmidt, Z. - 1994 \n
co\\_Kartoatmodjo\\_Schmidt(Py, Ty, Pb, API, Ge_gas) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Kartoatmodjo_Schmidt(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64

    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end

    rs = c1 * (Ge_gas^c2) * Pb^(1/c4) *10^(c3*API/(Ty + 460))

    co = (6.8257*rs^0.5002)*(API^0.3613)*(Ty^0.76606)*(Ge_gas^0.35505)/(Py*10^6)

    return round(co, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de McCain, W.D.Jr., Rollins, J.B. y Villena-Lanzi, A.J.
""" Correlación de McCain, W.D.Jr., Rollins, J.B. y Villena-Lanzi, A.J. - 1988 \n
co\\_McCain\\_Rollins\\_Villena\\_1(Py, Ty, Pb, Rsb, API) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
Rsb -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_McCain_Rollins_Villena_1(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, Rsb::Union{Int64, Float64}, API::Union{Int64, Float64})::Float64

    ln_co = -7.573 - 1.45*log(Py) - 0.383*log(Pb) + 1.402*log(Ty + 460) + 0.256*log(API) + 0.449*log(Rsb)
    co = exp(ln_co)

    return round(co, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de McCain, W.D.Jr., Rollins, J.B. y Villena-Lanzi, A.J.
""" Correlación de McCain, W.D.Jr., Rollins, J.B. y Villena-Lanzi, A.J. - 1988 \n
co\\_McCain\\_Rollins\\_Villena\\_2(Py, Ty, Rsb, API) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rsb -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_McCain_Rollins_Villena_2(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rsb::Union{Int64, Float64}, API::Union{Int64, Float64})::Float64

    ln_co = -7.663 - 1.497*log(Py) + 1.115*log(Ty + 460) + 0.533*log(API) + 0.184*log(Rsb)
    co = exp(ln_co)

    return round(co, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de McCain, W.D.Jr., Rollins, J.B. y Villena-Lanzi, A.J.
""" Correlación de McCain, W.D.Jr., Rollins, J.B. y Villena-Lanzi, A.J. - 1988 \n
co\\_McCain\\_Rollins\\_Villena\\_3(Py, Ty, API, Ge_gas) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_McCain_Rollins_Villena_3(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64

    ln_co = -7.114 - 1.394*log(Py) + 0.981*log(Ty + 460) + 0.77*log(API) + 0.446*log(Ge_gas) 
    co = exp(ln_co)

    return round(co, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Standing, M.B.
""" Correlación de Standing M.B. - 1977 \n
co\\_Standing(Py, Ty, Rs, API, z, Ge_gas, Bo) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rs -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @yac/bl @std) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Standing(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rs::Union{Int64, Float64}, API::Union{Int64, Float64}, z::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Bo::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)
    Psc = 756.8 - 131*Ge_gas - 3.6*Ge_gas^2
    Tsc = 169.2 + 349.5*Ge_gas - 74*Ge_gas^2
    Psr = Py/Psc
    Tsr = (Ty + 460)/Tsc

    Bg = 50.3e-4*(z*(Ty + 460)/Py)
    dBo = 14.4e-5*(Ge_gas/Ge_oil)^0.5 * (Rs*(Ge_gas/Ge_oil)^0.5 + 1.25*Ty)^0.2
    dRs = Rs/(0.83*Py + 21.15)
    co = (-1/Bo)*(dRs)*(dBo - Bg)

    return round(co, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de Vásquez, M.E. y Beggs, H.D.
""" Correlación de Vásquez, M.E. y Beggs, H.D. - 1980 \n
co\\_Vasquez\\_Beggs\\_sat\\_sp(Py, Ty, Psp, Tsp, Rs, API, z, Ge_gas, Bo) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Psp -> Presión en el separador - (psia) \n
Tsp -> Temperatura en el separador - (°F) \n
Rs -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @yac/bl @std) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Vasquez_Beggs_sat_sp(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64, Float64}, Rs::Union{Int64, Float64}, API::Union{Int64, Float64}, z::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Bo::Union{Int64, Float64})::Float64

    if API <= 30
        c1, c2, c3 = 4.677e-4, 1.751e-5, 1.8106e-8
    else
        c1, c2, c3 = 4.67e-4, 1.1e-5, 1.337e-9
    end
    
    Gec_gas = Ge_gas * (1 + 5.912e-5*API*Tsp*log10(Psp/114.7))

    Psc = 756.8 - 131*Gec_gas - 3.6*Gec_gas^2
    Tsc = 169.2 + 349.5*Gec_gas - 74*Gec_gas^2
    Psr = Py/Psc
    Tsr = (Ty + 460)/Tsc

    Bg = 50.3e-4*(z*(Ty + 460)/Py)
    dBo = c1 + c3*(Ty - 60)*(API/Gec_gas)
    dRs = c2*(Rs/Py)
    co = (-1/Bo)*(dRs)*(dBo - Bg)

    return round(co, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de Vásquez, M.E. y Beggs, H.D.
""" Correlación de Vásquez, M.E. y Beggs, H.D. - 1980 \n
co\\_Vasquez\\_Beggs\\_sat(Py, Ty, Rs, API, z, Ge_gas, Bo) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rs -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @yac/bl @std) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Vasquez_Beggs_sat(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rs::Union{Int64, Float64}, API::Union{Int64, Float64}, z::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Bo::Union{Int64, Float64})::Float64

    if API <= 30
        c1, c2, c3 = 4.677e-4, 1.751e-5, 1.8106e-8
    else
        c1, c2, c3 = 4.67e-4, 1.1e-5, 1.337e-9
    end

    Psc = 756.8 - 131*Ge_gas - 3.6*Ge_gas^2
    Tsc = 169.2 + 349.5*Ge_gas - 74*Ge_gas^2
    Psr = Py/Psc
    Tsr = (Ty + 460)/Tsc

    Bg = 50.3e-4*(z*(Ty + 460)/Py)
    dBo = c1 + c3*(Ty - 60)*(API/Ge_gas)
    dRs = c2*(Rs/Py)
    co = (-1/Bo)*(dRs)*(dBo - Bg)

    return round(co, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de Glaso, O.
""" Correlación de Glaso, O. - 1980 \n
co\\_Glaso(Py, Ty, Rs, API, z, Ge_gas, Bo) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rs -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @yac/bl @std) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Glaso(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rs::Union{Int64, Float64}, API::Union{Int64, Float64}, z::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Bo::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)
    F = Rs*(Ge_gas/Ge_oil)^0.526 + 0.968*Ty
    Psc = 756.8 - 131*Ge_gas - 3.6*Ge_gas^2
    Tsc = 169.2 + 349.5*Ge_gas - 74*Ge_gas^2
    Psr = Py/Psc
    Tsr = (Ty + 460)/Tsc

    Bg = 50.3e-4*(z*(Ty + 460)/Py)
    dBo = (2.91329 - 0.55366*log10(F))*10^(-6.58511 + 2.9132*log10(F) - 0.27683*log10(F)^2)*(1/F)*(Ge_gas/Ge_oil)^0.526
    dRs = 2.02777*(Rs/(Py*(14.1811 - 3.3093*log10(Py))^0.5))
    co = (-1/Bo)*(dRs)*(dBo - Bg)

    return round(co, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de la Total, C.F.P.
""" Correlación de la Total, C.F.P. - 1983 \n
co\\_Total(Py, Ty, Rs, API, z, Ge_gas, Bo) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rs -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @yac/bl @std) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Total(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rs::Union{Int64, Float64}, API::Union{Int64, Float64}, z::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Bo::Union{Int64, Float64})::Float64

    if API <= 10
        c1, c2, c3, c4 = 12.2651, 0.030405, 0, 0.9669
    elseif 10 < API <= 35
        c1, c2, c3, c4 = 15.0057, 0.0152, 4.484e-4, 1.095
    elseif 35 < API <= 45
        c1, c2, c3, c4 = 112.925, 0.0248, -1.469e-3, 1.129
    end

    Psc = 756.8 - 131*Ge_gas - 3.6*Ge_gas^2
    Tsc = 169.2 + 349.5*Ge_gas - 74*Ge_gas^2
    Psr = Py/Psc
    Tsr = (Ty + 460)/Tsc

    Bg = 50.3e-4*(z*(Ty + 460)/Py)
    dBo = 4.857e-4 + (17.569e-9)*(Ty - 60)*(API/Ge_gas) 
    dRs = c4*(Rs/Py)
    co = (-1/Bo)*(dRs)*(dBo - Bg)

    return round(co, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de Al-Marhoun, M.A.
""" Correlación de Al-Marhoun, M.A. - 1988 \n
co\\_Almarhoun(Py, Ty, Rs, API, z, Ge_gas, Bo) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rs -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @yac/bl @std) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_AlMarhoun(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rs::Union{Int64, Float64}, API::Union{Int64, Float64}, z::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Bo::Union{Int64, Float64})::Float64
    
    Ge_oil = 141.5/(API + 131.5)
    F = (Rs^0.74239)*(Ge_gas)^0.323294*(Ge_oil)^-1.20204

    Psc = 756.8 - 131*Ge_gas - 3.6*Ge_gas^2
    Tsc = 169.2 + 349.5*Ge_gas - 74*Ge_gas^2
    Psr = Py/Psc
    Tsr = (Ty + 460)/Tsc

    Bg = 50.3e-4*(z*(Ty + 460)/Py)
    dBo = (1.35556e-3*F + 4.723e-6*F^2)/Rs
    dRs = 1.3984*Rs/Py
    co = (-1/Bo)*(dRs)*(dBo - Bg)

    return round(co, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de Dokla, M.E. y Osman, M.E.
""" Correlación de Dokla, M.E. y Osman, M.E. - 1992 \n
co\\_Dokla\\_Osman(Py, Ty, Rs, API, z, Ge_gas, Bo) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rs -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @yac/bl @std) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Dokla_Osman(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rs::Union{Int64, Float64}, API::Union{Int64, Float64}, z::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Bo::Union{Int64, Float64})::Float64
    
    Ge_oil = 141.5/(API + 131.5)
    F = (Rs^0.773572)*(Ge_gas^0.40402)*(Ge_oil)^-0.882605
    Psc = 756.8 - 131*Ge_gas - 3.6*Ge_gas^2
    Tsc = 169.2 + 349.5*Ge_gas - 74*Ge_gas^2
    Psr = Py/Psc
    Tsr = (Ty + 460)/Tsc

    Bg = 50.3e-4*(z*(Ty + 460)/Py)
    dBo = (1.08126e-3*F + 5.887e-6*F^2)/Rs
    dRs = 1.3811*Rs/Py
    co = (-1/Bo)*(dRs)*(dBo - Bg)

    return round(co, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de Petrosky, G.E. Jr y Farshad, F.F.
""" Correlación de Petrosky, G.E. Jr y Farshad, F.F. - 1998 \n
co\\_Petrosky\\_Farshad\\_sat(Py, Ty, Rs, API, z, Ge_gas, Bo) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rs -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @yac/bl @std) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Petrosky_Farshad_sat(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rs::Union{Int64, Float64}, API::Union{Int64, Float64}, z::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Bo::Union{Int64, Float64})::Float64
    
    Ge_oil = 141.5/(API + 131.5)
    F = (Rs^0.3738)*(Ge_gas^0.2914/Ge_oil^0.6265) + 0.24626*Ty^0.5371
    Psc = 756.8 - 131*Ge_gas - 3.6*Ge_gas^2
    Tsc = 169.2 + 349.5*Ge_gas - 74*Ge_gas^2
    Psr = Py/Psc
    Tsr = (Ty + 460)/Tsc

    Bg = 50.3e-4*(z*(Ty + 460)/Py)
    dBo = 83.313e-6*(F^2.0936)*(Rs^-0.6262)*(Ge_gas^0.2914/Ge_oil^0.6265)
    dRs = 1.73184*Rs/(Py + 139.05)
    co = (-1/Bo)*(dRs)*(dBo - Bg)

    return round(co, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T. y Schmidt, Z.
""" Correlación de Kartoatmodjo, T. y Schmidt, Z. - 1994 \n
co\\_Kartoatmodjo\\_Schmidt\\_sat(Py, Ty, Rs, API, z, Ge_gas, Bo) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rs -> Relación de solubilidad a la presión de burbuja - (pie^3/bls) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @yac/bl @std) {Int64, Float64} \n
Co -> Compresibilidad del aceite - (psia^-1) {Float64} \n
"""
function co_Kartoatmodjo_Schmidt_sat(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rs::Union{Int64, Float64}, API::Union{Int64, Float64}, z::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Bo::Union{Int64, Float64})::Float64
    
    Ge_oil = 141.5/(API + 131.5)

    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end
    
    F = F = (Rs^0.755)*(Ge_gas^0.25)*(Ge_oil^-1.5) + 0.45*Ty
    Psc = 756.8 - 131*Ge_gas - 3.6*Ge_gas^2
    Tsc = 169.2 + 349.5*Ge_gas - 74*Ge_gas^2
    Psr = Py/Psc
    Tsr = (Ty + 460)/Tsc

    Bg = 50.3e-4*(z*(Ty + 460)/Py)
    dBo = 1.1325e-4*(F^1.5 - 0.45*Ty*F^0.5)/Rs
    dRs = Rs/(c4*Py)
    co = (-1/Bo)*(dRs)*(dBo - Bg)

    return round(co, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------