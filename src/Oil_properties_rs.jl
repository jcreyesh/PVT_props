# ---- Relación de Solubilidad - Rs ----
# Correlación de Standing, M.B.
"""
Correlación de Standing, M.B. - 1977   \n
rs\\_Standing(Py, Ty, Pb, API, Ge_gas) \n
\t Py -> Presión de yacimiento - (psia) \n
\t Ty -> Temperatura del yacimiento - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl)
"""
function rs_Standing(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})
    
    if Py >= Pb
        rs = Ge_gas * ((Pb/18.2 + 1.4) * 10^(12.5e-3*API - 91e-5*Ty))^1.2048
    else
        rs = Ge_gas * ((Py/18.2 + 1.4) * 10^(12.5e-3*API - 91e-5*Ty))^1.2048
    end
    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Lasater, J.A.
"""
Correlación de Lasater, J.A. - 1958   \n
rs\\_Lasater(Py, Ty, Pb, API, Ge_gas) \n
\t Py -> Presión de yacimiento - (psia) \n
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_Lasater(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    Ge_oil = 141.5/(API + 131.5)

    if API <= 40
        Mo = 630 - 10*API
    else
        Mo = 73110*API^-1.562
    end

    if Py >= Pb
        pf = Pb*Ge_gas/(Ty + 460)
        Py = Pb
    else 
        pf = Py*Ge_gas/(Ty + 460)
    end

    if pf < 3.29
        yg = 0.359*log(1.476*Py*Ge_gas/(Ty + 460) + 0.476)
    else
        yg = (0.121*Py*Ge_gas/(Ty + 460) - 0.236)^0.281
    end

    rs = 132755*Ge_oil*yg/(Mo*(1 - yg))

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Vásquez, M.E. y Beggs, H.D.
"""
Correlación de Vásquez, M.E. y Beggs, H.D. - 1980 \n
rs\\_Vasquez\\_Beggs\\_sp(Py, Ty, Pb, API, Ge\\_gas) \n
\t Py -> Presión de yacimiento - (psia)
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Psp -> Presión en el separador - (psia) \n
\t Tsp -> Temperatura en el separador - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_Vasquez_Beggs_sp(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    
    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else 
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end    

    Gec_gas = Ge_gas * (1 + 5.912e-5*API*Tsp*log10(Psp/114.7))

    if Py >= Pb
        Py = Pb
    end

    rs = c1*Gec_gas*Py^c2*exp(c3*API/(Ty + 460))

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Vásquez, M.E. y Beggs, H.D.
"""
Correlación de Vásquez, M.E. y Beggs, H.D. - 1980  \n
rs\\_Vasquez\\_Beggs(Py, Ty, Pb, API, Ge\\_gas) \n
\t Py -> Presión de yacimiento - (psia)
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_Vasquez_Beggs(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    
    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else 
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end    

    if Py >= Pb
        Py = Pb
    end

    rs = c1*Ge_gas*Py^c2*exp(c3*API/(Ty + 460))

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Glaso, O.
"""
Correlación de Glaso, O. - 1980     \n
rs\\_Glaso(Py, Ty, Pb, API, Ge\\_gas) \n
\t Py -> Presión de yacimiento - (psia)
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_Glaso(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    
    if Py >= Pb
        Py = Pb
    end

    F = 10^(2.8869 - (14.1811 - 3.3093*log10(Py))^0.5)
    rs = Ge_gas*(F*API^0.989/Ty^0.172)^1.2255

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de TOTAL, C.F.P.
"""
Correlación de TOTAL, C.F.P. - 1983        \n
rs\\_Total\\_CFP(Py, Ty, Pb, API, Ge\\_gas) \n
\t Py -> Presión de yacimiento - (psia)
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_Total_CFP(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    
    if API <= 10
        c1, c2, c3, c4 = 12.2651, 0.030405, 0, 0.9669
    elseif 10 < API <= 35
        c1, c2, c3, c4 = 15.0057, 0.0152, 4.484e-4, 1.0950
    elseif 35 < API <= 45
        c1, c2, c3, c4 = 112.925, 0.0248, -1.469e-3, 1.129
    end

    if Py >= Pb
        Py = Pb
    end

    rs = Ge_gas*(Py*10^(c2*API - c3*Ty)/c1)^c4

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Al-Marhoun, M.A.
"""
Correlación de Al-Marhoun, M.A. - 1988   \n
rs\\_Al-Marhoun(Py, Ty, Pb, API, Ge\\_gas) \n
\t Py -> Presión de yacimiento - (psia)
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_AlMarhoun(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    
    Ge_oil = 141.5/(API + 131.5)

    if Py >= Pb
        Py = Pb
    end

    rs = (185.84321*Py*(Ge_gas^1.87784)*(Ge_oil^-3.1437)*((Ty + 460)^-1.32657))^1.3984

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Dokla, M.E. y Osman,M.E.
"""
Correlación de Dokla, M.E. y Osman,M.E. - 1992    \n
rs\\_Dokla\\_Osman(Py, Ty, Pb, API, Ge\\_gas) \n
\t Py -> Presión de yacimiento - (psia)
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_Dokla_Osman(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    
    Ge_oil = 141.5/(API + 131.5)

    if Py >= Pb
        Py = Pb
    end

    rs = (0.11956e-3*Py*(Ge_gas^1.01049)*(Ge_oil^-0.107991)*((Ty + 460)^0.952584))^1.3811

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Petrosky, G.E.Jr. y Farshad, F.F.
"""
Correlación de Petrosky, G.E. y Farshad, F.F. - 1993   \n
rs\\_Petrosky\\_Farshad(Py, Ty, Pb, API, Ge\\_gas) \n
\t Py -> Presión de yacimiento - (psia)
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_Petrosky_Farshad(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    
    if Py >= Pb
        Py = Pb
    end

    rs = ((Ge_gas^0.8439) * (Py/112.727 + 12.34) * 10^(7.916e-4*API^1.541 - 4.561e-5*Ty^1.3911))^1.73184

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T. y Schmidt, Z.
"""
Correlación de Kartoatmodjo, T. y Schmidt, Z. - 1994    \n
rs\\_Kartoatmodjo\\_Schmidt\\_sp(Py, Ty, Pb, API, Ge\\_gas) \n
\t Py -> Presión de yacimiento - (psia) \n
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Psp -> Presión en el separador - (psia) \n
\t Tsp -> Temperatura en el separador - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
\t rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_Kartoatmodjo_Schmidt_sp(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    
    Gec_gas = Ge_gas * (1 + 0.1595*API^0.4078*Tsp^-0.2466*log10(Psp/114.7))

    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end

    if Py >= Pb
        Py = Pb
    end

    rs = c1 * (Gec_gas^c2) * Py^(1/c4) *10^(c3*API/(Ty + 460))

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T. y Schmidt, Z.
"""
Correlación de Kartoatmodjo, T. y Schmidt, Z. - 1994    \n
rs\\_Kartoatmodjo\\_Schmidt(Py, Ty, Pb, API, Ge\\_gas)    \n
\t Py -> Presión de yacimiento - (psia) \n
\t Ty -> Temperatura de yacimiento - (°F) \n
\t Pb -> Presión de burbuja - (psia) \n
\t API -> Gravedad API - (grados) \n
\t Ge_gas -> Gravedad específica del gas - (adim) \n
rs -> Relación de solubilidad - (pie3/bl) \n
"""
function rs_Kartoatmodjo_Schmidt(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64
    
    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end

    if Py >= Pb
        Py = Pb
    end

    rs = c1 * (Ge_gas^c2) * Py^(1/c4) *10^(c3*API/(Ty + 460))

    return round(rs, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------