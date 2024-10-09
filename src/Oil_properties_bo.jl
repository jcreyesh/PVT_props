# ---- Factor de volumen del aceite - Bo ----
# Correlación de Standing, M.B.
""" Correlación de Standing, M.B. para Bo - 1977 \n
bo\\_Standing(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
\t Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_Standing(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)

    if Py <= Pb
        rs = Ge_gas * ((Py/18.2 + 1.4) * 10^(12.5e-3*API - 91e-5*Ty))^1.2048
        F = rs * sqrt(Ge_gas/Ge_oil) + 1.25*Ty
        bo = 0.9759 + 12e-5*F^1.2
    
    elseif Py > Pb
        rsb = Ge_gas * ((Pb/18.2 + 1.4) * 10^(12.5e-3*API - 91e-5*Ty))^1.2048
        Fb = rsb * sqrt(Ge_gas/Ge_oil) + 1.25*Ty
        bob = 0.9759 + 12e-5*Fb^1.2
        bo = bob*exp(Co*(Pb - Py))
    
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Vasquez, M.E. y Beggs, H.D.
""" Correlación de Vásquez, M.E. y Beggs, H.D. para Bo - 1980 \n
bo\\_Vasquez\\_Beggs(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t Psp -> Presión en el separador - (psia) \n
\t Tsp -> Temperatura en el separador - (°F) \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
\t Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_Vasquez_Beggs_sp(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)

    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else 
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end

    Gec_gas = Ge_gas * (1 + 5.912e-5*API*Tsp*log10(Psp/114.7))

    if Py < Pb
        rs = c1*Gec_gas*Py^c2*exp(c3*API/(Ty + 460))
    elseif Py >= Pb
        rs = c1*Gec_gas*Pb^c2*exp(c3*API/(Ty + 460))
    end
    
    if API <= 30
        c1, c2, c3 = 4.677e-4, 1.751e-5, -1.8106e-8
    else 
        c1, c2, c3 = 4.67e-4, 1.1e-5, 1.3370e-9
    end

    if Py <= Pb
        bo = 1 + c1*rs + c2*(Ty - 60)*(API/Gec_gas) + c3*rs*(Ty - 60)*(API/Gec_gas)
    elseif Py > Pb
        bob = 1 + c1*rs + c2*(Ty - 60)*(API/Gec_gas) + c3*rs*(Ty - 60)*(API/Gec_gas)
        bo = bob*exp(Co*(Pb - Py))
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Vasquez, M.E. y Beggs, H.D.
""" Correlación de Vásquez, M.E. y Beggs, H.D. para Bo - 1980 \n
bo\\_Vasquez\\_Beggs(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
\t Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_Vasquez_Beggs(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})

    Ge_oil = 141.5/(API + 131.5)

    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else 
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end

    if Py < Pb
        rs = c1*Ge_gas*Py^c2*exp(c3*API/(Ty + 460))
    elseif Py >= Pb
        rs = c1*Ge_gas*Pb^c2*exp(c3*API/(Ty + 460))
    end
    
    if API <= 30
        c1, c2, c3 = 4.677e-4, 1.751e-5, -1.8106e-8
    else 
        c1, c2, c3 = 4.67e-4, 1.1e-5, 1.3370e-9
    end

    if Py <= Pb
        bo = 1 + c1*rs + c2*(Ty - 60)*(API/Ge_gas) + c3*rs*(Ty - 60)*(API/Ge_gas)
    elseif Py > Pb
        bob = 1 + c1*rs + c2*(Ty - 60)*(API/Ge_gas) + c3*rs*(Ty - 60)*(API/Ge_gas)
        bo = bob*exp(Co*(Pb - Py))
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Glaso, O.
""" Correlación de Glaso, O. para Bo - 1980 \n
bo\\_Glaso(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
\t Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_Glaso(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)

    if Py <= Pb
        E = 10^(2.8869 - (14.1811 - 3.3093*log10(Py))^0.5)
        rs = Ge_gas*(E*API^0.989/Ty^0.172)^1.2255
        F = rs*(Ge_gas/Ge_oil)^0.526 + 0.986*Ty
        bo = 1 + 10^(-6.58511 + 2.91329*log10(F) - 0.27683*log10(F)^2)
    
    elseif Py > Pb
        Eb = 10^(2.8869 - (14.1811 - 3.3093*log10(Pb))^0.5)
        rsb = Ge_gas*(Eb*API^0.989/Ty^0.172)^1.2255
        Fb = rsb*(Ge_gas/Ge_oil)^0.526 + 0.986*Ty
        bob = 1 + 10^(-6.58511 + 2.91329*log10(Fb) - 0.27683*log10(Fb)^2)
        bo = bob*exp(Co*(Pb - Py))
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de TOTAL, C.F.P.
""" Correlación de la Total, C.F.P. para Bo - 1983 \n
bo\\_Total\\_CFP(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n|
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
\t Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_Total_CFP(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)
    
    if API <= 10
        c1, c2, c3, c4 = 12.2651, 0.030405, 0, 0.9669
    elseif 10 < API <= 35
        c1, c2, c3, c4 = 15.0057, 0.0152, 4.484e-4, 1.0950
    elseif 35 < API <= 45
        c1, c2, c3, c4 = 112.925, 0.0248, -1.469e-3, 1.129
    end

    if Py <= Pb
        rs = Ge_gas*(Py*10^(c2*API - c3*Ty)/c1)^c4
        bo = 1.022 + 4.857e-4*rs - 2.009e-6*(Ty - 60)*(API/Ge_gas) + 17.569e-9*rs*(Ty - 60)*(API/Ge_gas)
    
    elseif Py > Pb
        rsb = Ge_gas*(Pb*10^(c2*API - c3*Ty)/c1)^c4
        bob = 1.022 + 4.857e-4*rsb - 2.009e-6*(Ty - 60)*(API/Ge_gas) + 17.569e-9*rsb*(Ty - 60)*(API/Ge_gas)
        bo = bob*exp(Co*(Pb - Py))
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Al-Marhoun, M.A.
""" Correlación de Al-Marhoun, M.A. para Bo - 1988 \n
bo\\_AlMarhoun(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
\t Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_AlMarhoun(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)
    
    if Py <= Pb
        rs = (185.84321*Py*(Ge_gas^1.87784)*(Ge_oil^-3.1437)*((Ty + 460)^-1.32657))^1.3984
        F = (rs^0.74239)*(Ge_gas)^0.323294*(Ge_oil)^-1.20204
        bo = 0.497069 + 0.862963e-3*(Ty + 460) + 0.182594e-2*F + 0.318099e-5*F^2
    
    elseif Py > Pb
        rsb = (185.84321*Pb*(Ge_gas^1.87784)*(Ge_oil^-3.1437)*((Ty + 460)^-1.32657))^1.3984
        Fb = (rsb^0.74239)*(Ge_gas)^0.323294*(Ge_oil)^-1.20204
        bob = 0.497069 + 0.862963e-3*(Ty + 460) + 0.182594e-2*Fb + 0.318099e-5*Fb^2
        bo = bob*exp(Co*(Pb - Py))
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Dokla, M.E. y Osman, M.E.
""" Correlación de Dokla, M.E. y Osman, M.E. para Bo - 1992 \n
bo\\_Dokla\\_Osman(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
\t Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_Dokla_Osman(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)
    
    if Py <= Pb
        rs = (0.11956e-3*Py*(Ge_gas^1.01049)*(Ge_oil^-0.107991)*((Ty + 460)^0.952584))^1.3811
        F = (rs^0.773572)*(Ge_gas^0.40402)*(Ge_oil)^-0.882605
        bo = 0.431935e-1 + 0.156667e-2*(Ty + 460) + 0.139775e-2*F + 0.380525e-5*F^2
    
    elseif Py > Pb
        rsb = (0.11956e-3*Pb*(Ge_gas^1.01049)*(Ge_oil^-0.107991)*((Ty + 460)^0.952584))^1.3811
        Fb = (rsb^0.773572)*(Ge_gas^0.40402)*(Ge_oil)^-0.882605
        bob = 0.431935e-1 + 0.156667e-2*(Ty + 460) + 0.139775e-2*Fb + 0.380525e-5*Fb^2
        bo = bob*exp(Co*(Pb - Py))
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
""" Correlación de Petrosky, G.E. y Farshad, F.F. - 1993 \n
bo\\_Petrosky\\_Farshad(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
\t Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_Petrosky_Farshad(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)
    
    if Py <= Pb
        rs = ((Ge_gas^0.8439) * (Py/112.727 + 12.34) * 10^(7.916e-4*API^1.541 - 4.561e-5*Ty^1.3911))^1.73184
        F = (rs^0.3738)*(Ge_gas^0.2914/Ge_oil^0.6265) + 0.24626*Ty^0.5371
        bo = 1.0113 + 7.2046e-5*F^3.0936
    
    elseif Py > Pb
        rsb = ((Ge_gas^0.8439) * (Pb/112.727 + 12.34) * 10^(7.916e-4*API^1.541 - 4.561e-5*Ty^1.3911))^1.73184
        Fb = (rsb^0.3738)*(Ge_gas^0.2914/Ge_oil^0.6265) + 0.24626*Ty^0.5371
        bob = 1.0113 + 7.2046e-5*Fb^3.0936
        bo = bob*exp(Co*(Pb - Py))
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
""" Correlación de Kartoatmodjo, T. y Schmidt, Z. - 1994 \n
bo\\_Kartoatmodjo\\_Schmidt\\_sp(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
\t Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_Kartoatmodjo_Schmidt_sp(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)

    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end
    
    Gec_gas = Ge_gas * (1 + 0.1595*API^0.4078*Tsp^-0.2466*log10(Psp/114.7))

    if Py <= Pb
        rs = c1 * (Gec_gas^c2) * Py^(1/c4) *10^(c3*API/(Ty + 460))
        F = (rs^0.755)*(Gec_gas^0.25)*(Ge_oil^-1.5) + 0.45*Ty
        bo = 0.98496 + 1e-4*F^1.5
    
    elseif Py > Pb
        rsb = c1 * (Gec_gas^c2) * Pb^(1/c4) *10^(c3*API/(Ty + 460))
        Fb = (rsb^0.755)*(Gec_gas^0.25)*(Ge_oil^-1.5) + 0.45*Ty
        bob = 0.98496 + 1e-4*Fb^1.5
        bo = bob*exp(Co*(Pb - Py))
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
""" Correlación de Kartoatmodjo, T. y Schmidt, Z. - 1994 \n
bo\\_Kartoatmodjo\\_Schmidt(Py, Ty, Pb, API, Ge_gas, Co) \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
API -> Gravedad API - (grados) {Int64, Float64} \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
Co -> Compresibilidad del aceite- (psi^-1) {Int64, Float64} \n
Bo -> Factor de volumen del aceite - (bl @std/bl @yac) {Float64} \n
"""
function bo_Kartoatmodjo_Schmidt(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Co::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)
    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end
    
    if Py <= Pb
        rs = c1 * (Ge_gas^c2) * Py^(1/c4) *10^(c3*API/(Ty + 460))
        F = (rs^0.755)*(Ge_gas^0.25)*(Ge_oil^-1.5) + 0.45*Ty
        bo = 0.98496 + 1e-4*F^1.5
    
    elseif Py > Pb
        rsb = c1 * (Ge_gas^c2) * Pb^(1/c4) *10^(c3*API/(Ty + 460))
        Fb = (rsb^0.755)*(Ge_gas^0.25)*(Ge_oil^-1.5) + 0.45*Ty
        bob = 0.98496 + 1e-4*Fb^1.5
        bo = bob*exp(Co*(Pb - Py))
    end

    return round(bo, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------