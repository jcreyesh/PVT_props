# ---- Relación de solubilidad del agua de formación - Rsw ----
# Correlación de Culberson, O.L. y Mcketta, J.J.Jr.
""" Correlación de Culberson, O.L. y Mcketta, J.J.Jr. - 1951 \n
rsw\\_Culberson\\_Mcketta(Py, Ty, s) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
Rsw -> Relación de solubilidad del agua - (pie^3/bls) {Float64}\n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function rsw_Culberson_Mcketta(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, s::Union{Int64, Float64})::Float64
    s = s/10
    A = 8.15839 - 6.12265e-2*Ty + 1.91663e-4*Ty^2 - 2.1654e-7*Ty^3
    B = 1.01021e-2 - 7.44241e-5*Ty + 3.05553e-7*Ty^2 - 2.94883e-10*Ty^3
    C = (-9.02505 + 0.130237*Ty - 8.53425e-4*Ty^2 + 2.34122e-6*Ty^3 - 2.37049e-9*Ty^4)*10^-7
    Rswp = A + B*Py + C*Py^2
    Rsw = 10^(-0.0840655*s*Ty^-0.285854)*Rswp

    return round(Rsw, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de McCoy, R.L.
""" Correlación de McCoy, R.L. - 1983 \n
rsw\\_McCoy(Py, Ty, s) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
Rsw -> Relación de solubilidad del agua - (pie^3/bls) {Float64}\n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function rsw_McCoy(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, s::Union{Int64, Float64})::Float64
    s = s/10
    A = 2.12 + 3.45e-3*Ty - 3.59e-5*Ty^2
    B = 0.0107 - 5.26e-5*Ty + 1.48e-7*Ty^2
    C = -8.75e-7 + 3.9e-9*Ty - 1.02e-11*Ty^2
    Rswp = A + B*Py + C*Py^2
    Rsw = (1-(0.0753 - 1.73e-4*Ty)*s)*Rswp

    return round(Rsw, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de McCain, W.D., Jr.
""" Correlación de McCain, W.D., Jr. - 1990 \n
bw\\_McCain(Py, Ty, s) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Bw -> Factor de volumen del agua - (bls @yac/bls @std) {Float64}\n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function bw_McCain(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64})::Float64
    # s = s/10
    dvwt = -1.0001e-2 + 1.33391e-4*Ty + 5.50654e-7*Ty^2
    dvwp = -1.95301e-9*Py*Ty - 1.72834e-13*Py^2*Ty - 3.58922e-7*Py - 2.25341e-10*Py^2
    bw = (1 + dvwp)*(1 + dvwt)

    return round(bw, digits=4)
end

#----------------------------------------------------------------------------------------------------------------------------
# Correlación de McCoy, R.L.
""" Correlación de McCoy, R.L. - 1983 \n
bw\\_McCoy(Py, Ty, s) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Pb -> Presión de burbija - (psia) {Int64, Float64} \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
Bw -> Factor de volumen del agua - (bls @yac/bls @std) {Float64}\n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function bw_McCoy(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Pb::Union{Int64, Float64}, s::Union{Int64, Float64})::Float64
    s = s/10

    if Py > Pb
        A = 0.9911 + 6.35e-5*Ty + 8.5e-7*Ty^2
        B = -1.093e-6 - 3.497e-9*Ty + 4.57e-12*Ty^2
        C = -5e-11 + 6.429e-13*Ty - 1.43e-15*Ty^2
    else
        A = 0.9947 + 5.8e-6*Ty + 1.02e-6*Ty^2
        B = -4.228e-6 + 1.8376e-8*Ty - 6.77e-11*Ty^2
        C = 1.3e-10 - 1.3855e-12*Ty + 4.285e-15*Ty^2
    end

    bwp = A + B*Py + C*Py^2
    bw = (1 + s*(5.1e-8*Py + (5.47e-6 - 1.95e-10*Py)*(Ty - 60) - (3.23e-8 - 8.5e-13*Py)*(Ty - 60)^2))*bwp

    return round(bw, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de McCoy, R.L.
""" Correlación de Dodson, C.R. y Standing, M.B. - 1944 \n
cwb\\_Dodson\\_Standing(Py, Ty, Rsw, s) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Rsw -> Relación de solubilidad del agua - (pie^3/bls) \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
cwb -> Compresibilidad del agua bajosaturada (bls @yac/bls @std) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function cwb_Dodson_Standing(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, Rsw::Union{Int64, Float64}, s::Union{Int64, Float64})
    s = s/10
    A = 3.8546 - 1.34e-4*Py
    B = -0.01052 + 4.77e-7*Py
    C = 3.9267e-5 - 8.8e-10*Py
    cwp = (A + B*Ty + C*Ty^2)/10^6
    cwb = (1 + 8.9e-3*Rsw)*cwp
    fc_sal = 1 + s^0.7*(-5.2e-2 + 2.7e-4*Ty - 1.14e-6*Ty^2 + 1.121e-9*Ty^3)
    cwb = fc_sal*3.43e-6

    return round(cwb, digits=10)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Osif, T.L.
""" Correlación de Osif, T.L. - 1988 \n
cwb\\_Osif(Py, Ty, s) \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
cwb -> Compresibilidad del agua bajosaturada (bls @yac/bls @std) {Float64}  \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                      \n
1000 < Presión (psia) < 20,000                                    \n
200 < Temperatura (°F) < 270                                      \n
0 < Salinidad del agua (g/lt) < 200                               \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function cwb_Osif(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, s::Union{Int64, Float64})
    
    s = s*1000
    cwb = 1/(7.033*Py + 541.5*s/58443 - 537*Ty + 403300)

    return round(cwb, digits=10)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Ramey, H.J., Jr
""" Correlación de Ramey. H.J., Jr. - 1964 \n
cws\\_Ramey(Py, Ty, cwb, Bg, Bw, s) \n
\\--------------------------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
cwb -> Compresibilidad del agua bajosaturada - (psia^-1) {Int64, Float64} \n
Bg -> Factor de volumen del gas - (pie^3 @yac/pie^3 @std) {Int64, Float64} \n
Bw -> Factor de volumen del agua - (bls^3 @yac/bls^3 @std) {Int64, Float64} \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
cws -> Compresibilidad del agua saturada - (psi-1) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------------------------------------------\n
"""
function cws_Ramey(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, cwb::Union{Int64, Float64}, Bg::Union{Int64, Float64}, Bw::Union{Int64, Float64}, s::Union{Int64, Float64})
    
    s = s/10
    A = 8.15839 - 6.12265e-2*Ty + 1.91663e-4*Ty^2 - 2.1654e-7*Ty^3
    B = 1.01021e-2 - 7.44241e-5*Ty + 3.05553e-7*Ty^2 - 2.94883e-10*Ty^3
    C = (-9.02505 + 0.130237*Ty - 8.53425e-4*Ty^2 + 2.34122e-6*Ty^3 - 2.37049e-9*Ty^4)*10^-7

    fc_sal = 1 + s^0.7*(-5.2e-2 + 2.7e-4*Ty - 1.14e-6*Ty^2 + 1.121e-9*Ty^3)
    dRsw = B + 2*C*Py
    dRsw_c = dRsw * fc_sal

    cws = cwb + Bg*dRsw_c/Bw

    return round(cws, digits=10)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Van Wingen, N.
""" Correlación de Van Wingen - 1950 \n
uw\\_Van_Wingen(Ty) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
uw -> Viscosidad del agua saturada - (cp) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function uw_Van_Wingen(Ty::Union{Int64, Float64})::Float64

    uw = exp(1.003 - 1.479e-2*Ty + 1.982e-5*Ty^2)

    return round(uw, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Matthews, C.S. y Russel, D.G.
""" Correlación de Matthews, C.S. y Russell, D.G. - 1967 \n
uw\\_Matthews_Russell(Py, Ty, s) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
uw -> Viscosidad del agua saturada - (cp) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function uw_Matthews_Russell(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, s::Union{Int64, Float64})::Float64

    s = s/10
    A = -0.04518 + 93.13e-4*s - 3.93e-4*s^2
    B = 70.634 + 95.76e-3*s^2
    uw1 = A + B/Ty
    fc = 1 + 3.5e-12*Py^2*(Ty - 40)
    uw = uw1*fc

    return round(uw, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de McCain, W.D., Jr.
""" Correlación de McCain, W.D., Jr. - 1990 \n
uw\\_McCain(Py, Ty, s) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
uw -> Viscosidad del agua saturada - (cp) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function uw_McCain(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, s::Union{Int64, Float64})::Float64

    s = s/10
    A = 109.574 - 8.40564*s + 0.313314*s^2 + 8.72213e-3*s^3
    B = -1.12166 + 2.63951e-2*s - 6.79461e-4*s^2 - 5.47119e-5*s^3 + 1.55586e-6*s^4
    uw1 = A*Ty^B
    uw = (99.94e-2 + 4.0295e-5*Py + 3.1062e-9*Py^2)*uw1

    return round(uw, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de McCoy, R.L.
""" Correlación de McCoy - 1983 \n
uw\\_McCoy(Py, Ty, s) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
uw -> Viscosidad del agua saturada - (cp) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function uw_McCoy(Ty::Union{Int64, Float64}, s::Union{Int64, Float64})::Float64
    
    Ty = 5/9*Ty + 255.37
    s = s/10
    uwp = 24.14e-3*10^(247.8/(Ty - 140))
    uw = (1 - 1.87e-3*s^0.5 + 2.18e-4*s^2.5 + (Ty^0.5 - 1.35e-2*Ty)*(2.76e-3*s - 3.44e-4*s^1.5))*uwp

    return round(uw, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de McCoy, R.L.
""" Correlación de McCain, W.D., Jr. - 1990 \n
dw\\_McCain(Bw, s) \n
\\---------------------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Bw -> Factor de volumen del agua - (bls^3 @yac/bls^3 @std) {Int64, Float64} \n
s -> Salinidad del agua - (Mppm) {Int64, Float64} \n
dw -> Densidad del agua - (lbs/pie^3) {Float64} \n
\\---------------------------------------------------------------------------------------------------------------------------------------------\n
"""
function dw_McCain(Bw::Union{Int64, Float64}, s::Union{Int64, Float64})::Float64
    
    s = s/10
    dw1 = 62.368 + 0.438603*s + 1.60074e-3*s^2
    dw = dw1/Bw

    return round(dw, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Jennings, H.Y.Jr. y Newman, G.H.
""" Correlación de Jennings, H.Y.Jr. y Newman, G.H. - 1971 \n
tigw\\_Jennings_Newman(Py, Ty) \n
\\--------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py -> Presión yacimiento - (psia) {Int64, Float64} \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
tigw -> Tensión interfacial gas-agua - (dinas/cm) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------------------\n
"""
function tigw_Jennings_Newman(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64})::Float64
    
    A = 79.1618 - 0.118978*Ty
    B = -5.28473e-3 + 9.87913e-6*Ty
    C = (2.33814 - 4.57194e-4*Ty - 7.52678e-6*Ty^2)*10^-7
    tigw = A + B*Py + C*Py^2

    return round(tigw, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------