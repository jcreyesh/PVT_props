# ---- Viscosidad del aceite - uo ----
# ---- Aceite bajosaturado Py > Pb ----
# Correlación de Beal, C.
""" Correlación de Beal, C. - 1946 \n
uob\\_Beal(Ty, API) \n
\\- Argumentos - \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t uob -> Viscosidad del aceite saturado - (cp) {Float64}
\n
\\- Rangos - \n
\t 98 < Temperatura (°F) < 250              \n
\t 10 < Gravedad del aceite (°API) < 52.5   \n
\t 0.865 < Viscosidad (cp) < 1.55           \n
"""
function uob_Beal(Ty::Union{Int64, Float64}, API::Union{Int64, Float64})::Float64
    
    a = exp10(0.43 + 8.33/API)
    uob = (0.32 + 1.8e7/API^4.53)*(360/(Ty + 200))^a

    return round(uob, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Beggs, H.D. y Robinson, J.R.
""" Correlación de Beggs, H.D. y Robinson, J.R. - 1975 \n
uob\\_Beggs\\_Robinson(Ty, API) \n
\\- Argumentos - \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t uob -> Viscosidad del aceite saturado - (cp) {Float64}
\\- Rangos -                           \n
\t 15 < Presión (psia) < 5265             \n
\t 70 < Temperatura (°F) < 295            \n
\t 20 < Rs (pie^3/bls) < 2070             \n
\t 16 < Gravedad del aceite (°API) < 58   \n
"""
function uob_Beggs_Robinson(Ty::Union{Int64, Float64}, API::Union{Int64, Float64})::Float64
    
    z = 3.0324 - 0.02023*API
    y = 10^z
    x = y*Ty^-1.163
    uob = 10^x - 1

    return round(uob, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Glaso, O.
""" Correlación de Glaso, O. - 1980 \n
uob\\_Glaso(Ty, API) \n
\\- Argumentos - \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t uob -> Viscosidad del aceite saturado - (cp) {Float64}
\\- Rangos -                               \n
\t 50 < Temperatura (°F) < 300                \n
\t 20.1 < Gravedad del aceite (°API) < 48.1   \n
\t 0.616 < Viscosidad (cp) < 39.1             \n
"""
function uob_Glaso(Ty::Union{Int64, Float64}, API::Union{Int64, Float64})::Float64

    uob = 3.141e10*Ty^-3.444*(log10(API)^(10.313*log10(Ty) - 36.447))

    return round(uob, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Egbogad, E.O.
""" Correlación de Egbogad, O. - 1983                         \n
uob\\_Egbogad(Ty, API)                                        \n
\\- Argumentos -                                              \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}   \n
\t API -> Gravedad API - (grados) {Int64, Float64}            \n
\t uob -> Viscosidad del aceite saturado - (cp) {Float64}
\\- Rangos -                                                  \n
\t 59 < Temperatura (°F) < 176                                \n
\t 5 < Gravedad del aceite (°API) < 58                        \n
"""
function uob_Egbogad(Ty::Union{Int64, Float64}, API::Union{Int64, Float64})::Float64

    a = exp10(exp10(1.8653 - 0.025086*API - 0.5644*log10(Ty)))
    uob = a - 1 

    return round(uob, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T., Schmidt, Z.
""" Correlación de Kartoatmodjo, T., Schmidt, Z. - 1994 \n
uob\\_Kartoatmodjo\\_Schmidt(Ty, API) \n
\\- Argumentos - \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64}          \n
\t uob -> Viscosidad del aceite saturado - (cp) {Float64}   \n
\\- Rangos -                                                \n
\t 75 < Temperatura (°F) < 320                              \n
\t 14.4 < Gravedad del aceite (°API) < 58.9                 \n
\t 0.5 < Viscosidad (cp) < 682                              \n
"""
function uob_Kartoatmodjo_Schmidt(Ty::Union{Int64, Float64}, API::Union{Int64, Float64})::Float64

    uob = 16e8*Ty^-2.8177*log10(API)^(5.7526*log10(Ty) - 26.9718)

    return round(uob, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# ---- Aceite saturado Py = Pb ----
# Correlación de Chew, J.N. y Conally, C.A.Jr.
""" Correlación de Chew, J.N. y Conally, C.A.Jr. - 1959 \n
uos\\_Chew\\_Conally(Rs, uob) \n
\\- Argumentos - \n
\t Rs -> Relación de solubilidad @ Pb - (pie^3/bls) {Int64, Float64} \n
\t uob -> Viscosidad del aceite bajosaturado - (cp) {Float64} \n
\t uob -> Viscosidad del aceite saturado - (cp) {Float64}
\t uos -> Viscosidad del aceite saturado - (cp) {Float64} \n
\\- Rangos -                                  \n
\t 132 < Presión de burbuja (psia) < 5645     \n
\t 72 < Temperatura (°F) < 292                \n
\t 51 < Rs (pie^3/bls) < 3544                 \n
\t 0.377 < Viscosidad del aceite (cp) < 50    \n
"""
function uos_Chew_Conally(Rs::Union{Int64, Float64}, uob::Union{Int64, Float64})::Float64
    a = exp10(Rs*(2.2e-7*Rs - 7.4e-4))
    b = 0.68/(10^(8.62e-5*Rs)) + 0.25/(10^(1.1e-3*Rs)) + 0.062/(10^(3.74e-3*Rs))
    uos = a*uob^b

    return round(uos, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Beggs, H.D. y Robinson, J.R.
""" Correlación de Beggs, H.D. y Robinson, J.R. - 1975 \n
uos\\_Beggs\\_Robinson(Rs, uob) \n
\\- Argumentos - \n
\t Rs -> Relación de solubilidad @ Pb - (pie^3/bls) {Int64, Float64} \n
\t uob -> Viscosidad del aceite bajosaturado - (cp) {Float64} \n
\t uos -> Viscosidad del aceite saturado - (cp) {Float64} \n
\\- Rangos -                               \n
\t 15 < Presión (psia) < 5265                 \n
\t 70 < Temperatura (°F) < 295                \n
\t 20 < Rs (pie^3/bls) < 2070                 \n
\t 16 < Gravedad del aceite (°API) < 58       \n
"""
function uos_Beggs_Robinson(Rs::Union{Int64, Float64}, uob::Union{Int64, Float64})::Float64
    a = 10.715*(Rs + 100)^-0.515
    b = 5.44*(Rs + 150)^-0.338
    uos = a*uob^b

    return round(uos, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T., Schmidt, Z.
""" Correlación de Kartoatmodjo, T., Schmidt, Z. - 1994 \n
uos\\_Kartoadmodjo\\_Schmidt(Rs, uob) \n
\\- Argumentos - \n
\t Rs -> Relación de solubilidad @ Pb - (pie^3/bls) {Int64, Float64} \n
\t uob -> Viscosidad del aceite bajosaturado - (cp) {Float64} \n
\t uob -> Viscosidad del aceite saturado - (cp) {Float64} \n
\t uos -> Viscosidad del aceite saturado - (cp) {Float64} \n
"""
function uos_Kartoadmodjo_Schmidt(Rs::Union{Int64, Float64}, uob::Union{Int64, Float64})::Float64
    
    b = 10^(-81e-5*Rs)
    a = (0.2001 + 0.8428*10^(-0.000845*Rs))*(uob)^(0.43 + 0.5165*b)
    uos = -0.06821 + 0.9824*a + 40.34e-5*a^2

    return round(uos, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# ---- Aceite subsaturado Py < Pb ----
# Correlación de Beal, C.
""" Correlación de Beal, C. - 1946 \n
uosb\\_Beal(Py, Pb, uob) \n
\\- Argumentos - \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t uob -> Viscosidad del aceite saturado - (cp) {Float64} \n
\t uosb -> Viscosidad del aceite subsaturado - (cp) {Float64} \n
\\- Rangos -                               \n
\t 140 < Presión (psia) < 4135                \n
\t 12 < Rs (pie^3/bls) < 1827                 \n
\t 0.142 < Viscosidad del aceite (cp) < 127    \n
"""
function uosb_Beal(Py::Union{Int64, Float64}, Pb::Union{Int64, Float64}, uob::Union{Int64, Float64})::Float64
    
    uosb = (0.024*uob^1.6 + 0.038*uob^0.56)*1e-3*(Py - Pb) + uob

    return round(uosb, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Beal, C.
""" Correlación de Vásquez, M.E. y Beggs, H.D. - 1980 \n
uosb\\_Vasquez\\_Beggs(Py, Pb, uob) \n
\\- Argumentos - \n
Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
uob -> Viscosidad del aceite saturado - (cp) {Float64} \n
uosb -> Viscosidad del aceite subsaturado - (cp) {Float64} \n
\\- Rangos -                                        \n
\t 141 < Presión (psia) < 9515                         \n
\t 9.3 < Rs (pie^3/bls) < 2199                         \n
\t 15.3 < Gravedad del aceite (°API) < 59.5            \n
\t 0.511 < Gravedad específica del gas (adim) < 1.351  \n
\t 0.117 < Viscosidad del aceite (cp) < 148            \n
"""
function uosb_Vasquez_Beggs(Py::Union{Int64, Float64}, Pb::Union{Int64, Float64}, uob::Union{Int64, Float64})
    
    m = (2.6*Py^1.187)*exp(-11.513 - 8.98e-5*Py)
    uosb = uob*(Py/Pb)^m

    return round(uosb, digits=4)
end
#---------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T., Schmidt, Z.
""" Correlación de Kartoatmodjo, T., Schmidt, Z. - 1994 \n
uosb\\_Kartoatmodjo\\_Schmidt(Py, Pb, uob) \n
\\- Argumentos - \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Pb -> Presión de burbuja - (psia) {Int64, Float64} \n
\t uob -> Viscosidad del aceite saturado - (cp) {Float64} \n
\t uosb -> Viscosidad del aceite subsaturado - (cp) {Float64} \n
"""
function uosb_Kartoatmodjo_Schmidt(Py::Union{Int64, Float64}, Pb::Union{Int64, Float64}, uob::Union{Int64, Float64})::Float64
    
    uosb = 1.00081*uob + 1.127e-3*(Py - Pb)*(-65.17e-4*uob^1.8148 + 0.038*uob^1.59)

    return round(uosb, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------