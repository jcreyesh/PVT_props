# ---- Presión de Burbuja ----
# Correlación de Standing, M.B.
""" Correlación de Standing, M.B. - 1947                               \n
pb\\_Standing(API, Ge\\_gas, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)         \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                       \n
API -> Gravedad °API - (grados) {Int64, Float64}                       \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}        \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}               \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}           \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                         \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                         \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                           \n
pb -> Presión de burbuja - (psia) {Float64}                            \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Standing(API, Ge\\_gas, Ty, Rs, fm\\_nhc)                             \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                       \n
API -> Gravedad °API - (grados) {Int64, Float64}                       \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}        \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}               \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}           \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}         \n
pb -> Presión de burbuja - (psia) {Float64}                            \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                           \n
130 < Pb - (psia) < 7,000                                              \n
100 < Ty - (°F) < 258                                                  \n
1.024 < Bo - (bls @yac/bls @std) < 2.15                                \n
20 < Rs - (pie^3/bls) < 1425                                           \n
16.5 < °API < 63.8                                                     \n
265 < Presión separador etapa-1 - (psia) < 465                         \n
Presión separador etapa-2 : 14.7 - (psia)                              \n
Temperatura del separador : 100 - (°F)                                 \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Standing(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2
    
    F = ((Rs/Ge_gas)^0.83) * (10^(91e-5*Ty - 12.5e-3*API))
    pb = 18.2*(F - 1.4)
    pbc = pb * C_co2 * C_h2s * C_n2

    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Standind, M.B.
function pb_Standing(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2
    
    F = ((Rs/Ge_gas)^0.83) * (10^(91e-5*Ty - 12.5e-3*API))
    pb = 18.2*(F - 1.4)
    pbc = pb * C_co2 * C_h2s * C_n2

    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Lasater, J.A.
""" Correlación de Lasater, J.A. - 1958                           \n
pb\\_Lasater(API, Ge\\_gas, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)     \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                  \n
API -> Gravedad °API - (grados) {Int64, Float64}                  \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}   \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}          \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}      \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                    \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                    \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                      \n
pb -> Presión de burbuja - (psia) {Float64}                       \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Lasater(API, Ge\\_gas, Ty, Rs, fm\\_nhc)                         \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                  \n
API -> Gravedad °API - (grados) {Int64, Float64}                  \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}   \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}          \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}      \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}    \n
pb -> Presión de burbuja - (psia) {Float64}                       \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                      \n
48 < Pb - (psia) < 5,780                                          \n
82 < Ty - (°F) < 272                                              \n
3 < Rs - (pie^3/bls) < 2905                                       \n
17.9 < °API < 51.1                                                \n
0.574 < Ge_gas < 1.233                                            \n
15 < Presión separador etapa-1 - (psia) < 605                     \n
34 < Temperatura del separador - (°F) < 106                       \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Lasater(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2
    
    Ge_oil = 141.5/(API + 131.5)
    Mo = 630 - 10*API
    yg = (Rs/379.3)/(Rs/379.3 + 350*Ge_oil/Mo)

    if yg <= 0.6
        pf =  0.679*exp(2.786*yg) - 0.323
    else
        pf = 8.26*yg^3.56 + 1.95
    end

    pb = pf * (Ty + 460)/Ge_gas
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Lasater, J.A.
function pb_Lasater(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3] 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2
    
    Ge_oil = 141.5/(API + 131.5)
    Mo = 630 - 10*API
    yg = (Rs/379.3)/(Rs/379.3 + 350*Ge_oil/Mo)

    if yg <= 0.6
        pf =  0.679*exp(2.786*yg) - 0.323
    else
        pf = 8.26*yg^3.56 + 1.95
    end

    pb = pf * (Ty + 460)/Ge_gas
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Vásquez, M.E. y Beggs, H.D.
""" Correlación de Vásquez, M.E. y Beggs, H.D. - 1980                                     \n
pb\\_Vasquez\\_Beggs\\_sp(API, Ge\\_gas, Psp, Tsp, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)      \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                                 \n
API -> Gravedad °API - (grados) {Int64, Float64}                                 \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}                  \n
Psp -> Presión en el separador - (psia) {Int64, Float64}                         \n
Tsp -> Temperatura en el separador - (°F) {Int64, Float64}                       \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                         \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                     \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                                   \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                                   \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                                     \n
pb -> Presión de burbuja - (psia) {Float64}                                      \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Vasquez\\Beggs\\_sp(API, Ge\\_gas, Psp, Tsp, Ty, Rs, fm\\_nhc)                  \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                                 \n
API -> Gravedad °API - (grados) {Int64, Float64}                                 \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}                  \n
Psp -> Presión en el separador - (psia) {Int64, Float64}                         \n
Tsp -> Temperatura en el separador - (°F) {Int64, Float64}                       \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                         \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                     \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}                   \n
pb -> Presión de burbuja - (psia) {Float64}                                      \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                                     \n
\\-------------------------------------------------------------------------------------------------------------------------------\n                                                                                 \n
API <= 30                                                                        \n
15 < Pb (psia) < 4,572                                                           \n
Temperatura promedio (°F) : 162                                                  \n
1.042 < Bo (bls @yac/bls @ std) < 1.545                                          \n
0 < Rs (pie^3/bls) < 831                                                         \n
5.3 < °API < 30                                                                  \n
0.511 < Ge_gas < 1.351                                                           \n
\\-------------------------------------------------------------------------------------------------------------------------------\n                                                                                 \n
API > 30                                                                         \n
15 < Pb (psia) < 6,055                                                           \n
Temperatura promedio (°F) : 180                                                  \n
1.028 < Bo (bls @yac/bls @ std) < 2.226                                          \n
0 < Rs (pie^3/bls) < 2,199                                                       \n
30.6 < °API < 59.5                                                               \n
0.53 < Ge_gas < 1.259                                                            \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Vasquez_Beggs_sp(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2
    
    Gec_gas = Ge_gas * (1 + 5.912e-5*API*Tsp*log10(Psp/114.7))

    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else 
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end

    pb = (Rs/(c1*Gec_gas*exp(c3*API/(Ty + 460))))^(1/c2)
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Vásquez, M.E. y Beggs, H.D.
function pb_Vasquez_Beggs_sp(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2
    
    Gec_gas = Ge_gas * (1 + 5.912e-5*API*Tsp*log10(Psp/114.7))

    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else 
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end

    pb = (Rs/(c1*Gec_gas*exp(c3*API/(Ty + 460))))^(1/c2)
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Vásquez, M.E. y Beggs, H.D.
""" Correlación de Vásquez, M.E. y Beggs, H.D. - 1980                                     \n
pb\\_Vasquez\\_Beggs\\(API, Ge\\_gas, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)      \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                                 \n
API -> Gravedad °API - (grados) {Int64, Float64}                                 \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}                  \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                         \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                     \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                                   \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                                   \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                                     \n
pb -> Presión de burbuja - (psia) {Float64}                                      \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Vasquez\\Beggs(API, Ge\\_gas, Ty, Rs, fm\\_nhc)                  \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                                 \n
API -> Gravedad °API - (grados) {Int64, Float64}                                 \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}                  \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                         \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                     \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}                   \n
pb -> Presión de burbuja - (psia) {Float64}                                      \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                                     \n
\\-------------------------------------------------------------------------------------------------------------------------------\n                                                                                 \n
API <= 30                                                                        \n
15 < Pb (psia) < 4,572                                                           \n
Temperatura promedio (°F) : 162                                                  \n
1.042 < Bo (bls @yac/bls @ std) < 1.545                                          \n
0 < Rs (pie^3/bls) < 831                                                         \n
5.3 < °API < 30                                                                  \n
0.511 < Ge_gas < 1.351                                                           \n
\\-------------------------------------------------------------------------------------------------------------------------------\n                                                                                 \n
API > 30                                                                         \n
15 < Pb (psia) < 6,055                                                           \n
Temperatura promedio (°F) : 180                                                  \n
1.028 < Bo (bls @yac/bls @ std) < 2.226                                          \n
0 < Rs (pie^3/bls) < 2,199                                                       \n
30.6 < °API < 59.5                                                               \n
0.53 < Ge_gas < 1.259                                                            \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Vasquez_Beggs(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2
    
    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else 
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end

    pb = (Rs/(c1*Ge_gas*exp(c3*API/(Ty + 460))))^(1/c2)
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Vásquez, M.E. y Beggs, H.D.
function pb_Vasquez_Beggs(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2
    
    if API <= 30
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else 
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    end

    pb = (Rs/(c1*Ge_gas*exp(c3*API/(Ty + 460))))^(1/c2)
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Glaso, O.
""" Correlación de Glaso, O. - 1980                                    \n
pb\\_Glaso(API, Ge\\_gas, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)            \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                       \n
API -> Gravedad °API - (grados) {Int64, Float64}                       \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}        \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}               \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}           \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                         \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                         \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                           \n
pb -> Presión de burbuja - (psia) {Float64}                            \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Glaso(API, Ge\\_gas, Ty, Rs, fm\\_nhc)                            \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                       \n
API -> Gravedad °API - (grados) {Int64, Float64}                       \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}        \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}               \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}           \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}         \n
pb -> Presión de burbuja - (psia) {Float64}                            \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                           \n
165 < Pb - (psia) < 7,142                                              \n
80 < Ty - (°F) < 280                                                   \n
1.025 < Bo - (bls @yac/bls @std) < 2.588                               \n
90 < Rs - (pie^3/bls) < 2637                                           \n
22.3 < °API < 48.1                                                     \n
0.65 < Ge_gas < 1.276                                                  \n
Presión separador:                                                     \n
etapa-1 - (psia) : 415                                                 \n
etapa-2 - (psia) : 15                                                  \n
Temperatura del separador - (°F) : 125                                 \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Glaso(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2

    F = (Rs/Ge_gas)^0.816 * ((Ty^0.172)/(API^0.989))
    pb = 10^(1.7669 + 1.7447*log10(F) - 0.30218*(log10(F))^2)
    pbc = pb * C_co2 * C_h2s * C_n2

    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Glaso, O.
function pb_Glaso(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2

    F = (Rs/Ge_gas)^0.816 * ((Ty^0.172)/(API^0.989))
    pb = 10^(1.7669 + 1.7447*log10(F) - 0.30218*(log10(F))^2)
    pbc = pb * C_co2 * C_h2s * C_n2

    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Total, C.F.P.
""" Correlación de Total, C.F.P. - 1980                                \n
pb\\_Total(API, Ge\\_gas, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)                    \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                       \n
API -> Gravedad °API - (grados) {Int64, Float64}                       \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}        \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}               \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}           \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                         \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                         \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                           \n
pb -> Presión de burbuja - (psia) {Float64}                            \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Total(API, Ge\\_gas, Ty, Rs, fm\\_nhc)                                \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                       \n
API -> Gravedad °API - (grados) {Int64, Float64}                       \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}        \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}               \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}           \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}         \n
pb -> Presión de burbuja - (psia) {Float64}                            \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Total_CFP(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2

    if API <= 10
        c1, c2, c3, c4 = 12.847, 0.9636, 9.93e-4, 34.17e-3
    elseif 10 < API <= 35
        c1, c2, c3, c4 = 25.2755, 0.7617, 8.35e-4, 11.292e-3
    elseif 35 < API <= 45
        c1, c2, c3, c4 = 216.4711, 0.6922, -4.27e-4, 23.14e-3
    end

    pb = c1 * (Rs/Ge_gas)^c2 * 10^(c3*Ty - c4*API)
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Total, C.F.P.
function pb_Total_CFP(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2

    if API <= 10
        c1, c2, c3, c4 = 12.847, 0.9636, 9.93e-4, 34.17e-3
    elseif 10 < API <= 35
        c1, c2, c3, c4 = 25.2755, 0.7617, 8.35e-4, 11.292e-3
    elseif 35 < API <= 45
        c1, c2, c3, c4 = 216.4711, 0.6922, -4.27e-4, 23.14e-3
    end

    pb = c1 * (Rs/Ge_gas)^c2 * 10^(c3*Ty - c4*API)
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Al-Marhun, M.A.
""" Correlación de Al-Marhoun, M.A. - 1947                               \n
pb\\_AlMarhoun(API, Ge\\_gas, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)          \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                         \n
API -> Gravedad °API - (grados) {Int64, Float64}                         \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}          \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                 \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}             \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                           \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                           \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                             \n
pb -> Presión de burbuja - (psia) {Float64}                              \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_AlMarhoun(API, Ge\\_gas, Ty, Rs, fm\\_nhc)                              \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                         \n
API -> Gravedad °API - (grados) {Int64, Float64}                         \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}          \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                 \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}             \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}           \n
pb -> Presión de burbuja - (psia) {Float64}                              \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                             \n
20 < Pb - (psia) < 3,573                                                 \n
74 < Ty - (°F) < 240                                                     \n
1.032 < Bo - (bls @yac/bls @std) < 1.997                                 \n
1.032 < Bt - (bls @yac/bls @std) < 6.982                                 \n
26 < Rs - (pie^3/bls) < 1602                                             \n
19.4 < °API < 44.6                                                       \n
0.752 < Ge_gas < 1.367                                                   \n
265 < Presión separador etapa-1 - (psia) < 465                           \n
0 < n2 en superficie - (% molar) < 3.89                                  \n
0 < co2 en superficie - (% molar) < 16.38                                \n
0 < h2s en superficie - (% molar) < 16.13                                \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_AlMarhoun(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2

    Ge_oil = 141.5/(API + 131.5)
    pb = (5.38088e-3*Rs^0.715082) * (Ge_gas)^-1.87787 * (Ge_oil)^3.1437 * (Ty + 460)^1.32657
    pbc = pb * C_co2 * C_h2s * C_n2

    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Marhoun, M.A.
function pb_AlMarhoun(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2

    Ge_oil = 141.5/(API + 131.5)
    pb = (5.38088e-3*Rs^0.715082) * (Ge_gas)^-1.87787 * (Ge_oil)^3.1437 * (Ty + 460)^1.32657
    pbc = pb * C_co2 * C_h2s * C_n2

    return pbc
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Dokla, M.E. y Osman,M.E.
""" Correlación de Dokla, M.E. y Osman,M.E. - 1992                       \n
pb\\_Dokla\\_Osman(API, Ge\\_gas, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)      \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                         \n
API -> Gravedad °API - (grados) {Int64, Float64}                         \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}          \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                 \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}             \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                           \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                           \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                             \n
pb -> Presión de burbuja - (psia) {Float64}                              \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Dokla\\_Osman(API, Ge\\_gas, Ty, Rs, fm\\_nhc)                          \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                         \n
API -> Gravedad °API - (grados) {Int64, Float64}                         \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}          \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                 \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}             \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}           \n
pb -> Presión de burbuja - (psia) {Float64}                              \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                             \n
590 < Pb - (psia) < 4,640                                                \n
190 < Ty - (°F) < 275                                                    \n
1.216 < Bo - (bls @yac/bls @std) < 2.493                                 \n
1.032 < Bt - (bls @yac/bls @std) < 6.982                                 \n
81 < Rs - (pie^3/bls) < 2,266                                            \n
0.8236 < Ge_oil < 0.886                                                  \n
0.789 < Ge_gas < 1.29                                                    \n
0.1 < n2 en superficie - (% molar) < 1.85                                \n
0.37 < co2 en superficie - (% molar) < 8.9                               \n
0 < h2s en superficie - (% molar) < 6.02                                 \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Dokla_Osman(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2

    Ge_oil = 141.5/(API + 131.5)
    pb = (0.836386e4*Rs^0.724047) * (Ge_gas)^-1.01049 * (Ge_oil)^0.107991 * (Ty + 460)^-0.952584
    pbc = pb * C_co2 * C_h2s * C_n2

    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Dokla, M.E. y Osman,M.E.
function pb_Dokla_Osman(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2

    Ge_oil = 141.5/(API + 131.5)
    pb = (0.836386e4*Rs^0.724047) * (Ge_gas)^-1.01049 * (Ge_oil)^0.107991 * (Ty + 460)^-0.952584
    pbc = pb * C_co2 * C_h2s * C_n2

    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Petrosky, G.E.Jr. y Farshad, F.F.
""" Correlación de Petrosky, G.E.Jr. y Farshad, F.F. - 1993                   \n
pb\\_Petrosky\\_Farshad(API, Ge\\_gas, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)      \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                              \n
API -> Gravedad °API - (grados) {Int64, Float64}                              \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}               \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                      \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                  \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                                \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                                \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                                  \n
pb -> Presión de burbuja - (psia) {Float64}                                   \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Petrosky\\_Farshad(API, Ge\\_gas, Ty, Rs, fm\\_nhc)                      \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                              \n
API -> Gravedad °API - (grados) {Int64, Float64}                              \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}               \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                      \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                  \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}                \n
pb -> Presión de burbuja - (psia) {Float64}                                   \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                                  \n
1,574 < Pb - (psia) < 6,523                                                   \n
1,700 < Py - (psia) < 10,692                                                  \n
114 < Ty - (°F) < 288                                                         \n
1.1178 < Bo - (bls @yac/bls @std) < 1.6229                                    \n
217 < Rs - (pie^3/bls) < 1,406                                                \n
16.3 < °API < 45                                                              \n
0.5781 < Ge_gas < 0.8519                                                      \n
3.507 < Co x 10^-6 (psi^-1) < 24.64                                           \n
0 < n2 en superficie - (% molar) < 3.72                                       \n
0 < co2 en superficie - (% molar) < 0.79                                      \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Petrosky_Farshad(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2

    Ge_oil = 141.5/(API + 131.5)
    F = (Rs^0.5774/Ge_gas^0.8439) * 10^(4.561e-5*Ty^1.3911 - 7.916e-4*API^1.541)
    pb = 112.727*(F - 12.34)
    pbc = pb * C_co2 * C_h2s * C_n2

    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Petrosky, G.E.Jr. y Farshad, F.F.
function pb_Petrosky_Farshad(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2

    Ge_oil = 141.5/(API + 131.5)
    F = (Rs^0.5774/Ge_gas^0.8439) * 10^(4.561e-5*Ty^1.3911 - 7.916e-4*API^1.541)
    pb = 112.727*(F - 12.34)
    pbc = pb * C_co2 * C_h2s * C_n2

    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T. y Schmidt, Z.
""" Correlación de Kartoatmodjo, T. y Schmidt, Z. - 1994                                         \n
pb\\_Kartoatmodjo\\_Schmidt\\_sp(API, Ge\\_gas, Psp, Tsp, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)      \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                                 \n
API -> Gravedad °API - (grados) {Int64, Float64}                                 \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}                  \n
Psp -> Presión en el separador - (psia) {Int64, Float64}                         \n
Tsp -> Temperatura en el separador - (°F) {Int64, Float64}                       \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                         \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                     \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                                   \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                                   \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                                     \n
pb -> Presión de burbuja - (psia) {Float64}                                      \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Kartoatmodjo\\_Schmidt\\_sp(API, Ge\\_gas, Psp, Tsp, Ty, Rs, fm\\_nhc)      \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                                 \n
API -> Gravedad °API - (grados) {Int64, Float64}                                 \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}                  \n
Psp -> Presión en el separador - (psia) {Int64, Float64}                         \n
Tsp -> Temperatura en el separador - (°F) {Int64, Float64}                       \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                         \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                     \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}                   \n
pb -> Presión de burbuja - (psia) {Float64}                                      \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                                     \n
14.7 < Pb (psia) < 6,054.7                                                       \n
75 < Ty - (°F) < 320                                                             \n
1.007 < Bo (bls @yac/bls @ std) < 2.144                                          \n
0 < Rs (pie^3/bls) < 2,890                                                       \n
14.4 < °API < 58.9                                                               \n
0.379 < Ge_gas < 1.709                                                           \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Kartoatmodjo_Schmidt_sp(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2
           
    
    Gec_gas = Ge_gas * (1 + 0.1595*API^0.4078*Tsp^-0.2466*log10(Psp/114.7))

    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end

    pb = (Rs/(c1 * Gec_gas^c2 * 10^(c3*API/(Ty + 460))))^c4
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T. y Schmidt, Z.
function pb_Kartoatmodjo_Schmidt_sp(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Psp::Union{Int64, Float64}, Tsp::Union{Int64,Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2
    
    Gec_gas = Ge_gas * (1 + 0.1595*API^0.4078*Tsp^-0.2466*log10(Psp/114.7))

    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end

    pb = (Rs/(c1 * Gec_gas^c2 * 10^(c3*API/(Ty + 460))))^c4
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T. y Schmidt, Z.
""" Correlación de Kartoatmodjo, T. y Schmidt, Z. - 1994                         \n
pb\\_Kartoatmodjo\\_Schmidt(API, Ge\\_gas, Ty, Rs, y\\_co2, y\\_h2s, y\\_n2)     \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                                 \n
API -> Gravedad °API - (grados) {Int64, Float64}                                 \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}                  \n
Psp -> Presión en el separador - (psia) {Int64, Float64}                         \n
Tsp -> Temperatura en el separador - (°F) {Int64, Float64}                       \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                         \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                     \n
y_co2 -> Fracción molar co2 - (adim) {Float64}                                   \n
y_h2s -> Fracción molar h2s - (adim) {Float64}                                   \n
y_n2 -> Fracción molar n2  -(adim) {Float64}                                     \n
pb -> Presión de burbuja - (psia) {Float64}                                      \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
pb\\_Kartoatmodjo\\_Schmidt(API, Ge\\_gas, Psp, Tsp, Ty, Rs, fm\\_nhc)      \n
\\--------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos -                                                                 \n
API -> Gravedad °API - (grados) {Int64, Float64}                                 \n
Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}                  \n
Psp -> Presión en el separador - (psia) {Int64, Float64}                         \n
Tsp -> Temperatura en el separador - (°F) {Int64, Float64}                       \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                         \n
Rs -> Relación de solubilidad - (pie^3/bls) {Int64, Float64}                     \n
fm_nhc -> Fracciones molares no-hc's - (adim) {AbstractVector}                   \n
pb -> Presión de burbuja - (psia) {Float64}                                      \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
\\- Rangos -                                                                     \n
14.7 < Pb (psia) < 6,054.7                                                       \n
75 < Ty - (°F) < 320                                                             \n
1.007 < Bo (bls @yac/bls @ std) < 2.144                                          \n
0 < Rs (pie^3/bls) < 2,890                                                       \n
14.4 < °API < 58.9                                                               \n
0.379 < Ge_gas < 1.709                                                           \n
\\-------------------------------------------------------------------------------------------------------------------------------\n
"""
function pb_Kartoatmodjo_Schmidt(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, y_co2::Union{Int64,Float64}, y_h2s::Union{Int64,Float64}, y_n2::Union{Int64,Float64})::Float64
    C_co2 = 1 - 693.8*y_co2*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*y_h2s + 0.019*(45 - API)*y_h2s^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*y_n2 
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*y_n2^2
           
    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end

    pb = (Rs/(c1 * Ge_gas^c2 * 10^(c3*API/(Ty + 460))))^c4
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kartoatmodjo, T. y Schmidt, Z.
function pb_Kartoatmodjo_Schmidt(API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64}, Ty::Union{Int64,Float64}, Rs::Union{Int64,Float64}, fm_nhc::AbstractVector)::Float64
    C_co2 = 1 - 693.8*fm_nhc[1]*Ty^-1.553
    C_h2s = 1 - (0.9035 + 0.0015*API)*fm_nhc[2] + 0.019*(45 - API)*fm_nhc[2]^2
    C_n2 = 1  + ((-2.65e-4*API + 5.5e-3)*Ty + (0.0931*API - 0.8295))*fm_nhc[3]
           + ((1.945e-11*API^4.699)*Ty + (0.027*API - 2.366))*fm_nhc[3]^2
           
    if API <= 30
        c1, c2, c3, c4 = 0.05958, 0.7972, 13.1405, 0.9986
    else 
        c1, c2, c3, c4 = 0.03150, 0.7587, 11.2895, 0.9143
    end

    pb = (Rs/(c1 * Ge_gas^c2 * 10^(c3*API/(Ty + 460))))^c4
    pbc = pb * C_co2 * C_h2s * C_n2
    
    return round(pbc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------