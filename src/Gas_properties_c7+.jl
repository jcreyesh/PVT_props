# ---- Psc y Tsc de la fracción C7+ ----
#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kessler y Lee para la fracción C7+
""" Correlación de Kessler, M.G. y Lee, B.I. - 1976 \n
Kessler\\_Lee(M\\_c7m, Ge\\_c7m) \n
\\--------------------------------------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
M_c7m -> Peso Molecular de la fracción C7+ - (lbm/lbm-mol) {Int64, Float64} \n
Ge_c7m -> Gravedad específica de la fracción C7+ - (adim) {Int64, Float64} \n
Psc_c7m -> Presión pseducrítica de la fracción C7+ - (psia) {Float64} \n
Tsc_c7m -> Temperatura pseducrítica de la fracción C7+ - (°F) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------------------------------------------------------\n
"""
function PTsc_Kessler_Lee(M_c7m::Union{Int64,Float64}, Ge_c7m::Union{Int64, Float64})::Tuple{Float64, Float64}
    Tb = (4.5579*M_c7m^0.15178 * Ge_c7m^0.15427)^3

    Psc_c7m = exp(8.3634 - (0.0566/Ge_c7m) - (0.24244 + (2.2898/Ge_c7m) + (0.11857/Ge_c7m^2))*Tb*10^-3 +
    (1.4685 + (3.648/Ge_c7m) + (0.47227/Ge_c7m^2))*Tb^2*10^-7 -(0.42019 + (1.6977/Ge_c7m^2))*(Tb^3*10^-10))

    Tsc_c7m = 341.7 + 811*Ge_c7m + (0.4244 + 0.1174*Ge_c7m)*Tb + (0.4669 - 3.2623*Ge_c7m)*(10^5/Tb)
    return round(Psc_c7m, digits=4), round(Tsc_c7m, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Matthews, T.A., Roland C.H. y Katz, D.L. para la fracción C7+
""" Correlación de Matthews, T.A., Roland C.H. y Katz, D.L. - 1942 \n
Matthews\\_Roland(M\\_c7m, Ge\\_c7m) \n
\\-------------------------------------------------------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
M_c7m -> Peso Molecular de la fracción C7+ - (lbm/lbm-mol) {Int64, Float64} \n
Ge_c7m -> Gravedad específica de la fracción C7+ - (adim) {Int64, Float64} \n
Psc_c7m -> Presión pseducrítica de la fracción C7+ - (psia) {Float64} \n
Tsc_c7m -> Temperatura pseducrítica de la fracción C7+ - (°F) {Float64} \n
\\-------------------------------------------------------------------------------------------------------------------------------------------------------------\n
"""
function PTsc_Matthews_Roland(M_c7m::Union{Int64,Float64}, Ge_c7m::Union{Int64, Float64})::Tuple{Float64, Float64}
    Psc_c7m = 1188 - 431*log10(M_c7m - 61.1) + (2319 - 852*log10(M_c7m - 53.71))*(Ge_c7m - 0.8)
    Tsc_c7m = 608 + 364*log10(M_c7m - 71.2) + (2450*log10(M_c7m) - 3800)*log10(Ge_c7m)
    return round(Psc_c7m, digits=4), round(Tsc_c7m, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------